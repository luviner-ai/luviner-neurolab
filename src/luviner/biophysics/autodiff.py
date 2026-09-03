"""
Reverse-mode automatic differentiation over NumPy arrays.

Simulators of biophysical neurons run forwards. They tell you what a set
of conductances does, and if the answer is wrong you change the numbers
and run again. Every parameter in this package so far was found that
way -- by hand, one sweep at a time.

Differentiating the simulation removes that loop. With a gradient you
can ask the inverse question: what conductances would produce *this*
behaviour? Fit a cell to a recording, tune a circuit until it computes
something, or put a biophysical module inside a neural network and train
the whole thing end to end.

Conductance-based neurons are unusually well suited to this. A
leaky-integrate-and-fire neuron resets its voltage discontinuously when
it spikes, which breaks the gradient and forces spiking networks to use
approximate surrogate derivatives. Hodgkin-Huxley has no reset: the
voltage is a smooth solution of an ODE, every operation in the update is
differentiable, and the gradient that comes out is exact.

This is a tape-based reverse-mode implementation. Each `Var` records the
operation that produced it; `backward()` walks the tape in reverse
applying the chain rule. Arrays stay NumPy throughout, so the numerical
code reads the same as the forward-only version.

Only what the biophysics needs is implemented. It is not a general
autodiff framework and does not try to be.
"""

import numpy as np


class Var:
    """An array with a derivative.

    `value` is a NumPy array; `grad` accumulates d(output)/d(this) once
    `backward()` has run on a scalar downstream of it.
    """

    __slots__ = ('value', 'grad', '_backward', '_parents', 'requires_grad')

    def __init__(self, value, requires_grad=False, _parents=(), _backward=None):
        self.value = np.asarray(value, dtype=float)
        self.requires_grad = requires_grad or any(
            p.requires_grad for p in _parents)
        self.grad = None
        self._parents = _parents
        self._backward = _backward

    # ---- plumbing --------------------------------------------------------

    @property
    def shape(self):
        return self.value.shape

    def __repr__(self):
        return f"Var(shape={self.value.shape}, requires_grad={self.requires_grad})"

    def _accumulate(self, g):
        g = _unbroadcast(g, self.value.shape)
        self.grad = g if self.grad is None else self.grad + g

    def backward(self):
        """Propagate gradients back from this (scalar) node."""
        if self.value.size != 1:
            raise ValueError("backward() requires a scalar output")
        order, seen = [], set()

        def visit(v):
            if id(v) in seen:
                return
            seen.add(id(v))
            for p in v._parents:
                visit(p)
            order.append(v)

        visit(self)
        for v in order:
            v.grad = None
        self.grad = np.ones_like(self.value)
        for v in reversed(order):
            if v._backward is not None and v.grad is not None:
                v._backward(v.grad)
        return self

    def zero_grad(self):
        self.grad = None
        return self

    # ---- arithmetic ------------------------------------------------------

    def __add__(self, other):
        other = _as_var(other)
        out = Var(self.value + other.value, _parents=(self, other))
        out._backward = lambda g: (self._accumulate(g), other._accumulate(g))
        return out

    __radd__ = __add__

    def __neg__(self):
        out = Var(-self.value, _parents=(self,))
        out._backward = lambda g: self._accumulate(-g)
        return out

    def __sub__(self, other):
        return self + (-_as_var(other))

    def __rsub__(self, other):
        return _as_var(other) + (-self)

    def __mul__(self, other):
        other = _as_var(other)
        out = Var(self.value * other.value, _parents=(self, other))
        out._backward = lambda g: (self._accumulate(g * other.value),
                                   other._accumulate(g * self.value))
        return out

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = _as_var(other)
        out = Var(self.value / other.value, _parents=(self, other))

        def back(g):
            self._accumulate(g / other.value)
            other._accumulate(-g * self.value / (other.value ** 2))
        out._backward = back
        return out

    def __rtruediv__(self, other):
        return _as_var(other) / self

    def __pow__(self, power):
        p = float(power)
        out = Var(self.value ** p, _parents=(self,))
        out._backward = lambda g: self._accumulate(
            g * p * self.value ** (p - 1.0))
        return out

    def __getitem__(self, idx):
        out = Var(self.value[idx], _parents=(self,))

        def back(g):
            full = np.zeros_like(self.value)
            np.add.at(full, idx, g)
            self._accumulate(full)
        out._backward = back
        return out

    # ---- reductions ------------------------------------------------------

    def sum(self, axis=None):
        out = Var(self.value.sum(axis=axis), _parents=(self,))

        def back(g):
            g = np.asarray(g)
            if axis is not None:
                g = np.expand_dims(g, axis)
            self._accumulate(np.broadcast_to(g, self.value.shape).copy())
        out._backward = back
        return out

    def mean(self, axis=None):
        n = (self.value.size if axis is None else self.value.shape[axis])
        return self.sum(axis) / float(n)


# ---- helpers -------------------------------------------------------------

def _as_var(x):
    return x if isinstance(x, Var) else Var(x)


def _unbroadcast(g, shape):
    """Sum a gradient back down to `shape` after NumPy broadcasting."""
    g = np.asarray(g, dtype=float)
    while g.ndim > len(shape):
        g = g.sum(axis=0)
    for i, s in enumerate(shape):
        if s == 1 and g.shape[i] != 1:
            g = g.sum(axis=i, keepdims=True)
    return g.reshape(shape) if g.shape != shape else g


# ---- elementwise functions ----------------------------------------------

def exp(x):
    x = _as_var(x)
    out = Var(np.exp(x.value), _parents=(x,))
    out._backward = lambda g: x._accumulate(g * out.value)
    return out


def log(x):
    x = _as_var(x)
    out = Var(np.log(x.value), _parents=(x,))
    out._backward = lambda g: x._accumulate(g / x.value)
    return out


def tanh(x):
    x = _as_var(x)
    t = np.tanh(x.value)
    out = Var(t, _parents=(x,))
    out._backward = lambda g: x._accumulate(g * (1.0 - t * t))
    return out


def sigmoid(x):
    x = _as_var(x)
    v = x.value
    s = np.where(v >= 0, 1.0 / (1.0 + np.exp(-np.abs(v))),
                 np.exp(-np.abs(v)) / (1.0 + np.exp(-np.abs(v))))
    out = Var(s, _parents=(x,))
    out._backward = lambda g: x._accumulate(g * s * (1.0 - s))
    return out


def softplus(x):
    """log(1+exp(x)), evaluated stably. Useful to keep a parameter positive."""
    x = _as_var(x)
    v = x.value
    out = Var(np.log1p(np.exp(-np.abs(v))) + np.maximum(v, 0.0), _parents=(x,))
    s = np.where(v >= 0, 1.0 / (1.0 + np.exp(-np.abs(v))),
                 np.exp(-np.abs(v)) / (1.0 + np.exp(-np.abs(v))))
    out._backward = lambda g: x._accumulate(g * s)
    return out


def maximum(x, other):
    x = _as_var(x)
    other = _as_var(other)
    mask = x.value >= other.value
    out = Var(np.maximum(x.value, other.value), _parents=(x, other))

    def back(g):
        x._accumulate(g * mask)
        other._accumulate(g * (~mask))
    out._backward = back
    return out


def minimum(x, other):
    x = _as_var(x)
    other = _as_var(other)
    mask = x.value <= other.value
    out = Var(np.minimum(x.value, other.value), _parents=(x, other))

    def back(g):
        x._accumulate(g * mask)
        other._accumulate(g * (~mask))
    out._backward = back
    return out


def clip(x, lo, hi):
    return minimum(maximum(x, lo), hi)


def stack(vars_):
    """Stack a list of Vars along a new leading axis."""
    vars_ = [_as_var(v) for v in vars_]
    out = Var(np.stack([v.value for v in vars_]), _parents=tuple(vars_))

    def back(g):
        for i, v in enumerate(vars_):
            v._accumulate(g[i])
    out._backward = back
    return out


def where(condition, a, b):
    """Elementwise select. `condition` is a plain boolean array."""
    a, b = _as_var(a), _as_var(b)
    condition = np.asarray(condition, dtype=bool)
    out = Var(np.where(condition, a.value, b.value), _parents=(a, b))

    def back(g):
        a._accumulate(g * condition)
        b._accumulate(g * (~condition))
    out._backward = back
    return out


def x_over_1_minus_exp_neg_x(x):
    """The HH rate-function helper, with its removable singularity.

    Differentiable at x = 0: the Taylor branch has a derivative too, so
    the gradient does not blow up where the value is finite.
    """
    x = _as_var(x)
    v = x.value
    small = np.abs(v) < 1e-4
    safe = np.where(small, 1.0, v)
    val = np.where(small, 1.0 + v / 2.0 + v * v / 12.0,
                   safe / (1.0 - np.exp(-safe)))
    out = Var(val, _parents=(x,))

    def back(g):
        e = np.exp(-safe)
        denom = 1.0 - e
        d_exact = (denom - safe * e) / (denom ** 2)
        d_small = 0.5 + v / 6.0
        out_grad = np.where(small, d_small, d_exact)
        x._accumulate(g * out_grad)
    out._backward = back
    return out


# ---- optimisation --------------------------------------------------------

class Adam:
    """Adam optimizer over a list of Var parameters."""

    def __init__(self, params, lr=0.01, beta1=0.9, beta2=0.999, eps=1e-8):
        self.params = list(params)
        self.lr = lr
        self.beta1, self.beta2, self.eps = beta1, beta2, eps
        self.m = [np.zeros_like(p.value) for p in self.params]
        self.v = [np.zeros_like(p.value) for p in self.params]
        self.t = 0

    def step(self):
        self.t += 1
        for i, p in enumerate(self.params):
            if p.grad is None:
                continue
            g = p.grad
            if not np.all(np.isfinite(g)):
                # A diverged forward pass produces NaN gradients. Writing
                # them into the parameter poisons it permanently, since
                # NaN survives every subsequent clip and update. Skip.
                continue
            self.m[i] = self.beta1 * self.m[i] + (1 - self.beta1) * g
            self.v[i] = self.beta2 * self.v[i] + (1 - self.beta2) * g * g
            mhat = self.m[i] / (1 - self.beta1 ** self.t)
            vhat = self.v[i] / (1 - self.beta2 ** self.t)
            p.value = p.value - self.lr * mhat / (np.sqrt(vhat) + self.eps)

    def zero_grad(self):
        for p in self.params:
            p.grad = None


def numerical_gradient(f, params, eps=1e-6):
    """Central-difference gradient, for checking the analytic one."""
    grads = []
    for p in params:
        g = np.zeros_like(p.value)
        flat = p.value.reshape(-1)
        gflat = g.reshape(-1)
        for i in range(flat.size):
            old = flat[i]
            flat[i] = old + eps
            hi = float(np.asarray(f().value))
            flat[i] = old - eps
            lo = float(np.asarray(f().value))
            flat[i] = old
            gflat[i] = (hi - lo) / (2 * eps)
        grads.append(g)
    return grads
