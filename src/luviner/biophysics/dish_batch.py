"""
Ten dishes integrated as one population.

`SiliconDish` is vectorised over cells and not over dishes, and at 60
cells that is the wrong axis: a profile of one run spends 0.37 s of every
0.53 s inside `_rates` and `_derivatives`, which are pure element-wise
numpy on 60-element arrays. The cost is ufunc dispatch, not arithmetic.
N10's autopsy put the number on it -- a 17,000-stimulus ladder at ~1.5 s
per stimulus is 45 h of wall clock for an experiment whose arithmetic
would fit in minutes -- and named the fix: integrate the seeds together.

`DishBatch` holds K dishes as one `NeuronPopulation` of K x 60 cells with
block-diagonal wiring, one `PlasticConnections` over the union of their
plastic synapses, and one homeostatic budget whose per-neuron totals are
independent because the blocks are. Nothing about the equations changes.
Every array operation that used to run on 60 elements now runs on 600,
and the interpreter overhead is paid once instead of ten times.

Why it is exact, and how that is checked
----------------------------------------
The claim is not "close enough". Batched and separate must agree to
**1e-9**, and the reason that is the right acceptance test rather than a
courtesy is the one this document keeps arriving at: build the experiment
so a failure produces an impossible number. Two block-diagonal dishes
that share nothing can only disagree by a bug, and a bug in the offsets
or the wiring produces gross disagreement, never a plausible drift.

Three things have to hold for the equality to be bit-for-bit, and each is
a design constraint rather than a hope:

- **Element-wise arithmetic is position-independent.** Checked on this
  numpy: `exp`, `log`, multiplication and powers on the first 60 elements
  of a 600-array are bit-identical to the same operation on those 60
  alone.
- **Accumulation order per postsynaptic cell is preserved.** The
  synapses are concatenated *per `connect()` call* -- every dish's AMPA
  block, then every dish's NMDA block, and so on -- so the synapses
  reaching a given cell appear in the same relative order as they do in
  its own dish, and `np.add.at` sums them in that order.
- **Each dish keeps its own noise stream.** The background noise is
  drawn from the dish's own `RandomState`, in blocks of whole timesteps
  (`RandomState.randn(m, n)` is the same stream as `m` successive
  `randn(n)` calls), so the sequence a dish sees does not depend on how
  the batch chunks time.

What the batch may not do
-------------------------
All the dishes in a batch share one rule instance, one budget and one
`plastic` flag, because those are scalars applied to the whole weight
vector. A ladder whose arms differ in the *rule* batches over its seeds,
one arm at a time -- which is where the ten-fold win is anyway.

Divergent protocols
-------------------
Time is shared but control flow need not be. `run_programs` takes one
generator per dish, each yielding `(duration_ms, current)` segments, and
advances the batch in chunks bounded by the next segment boundary of any
dish. A protocol that stops when the network complies -- Shahaf-Marom,
where every dish reaches its criterion at a different stimulus -- runs
each dish through *its own* sequence of segments, of the same durations
and in the same order as it would alone.

A dish whose program ends before the others keeps being integrated,
because it is a block of a population that is still stepping. Everything
it recorded was recorded before that point, so nothing measured is
contaminated; but a quantity read off the batch *after* a dish has
finished is not that dish's final state. Programs snapshot what they
need when they end.
"""

import numpy as np

from .population import NeuronPopulation
from .plasticity import PlasticConnections, WeightNormalization

# Noise is drawn a block of timesteps at a time. Long enough that the
# per-call overhead disappears, short enough that a 2-second rest does
# not allocate a 200 MB array.
_NOISE_BLOCK = 1024


class DishBatch:
    """K dishes stepped as one population.

    Parameters
    ----------
    dishes : the dishes to batch. They must agree on everything that is
        a scalar of the whole population -- `dt`, the noise and drive
        levels, the NMDA ratio, the plasticity rule, the budget and the
        `plastic` flag -- and may differ in everything that is per-cell
        or per-synapse, which is to say in their wiring and their seed.

    The dishes are used for their wiring and their noise streams and are
    then no longer stepped themselves; their own `pop` and `conn` are
    left untouched, so a dish can still be run alone as a reference.
    """

    def __init__(self, dishes):
        self.dishes = list(dishes)
        if not self.dishes:
            raise ValueError("a batch needs at least one dish")
        d0 = self.dishes[0]
        for d in self.dishes[1:]:
            # `drive_level` is deliberately absent: the background drive is
            # a per-neuron vector that the batch concatenates per dish, so
            # dishes at different drives are as separable as dishes with
            # different wiring. `noise` is not -- it multiplies the whole
            # noise array as a scalar -- and the difference is exactly the
            # difference between a per-dish quantity and a batch one.
            for attr in ('n', 'n_exc', 'n_inh', 'n_sites', 'dt', 'noise',
                         'nmda_ratio', 'plastic', 'budget', 'tau_norm'):
                if getattr(d, attr) != getattr(d0, attr):
                    raise ValueError(
                        f"dishes disagree on {attr!r}: batching would apply "
                        f"one value to all of them")
            for attr in ('A_plus', 'A_minus', 'tau_plus', 'tau_minus',
                         'w_min', 'w_max', 'mu'):
                if getattr(d.rule, attr) != getattr(d0.rule, attr):
                    raise ValueError(
                        f"dishes disagree on rule.{attr}: one arm per batch")

        self.K = len(self.dishes)
        self.n = d0.n
        self.dt = d0.dt
        self.noise = d0.noise
        self.drives = [d.drive_level for d in self.dishes]
        self.nmda_ratio = d0.nmda_ratio
        self.plastic = d0.plastic
        self.rule = d0.rule

        self._build()

    # ---- construction ----------------------------------------------------

    def _build(self):
        """One population, block-diagonal, assembled call by call.

        The per-`connect()` arrays are read off the dishes' own
        populations rather than rebuilt from a synapse prototype, so the
        batch cannot drift from whatever the dish class wired -- which
        matters because `StructuredDish` wires its inhibition differently
        and the batch has to be blind to that.
        """
        n, K = self.n, self.K
        calls = len(self.dishes[0].pop._pre)
        for d in self.dishes:
            if len(d.pop._pre) != calls:
                raise ValueError("dishes were wired with a different number "
                                 "of connect() calls; the batch cannot align "
                                 "them")

        # g_M is per-neuron in the population, so dishes with different
        # slow-current conductances batch as naturally as different drives.
        g_M = np.concatenate([np.full(n, getattr(d, 'g_M', 0.0))
                              for d in self.dishes])
        pop = NeuronPopulation.cortical(K * n, g_M=g_M)
        start = 0
        for j in range(calls):
            cat = lambda name: np.concatenate(
                [getattr(d.pop, name)[j] for d in self.dishes])
            pre = np.concatenate([d.pop._pre[j] + k * n
                                  for k, d in enumerate(self.dishes)])
            post = np.concatenate([d.pop._post[j] + k * n
                                   for k, d in enumerate(self.dishes)])
            pop._pre.append(pre)
            pop._post.append(post)
            pop._g_max.append(cat('_g_max'))
            pop._E_rev.append(cat('_E_rev'))
            pop._tau_r.append(cat('_tau_r'))
            pop._tau_d.append(cat('_tau_d'))
            pop._delay.append(cat('_delay'))
            pop._mg.append(cat('_mg'))
            pop._slices.append(slice(start, start + pre.size))
            start += pre.size
        pop.reset()
        pop._build(self.dt)
        self.pop = pop
        self._slice_ampa = pop.connection_slice(0)
        self._slice_nmda = pop.connection_slice(1)

        # The plastic set is the E->E block of every dish, in dish order,
        # which is exactly how connect() call 0 was concatenated -- so the
        # writeback into g_max is a straight assignment.
        self.pre = np.concatenate([d.pre + k * n
                                   for k, d in enumerate(self.dishes)])
        self.post = np.concatenate([d.post + k * n
                                    for k, d in enumerate(self.dishes)])
        w0 = np.concatenate([d.conn.weights for d in self.dishes])
        sizes = [d.pre.size for d in self.dishes]
        edges = np.cumsum([0] + sizes)
        self.syn_slice = [slice(int(a), int(b))
                          for a, b in zip(edges[:-1], edges[1:])]

        d0 = self.dishes[0]
        homeo = (WeightNormalization(budget=d0.budget, tau=d0.tau_norm,
                                     w_max=self.rule.w_max)
                 if d0.budget is not None else None)
        self.conn = PlasticConnections(self.pre, self.post, w0,
                                       rule=self.rule, n_neurons=K * n,
                                       homeostasis=homeo)
        # The dishes are read for their wiring at whatever weights they
        # currently hold, so the conductances the population steps on are
        # taken from `conn` rather than from the connect()-time arrays.
        # For a fresh dish this is a no-op; for one that has been run it
        # is the difference between its weights and its initial ones.
        pop.g_max[self._slice_ampa] = self.conn.weights
        pop.g_max[self._slice_nmda] = self.conn.weights * self.nmda_ratio
        self._drive = np.concatenate([d._drive for d in self.dishes])

    # ---- state -----------------------------------------------------------

    @property
    def t(self):
        return self.pop.t

    def cells(self, k, cells=None):
        """Population indices of dish `k`'s cells (all of them by default)."""
        base = k * self.n
        if cells is None:
            return base + np.arange(self.n)
        return base + np.atleast_1d(np.asarray(cells, dtype=int))

    def weights(self, k):
        """Dish `k`'s plastic conductances, in its own synapse order."""
        return self.conn.weights[self.syn_slice[k]]

    def spike_count(self, k, cells, t_lo, t_hi):
        total = 0
        for c in self.cells(k, cells):
            s = self.pop.spike_times[int(c)]
            if s:
                a = np.asarray(s)
                total += int(np.count_nonzero((a > t_lo) & (a <= t_hi)))
        return total

    def response(self, k, site, t0, lo=3.0, hi=25.0):
        """As `SiliconDish.response`, for one dish of the batch."""
        return self.spike_count(k, self.dishes[k].sites[int(site)],
                                t0 + lo, t0 + hi)

    def pathway_weight(self, k, pre_site, post_site):
        """As `SiliconDish.pathway_weight`, for one dish of the batch."""
        d = self.dishes[k]
        m = ((d.site_of[d.pre] == int(pre_site))
             & (d.site_of[d.post] == int(post_site)))
        w = self.weights(k)
        return float(w[m].mean()) if m.any() else 0.0

    def pathway_matrix(self, k):
        n = self.dishes[k].n_sites
        return np.array([[self.pathway_weight(k, a, b) for b in range(n)]
                         for a in range(n)])

    # ---- running ---------------------------------------------------------

    def _advance(self, n_steps, extra, learn):
        """`n_steps` timesteps with a constant injected current.

        The body is `SiliconDish.run`'s, term for term and in the same
        order -- `drive + randn * noise`, then `+ extra` -- because the
        equality with a separate run is checked to 1e-9 and float
        addition is not associative.
        """
        if n_steps <= 0:
            return self
        pop, conn, dt = self.pop, self.conn, self.dt
        sa, sn, ratio = self._slice_ampa, self._slice_nmda, self.nmda_ratio
        write = learn or conn.homeostasis is not None
        scale, drive = self.noise, self._drive
        done = 0
        while done < n_steps:
            m = min(_NOISE_BLOCK, n_steps - done)
            noise = np.concatenate(
                [d._noise_rng.randn(m, self.n) for d in self.dishes], axis=1)
            for i in range(m):
                I = drive + noise[i] * scale
                if extra is not None:
                    I = I + extra
                pop.step(I, dt)
                conn.step(pop.last_fired, dt, learn=learn)
                if write:
                    pop.g_max[sa] = conn.weights
                    pop.g_max[sn] = conn.weights * ratio
            done += m
        return self

    def run(self, duration, extra=None, learn=True):
        """Advance every dish by `duration` ms with one shared schedule.

        `extra` is `None`, a `(K, n)` array, or a flat `(K*n,)` one.
        """
        learn = bool(learn) and self.plastic
        if extra is not None:
            extra = np.asarray(extra, dtype=float).reshape(self.K * self.n)
        return self._advance(int(round(duration / self.dt)), extra, learn)

    def pulse(self, sites, amplitude=40.0, width=1.0, follow=200.0,
              learn=True):
        """One stimulus per dish, then `follow` ms of free running.

        `sites` is one site index per dish; `None` in a slot means that
        dish receives no stimulus while the others do. Returns the onset
        time, which is shared -- that is what makes the batch a batch.
        """
        t0 = self.pop.t
        I = self.stimulus(sites, amplitude)
        self.run(width, extra=I, learn=learn)
        self.run(follow, learn=learn)
        return t0

    def stimulus(self, sites, amplitude=40.0):
        """The flat injected-current vector for one site per dish."""
        I = np.zeros(self.K * self.n)
        for k, s in enumerate(np.atleast_1d(sites)):
            if s is None:
                continue
            I[self.cells(k, self.dishes[k].sites[int(s)])] = float(amplitude)
        return I

    def rest(self, duration, learn=True):
        return self.run(duration, learn=learn)

    # ---- divergent protocols ---------------------------------------------

    def run_programs(self, programs, learn=True):
        """Advance the batch while each dish follows its own segments.

        `programs` is one generator per dish, each yielding
        `(duration_ms, current)` where `current` is `None` or an `(n,)`
        vector for that dish's cells. The batch advances in chunks
        bounded by the next segment boundary of any dish, so every dish
        sees exactly the segment durations its own protocol asked for.

        Durations are rounded to whole timesteps the way
        `SiliconDish.run` rounds them, and a dish whose program has ended
        keeps being integrated with no injected current -- see the module
        docstring on what that does and does not contaminate.
        """
        learn = bool(learn) and self.plastic
        n, K = self.n, self.K
        pending = [None] * K          # [steps remaining, current or None]
        alive = [True] * K
        while True:
            for k in range(K):
                while alive[k] and pending[k] is None:
                    try:
                        dur, cur = next(programs[k])
                    except StopIteration:
                        alive[k] = False
                        break
                    steps = int(round(dur / self.dt))
                    if steps > 0:
                        pending[k] = [steps, cur]
            live = [p[0] for p in pending if p is not None]
            if not live:
                return self
            m = min(live)
            extra = None
            for k in range(K):
                if pending[k] is not None and pending[k][1] is not None:
                    if extra is None:
                        extra = np.zeros(K * n)
                    extra[k * n:(k + 1) * n] = pending[k][1]
            self._advance(m, extra, learn)
            for k in range(K):
                if pending[k] is not None:
                    pending[k][0] -= m
                    if pending[k][0] <= 0:
                        pending[k] = None


def shahaf_marom_batch(batch, stim_site=0, target_site=3, window=(3.0, 25.0),
                       isi=200.0, rest=2000.0, max_stim=30, n_trials=8,
                       contingent=True, yoked=None, learn=True, overshoot=0,
                       amplitude=40.0):
    """`shahaf_marom` on a batch, each dish stopping when *it* complies.

    One result dict per dish, with the same keys and the same meaning as
    the single-dish protocol. The per-stimulus checkpoints are taken by
    each dish's own program at its own stimuli, so they are that dish's
    trajectory and not a batch average.

    The one thing that is *not* shared is where a trial ends, which is
    the whole content of the experiment: the criterion is met at a
    different stimulus in every dish, and the batch exists precisely so
    that difference costs nothing.
    """
    K = batch.K
    if not contingent and yoked is None:
        raise ValueError("the yoked control needs the counts to replay")
    n_sites = batch.dishes[0].n_sites
    others = [s for s in range(n_sites)
              if s not in (stim_site, target_site)]
    out = [{'counts': [], 'reached': [], 'target_response': [],
            'other_response': [], 'pathway': [], 'pathway_others': [],
            'final_t': None, 'final_matrix': None, 'final_weights': None}
           for _ in range(K)]

    def program(k):
        rec = out[k]
        I = np.zeros(batch.n)
        I[batch.dishes[k].sites[int(stim_site)]] = float(amplitude)
        y = None if contingent else yoked[k]
        for trial in range(n_trials):
            quota = max_stim if contingent else int(y[trial % len(y)])
            count, reached = 0, None
            while count < max_stim:
                t0 = batch.t
                yield 1.0, I
                yield isi - 1.0, None
                count += 1
                r = batch.response(k, target_site, t0, *window)
                o = np.mean([batch.response(k, s, t0, *window)
                             for s in others])
                rec['target_response'].append(r)
                rec['other_response'].append(float(o))
                rec['pathway'].append(
                    batch.pathway_weight(k, stim_site, target_site))
                rec['pathway_others'].append(float(np.mean(
                    [batch.pathway_weight(k, stim_site, s) for s in others])))
                if r >= 1 and reached is None:
                    reached = count
                if contingent:
                    if reached is not None and count >= reached + overshoot:
                        break
                elif count >= quota:
                    break
            rec['counts'].append(count)
            rec['reached'].append(reached)
            yield rest, None
        # THE SNAPSHOT, and it belongs here rather than after
        # `run_programs` returns. This code runs when *this* dish's
        # program ends, which is the only moment its weights mean what
        # the protocol says they mean: a dish that finishes before the
        # others keeps being integrated, and reading its conductances
        # afterwards measures the extra run as well. Doing exactly that
        # cost a re-run of two experiments -- the module docstring warned
        # about it and the caller did it anyway, so the snapshot is now
        # taken here where a caller cannot get it wrong.
        rec['final_t'] = batch.t
        rec['final_matrix'] = batch.pathway_matrix(k).tolist()
        rec['final_weights'] = batch.weights(k).tolist()

    batch.run_programs([program(k) for k in range(K)], learn=learn)
    return [{'counts': np.array(r['counts'], dtype=float),
             'reached': r['reached'],
             'target_response': np.array(r['target_response'], dtype=float),
             'other_response': np.array(r['other_response'], dtype=float),
             'pathway': np.array(r['pathway'], dtype=float),
             'pathway_others': np.array(r['pathway_others'], dtype=float),
             # Taken when this dish's own program ended, not when the
             # batch did. See the note in `program`.
             'final_t': r['final_t'],
             'final_matrix': np.array(r['final_matrix'], dtype=float),
             'final_weights': np.array(r['final_weights'], dtype=float)}
            for r in out]


def site_specificity_batch(batch, stim=0, target=None, amplitude=40.0,
                           follow=200.0, window=(25.0, 60.0), n_trials=10,
                           settle=100.0):
    """`site_specificity` on a batch, one result dict per dish.

    Lockstep: every dish is settled together, then stimulated at the same
    site at the same times, because nothing here depends on what any dish
    did. The dishes may still differ in everything that is wiring -- the
    sweep puts the flat control, six length scales and two inhibition
    modes in one batch, which is where the throughput is.

    Scored with `score_sites`, the same function the single-dish version
    calls, rather than a copy of it. The locality fraction has been
    redefined once already; two implementations of it is how the sweep and
    the check that guards it end up disagreeing with neither being wrong.
    """
    from .dish_structured import score_sites

    n = batch.dishes[0].n_sites
    if target is None:
        target = (stim + 1) % n
    if target == stim:
        raise ValueError("the stimulated site cannot also be the target")

    batch.run(settle, learn=False)
    per_site = [{s: [] for s in range(n)} for _ in range(batch.K)]
    for _ in range(int(n_trials)):
        t0 = batch.pulse([stim] * batch.K, amplitude=amplitude,
                         follow=follow, learn=False)
        for k in range(batch.K):
            for s in range(n):
                per_site[k][s].append(batch.response(k, s, t0, *window))
    return [score_sites({s: float(np.mean(v)) for s, v in p.items()},
                        stim, target, n, window)
            for p in per_site]


def observed_lags_of(batch, k, t_lo, t_hi, max_lag=250.0):
    """`observed_lags` for one dish of a batch.

    The dish's cells are a block of the population, so its plastic
    synapses are its own `pre`/`post` shifted onto that block. Everything
    else -- including why `max_lag` defaults past the knee at 250 ms --
    is `lags_between`'s, because it is the same function.
    """
    from .dish import lags_between

    d = batch.dishes[k]
    off = k * batch.n
    return lags_between(batch.pop.spike_times, d.pre + off, d.post + off,
                        t_lo, t_hi, max_lag)
