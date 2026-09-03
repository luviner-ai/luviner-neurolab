"""
Short-term synaptic depression: the mechanism that terminates a burst.

The dish as built has no bursting regime. Measured on `dish-2` at global
inhibition, sigma 25, across the drive from 0.90 to 1.45: it fires
*exactly zero* spikes at 0.9000 and is 75-98% active at 0.9056, the first
drive above it. There is no graded middle, and the reason is structural
rather than a matter of finding the right drive. **A constant background
current at threshold has nothing that makes a burst stop.** Recurrent
excitation starts one, and nothing runs out, so it latches.

Cultures do not work that way, and the standard account of why is a
synapse whose transmitter is a finite pool: Tsodyks & Markram's
frequency-dependent synapse, in which a burst consumes the resource that
sustains it and therefore ends itself. Tsodyks, Uziel & Markram (2000)
show recurrent excitation plus depressing synapses producing
self-terminating population bursts at culture-like rates over a *wide*
range of drive, which is the property the constant-drive dish lacks.

The model, depressing case
--------------------------
Each synapse carries a fraction of available resource `x`, recovering to
1 with time constant `tau_rec`. A presynaptic spike releases a fraction
`U` of what is there and removes it:

    between spikes    dx/dt = (1 - x) / tau_rec
    on a spike        release = U * x,   then  x <- x - U * x

so the postsynaptic kick is `U * x` rather than 1. Facilitation is not
implemented: cortical pyramidal-to-pyramidal synapses of the type this
dish models are depressing, and adding a second process would be a
parameter nobody asked for.

Why this is checkable rather than merely plausible
--------------------------------------------------
The steady state under a regular train at frequency `f` has a closed
form, which is what makes the implementation falsifiable before any
network runs on it. With `dt = 1/f`:

    x_ss = (1 - E) / (1 - (1 - U) * E),      E = exp(-dt / tau_rec)
    A_ss / A_1 = x_ss

and in the high-frequency limit the *rate* of transmission saturates:

    A_ss * f  ->  1 / tau_rec

which is Tsodyks & Markram's limiting frequency result and depends on
`tau_rec` alone. `tests/test_shortterm.py` checks both against the
simulation, and the second is the one that matters: a depressing synapse
that reproduces its own steady-state formula but not the 1/f law has the
recovery in the wrong place.

References:
    Tsodyks MV, Markram H (1997). "The neural code between neocortical
        pyramidal neurons depends on neurotransmitter release
        probability." PNAS 94:719-723.
    Markram H, Wang Y, Tsodyks M (1998). "Differential signaling via the
        same axon of neocortical pyramidal neurons." PNAS 95:5323-5328.
    Tsodyks M, Uziel A, Markram H (2000). "Synchrony generation in
        recurrent networks with frequency-dependent synapses."
        J Neurosci 20:RC50.
"""

import numpy as np


class ShortTermDepression:
    """Tsodyks-Markram depression over a set of synapses.

    Parameters
    ----------
    U : release fraction per spike. 0.5 is the depressing
        pyramidal-to-pyramidal value.
    tau_rec : ms, recovery of the resource pool. 800 ms is the value
        Tsodyks-Markram fit to depressing cortical connections.
    mask : boolean over the population's synapse array, True where the
        synapse depresses. Synapses outside the mask are untouched and
        transmit at full strength, which is how E->E can depress while
        the inhibitory wiring does not.
    """

    def __init__(self, U=0.5, tau_rec=800.0, mask=None):
        """`U` and `tau_rec` may be scalars or per-synapse arrays.

        Per-synapse is not a generalisation for its own sake. A batched
        drive sweep can hold dishes at different drives because the drive
        is a per-neuron vector; it could not hold dishes at different
        depression parameters only because those were stored as scalars,
        which is an implementation choice and not physics. With arrays,
        a scan over `(U, tau_rec)` batches like everything else -- the
        alternative was four batches of four dishes paying the fixed
        interpreter cost four times, which the wall-clock guard refused
        at 16.2 minutes.
        """
        U_arr, tau_arr = np.asarray(U, dtype=float), np.asarray(tau_rec, dtype=float)
        if np.any(U_arr <= 0.0) or np.any(U_arr > 1.0):
            raise ValueError("U must be in (0, 1]")
        if np.any(tau_arr <= 0):
            raise ValueError("tau_rec must be positive")
        self.U = float(U) if U_arr.ndim == 0 else U_arr
        self.tau_rec = float(tau_rec) if tau_arr.ndim == 0 else tau_arr
        self.mask = None if mask is None else np.asarray(mask, dtype=bool)
        self.x = None

    def reset(self, n_syn, dt=0.05):
        """Full resources everywhere, which is the rested synapse.

        `dt` is taken here rather than in `release` so the recovery factor
        is computed once instead of per timestep; it is exact either way.
        """
        self.x = np.ones(int(n_syn))
        for name in ('U', 'tau_rec'):
            v = getattr(self, name)
            if not np.isscalar(v) and v.size not in (1, self.x.size):
                raise ValueError(f"{name} has {v.size} entries for "
                                 f"{self.x.size} synapses")
        self._recover = 1.0 - np.exp(-dt / self.tau_rec)
        if self.mask is not None and self.mask.size != self.x.size:
            raise ValueError(f"mask has {self.mask.size} entries for "
                             f"{self.x.size} synapses")
        return self

    def release(self, arrived, dt):
        """The kick per synapse for this timestep, and deplete.

        `arrived` is the population's boolean vector of spikes arriving
        *now*, after their axonal delay. Returns what would have been
        `arrived` in a model without depression: 1 where a spike arrives
        at a non-depressing synapse, `U * x` where it arrives at a
        depressing one, 0 elsewhere.

        Recovery is applied first and exactly -- `x <- 1 + (x-1)e^{-dt/tau}`
        rather than a forward Euler step -- because `tau_rec` is 800 ms
        against a 0.05 ms timestep and an approximation there would be a
        slow drift in the one quantity the burst duration depends on.
        """
        if self.x is None:
            self.reset(len(arrived))
        self.x += (1.0 - self.x) * self._recover
        kick = np.where(arrived, 1.0, 0.0)
        if self.mask is None:
            fire = arrived
        else:
            fire = arrived & self.mask
        if fire.any():
            u = self.U if np.isscalar(self.U) else self.U[fire]
            released = u * self.x[fire]
            kick[fire] = released
            self.x[fire] -= released
        return kick

    # ---- the closed forms the tests check against ------------------------

    def steady_state(self, freq_hz):
        """`x` at steady state under a regular train at `freq_hz`."""
        e = np.exp(-(1000.0 / np.asarray(freq_hz, dtype=float)) / self.tau_rec)
        return (1.0 - e) / (1.0 - (1.0 - self.U) * e)

    def limiting_rate(self):
        """`A_ss * f` as f -> infinity: 1 / tau_rec, in kHz-equivalent
        units of release per millisecond."""
        return 1.0 / self.tau_rec


def depressing(obj, U=0.5, tau_rec=800.0):
    """Attach depression to the E->E synapses of a dish or a batch.

    Works on anything carrying `pop`, `_slice_ampa` and `_slice_nmda` --
    which is `SiliconDish`, `StructuredDish` and `DishBatch` alike, so a
    batched drive sweep and a single dish get the identical preparation
    from the identical call.

    **Only the excitatory recurrence depresses.** The E->I, I->E and I->I
    conductances are fixed in this model and stay fixed: the mask is
    built from the two E->E `connect()` calls and nothing else, so the
    inhibitory wiring transmits at full strength however hard the network
    is firing. That asymmetry is the point -- a burst ends because the
    excitation that sustains it runs out, not because the brake gets
    stronger.
    """
    pop = obj.pop
    mask = np.zeros(pop.n_syn, dtype=bool)
    mask[obj._slice_ampa] = True
    mask[obj._slice_nmda] = True
    pop.stp = ShortTermDepression(U=U, tau_rec=tau_rec, mask=mask)
    pop.stp.reset(pop.n_syn)
    return obj
