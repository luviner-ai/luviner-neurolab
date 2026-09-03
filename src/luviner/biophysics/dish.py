"""
The silicon dish: a cultured network on a multi-electrode array, simulated.

A living culture on an MEA is offered as compute: electrodes in, spikes
out, the whole thing embodied in a task and left to learn. Shahaf & Marom
(2001) trained one by *stopping* the stimulation when it did what they
wanted; Kagan et al. (2022) put ~800k cells in Pong and reported longer
rallies within minutes. Both have learning curves. Neither has a
mechanism -- in a dish no synapse is readable and no run is repeatable,
so the question "which synapses changed, and under which rule?" cannot be
asked of the preparation that produced the curve.

This module is the same object with the same interface, built from
conductance-based cells with a measured plasticity rule, where every
synapse is readable and every run is repeatable. The point is not to
claim that silicon beats the wet rack. It is that a digital twin can be
asked the one question the original cannot answer about itself.

The preparation
---------------
`SiliconDish` holds an excitatory/inhibitory network of cortical cells
whose recurrent E->E synapses follow STDP under a total-input budget.
The excitatory cells are divided into `n_sites` groups, and a group is
what an electrode is: current injected into it is a stimulus, spikes
recorded from it are that electrode's channel. There is no spatial
structure beyond that, because there is none in the connectivity either;
sites are interchangeable by construction, which is what makes a
site-specific result mean something.

The dish is run *silent between stimuli*: the background drive sits below
rheobase, so the cells fire only when something makes them. That is a
deliberate choice of preparation, in the same spirit as blocking sodium
to measure an EPSP. It costs the spontaneous network bursts real cultures
show, and it buys a baseline of exactly zero spikes, against which an
evoked response cannot be confused with ongoing activity. Every
measurement here is anchored to the stimulus onset, never to the length
of the run.

The knob that decides the mechanism
-----------------------------------
`rule_with_integral` builds an STDP rule with a prescribed window
integral, holding `A_plus` and both time constants fixed and solving for
`A_minus`. The integral -- not the peak -- is what governs the change
when pre- and post-spike times land on both sides of zero, which is the
situation in a network driven by a stimulus rather than by a paired
protocol. Setting it to zero disables the *net* effect of the rule while
leaving individual pairs as strong as they were, so it is a mechanism
ablation rather than a parameter sweep: if an effect survives a zero
integral, STDP is not what produced it.

References:
    Shahaf G, Marom S (2001). "Learning in networks of cortical neurons."
        J Neurosci 21:8782-8788.
    Kagan BJ, Kitchen AC, Tran NT, et al. (2022). "In vitro neurons learn
        and exhibit sentience when embodied in a simulated game-world."
        Neuron 110:3952-3969.
    Cai H, Ao Z, Tian C, et al. (2023). "Brain organoid reservoir
        computing for artificial intelligence." Nat Electron 6:1032-1039.
    Marom S, Shahaf G (2002). "Development, learning and memory in large
        random networks of cortical neurons." Q Rev Biophys 35:63-87.
"""

import numpy as np

from .population import NeuronPopulation
from .synapse import Synapse
from .plasticity import STDP, PlasticConnections, WeightNormalization


def rule_with_integral(integral, A_plus=0.012, tau_plus=16.8,
                       tau_minus=33.7, w_min=0.0, w_max=0.06, mu=0.0):
    """An STDP rule whose window integral is exactly `integral`.

    `A_minus` is solved for; everything else is held. Three settings
    matter and they are a mechanism ablation, not a sweep:

        integral < 0   the rule dissolves co-active groups (the standard
                       choice, Song-Miller-Abbott 2000)
        integral = 0   individual pairs still potentiate and depress at
                       full strength, but nothing accumulates
        integral > 0   co-active groups bind

    Raises if the requested integral needs a negative `A_minus`, which
    would be a rule that potentiates in both directions.
    """
    A_minus = (A_plus * tau_plus - float(integral)) / tau_minus
    if A_minus < 0:
        raise ValueError(
            f"integral {integral} needs A_minus={A_minus:.4f} < 0 at "
            f"A_plus={A_plus}; raise tau_minus or lower the integral")
    return STDP(A_plus=A_plus, A_minus=A_minus, tau_plus=tau_plus,
                tau_minus=tau_minus, w_min=w_min, w_max=w_max, mu=mu)


def lag_kernels(lags, tau_plus=16.8, tau_minus=33.7):
    """The two numbers an observed lag distribution reduces to.

    `STDP.expected_change` on a set of lags is `A_plus * P - A_minus * M`,
    where P and M are the means of the two exponentials over the lags that
    fall on their own side of zero. Everything a rule can do to a given
    network state is contained in that pair.
    """
    lags = np.asarray(lags, dtype=float)
    if lags.size == 0:
        return 0.0, 0.0
    pos = lags >= 0
    P = float(np.sum(np.exp(-lags[pos] / tau_plus)) / lags.size)
    M = float(np.sum(np.exp(lags[~pos] / tau_minus)) / lags.size)
    return P, M


def rule_with_expected_change(target, lags, A_plus=0.012, tau_plus=16.8,
                              tau_minus=33.7, w_min=0.0, w_max=0.06, mu=0.0):
    """An STDP rule whose mean change per pair, on *these* lags, is `target`.

    `rule_with_integral(0)` is the obvious ablation and it is the wrong
    one. The integral is the net change per pair when the lags are spread
    evenly over the window; the lags a stimulated network actually
    produces are concentrated near zero, where the window is
    potentiation-dominated, so a rule with zero area still potentiates.
    Zeroing the *area* is a statement about the rule; zeroing the
    *expected change* is a statement about the rule applied to this
    network, and only the second is a mechanism ablation.

    Neutrality has a scale-free form. `A_plus * P - A_minus * M = 0` is
    `A_minus / A_plus = P / M`: the *ratio* of the two learning rates is
    what the network fixes, not either rate and not the area. Written as
    an integral it is `A_plus * (tau_plus - tau_minus * P / M)`, which is
    why the neutral integral moves when `A_plus` does while the ratio
    does not.

    The integral of the rule this returns is a prediction in its own
    right: it says how far from zero the neutral point sits for a network
    whose lags look like these.
    """
    P, M = lag_kernels(lags, tau_plus, tau_minus)
    if M <= 0:
        raise ValueError("no anticausal lags: the depression side is "
                         "unconstrained on this distribution")
    A_minus = (A_plus * P - float(target)) / M
    if A_minus < 0:
        raise ValueError(
            f"target {target} needs A_minus={A_minus:.4f} < 0")
    return STDP(A_plus=A_plus, A_minus=A_minus, tau_plus=tau_plus,
                tau_minus=tau_minus, w_min=w_min, w_max=w_max, mu=mu)


class SiliconDish:
    """An E/I culture on a multi-electrode array.

    Parameters
    ----------
    n_exc, n_inh : cell counts. The excitatory cells are split evenly
        into `n_sites` electrode groups.
    connectivity : probability of an E->E connection. E<->I is all-to-all
        and fixed; only the excitatory recurrence is plastic, which is
        where the measured plasticity lives.
    w_rec : initial E->E conductance, drawn with 20% jitter.
    drive, noise : background current, uA/cm^2. `drive` is below rheobase
        (1.0 for this cell) on purpose -- see the module docstring.
    rule : an STDP instance. Defaults to a *net-depressing* window, which
        is the standard setting and the one the DishBrain prediction says
        the effect requires.
    budget : total incoming E->E conductance held per cell by
        `WeightNormalization`. None disables homeostasis.
    plastic : False switches off the *Hebbian* rule. The homeostatic
        budget, if one is set, keeps running -- it is a property of the
        cell, not of the training protocol, and leaving it in both arms
        is what makes the control differ from the experiment in exactly
        one thing. With `budget=None` the weights are then frozen
        outright.
    """

    RHEOBASE = 1.0          # measured for this cell in test_cortical

    def __init__(self, n_exc=48, n_inh=12, n_sites=6, connectivity=0.15,
                 w_rec=0.015, nmda_ratio=2.0, drive=0.90, noise=0.25,
                 rule=None, budget=None, tau_norm=200.0, plastic=True,
                 dt=0.05, seed=0, g_M=0.0):
        self.n_exc = int(n_exc)
        self.n_inh = int(n_inh)
        self.n = self.n_exc + self.n_inh
        if self.n_exc % n_sites:
            raise ValueError("n_exc must divide evenly into n_sites")
        self.n_sites = int(n_sites)
        self.per_site = self.n_exc // self.n_sites
        self.sites = [np.arange(s * self.per_site, (s + 1) * self.per_site)
                      for s in range(self.n_sites)]
        self.site_of = np.repeat(np.arange(self.n_sites), self.per_site)
        self.dt = float(dt)
        self.drive_level = float(drive)
        self.noise = float(noise)
        self.plastic = bool(plastic)
        self.seed = int(seed)
        # Muscarinic slow K+ conductance, per cell. Zero is the dish as it
        # was. NOTE THE SCALE: these cells run the CORTICAL set, g_L = 0.1,
        # and `channels.PotassiumM` defaults to g_max = 0.5 -- five times
        # the leak. Measured on single dish cells, that default drops the
        # firing rate to 0.7 Hz at I = 1.5 and 6.7 Hz at I = 3.0, i.e. it
        # silences them. The usable range here is 0.02-0.2.
        self.g_M = float(g_M)

        self.exc = np.arange(self.n_exc)
        self.inh = self.n_exc + np.arange(self.n_inh)

        rng = np.random.RandomState(seed)
        mask = rng.rand(self.n_exc, self.n_exc) < connectivity
        np.fill_diagonal(mask, False)
        self.pre, self.post = (a.copy() for a in np.nonzero(mask))
        w0 = w_rec * (1.0 + 0.2 * rng.randn(self.pre.size))
        self.w0 = np.clip(w0, 0.0, None)

        self.rule = rule if rule is not None else rule_with_integral(-0.134)
        self.nmda_ratio = float(nmda_ratio)
        self.budget = budget
        self.tau_norm = float(tau_norm)
        self._noise_rng = np.random.RandomState(seed + 1000)

        self._build(self.w0)

    # ---- construction ----------------------------------------------------

    def _build(self, weights):
        pop = NeuronPopulation.cortical(
            self.n, g_M=np.full(self.n, getattr(self, 'g_M', 0.0)))
        pop.connect(self.pre, self.post,
                    Synapse.ampa(g_max=1.0, delay=1.0), weights=weights)
        pop.connect(self.pre, self.post,
                    Synapse.nmda(g_max=1.0, delay=1.0),
                    weights=weights * self.nmda_ratio)
        for a, b, g, kind in ((self.exc, self.inh, 0.02, 'ampa'),
                              (self.inh, self.exc, 0.25, 'gaba'),
                              (self.inh, self.inh, 0.02, 'gaba')):
            p, q = (x.ravel() for x in np.meshgrid(a, b, indexing='ij'))
            syn = (Synapse.ampa(g_max=1.0, delay=1.0) if kind == 'ampa'
                   else Synapse.gaba_a(g_max=1.0, delay=1.0, E_rev=-75.0))
            pop.connect(p, q, syn, weights=np.full(p.size, g))
        pop.reset()
        pop._build(self.dt)
        self.pop = pop
        self._slice_ampa = pop.connection_slice(0)
        self._slice_nmda = pop.connection_slice(1)

        homeo = (WeightNormalization(budget=self.budget, tau=self.tau_norm,
                                     w_max=self.rule.w_max)
                 if self.budget is not None else None)
        self.conn = PlasticConnections(self.pre, self.post,
                                       np.asarray(weights, dtype=float),
                                       rule=self.rule, n_neurons=self.n,
                                       homeostasis=homeo)
        self._drive = np.zeros(self.n)
        self._drive[:self.n_exc] = self.drive_level
        self._drive[self.n_exc:] = self.drive_level * 0.75

    def fresh_copy(self, weights=None, seed=None):
        """A new dish with the same wiring and the given weights.

        Probing the network that has just been trained scores its
        leftover activity rather than its weights. Every claim about what
        the weights hold is made on one of these.
        """
        twin = SiliconDish.__new__(SiliconDish)
        twin.__dict__.update({k: v for k, v in self.__dict__.items()
                              if k not in ('pop', 'conn', '_drive',
                                           '_slice_ampa', '_slice_nmda',
                                           '_noise_rng')})
        twin._noise_rng = np.random.RandomState(
            (self.seed if seed is None else seed) + 2000)
        twin._build(self.w0 if weights is None else np.asarray(weights))
        return twin

    # ---- state -----------------------------------------------------------

    @property
    def t(self):
        return self.pop.t

    @property
    def weights(self):
        return self.conn.weights

    def spike_count(self, cells, t_lo, t_hi):
        """Spikes from `cells` in the absolute interval (t_lo, t_hi]."""
        total = 0
        for c in np.atleast_1d(cells):
            s = self.pop.spike_times[int(c)]
            if s:
                a = np.asarray(s)
                total += int(np.count_nonzero((a > t_lo) & (a <= t_hi)))
        return total

    # ---- running ---------------------------------------------------------

    def run(self, duration, extra=None, learn=True):
        """Advance the dish by `duration` ms.

        `extra` is an injected current on top of the background, either a
        constant vector or a callable of the step index.
        """
        learn = bool(learn) and self.plastic
        n_steps = int(round(duration / self.dt))
        at = extra if callable(extra) else (lambda i: extra)
        pop, conn, dt = self.pop, self.conn, self.dt
        sa, sn, ratio = self._slice_ampa, self._slice_nmda, self.nmda_ratio
        for i in range(n_steps):
            I = self._drive + self._noise_rng.randn(self.n) * self.noise
            e = at(i)
            if e is not None:
                I = I + e
            pop.step(I, dt)
            conn.step(pop.last_fired, dt, learn=learn)
            if learn or conn.homeostasis is not None:
                pop.g_max[sa] = conn.weights
                pop.g_max[sn] = conn.weights * ratio
        return self

    def pulse(self, site, amplitude=40.0, width=1.0, follow=200.0,
              learn=True):
        """One stimulus at `site`, then `follow` ms of free running.

        Returns the absolute time of stimulus onset. Responses are
        measured against that time, never against the run's length: a
        window anchored to the run is an average in disguise and its value
        moves with the run's duration.
        """
        t0 = self.pop.t
        I = np.zeros(self.n)
        I[self.sites[int(site)]] = float(amplitude)
        self.run(width, extra=I, learn=learn)
        self.run(follow, learn=learn)
        return t0

    def rest(self, duration, learn=True):
        """No stimulus. The dish is silent, so this is where nothing
        happens -- which is the point of the Shahaf-Marom reward."""
        return self.run(duration, learn=learn)

    # ---- readout ---------------------------------------------------------

    def response(self, site, t0, lo=3.0, hi=25.0):
        """Spikes at `site` in the window [t0+lo, t0+hi] after a stimulus.

        The default window is the short-latency, pathway-specific one. The
        delayed network burst that follows any stimulus at ~50-150 ms is
        not site-specific and saturates at 100% for every site, so a
        criterion placed there measures nothing.
        """
        return self.spike_count(self.sites[int(site)], t0 + lo, t0 + hi)

    def responded(self, site, t0, lo=3.0, hi=25.0, threshold=1):
        return self.response(site, t0, lo, hi) >= threshold

    # ---- what the wet dish cannot read -----------------------------------

    def pathway_weight(self, pre_site, post_site):
        """Mean E->E conductance from one site's cells to another's.

        No electrode measures this. It is the quantity the whole exercise
        exists to expose.
        """
        m = ((self.site_of[self.pre] == int(pre_site))
             & (self.site_of[self.post] == int(post_site)))
        return float(self.conn.weights[m].mean()) if m.any() else 0.0

    def pathway_matrix(self):
        """`pathway_weight` for every ordered pair of sites."""
        return np.array([[self.pathway_weight(a, b)
                          for b in range(self.n_sites)]
                         for a in range(self.n_sites)])


# ---- protocol: Shahaf & Marom (2001) ------------------------------------

def shahaf_marom(dish, stim_site=0, target_site=3, window=(3.0, 25.0),
                 isi=200.0, rest=2000.0, max_stim=30, n_trials=8,
                 contingent=True, yoked=None, learn=True, overshoot=0):
    """Train a dish by withdrawing the stimulation when it complies.

    The original: a culture is stimulated at low rate through one
    electrode until a chosen *other* electrode responds within a target
    latency; stimulation then stops for a while, and the cycle repeats.
    The number of stimuli needed falls across cycles. The reward is the
    cessation, which is why the result is interesting -- there is no
    reward pathway in a dish to carry anything else.

    Parameters
    ----------
    contingent : True is the experiment. False is the control that gives
        it its meaning: the rests arrive after the same number of stimuli
        as in a paired contingent run (`yoked`), so the stimulation and
        the rest are statistically identical, but the rest no longer
        depends on what the network did. If the count still falls, the
        fall is not selection by stopping.
    yoked : the per-trial stimulus counts to replay when `contingent` is
        False. Required in that case.
    learn : False freezes the synapses. Any fall in the count then comes
        from the cells' own state, not from plasticity.
    overshoot : how many further stimuli to deliver after the criterion
        is met, before the rest. 0 is the experiment as published --
        stop the moment it complies. Larger values are the *overshoot*
        arm of N10: the pairing continues at the response's own lag,
        which is where the plasticity window may already be depressing.
        The search allowance is unchanged, so `max_stim` must leave room
        for both: a trial delivers at most `max_stim` stimuli whatever
        the overshoot asks for.

    Returns a dict with the per-trial counts, and -- checkpointed once per
    stimulus, because an endpoint cannot tell a search that never found
    the structure from one that found it and lost it -- the target
    response, the mean response of the other sites, and the mean
    conductance of the stimulated->target pathway.
    """
    if not contingent and yoked is None:
        raise ValueError("the yoked control needs the counts to replay")
    others = [s for s in range(dish.n_sites)
              if s not in (stim_site, target_site)]

    counts, hit_at, resp_t, resp_o, path, path_o = [], [], [], [], [], []
    for trial in range(n_trials):
        quota = max_stim if contingent else int(yoked[trial % len(yoked)])
        n, reached = 0, None
        while n < max_stim:
            t0 = dish.pulse(stim_site, follow=isi - 1.0, learn=learn)
            n += 1
            r = dish.response(target_site, t0, *window)
            o = np.mean([dish.response(s, t0, *window) for s in others])
            resp_t.append(r)
            resp_o.append(float(o))
            path.append(dish.pathway_weight(stim_site, target_site))
            path_o.append(float(np.mean([dish.pathway_weight(stim_site, s)
                                         for s in others])))
            if r >= 1 and reached is None:
                reached = n
            if contingent:
                if reached is not None and n >= reached + overshoot:
                    break
            elif n >= quota:
                break
        counts.append(n)
        hit_at.append(reached)
        dish.rest(rest, learn=learn)
    return {'counts': np.array(counts, dtype=float),
            'reached': hit_at,
            'target_response': np.array(resp_t, dtype=float),
            'other_response': np.array(resp_o, dtype=float),
            'pathway': np.array(path, dtype=float),
            'pathway_others': np.array(path_o, dtype=float)}


# ---- protocol: Kagan et al. (2022), the DishBrain contingency ------------

def stimulation_train(dish, kind, sites, n_pulses=4, interval=100.0,
                      jitter=75.0, rng=None, amplitude=40.0, learn=True):
    """The feedback the dish receives after it acts.

    DishBrain's contingency is not reward and punishment in the usual
    sense; it is *predictable* versus *unpredictable* stimulation. A hit
    is followed by a regular train at a fixed site, a miss by pulses at
    random sites and random times. Here the two are matched in pulse
    count, amplitude and mean rate, so predictability is the only
    difference between them -- which the original protocol did not do,
    and which is what makes an effect attributable.

    kind : 'predictable' or 'unpredictable'.
    sites : the pool the unpredictable train draws from; the predictable
        train uses its first element.

    Returns the (site, onset) pairs delivered, so a caller can anchor its
    windows to them.
    """
    rng = rng if rng is not None else np.random.RandomState(0)
    sites = list(sites)
    delivered = []
    for _ in range(int(n_pulses)):
        if kind == 'predictable':
            site, gap = sites[0], interval
        elif kind == 'unpredictable':
            site = sites[rng.randint(len(sites))]
            gap = interval + jitter * (2.0 * rng.rand() - 1.0)
        else:
            raise ValueError("kind must be 'predictable' or 'unpredictable'")
        delivered.append((site, dish.pulse(site, amplitude=amplitude,
                                           follow=gap - 1.0, learn=learn)))
    return delivered


def observed_lags(dish, t_lo, t_hi, max_lag=250.0):
    """Every pre/post spike lag across the plastic synapses in a window.

    Returned as `t_post - t_pre`, the sign convention of `STDP.window`.
    This is the quantity that decides what the rule does, and it is not
    observable in a dish: it needs the wiring and every spike time on
    both sides of every synapse.

    Feed it to `STDP.expected_change` to get the mean weight change per
    pair the rule will actually produce -- which equals the window
    integral only when the lags are spread evenly, and departs from it
    exactly when the stimulation imposes an order.

    `max_lag` is not a matter of taste. The two kernels P and M are
    normalised by the *number of lags collected*, so a cutoff that
    truncates the depression side -- the wider of the two, tau_minus =
    33.7 ms -- inflates P/M and makes the neutral point look far more
    negative than it is. Measured on the stimulation protocol of
    experiments/dish/neutral_point.py, P/M runs 0.734, 0.625, 0.573,
    0.557, 0.551, 0.551 at cutoffs of 40, 60, 100, 150, 250 and 400 ms:
    converged by 250, and wrong by a factor of 2.4 in the neutral
    integral at 60. The default is past the knee.
    """
    return lags_between(dish.pop.spike_times, dish.pre, dish.post,
                        t_lo, t_hi, max_lag)


def lags_between(spike_times, pre, post, t_lo, t_hi, max_lag=250.0):
    """`observed_lags` on raw spike trains and index arrays.

    Split out so a batched dish -- whose cells are a block of a larger
    population -- can collect the same quantity through the same code
    with offset indices, instead of a second copy of this loop. The
    `max_lag` reasoning in `observed_lags` applies here unchanged and is
    the reason this is one function rather than two.
    """
    spikes = [np.asarray([s for s in st if t_lo < s <= t_hi])
              for st in spike_times]
    out = []
    for p, q in zip(pre, post):
        a, b = spikes[p], spikes[q]
        if a.size == 0 or b.size == 0:
            continue
        d = b[None, :] - a[:, None]
        out.append(d[np.abs(d) <= max_lag])
    return np.concatenate(out) if out else np.zeros(0)


def two_alternative(dish, cues=(0, 1), motors=(4, 5), feedback_sites=(2, 3),
                    n_trials=60, window=(3.0, 60.0), rest=150.0,
                    n_pulses=4, interval=100.0, seed=0, learn=True,
                    record_lags=False):
    """DishBrain's contingency on a two-alternative task.

    Not Pong. Controlling a paddle needs a graded motor readout the
    preparation does not have -- both motor sites are swept by the same
    network burst, so their difference is noise until something makes it
    otherwise. The contingency, which is what the DishBrain result turns
    on, transfers without the paddle: one of two cue sites is stimulated,
    the dish 'chooses' whichever motor site responds more, and the
    feedback is predictable stimulation when the choice matches the cue
    and unpredictable stimulation when it does not.

    Returns the per-trial correctness, the margin between the two motor
    sites, and -- when `record_lags` is on -- the pre/post lags observed
    separately under each kind of feedback, which is where the mechanism
    is visible.
    """
    rng = np.random.RandomState(seed)
    correct, margin, chosen, cue_log = [], [], [], []
    lags = {'predictable': [], 'unpredictable': []}
    for trial in range(int(n_trials)):
        c = rng.randint(2)
        t0 = dish.pulse(cues[c], follow=window[1] + rest, learn=learn)
        r = [dish.response(m, t0, *window) for m in motors]
        pick = int(np.argmax(r)) if r[0] != r[1] else int(rng.randint(2))
        ok = (pick == c)
        correct.append(ok)
        chosen.append(pick)
        cue_log.append(c)
        margin.append(float(r[c] - r[1 - c]))

        kind = 'predictable' if ok else 'unpredictable'
        t_fb = dish.t
        stimulation_train(dish, kind, feedback_sites, n_pulses=n_pulses,
                          interval=interval, rng=rng, learn=learn)
        if record_lags:
            lags[kind].append(observed_lags(dish, t_fb, dish.t))
    out = {'correct': np.array(correct, dtype=float),
           'margin': np.array(margin, dtype=float),
           'cue': np.array(cue_log), 'choice': np.array(chosen)}
    if record_lags:
        out['lags'] = {k: (np.concatenate(v) if v else np.zeros(0))
                       for k, v in lags.items()}
    return out
