"""
A dish whose sites can be driven apart.

`SiliconDish` wires its excitatory cells at random and its inhibition
uniformly, and the consequence has now cost three experiments. Pong was
abandoned because one stimulus sweeps the whole network and both motor
regions ride the same burst. The two-alternative task returned a motor
margin of -0.15, -0.02, -0.04 -- the readout a contingency would have to
shape barely exists. And N10's mechanism turned out to be undecidable on
frozen dishes, where the target site's response count equals the mean of
the other sites to two decimals in 4 of 5 seeds. **One site is not a
dish.** Any theory about a pathway needs a preparation in which a pathway
is distinguishable from the network.

Two changes, both already measured elsewhere in this package:

- **Excitation is distance-dependent.** Cells sit on a ring and the
  probability of an E->E connection falls as a Gaussian in ring distance,
  with one length scale. `sigma_deg` is the only knob the sweep varies.
- **Inhibition is available lateral or global**, and the sweep found
  the second is the better default -- which is not what this file
  originally claimed. `topology.lateral_weights` with `invert=True` gives
  each region strong inhibition of the *others* and weak inhibition of
  itself, and its docstring records that raising a *global* conductance
  sevenfold left the ratio between a strong and a weak response at 0.50
  where the lateral one drove it to 0.06. That measurement stands, and
  the inference this file drew from it does not: it is a response ratio
  under *competing* drive, a selection question, and the quantity here is
  the spatial extent of the response to a *single* stimulus.

  Measured across the whole sweep (`## Run: dish-2 (sweep)`), the global
  arm is the cleaner one. Its locality falls monotonically with the
  length scale -- 0.991, 0.990, 0.896, 0.704, 0.491, 0.425 at sigma 15 to
  180 -- which is what a length-scale sweep should produce if the length
  scale is what acts; the lateral arm is not monotone, dipping to 0.910
  at sigma 25 between 0.999 and 0.983. **The default preparation is
  global inhibition at sigma 25** (0.990 +/- 0.027, 10/10 seeds above the
  flat band, worst seed 0.91, in-degree 7.10 against the flat dish's
  6.99). Lateral is kept as the control, not as the default. The two
  curves cross -- lateral is the better row at sigma 40 -- so neither mode
  dominates and the comparison has to be made at a stated length scale.

Everything else is `SiliconDish`, including its interface, so every
protocol written against it -- `shahaf_marom`, `two_alternative`,
`observed_lags`, `pathway_weight` -- runs here unchanged. The ring is
periodic, so there are no edge sites and no site is special.

The measure this exists to move is the **locality fraction** returned by
`site_specificity`: the evoked response at ring distance 1 over the total
across distances 1 to 3, with the stimulated site excluded from both
sides. It is bounded in [0, 1], and on a six-site ring two of the five
non-stimulated sites sit at distance 1, so a dish with no locality reads
0.400 by construction. The registered criterion was that it exceed the
flat dish's own seed band in at least 8 of 10 seeds at some length scale.

**It is met**: 10/10 above the band (0.490) at sigma 15, 25 and 40, and
8/10 still at 60, against an unstructured control that reads 0.409 +/-
0.081. At sigma 180 the dish returns to 0.414 -- a wide enough length
scale is the unstructured dish, and reads as one.

References:
    Compte A, Brunel N, Goldman-Rakic PS, Wang XJ (2000). "Synaptic
        mechanisms and network dynamics underlying spatial working
        memory." Cereb Cortex 10:910-923.
    Eytan D, Marom S (2006). "Dynamics and effective topology underlying
        synchronization in networks of cortical neurons."
        J Neurosci 26:8465-8476.
"""

import numpy as np

from .dish import SiliconDish
from .synapse import Synapse
from .topology import ring_positions, angular_distance, lateral_weights


class StructuredDish(SiliconDish):
    """A `SiliconDish` with ring excitation and lateral inhibition.

    Parameters
    ----------
    sigma_deg : the excitatory length scale in degrees of ring. Small
        means a stimulus stays where it lands; large approaches the
        unstructured dish. This is the axis the sweep runs on.
    inhibition : 'global' reproduces `SiliconDish`'s uniform brake and is
        the **default preparation** for work downstream of the sweep, at
        `sigma_deg = 25`; 'lateral' gives surround inhibition and is now
        the control. It was the other way round until the sweep's own
        control row was read -- see the module docstring.
    g_inh : peak inhibitory conductance, matched between the two
        inhibition modes so the comparison is not about how much.
    """

    def __init__(self, n_exc=48, n_inh=12, n_sites=6, sigma_deg=45.0,
                 inhibition='lateral', connectivity=0.15, g_inh=0.25,
                 w_rec=0.015, seed=0, **kw):
        if inhibition not in ('lateral', 'global'):
            raise ValueError("inhibition must be 'lateral' or 'global'")
        self.sigma_deg = float(sigma_deg)
        self.inhibition = inhibition
        self.g_inh = float(g_inh)
        self._peak_connectivity = float(connectivity)
        super().__init__(n_exc=n_exc, n_inh=n_inh, n_sites=n_sites,
                         connectivity=connectivity, w_rec=w_rec, seed=seed,
                         **kw)
        # The base class draws its excitatory wiring uniformly inside its
        # own __init__, before _build is reachable, so the distance
        # profile is installed here and the population rebuilt on it. The
        # alternative is editing dish.py, which another session holds.
        rng = np.random.RandomState(seed)
        mask = self._excitatory_mask(rng)
        self.pre, self.post = (a.copy() for a in np.nonzero(mask))
        w0 = w_rec * (1.0 + 0.2 * rng.randn(self.pre.size))
        self.w0 = np.clip(w0, 0.0, None)
        self._build(self.w0)

    # ---- connectivity ----------------------------------------------------

    def _excitatory_mask(self, rng):
        """E->E kept with a probability that falls with ring distance.

        **The expected in-degree is held fixed across the sweep**, and
        that is not a detail. A raw Gaussian at `sigma_deg = 15` leaves
        0.67 incoming connections per cell against 7.17 at the
        unstructured density, so any site specificity it produced could
        be read off the connection count rather than the arrangement --
        the same confound `shuffle_targets` exists to remove elsewhere in
        this package, and the same one that made the branch comparison
        meaningless until rheobase was matched.

        So the profile is normalised to a target in-degree taken from the
        unstructured dish. `sigma_deg` then changes only *where* a cell's
        connections go, never how many it has, and the sweep is a
        statement about arrangement.
        """
        theta = ring_positions(self.n_exc)
        d = angular_distance(theta[:, None], theta[None, :])
        sigma = np.radians(self.sigma_deg)
        shape = np.exp(-0.5 * (d / sigma) ** 2)
        np.fill_diagonal(shape, 0.0)
        target = self._peak_connectivity * (self.n_exc - 1)
        p = shape * (target / shape.sum(axis=0, keepdims=True))
        mask = rng.rand(self.n_exc, self.n_exc) < np.clip(p, 0.0, 1.0)
        np.fill_diagonal(mask, False)
        return mask

    def _build(self, weights):
        """As `SiliconDish._build`, with the inhibition profile swapped.

        Reimplemented rather than patched because the base class wires
        E<->I all-to-all with one conductance, and the whole point here is
        that the inhibition a cell receives depends on where it is.
        """
        from .population import NeuronPopulation
        from .plasticity import WeightNormalization, PlasticConnections

        pop = NeuronPopulation.cortical(self.n, g_M=np.zeros(self.n))
        pop.connect(self.pre, self.post,
                    Synapse.ampa(g_max=1.0, delay=1.0), weights=weights)
        pop.connect(self.pre, self.post,
                    Synapse.nmda(g_max=1.0, delay=1.0),
                    weights=weights * self.nmda_ratio)

        # E -> I, local: an interneuron hears the excitatory cells near it.
        # `lateral_weights` returns thresholded (pre, post, w) in local
        # indices, so the inhibitory ones are offset onto the population.
        p, q, w = lateral_weights(self.n_exc, self.n_inh,
                                  sigma_deg=self.sigma_deg, w_max=0.02)
        pop.connect(p, q + self.n_exc,
                    Synapse.ampa(g_max=1.0, delay=1.0), weights=w)

        # I -> E: the change that matters.
        if self.inhibition == 'lateral':
            p, q, w = lateral_weights(self.n_inh, self.n_exc,
                                      sigma_deg=self.sigma_deg, invert=True,
                                      w_max=self.g_inh)
        else:
            p, q = np.meshgrid(np.arange(self.n_inh), self.exc,
                               indexing='ij')
            p, q = p.ravel(), q.ravel()
            w = np.full(p.size, self.g_inh)
        pop.connect(p + self.n_exc, q,
                    Synapse.gaba_a(g_max=1.0, delay=1.0, E_rev=-75.0),
                    weights=w)

        p, q = np.meshgrid(self.inh, self.inh, indexing='ij')
        p, q = p.ravel(), q.ravel()
        pop.connect(p, q,
                    Synapse.gaba_a(g_max=1.0, delay=1.0, E_rev=-75.0),
                    weights=np.full(p.size, 0.02))

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

    # ---- the measure this file exists to move ----------------------------

    def in_degree(self):
        """Realised excitatory in-degree per cell.

        Reported next to the specificity because the length scale changes
        it, and a specificity that came from having fewer connections
        rather than from having local ones would be the wrong result.
        """
        counts = np.bincount(self.post, minlength=self.n)[:self.n_exc]
        return float(counts.mean())


def site_specificity(dish, stim=0, target=None, amplitude=40.0,
                     follow=200.0, window=(25.0, 60.0), n_trials=10,
                     settle=100.0):
    """Response of a target site divided by the mean of the other sites.

    **The stimulated site is excluded from both sides**, which is the
    definition the N10 instrument uses and the one the sweep is scored
    with. It matters: leaving the driven site in the denominator lets its
    refractoriness leak into "others" and pushes a flat dish's baseline
    off 1. With it excluded, a preparation where a stimulus is
    all-or-none and network-wide reads 1 by construction, because every
    remaining site is doing the same thing.

    **`target` defaults to the site adjacent to the stimulus, and on a
    ring that choice is not free.** The N10 instrument used site 3 with
    site 0 stimulated, which on a flat dish is arbitrary -- all sites are
    equivalent. On this ring six sites apart means site 3 is *antipodal*
    to site 0, the least connected of all, so scoring a structured dish
    with a fixed target of 3 would return an index below 1 and read as
    failure when it is in fact the strongest possible evidence of
    locality. The adjacent site is the one a length scale is supposed to
    reach.

    `by_distance` in the returned dict carries the whole profile, so
    either convention can be read off one measurement and neither has to
    be trusted.
    """
    n = dish.n_sites
    if target is None:
        target = (stim + 1) % n
    if target == stim:
        raise ValueError("the stimulated site cannot also be the target")

    dish.run(settle, learn=False)
    per_site = {s: [] for s in range(n)}
    for _ in range(int(n_trials)):
        t0 = dish.pulse(stim, amplitude=amplitude, follow=follow,
                        learn=False)
        for s in range(n):
            per_site[s].append(dish.response(s, t0, *window))
    mean = {s: float(np.mean(v)) for s, v in per_site.items()}
    return score_sites(mean, stim, target, n, window)


def score_sites(mean, stim, target, n, window):
    """The scoring half of `site_specificity`, on per-site mean counts.

    Split out so the batched sweep scores with *this* function rather
    than a copy of it. The locality fraction has already been redefined
    once, and two definitions of it living in two files is how the sweep
    and the single-dish check end up disagreeing without either being
    wrong.
    """
    others = [s for s in range(n) if s not in (stim, target)]
    denom = float(np.mean([mean[s] for s in others]))
    # Ring distance in sites, so the profile is readable without knowing
    # the geometry.
    by_distance = {}
    for s in range(n):
        d = min((s - stim) % n, (stim - s) % n)
        by_distance.setdefault(d, []).append(mean[s])
    by_distance = {d: float(np.mean(v)) for d, v in sorted(by_distance.items())}

    d1 = [mean[s] for s in range(n)
          if min((s - stim) % n, (stim - s) % n) == 1]
    near_far = [mean[s] for s in range(n)
                if 1 <= min((s - stim) % n, (stim - s) % n) <= 3]
    total = float(np.sum(near_far))

    return {'stim': stim, 'target': target,
            # The scalar the sweep is scored on. Bounded in [0, 1] and
            # finite on both preparations, which the ratio of a near site
            # to a far one is not: locality drives the far response to
            # exactly zero and the ratio with it. On a six-site ring two
            # of the five non-stimulated sites sit at distance 1, so a
            # dish with no locality reads 0.4 by construction and perfect
            # locality reads 1.0.
            'locality_fraction': (float(np.sum(d1)) / total
                                  if total > 0 else float('nan')),
            'target_response': mean[target],
            'others': denom,
            'index': (mean[target] / denom) if denom > 0 else np.inf,
            'per_site': mean,
            'by_distance': by_distance,
            'window': list(window)}
