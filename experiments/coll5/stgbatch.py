"""COLL-5: K STG dishes stepped as one population.

`dish-batch` did this for the cortical `SiliconDish`; the STG line never got
it, and `## Run: COLL-4`'s registration recorded the absence as a fact. The
measurement that makes it worth building, taken before any of COLL-5 ran:

    bare STGNetwork step, AB/PD 2 cells
      K = 1  (60 cells)   121.2 us/step   121.2 per dish
      K = 4  (240)        165.1           41.3
      K = 12 (720)        265.8           22.1      -> 5.5x per dish

`STGDish.run` is dish-agnostic: it steps `self.syn`, `self.driver`,
`self.conn` and `self.spike_times`, all sized by `self.n`. So the batch
subclasses `STGDish` and overrides ONLY `__init__`, and the integration loop
that produces every number is the shipped one, unchanged. `gate()` holds it
to `dish-batch`'s own acceptance criterion: batched and separate must agree
bit for bit, not approximately.
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'src'))

from luviner.biophysics import stg
from luviner.biophysics.dish import rule_with_integral
from luviner.biophysics.plasticity import PlasticConnections, WeightNormalization
from luviner.biophysics.shortterm import ShortTermDepression
from luviner.biophysics.stg_dish import STGDish, STGDriver
from luviner.biophysics.synapse import Synapse


class BatchedSTGDish(STGDish):
    """K dishes as one network. Built from dishes, for their wiring only.

    The dishes must agree on everything that is a scalar of the whole
    population -- dt, the depression parameters, the budget, tau_norm and
    the rule -- and may differ in wiring, seed and `g_ie`, because those
    are per-synapse.
    """

    def __init__(self, dishes, g_e, g_i, g_ie, nmda_ratio=2.0, g_ei=0.3,
                 budget=None, tau_norm=200.0, rule=None, dt=0.05,
                 depress=(0.20, 400.0)):
        self.dishes = list(dishes)
        K = self.K = len(self.dishes)
        d0 = self.dishes[0]
        for d in self.dishes[1:]:
            for a in ('n', 'n_exc', 'n_inh', 'dt'):
                if getattr(d, a) != getattr(d0, a):
                    raise ValueError(f"dishes disagree on {a!r}")
        n1 = d0.n
        self.n_exc, self.n_inh = d0.n_exc * K, d0.n_inh * K
        self.n = n1 * K
        self.dt = float(dt)
        self.per = n1
        self.n_exc1, self.n_inh1 = d0.n_exc, d0.n_inh

        cells = []
        for _ in range(K):
            cells += [stg.STGNeuron(list(g_e), name=f'E{i}')
                      for i in range(self.n_exc1)]
            cells += [stg.STGNeuron(list(g_i), name=f'I{i}')
                      for i in range(self.n_inh1)]
        self.driver = STGDriver(stg.STGNetwork(cells, []))

        pre_e = np.concatenate([d.pre_e + k * n1
                                for k, d in enumerate(self.dishes)])
        post_e = np.concatenate([d.post_e + k * n1
                                 for k, d in enumerate(self.dishes)])
        w0 = np.concatenate([d.w0 for d in self.dishes])
        self.edge_slice = []
        off = 0
        for d in self.dishes:
            self.edge_slice.append(slice(off, off + d.pre_e.size))
            off += d.pre_e.size

        pre_ei, post_ei, pre_ie, post_ie = [], [], [], []
        g_ie_all = []
        for k, d in enumerate(self.dishes):
            exc = np.arange(self.n_exc1) + k * n1
            inh = self.n_exc1 + np.arange(self.n_inh1) + k * n1
            a, b = (x.ravel() for x in np.meshgrid(exc, inh, indexing='ij'))
            pre_ei.append(a); post_ei.append(b)
            c, e = (x.ravel() for x in np.meshgrid(inh, exc, indexing='ij'))
            pre_ie.append(c); post_ie.append(e)
            g_ie_all.append(np.full(c.size, float(np.atleast_1d(g_ie)[k if np.ndim(g_ie) else 0])))
        pre_ei = np.concatenate(pre_ei); post_ei = np.concatenate(post_ei)
        pre_ie = np.concatenate(pre_ie); post_ie = np.concatenate(post_ie)
        g_ie_all = np.concatenate(g_ie_all)

        pre = np.concatenate([pre_e, pre_e, pre_ei, pre_ie])
        post = np.concatenate([post_e, post_e, post_ei, post_ie])
        kinds = ([Synapse.ampa(g_max=1.0, delay=1.0)] * pre_e.size
                 + [Synapse.nmda(g_max=1.0, delay=1.0)] * pre_e.size
                 + [Synapse.ampa(g_max=1.0, delay=1.0)] * pre_ei.size
                 + [Synapse.gaba_a(g_max=1.0, delay=1.0, E_rev=-75.0)]
                 * pre_ie.size)
        weights = np.concatenate([w0, w0 * nmda_ratio,
                                  np.full(pre_ei.size, g_ei), g_ie_all])
        self._slice_e = slice(0, pre_e.size)
        self._slice_nmda = slice(pre_e.size, 2 * pre_e.size)

        stp = None
        if depress is not None:
            m = np.zeros(pre.size, dtype=bool)
            m[self._slice_e] = True
            m[self._slice_nmda] = True
            stp = ShortTermDepression(U=depress[0], tau_rec=depress[1], mask=m)
        from luviner.biophysics.stg_dish import SpikeSynapses
        self.syn = SpikeSynapses(pre, post, kinds, weights, self.dt, stp=stp)

        r = rule if rule is not None else rule_with_integral(-0.134)
        homeo = (WeightNormalization(budget=budget, tau=tau_norm,
                                     w_max=r.w_max)
                 if budget is not None else None)
        self.conn = PlasticConnections(pre_e, post_e, w0.copy(), rule=r,
                                       n_neurons=self.n, homeostasis=homeo)
        self.t = 0.0
        self.spike_times = [[] for _ in range(self.n)]

    # ---- per-dish readout -------------------------------------------------
    def cell(self, k, i):
        """Global index of dish k's cell i."""
        return k * self.per + i

    def spikes_of(self, k, i, t_lo=0.0, t_hi=None):
        return self.spikes(self.cell(k, i), t_lo, t_hi)

    def weights_of(self, k):
        return self.conn.weights[self.edge_slice[k]]
