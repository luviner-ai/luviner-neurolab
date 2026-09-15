"""COLL-3: the same dish on two cell models.

`SiliconDish._build` hardcodes `NeuronPopulation.cortical`. The second cell
model the row asks for -- plain Hodgkin-Huxley -- turns out to be the SAME
class with `kinetics='squid'` and the textbook conductances (g_Na 120,
g_K 36, g_L 0.3, E_L -54.387), so the two arms differ in one constructor
call and in nothing else.

`_build`'s body is duplicated here rather than edited in `src/`, because the
shipped dish must stay the reference for every result already on the board.
The duplication is the risk, so `stage0` gates it: the cortical path through
this subclass must reproduce `SiliconDish` bit for bit.
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'src'))

from luviner.biophysics import SiliconDish
from luviner.biophysics.plasticity import PlasticConnections, WeightNormalization
from luviner.biophysics.population import NeuronPopulation
from luviner.biophysics.synapse import Synapse

N_EXC, N_INH, CONN, W_REC, BUDGET = 48, 12, 0.15, 0.015, 0.107


class KineticsDish(SiliconDish):
    """`SiliconDish` with the cell model as a parameter."""

    KINETICS = 'cortical'

    def _population(self):
        if self.KINETICS == 'cortical':
            return NeuronPopulation.cortical(
                self.n, g_M=np.full(self.n, getattr(self, 'g_M', 0.0)))
        return NeuronPopulation(self.n, kinetics='squid')

    def _build(self, weights):
        # ---- copied from SiliconDish._build, with one line replaced ----
        pop = self._population()
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
        # ---- end of the copied body ----


class CorticalDish(KineticsDish):
    KINETICS = 'cortical'


class HHDish(KineticsDish):
    KINETICS = 'squid'


def build(cls, seed, drive=0.90, plastic=True, budget=BUDGET):
    return cls(n_exc=N_EXC, n_inh=N_INH, connectivity=CONN, w_rec=W_REC,
               drive=drive, budget=budget, plastic=plastic, seed=seed)
