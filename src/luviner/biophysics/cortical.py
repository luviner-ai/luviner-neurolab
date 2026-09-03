"""
Cortical neuron: a Class I regular-spiking cell.

The 1952 model describes a squid axon, and its excitability class is a
structural property, not a parameter: it is Class II, so it cannot fire
below ~55 Hz. That floor is what limits spike frequency adaptation in
the base model -- adding I_M or I_KCa makes the cell slow down only until
it hits the floor, then it stops firing altogether rather than slowing
further. Adding I_A does not move the floor either; Connor and Stevens
obtained Class I with their *whole* model, not by adding one channel.

This module changes the base instead. `CorticalNeuron` uses the
Traub-Miles kinetics with an adjustable spike threshold V_T and the
cortical conductance set from Pospischil et al. (2008): a much smaller
delayed rectifier (g_Kd = 5 vs 36) and a smaller leak, which moves the
bifurcation and gives Class I excitability -- firing rate rises
continuously from near zero.

Everything else carries over: it is a subclass, so the extra channels,
voltage clamp, channel blockers, synapses and networks all work on it
unchanged. Put I_M on this cell and adaptation behaves the way it does
in a real pyramidal neuron.

Rate functions are written relative to V_T, so raising or lowering the
threshold shifts the whole spike-generating machinery together, which is
how Pospischil et al. fit different cell classes with one formalism.

References:
    Traub RD, Miles R (1991). Neuronal Networks of the Hippocampus.
        Cambridge University Press.
    Pospischil M, Toledo-Rodriguez M, Monier C, Piwkowska Z, Bal T,
        Fregnac Y, Markram H, Destexhe A (2008). "Minimal Hodgkin-Huxley
        type models for different classes of cortical and thalamic
        neurons." Biol Cybern 99:427-441.
"""

import numpy as np

from .hodgkin_huxley import HodgkinHuxleyNeuron, _x_over_1_minus_exp_neg_x


class CorticalNeuron(HodgkinHuxleyNeuron):
    """Regular-spiking cortical cell (Traub-Miles / Pospischil).

    Differs from HodgkinHuxleyNeuron in its conductances, reversal
    potentials and rate functions. The API is identical.

    Parameters
    ----------
    V_T : spike threshold, mV. The rate functions are all expressed
          relative to it, so it shifts excitability as a unit.
    """

    def __init__(self,
                 C_m=1.0,
                 g_Na=50.0, g_K=5.0, g_L=0.1,
                 E_Na=50.0, E_K=-90.0, E_L=-70.0,
                 V_T=-63.0, temperature=6.3, channels=(), calcium=None):
        # V_T must exist before reset() evaluates the rate functions.
        self.V_T = V_T
        super().__init__(C_m=C_m, g_Na=g_Na, g_K=g_K, g_L=g_L,
                         E_Na=E_Na, E_K=E_K, E_L=E_L,
                         temperature=temperature,
                         channels=channels, calcium=calcium)

    # ---- Traub-Miles kinetics, referred to V_T ---------------------------
    #
    # Overriding the base class's staticmethods with instance methods:
    # these depend on V_T, which is per-neuron. Calls go through
    # self.alpha_m(V) either way, so the base class machinery is unchanged.

    def alpha_m(self, V):
        return 1.28 * _x_over_1_minus_exp_neg_x((V - self.V_T - 13.0) / 4.0)

    def beta_m(self, V):
        return 1.4 * _x_over_1_minus_exp_neg_x(-(V - self.V_T - 40.0) / 5.0)

    def alpha_h(self, V):
        return 0.128 * np.exp(-(V - self.V_T - 17.0) / 18.0)

    def beta_h(self, V):
        return 4.0 / (1.0 + np.exp(-(V - self.V_T - 40.0) / 5.0))

    def alpha_n(self, V):
        return 0.16 * _x_over_1_minus_exp_neg_x((V - self.V_T - 15.0) / 5.0)

    def beta_n(self, V):
        return 0.5 * np.exp(-(V - self.V_T - 10.0) / 40.0)

    def reset(self, V0=None):
        """Rest the cell. Defaults to E_L rather than the squid's -65 mV."""
        return super().reset(self.E_L if V0 is None else V0)

    def block_sodium(self, fraction=1.0):
        self.g_Na = 50.0 * (1.0 - fraction)
        return self

    def block_potassium(self, fraction=1.0):
        self.g_K = 5.0 * (1.0 - fraction)
        return self
