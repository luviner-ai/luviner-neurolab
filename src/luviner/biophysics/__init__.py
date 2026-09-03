"""
Biophysics — real neuron simulation, not neural network engines.

Everything under luviner.biophysics integrates measured physical
quantities (ionic conductances, reversal potentials, gating kinetics).
Nothing here is trained, and nothing here is a machine-learning model:
this package simulates neurons rather than abstracting them.
"""

from .hodgkin_huxley import (
    HodgkinHuxleyNeuron,
    detect_spikes,
    firing_rate,
)
from .cortical import CorticalNeuron
from .channels import DendriticCaAP
from .dendrite import (Compartment, MultiCompartmentNeuron, build_star)
from .synapse import Synapse
from .network import Network, psp_amplitude, transmission_delay
from .population import NeuronPopulation
from .analysis import (population_rate, rate_spectrum, spectral_peak,
                       is_oscillatory, synchrony_index)
from .adaptive import AdaptiveIntegrator, hh_derivatives
from .adaptive_network import AdaptiveNetwork
from .blockstep import BlockTimestepNetwork
from .plasticity import (STDP, SynapticScaling, WeightNormalization,
                         PlasticConnections, pair_protocol)
from .shortterm import ShortTermDepression, depressing
from .dish import (SiliconDish, rule_with_integral, shahaf_marom,
                   lag_kernels, rule_with_expected_change,
                   stimulation_train, observed_lags, lags_between,
                   two_alternative)
from .branched import (BranchedPopulation, BranchGatedConnections,
                       BranchNormalization, GateHomeostasis)
from .dish_structured import StructuredDish, site_specificity
from .dish_batch import (DishBatch, shahaf_marom_batch,
                         site_specificity_batch, observed_lags_of)
from .modulation import ModulatedCircuit, modulated_fully_connected
from .topology import (ring_weights, uniform_weights, lateral_weights, shuffle_targets,
                       population_vector, bump_width, rates_in_window,
                       ring_positions, angular_distance)

__all__ = [
    'HodgkinHuxleyNeuron', 'CorticalNeuron', 'Compartment',
    'MultiCompartmentNeuron', 'build_star', 'DendriticCaAP', 'detect_spikes', 'firing_rate',
    'Synapse', 'Network', 'psp_amplitude', 'transmission_delay',
    'NeuronPopulation', 'population_rate', 'rate_spectrum', 'spectral_peak',
    'is_oscillatory', 'synchrony_index',
    'AdaptiveIntegrator', 'hh_derivatives', 'AdaptiveNetwork', 'BlockTimestepNetwork',
    'STDP', 'SynapticScaling', 'WeightNormalization', 'PlasticConnections',
    'pair_protocol', 'ShortTermDepression', 'depressing',
    'SiliconDish', 'rule_with_integral', 'shahaf_marom',
    'stimulation_train', 'observed_lags', 'lags_between',
    'two_alternative',
    'lag_kernels', 'rule_with_expected_change',
    'BranchedPopulation', 'BranchGatedConnections', 'BranchNormalization',
    'GateHomeostasis',
    'ModulatedCircuit', 'modulated_fully_connected',
    'StructuredDish', 'site_specificity',
    'DishBatch', 'shahaf_marom_batch', 'site_specificity_batch',
    'observed_lags_of',
    'ring_weights', 'uniform_weights', 'lateral_weights', 'shuffle_targets', 'population_vector',
    'bump_width', 'rates_in_window', 'ring_positions', 'angular_distance',
]
