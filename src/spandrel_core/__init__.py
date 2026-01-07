"""
Spandrel Core: Shared utilities for Type Ia supernova cosmology analysis.

This package provides the foundational components for the Spandrel hypothesis
framework, which tests for population evolution in Type Ia supernovae.

Modules:
    constants: Physical constants in CGS and cosmological units
    data: Data loaders for Pantheon+ and simulated SN datasets
    likelihood: Bayesian likelihood for Spandrel hypothesis testing

Usage:
    from spandrel_core import SNDataset, PantheonPlusLoader
    from spandrel_core.constants import C_LIGHT_KMS, H0_FIDUCIAL
    from spandrel_core.likelihood import SpandrelLikelihood
"""

__version__ = "0.1.0"
__author__ = "Deirikr Jaiusadastra Afrauthihinngreygaard"

# Cosmology utilities
from spandrel_core import constants
from spandrel_core.cosmology import (
    FlatLambdaCDM,
    angular_diameter_distance,
    comoving_distance,
    luminosity_distance,
)
from spandrel_core.data import (
    PantheonPlusLoader,
    SimulatedDataLoader,
    SNDataset,
    split_by_host_mass,
    split_by_redshift,
)

# Likelihood classes
from spandrel_core.likelihood import (
    CosmologyParams,
    EvolutionParams,
    SpandrelLikelihood,
    StandardizationParams,
    compute_model_comparison,
)

__all__ = [
    # Version
    "__version__",
    # Data structures
    "SNDataset",
    "PantheonPlusLoader",
    "SimulatedDataLoader",
    "split_by_host_mass",
    "split_by_redshift",
    # Likelihood
    "CosmologyParams",
    "StandardizationParams",
    "EvolutionParams",
    "SpandrelLikelihood",
    "compute_model_comparison",
    # Cosmology
    "FlatLambdaCDM",
    "luminosity_distance",
    "comoving_distance",
    "angular_diameter_distance",
    # Submodules
    "constants",
]
