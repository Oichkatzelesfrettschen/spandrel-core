# spandrel-core

Core utilities for Type Ia supernova cosmology analysis and Spandrel hypothesis testing.

## Overview

The Spandrel hypothesis proposes that Type Ia supernovae exhibit population evolution with redshift, potentially biasing dark energy constraints. This package provides the foundational tools for testing this hypothesis.

## Installation

```bash
pip install spandrel-core

# With FITS support (for astropy tables)
pip install spandrel-core[fits]

# With MCMC sampling support
pip install spandrel-core[sampling]
```

## Quick Start

```python
from spandrel_core import SNDataset, SimulatedDataLoader, SpandrelLikelihood
from spandrel_core.constants import H0_PLANCK, OMEGA_M_FIDUCIAL

# Generate test data
loader = SimulatedDataLoader(seed=42)
data = loader.generate(n_sn=1000, include_evolution=True, dM_dz=0.1)

# Create likelihood
likelihood = SpandrelLikelihood(
    z=data.z,
    m_obs=data.m_b,
    m_err=data.m_b_err,
    x1=data.x1,
    x1_err=data.x1_err,
    c=data.c,
    c_err=data.c_err,
    host_mass=data.host_mass,
    host_mass_err=data.host_mass_err,
)

# Evaluate log-likelihood
theta = [-19.3, 0.14, 3.1, 0.05, 0.1, 70.0, 0.3, -1.0, 0.0]
loglik = likelihood.log_likelihood_baseline(theta)
```

## Modules

- **constants**: Physical constants in CGS and cosmological units
- **data**: Data loaders for Pantheon+ and simulated SN datasets
- **likelihood**: Bayesian likelihood for Spandrel hypothesis testing

## License

GPL-3.0-or-later
