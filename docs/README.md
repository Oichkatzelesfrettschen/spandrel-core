# spandrel-core Documentation

Core utilities for Type Ia supernova cosmology and Spandrel hypothesis testing.

## Modules

- `spandrel_core.cosmology`: Cosmological distance calculators.
- `spandrel_core.data`: Data loaders for Pantheon+ and simulations.
- `spandrel_core.likelihood`: Bayesian likelihood for hypothesis testing.
- `spandrel_core.constants`: Physical and cosmological constants.

## Quick Start

```python
from spandrel_core.cosmology import luminosity_distance

dL = luminosity_distance(z=0.5)
print(f"Luminosity distance at z=0.5: {dL:.2f} Mpc")
```
