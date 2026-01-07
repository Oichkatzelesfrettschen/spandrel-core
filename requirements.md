# Installation Requirements (spandrel-core)

`spandrel-core` is a Python library submodule (`src/spandrel_core`) containing
shared primitives used by cosmology applications (data loading, likelihood,
distance-modulus utilities).

## Prerequisites

- Python `>=3.9`

## Install (dev)

```bash
cd spandrel-core
python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -e '.[dev]'
```

## Tests and gates

- Tests: `pytest`
- Strict gates (warnings-as-errors on the contract surface): `scripts/audit/run_tiers.sh`

