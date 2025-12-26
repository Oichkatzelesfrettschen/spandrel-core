"""Basic tests for spandrel-core package."""

import numpy as np
import pytest


def test_constants_import():
    """Verify constants module loads correctly."""
    from spandrel_core.constants import C_LIGHT_CGS, H0_FIDUCIAL, M_SUN

    assert C_LIGHT_CGS == pytest.approx(2.99792458e10)
    assert H0_FIDUCIAL == 70.0
    assert M_SUN == pytest.approx(1.989e33)


def test_sn_dataset():
    """Verify SNDataset dataclass works."""
    from spandrel_core import SNDataset

    data = SNDataset(
        name="test",
        z=np.array([0.1, 0.2]),
        z_cmb=np.array([0.1, 0.2]),
        m_b=np.array([20.0, 21.0]),
        m_b_err=np.array([0.1, 0.1]),
        x1=np.array([0.0, 0.5]),
        x1_err=np.array([0.1, 0.1]),
        c=np.array([0.0, 0.05]),
        c_err=np.array([0.02, 0.02]),
        host_mass=np.array([10.0, 11.0]),
        host_mass_err=np.array([0.1, 0.1]),
    )

    assert len(data) == 2
    assert data.name == "test"


def test_simulated_loader():
    """Verify SimulatedDataLoader generates valid data."""
    from spandrel_core import SimulatedDataLoader

    loader = SimulatedDataLoader(seed=42)
    data = loader.generate(n_sn=100)

    assert len(data) == 100
    assert data.z.min() >= 0.01
    assert data.z.max() <= 1.5
    assert np.all(np.isfinite(data.m_b))


def test_cosmology_params():
    """Verify CosmologyParams dataclass."""
    from spandrel_core import CosmologyParams

    cosmo = CosmologyParams(H0=70.0, Om=0.3, w0=-1.0, wa=0.0)

    assert cosmo.H0 == 70.0
    assert cosmo.Om == 0.3
    assert cosmo.w0 == -1.0


def test_host_mass_split():
    """Verify host mass split function."""
    from spandrel_core import SimulatedDataLoader, split_by_host_mass

    loader = SimulatedDataLoader(seed=42)
    data = loader.generate(n_sn=100)

    low, high = split_by_host_mass(data, threshold=10.0)

    assert len(low) + len(high) == len(data)
    assert np.all(low.host_mass < 10.0)
    assert np.all(high.host_mass >= 10.0)
