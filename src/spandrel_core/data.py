#!/usr/bin/env python3
"""
Data Loader for Spandrel Analysis

Loads and preprocesses Type Ia supernova data from various compilations:
- Pantheon+ (Scolnic et al. 2022)
- DES-SN5YR (DES Collaboration)
- Union3 (Rubin et al. 2023)

Also loads host galaxy properties from associated catalogs.
"""

import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple, Union

import numpy as np
import pandas as pd

# Try to import astropy for FITS support
try:
    from astropy.io import fits  # noqa: F401
    from astropy.table import Table  # noqa: F401

    HAS_ASTROPY = True
except ImportError:
    HAS_ASTROPY = False
    # Silent fallback - astropy optional for CSV-based workflows


@dataclass
class SNDataset:
    """Container for a Type Ia supernova dataset."""

    name: str
    z: np.ndarray  # Heliocentric redshift
    z_cmb: np.ndarray  # CMB-frame redshift
    m_b: np.ndarray  # Peak B-band magnitude (after MW correction)
    m_b_err: np.ndarray  # Magnitude uncertainty
    x1: np.ndarray  # SALT2 stretch
    x1_err: np.ndarray  # Stretch uncertainty
    c: np.ndarray  # SALT2 color
    c_err: np.ndarray  # Color uncertainty
    host_mass: np.ndarray  # log10(M_host/M_sun)
    host_mass_err: np.ndarray  # Host mass uncertainty
    cov_sys: Optional[np.ndarray] = None  # Systematic covariance matrix
    sn_names: Optional[np.ndarray] = None  # SN identifiers

    # Optional additional host properties
    host_sfr: Optional[np.ndarray] = None  # Star formation rate
    host_ssfr: Optional[np.ndarray] = None  # Specific SFR
    host_metallicity: Optional[np.ndarray] = None  # Metallicity proxy

    def __len__(self):
        return len(self.z)

    def __repr__(self):
        return f"SNDataset(name='{self.name}', n_sn={len(self)})"

    def select(self, mask: np.ndarray) -> "SNDataset":
        """Return a subset of the dataset based on boolean mask."""
        return SNDataset(
            name=f"{self.name}_subset",
            z=self.z[mask],
            z_cmb=self.z_cmb[mask],
            m_b=self.m_b[mask],
            m_b_err=self.m_b_err[mask],
            x1=self.x1[mask],
            x1_err=self.x1_err[mask],
            c=self.c[mask],
            c_err=self.c_err[mask],
            host_mass=self.host_mass[mask],
            host_mass_err=self.host_mass_err[mask],
            sn_names=self.sn_names[mask] if self.sn_names is not None else None,
            host_sfr=self.host_sfr[mask] if self.host_sfr is not None else None,
            host_ssfr=self.host_ssfr[mask] if self.host_ssfr is not None else None,
            host_metallicity=(
                self.host_metallicity[mask] if self.host_metallicity is not None else None
            ),
        )


class PantheonPlusLoader:
    """
    Loader for Pantheon+ sample (Scolnic et al. 2022).

    Expected data files:
    - Pantheon+SH0ES.dat: Main SN data
    - Pantheon+SH0ES_STAT+SYS.cov: Covariance matrix
    """

    DATA_URL = "https://github.com/PantheonPlusSH0ES/DataRelease"

    def __init__(self, data_dir: Union[str, Path]):
        """
        Initialize loader with data directory.

        Args:
            data_dir: Path to directory containing Pantheon+ files
        """
        self.data_dir = Path(data_dir)

    def load(
        self, z_min: float = 0.01, z_max: float = 2.3, require_host_mass: bool = True
    ) -> SNDataset:
        """
        Load Pantheon+ dataset.

        Args:
            z_min: Minimum redshift cut
            z_max: Maximum redshift cut
            require_host_mass: Only include SNe with host mass measurements

        Returns:
            SNDataset with Pantheon+ data
        """
        # Main data file
        data_file = self.data_dir / "Pantheon+SH0ES.dat"

        if not data_file.exists():
            raise FileNotFoundError(
                f"Pantheon+ data file not found: {data_file}\n" f"Download from: {self.DATA_URL}"
            )

        # Read data (space-separated, with header)
        df = pd.read_csv(data_file, sep=r"\s+", comment="#")

        # Expected columns (may vary by version)
        required_cols = ["zHD", "zHEL", "mB", "mBERR", "x1", "x1ERR", "c", "cERR"]

        # Check for required columns
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            raise ValueError(f"Missing columns in Pantheon+ file: {missing}")

        # Apply redshift cuts
        mask = (df["zHD"] >= z_min) & (df["zHD"] <= z_max)

        # Host mass (column name may vary)
        host_mass_col = None
        for col in ["HOST_LOGMASS", "HOSTGAL_LOGMASS", "host_logmass"]:
            if col in df.columns:
                host_mass_col = col
                break

        if host_mass_col is None:
            warnings.warn("No host mass column found; using placeholder values")
            host_mass = np.full(len(df), 10.5)
            host_mass_err = np.full(len(df), 0.5)
        else:
            host_mass = df[host_mass_col].values
            # Error column
            err_col = host_mass_col + "ERR"
            if err_col in df.columns:
                host_mass_err = df[err_col].values
            else:
                host_mass_err = np.full(len(df), 0.1)

            if require_host_mass:
                mask &= np.isfinite(host_mass) & (host_mass > 0)

        df_cut = df[mask]

        # Create dataset
        dataset = SNDataset(
            name="Pantheon+",
            z=df_cut["zHEL"].values,
            z_cmb=df_cut["zHD"].values,
            m_b=df_cut["mB"].values,
            m_b_err=df_cut["mBERR"].values,
            x1=df_cut["x1"].values,
            x1_err=df_cut["x1ERR"].values,
            c=df_cut["c"].values,
            c_err=df_cut["cERR"].values,
            host_mass=host_mass[mask],
            host_mass_err=host_mass_err[mask],
            sn_names=df_cut["CID"].values if "CID" in df_cut.columns else None,
        )

        # Try to load covariance matrix
        cov_file = self.data_dir / "Pantheon+SH0ES_STAT+SYS.cov"
        if cov_file.exists():
            try:
                cov = self._load_covariance(cov_file, len(dataset))
                dataset.cov_sys = cov
            except Exception as e:
                warnings.warn(f"Could not load covariance: {e}")

        return dataset

    def _load_covariance(self, cov_file: Path, n_sn: int) -> np.ndarray:
        """Load covariance matrix from file."""
        # Pantheon+ cov format: first line is N, then N*N values
        with open(cov_file) as f:
            n = int(f.readline().strip())
            values = np.array([float(x) for x in f.read().split()])

        cov = values.reshape(n, n)

        # Subset if needed (after redshift cuts)
        if n != n_sn:
            warnings.warn(f"Covariance matrix size ({n}) doesn't match dataset ({n_sn})")
            return None

        return cov


class SimulatedDataLoader:
    """
    Generate simulated SN data for testing.

    Useful for:
    - Testing analysis pipeline
    - Injection-recovery tests
    - Validation of Spandrel detection sensitivity
    """

    def __init__(self, seed: int = 42):
        """Initialize with random seed."""
        self.rng = np.random.default_rng(seed)

    def generate(
        self,
        n_sn: int = 1000,
        z_range: Tuple[float, float] = (0.01, 1.5),
        include_evolution: bool = False,
        dM_dz: float = 0.0,
        dx1_dz: float = 0.0,
        dc_dz: float = 0.0,
    ) -> SNDataset:
        """
        Generate simulated SN dataset.

        Args:
            n_sn: Number of SNe
            z_range: Redshift range (z_min, z_max)
            include_evolution: Whether to include population evolution
            dM_dz: Luminosity evolution (mag per unit z)
            dx1_dz: Stretch evolution (per unit z)
            dc_dz: Color evolution (per unit z)

        Returns:
            SNDataset with simulated data
        """
        # Redshift distribution (uniform in comoving volume, simplified)
        z = self.rng.uniform(z_range[0], z_range[1], n_sn)

        # Population parameters
        x1_mean = 0.0 + dx1_dz * z if include_evolution else np.zeros(n_sn)
        c_mean = 0.0 + dc_dz * z if include_evolution else np.zeros(n_sn)

        x1 = self.rng.normal(x1_mean, 1.0)
        c = self.rng.normal(c_mean, 0.1)

        # Host masses (bimodal distribution)
        host_mass = self.rng.choice(
            [self.rng.normal(10.0, 0.3), self.rng.normal(11.0, 0.3)], size=n_sn, p=[0.4, 0.6]
        )
        host_mass = np.array(
            [
                (
                    self.rng.normal(10.0, 0.3)
                    if self.rng.random() < 0.4
                    else self.rng.normal(11.0, 0.3)
                )
                for _ in range(n_sn)
            ]
        )

        # True cosmology
        Om = 0.3
        w0 = -1.0
        wa = 0.0
        H0 = 70.0

        # Distance modulus (flat ΛCDM approximation)
        # Using simple approximation for speed
        c_light = 299792.458
        dc = c_light / H0 * z * (1 - 0.25 * Om * z)  # Very rough approximation
        dl = dc * (1 + z)
        mu_true = 5 * np.log10(dl) + 25

        # Standardization
        M0 = -19.3
        alpha = 0.14
        beta = 3.1
        gamma = 0.05
        sigma_int = 0.1

        # Evolution term
        M_evol = dM_dz * z if include_evolution else 0

        # True absolute magnitude
        M_std = M0 + alpha * x1 - beta * c + gamma * np.maximum(host_mass - 10, 0) + M_evol

        # Observed magnitude
        m_true = mu_true - M_std
        obs_scatter = np.sqrt(0.15**2 + sigma_int**2)
        m_b = m_true + self.rng.normal(0, obs_scatter, n_sn)

        return SNDataset(
            name=f"Simulated(n={n_sn}, evol={include_evolution})",
            z=z,
            z_cmb=z,  # Same for simulation
            m_b=m_b,
            m_b_err=np.full(n_sn, 0.15),
            x1=x1 + self.rng.normal(0, 0.1, n_sn),  # Add measurement noise
            x1_err=np.full(n_sn, 0.1),
            c=c + self.rng.normal(0, 0.02, n_sn),
            c_err=np.full(n_sn, 0.02),
            host_mass=host_mass + self.rng.normal(0, 0.1, n_sn),
            host_mass_err=np.full(n_sn, 0.1),
            sn_names=np.array([f"SIM{i:05d}" for i in range(n_sn)]),
        )


def split_by_host_mass(dataset: SNDataset, threshold: float = 10.0) -> Tuple[SNDataset, SNDataset]:
    """
    Split dataset by host galaxy stellar mass.

    This is the primary Spandrel diagnostic: if population evolution
    correlates with host mass (as expected from metallicity effects),
    the two subsamples should show different Hubble residual trends.

    Args:
        dataset: Input SNDataset
        threshold: log10(M/M_sun) threshold for split

    Returns:
        (low_mass_dataset, high_mass_dataset)
    """
    high_mass_mask = dataset.host_mass >= threshold
    low_mass_mask = ~high_mass_mask

    low_mass = dataset.select(low_mass_mask)
    high_mass = dataset.select(high_mass_mask)

    low_mass.name = f"{dataset.name}_low_mass"
    high_mass.name = f"{dataset.name}_high_mass"

    return low_mass, high_mass


def split_by_redshift(dataset: SNDataset, z_split: float = 0.5) -> Tuple[SNDataset, SNDataset]:
    """
    Split dataset by redshift.

    Useful for testing whether evolution effects strengthen at high-z.

    Args:
        dataset: Input SNDataset
        z_split: Redshift threshold

    Returns:
        (low_z_dataset, high_z_dataset)
    """
    high_z_mask = dataset.z_cmb >= z_split
    low_z_mask = ~high_z_mask

    low_z = dataset.select(low_z_mask)
    high_z = dataset.select(high_z_mask)

    low_z.name = f"{dataset.name}_z<{z_split}"
    high_z.name = f"{dataset.name}_z>={z_split}"

    return low_z, high_z


if __name__ == "__main__":
    # Test with simulated data
    loader = SimulatedDataLoader(seed=42)

    # Generate without evolution
    data_no_evol = loader.generate(n_sn=1000, include_evolution=False)
    print(f"No evolution: {data_no_evol}")
    print(f"  z range: [{data_no_evol.z.min():.3f}, {data_no_evol.z.max():.3f}]")
    print(f"  <x1>: {data_no_evol.x1.mean():.3f} ± {data_no_evol.x1.std():.3f}")

    # Generate with evolution (Spandrel scenario)
    data_evol = loader.generate(
        n_sn=1000,
        include_evolution=True,
        dM_dz=0.1,  # 0.1 mag per unit z
        dx1_dz=0.5,  # Stretch increases with z
        dc_dz=0.02,  # Color slightly redder at high z
    )
    print(f"\nWith evolution: {data_evol}")
    print(f"  <x1> (z<0.5): {data_evol.x1[data_evol.z < 0.5].mean():.3f}")
    print(f"  <x1> (z>0.5): {data_evol.x1[data_evol.z > 0.5].mean():.3f}")

    # Host mass split test
    low_mass, high_mass = split_by_host_mass(data_evol)
    print("\nHost mass split:")
    print(f"  Low mass (N={len(low_mass)}): <m_b>={low_mass.m_b.mean():.3f}")
    print(f"  High mass (N={len(high_mass)}): <m_b>={high_mass.m_b.mean():.3f}")
