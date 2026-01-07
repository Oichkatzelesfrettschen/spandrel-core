"""Pantheon+SH0ES dataset interface.

This module lives in `spandrel-core` so other consumers can reuse a consistent,
typed interface. The caller is responsible for providing the dataset file path.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from spandrel_core.data import PantheonPlusLoader


@dataclass(frozen=True)
class DataStats:
    total_raw: int
    total_valid: int
    z_min: float
    z_max: float
    mu_min: float
    mu_max: float
    n_surveys: int
    n_calibrators: int


class PantheonData:
    """Validated interface around the Pantheon+SH0ES distance-modulus data."""

    def __init__(
        self,
        filepath: Path,
        *,
        z_min: float = 0.001,
        z_max: float = 2.5,
    ) -> None:
        self.filepath = Path(filepath)
        self.z_min = z_min
        self.z_max = z_max

        self._load_data()
        self._validate_data()
        self._compute_stats()

    def _load_data(self) -> None:
        if not self.filepath.exists():
            raise FileNotFoundError(f"Data file not found: {self.filepath}")

        self._raw_df = pd.read_csv(self.filepath, sep=r"\s+", comment="#")
        self._total_raw = int(len(self._raw_df))

    def _validate_data(self) -> None:
        loader = PantheonPlusLoader(data_dir=self.filepath.parent)
        z, mu, mu_err, df_cut = loader.load_distance_modulus(z_min=self.z_min, z_max=self.z_max)

        self.dataframe = df_cut
        self._z = z
        self._mu = mu
        self._mu_err = mu_err
        self._n_rejected = int(len(self._raw_df) - len(self.dataframe))

    def _compute_stats(self) -> None:
        df = self.dataframe
        self.stats = DataStats(
            total_raw=self._total_raw,
            total_valid=int(len(df)),
            z_min=float(df["zHD"].min()),
            z_max=float(df["zHD"].max()),
            mu_min=float(df["MU_SH0ES"].min()),
            mu_max=float(df["MU_SH0ES"].max()),
            n_surveys=int(df["IDSURVEY"].nunique()),
            n_calibrators=int((df["IS_CALIBRATOR"] == 1).sum()),
        )

    def get_cosmology_data(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return self._z, self._mu, self._mu_err

    def get_calibrator_subset(self) -> pd.DataFrame:
        return self.dataframe[self.dataframe["IS_CALIBRATOR"] == 1].copy()

    def get_hubble_flow_subset(self, z_cut: float = 0.01) -> pd.DataFrame:
        return self.dataframe[self.dataframe["zHD"] > z_cut].copy()

    def validate(self) -> dict[str, object]:
        z, mu, mu_err = self.get_cosmology_data()
        df = self.dataframe

        return {
            "total_entries": int(len(df)),
            "rejected_entries": int(self._n_rejected),
            "redshift_range": (float(z.min()), float(z.max())),
            "mu_range": (float(mu.min()), float(mu.max())),
            "has_nan_z": bool(np.isnan(z).any()),
            "has_nan_mu": bool(np.isnan(mu).any()),
            "has_negative_errors": bool((mu_err <= 0).any()),
            "surveys_present": df["IDSURVEY"].unique().tolist(),
            "calibrator_count": int((df["IS_CALIBRATOR"] == 1).sum()),
            "median_z": float(np.median(z)),
            "median_mu_err": float(np.median(mu_err)),
        }

    def __len__(self) -> int:
        return int(len(self.dataframe))

    def __repr__(self) -> str:
        return (
            f"PantheonData(n={self.stats.total_valid}, "
            f"z=[{self.stats.z_min:.4f}, {self.stats.z_max:.4f}], "
            f"surveys={self.stats.n_surveys})"
        )


def load_pantheon(
    filepath: Path,
    *,
    z_min: float = 0.001,
    z_max: float = 2.5,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = PantheonData(filepath=filepath, z_min=z_min, z_max=z_max)
    return data.get_cosmology_data()


def try_default_pantheon_file() -> Optional[Path]:
    """Best-effort discovery for monorepo checkouts; returns None if not found."""
    repo_root = Path(__file__).resolve().parents[3]
    cand = repo_root / "pantheon" / "data" / "Pantheon+SH0ES.dat"
    return cand if cand.exists() else None

