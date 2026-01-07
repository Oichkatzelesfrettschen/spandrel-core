from __future__ import annotations

from typing import Callable, cast

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad
from scipy.interpolate import interp1d

from spandrel_core.constants import C_LIGHT_KMS, H0_FIDUCIAL


class FlatLambdaCDM:
    """
    Simple flat LambdaCDM cosmology calculator.
    """
    def __init__(self, H0: float = H0_FIDUCIAL, Om0: float = 0.3):
        self.H0 = H0
        self.Om0 = Om0
        self.Ode0 = 1.0 - Om0
        self._interp_zmax = -1.0
        self._dc_interp: Callable[[np.ndarray], np.ndarray] | None = None

    def hubble_parameter(self, z: float | np.ndarray) -> float | np.ndarray:
        """E(z) = H(z)/H0"""
        hz = np.sqrt(self.Om0 * (1 + z) ** 3 + self.Ode0)
        if np.isscalar(z):
            return float(hz)
        return cast(np.ndarray, hz)

    def _comoving_integral(self, z: float) -> float:
        return 1.0 / float(self.hubble_parameter(z))

    def comoving_distance(self, z: float | np.ndarray) -> float | np.ndarray:
        """
        Comoving distance in Mpc.
        Handles scalar or array inputs.
        """
        if np.isscalar(z):
            z_float = float(z)
            inv_h, _ = quad(self._comoving_integral, 0, z_float)
            return float((C_LIGHT_KMS / self.H0) * float(inv_h))

        # Array input: Use interpolation strategy for speed
        z_arr = cast(np.ndarray, np.asarray(z, dtype=float))
        z_max = float(np.max(z_arr))

        # Check if we need to rebuild interpolator
        if z_max > self._interp_zmax or self._dc_interp is None:
            self._build_interpolator(max(z_max * 1.1, 2.0))  # At least z=2.0

        interp = self._dc_interp
        assert interp is not None
        scaled = (C_LIGHT_KMS / self.H0) * cast(np.ndarray, interp(z_arr))
        return cast(np.ndarray, scaled)

    def _build_interpolator(self, z_max: float, n_points: int = 1000):
        """Pre-compute comoving distance grid."""
        # Log-spaced grid + linear near 0
        z_grid = np.concatenate(
            [
                np.linspace(0, 0.1, 100),
                np.logspace(np.log10(0.1), np.log10(z_max), n_points - 100),
            ]
        )
        z_grid = cast(np.ndarray, np.unique(z_grid))
        z_grid[0] = 0.0

        inv_E = 1.0 / self.hubble_parameter(z_grid)
        dc_grid = cast(np.ndarray, cumulative_trapezoid(inv_E, z_grid, initial=0))

        interp = interp1d(
            z_grid, dc_grid, kind="cubic", bounds_error=False, fill_value="extrapolate"
        )
        self._dc_interp = cast(Callable[[np.ndarray], np.ndarray], interp)
        self._interp_zmax = z_max

    def luminosity_distance(self, z: float | np.ndarray) -> float | np.ndarray:
        """Luminosity distance in Mpc."""
        return (1 + z) * self.comoving_distance(z)

    def angular_diameter_distance(self, z: float | np.ndarray) -> float | np.ndarray:
        """Angular diameter distance in Mpc."""
        return self.comoving_distance(z) / (1 + z)


def luminosity_distance(
    z: float | np.ndarray, Om0: float = 0.3, H0: float = H0_FIDUCIAL
) -> float | np.ndarray:
    """Convenience function for luminosity distance."""
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    return cosmo.luminosity_distance(z)


def comoving_distance(
    z: float | np.ndarray, Om0: float = 0.3, H0: float = H0_FIDUCIAL
) -> float | np.ndarray:
    """Convenience function for comoving distance."""
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    return cosmo.comoving_distance(z)


def angular_diameter_distance(
    z: float | np.ndarray, Om0: float = 0.3, H0: float = H0_FIDUCIAL
) -> float | np.ndarray:
    """Convenience function for angular diameter distance."""
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    return cosmo.angular_diameter_distance(z)


def distance_modulus_from_luminosity_distance(dl_mpc: float | np.ndarray) -> float | np.ndarray:
    """Distance modulus from luminosity distance in Mpc.

    μ = 5 log10(D_L / 10 pc) where D_L is in parsecs.
    """
    mu = 5.0 * np.log10(np.asarray(dl_mpc, dtype=float) * 1e6 / 10.0)
    if np.isscalar(dl_mpc):
        return float(mu)
    return cast(np.ndarray, mu)


def distance_modulus_lcdm(
    z: float | np.ndarray, *, Om0: float = 0.3, H0: float = H0_FIDUCIAL
) -> float | np.ndarray:
    """LCDM distance modulus μ(z) using FlatLambdaCDM."""
    dl = luminosity_distance(z, Om0=Om0, H0=H0)
    return distance_modulus_from_luminosity_distance(dl)
