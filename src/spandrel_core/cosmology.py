import numpy as np
from scipy.integrate import quad, cumulative_trapezoid
from scipy.interpolate import interp1d
from typing import Union, Optional
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
        self._dc_interp = None
        
    def hubble_parameter(self, z: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """E(z) = H(z)/H0"""
        return np.sqrt(self.Om0 * (1 + z)**3 + self.Ode0)
        
    def _comoving_integral(self, z):
        return 1.0 / self.hubble_parameter(z)
        
    def comoving_distance(self, z: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Comoving distance in Mpc.
        Handles scalar or array inputs.
        """
        if np.isscalar(z):
            inv_h, _ = quad(self._comoving_integral, 0, z)
            return (C_LIGHT_KMS / self.H0) * inv_h
        
        # Array input: Use interpolation strategy for speed
        z_arr = np.asarray(z)
        z_max = np.max(z_arr)
        
        # Check if we need to rebuild interpolator
        if z_max > self._interp_zmax or self._dc_interp is None:
            self._build_interpolator(max(z_max * 1.1, 2.0)) # At least z=2.0
            
        return (C_LIGHT_KMS / self.H0) * self._dc_interp(z_arr)
    
    def _build_interpolator(self, z_max: float, n_points: int = 1000):
        """Pre-compute comoving distance grid."""
        # Log-spaced grid + linear near 0
        z_grid = np.concatenate([
            np.linspace(0, 0.1, 100),
            np.logspace(np.log10(0.1), np.log10(z_max), n_points - 100)
        ])
        z_grid = np.unique(z_grid)
        z_grid[0] = 0.0
        
        inv_E = 1.0 / self.hubble_parameter(z_grid)
        dc_grid = cumulative_trapezoid(inv_E, z_grid, initial=0)
        
        self._dc_interp = interp1d(z_grid, dc_grid, kind='cubic', bounds_error=False, fill_value="extrapolate")
        self._interp_zmax = z_max
        
    def luminosity_distance(self, z: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Luminosity distance in Mpc."""
        return (1 + z) * self.comoving_distance(z)
        
    def angular_diameter_distance(self, z: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Angular diameter distance in Mpc."""
        return self.comoving_distance(z) / (1 + z)

def luminosity_distance(z: Union[float, np.ndarray], Om0: float = 0.3, H0: float = H0_FIDUCIAL) -> Union[float, np.ndarray]:
    """Convenience function for luminosity distance."""
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    return cosmo.luminosity_distance(z)

def comoving_distance(z: Union[float, np.ndarray], Om0: float = 0.3, H0: float = H0_FIDUCIAL) -> Union[float, np.ndarray]:
    """Convenience function for comoving distance."""
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    return cosmo.comoving_distance(z)

def angular_diameter_distance(z: Union[float, np.ndarray], Om0: float = 0.3, H0: float = H0_FIDUCIAL) -> Union[float, np.ndarray]:
    """Convenience function for angular diameter distance."""
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    return cosmo.angular_diameter_distance(z)
