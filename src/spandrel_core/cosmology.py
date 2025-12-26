import numpy as np
from scipy.integrate import quad
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
        
    def hubble_parameter(self, z: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """E(z) = H(z)/H0"""
        return np.sqrt(self.Om0 * (1 + z)**3 + self.Ode0)
        
    def _comoving_integral(self, z):
        return 1.0 / self.hubble_parameter(z)
        
    def comoving_distance(self, z: float) -> float:
        """Comoving distance in Mpc/h if H0=100, or Mpc if H0 is provided."""
        inv_h, _ = quad(self._comoving_integral, 0, z)
        return (C_LIGHT_KMS / self.H0) * inv_h
        
    def luminosity_distance(self, z: float) -> float:
        """Luminosity distance in Mpc."""
        return (1 + z) * self.comoving_distance(z)
        
    def angular_diameter_distance(self, z: float) -> float:
        """Angular diameter distance in Mpc."""
        return self.comoving_distance(z) / (1 + z)

def luminosity_distance(z: float, Om0: float = 0.3, H0: float = H0_FIDUCIAL) -> float:
    """Convenience function for luminosity distance."""
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    return cosmo.luminosity_distance(z)

def comoving_distance(z: float, Om0: float = 0.3, H0: float = H0_FIDUCIAL) -> float:
    """Convenience function for comoving distance."""
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    return cosmo.comoving_distance(z)

def angular_diameter_distance(z: float, Om0: float = 0.3, H0: float = H0_FIDUCIAL) -> float:
    """Convenience function for angular diameter distance."""
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    return cosmo.angular_diameter_distance(z)
