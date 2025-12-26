#!/usr/bin/env python3
"""
Spandrel Likelihood Module

Implements the full Bayesian likelihood for the Spandrel hypothesis test:
comparing baseline SN standardization vs. population evolution models.

This module is designed for use with emcee or cobaya samplers.
"""

from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np
from scipy.integrate import quad
from scipy.stats import norm


@dataclass
class CosmologyParams:
    """Cosmological parameters for distance modulus calculation."""

    H0: float = 70.0  # km/s/Mpc
    Om: float = 0.3  # Matter density
    w0: float = -1.0  # DE equation of state at z=0
    wa: float = 0.0  # DE evolution parameter


@dataclass
class StandardizationParams:
    """SN standardization parameters (SALT2-like)."""

    M0: float = -19.3  # Baseline absolute magnitude
    alpha: float = 0.14  # Stretch coefficient
    beta: float = 3.1  # Color coefficient
    gamma: float = 0.05  # Host mass step coefficient
    sigma_int: float = 0.1  # Intrinsic scatter


@dataclass
class EvolutionParams:
    """Population evolution parameters (Spandrel hypothesis)."""

    dM_dz: float = 0.0  # Luminosity evolution: M(z) = M0 + dM_dz * z
    dx1_dz: float = 0.0  # Stretch drift: <x1>(z) = x1_0 + dx1_dz * z
    dc_dz: float = 0.0  # Color drift: <c>(z) = c_0 + dc_dz * z

    # Population scatter parameters
    sigma_x1: float = 1.0  # Intrinsic stretch scatter
    sigma_c: float = 0.1  # Intrinsic color scatter


class SpandrelLikelihood:
    """
    Full Bayesian likelihood for Spandrel hypothesis testing.

    Compares two models:
    1. Baseline: Standard SALT2 standardization + cosmology
    2. Spandrel: Baseline + redshift-dependent population evolution

    The key diagnostic is whether including evolution terms changes
    the inferred dark energy parameters (w0, wa).
    """

    def __init__(
        self,
        z: np.ndarray,
        m_obs: np.ndarray,
        m_err: np.ndarray,
        x1: np.ndarray,
        x1_err: np.ndarray,
        c: np.ndarray,
        c_err: np.ndarray,
        host_mass: np.ndarray,
        host_mass_err: np.ndarray,
        cov_matrix: Optional[np.ndarray] = None,
    ):
        """
        Initialize likelihood with SN data.

        Args:
            z: Redshifts (N,)
            m_obs: Observed peak B-band magnitudes (N,)
            m_err: Magnitude uncertainties (N,)
            x1: SALT2 stretch parameters (N,)
            x1_err: Stretch uncertainties (N,)
            c: SALT2 color parameters (N,)
            c_err: Color uncertainties (N,)
            host_mass: log10(M_host/M_sun) (N,)
            host_mass_err: Host mass uncertainties (N,)
            cov_matrix: Full covariance matrix (N, N) or None for diagonal
        """
        self.z = np.asarray(z)
        self.m_obs = np.asarray(m_obs)
        self.m_err = np.asarray(m_err)
        self.x1 = np.asarray(x1)
        self.x1_err = np.asarray(x1_err)
        self.c = np.asarray(c)
        self.c_err = np.asarray(c_err)
        self.host_mass = np.asarray(host_mass)
        self.host_mass_err = np.asarray(host_mass_err)

        self.n_sn = len(z)
        self.cov_matrix = cov_matrix

        # Precompute things that don't change
        self._validate_data()

    def _validate_data(self):
        """Validate data arrays have consistent shapes."""
        arrays = [
            self.z,
            self.m_obs,
            self.m_err,
            self.x1,
            self.x1_err,
            self.c,
            self.c_err,
            self.host_mass,
            self.host_mass_err,
        ]
        lengths = [len(a) for a in arrays]
        if len(set(lengths)) != 1:
            raise ValueError(f"Inconsistent array lengths: {lengths}")

    def distance_modulus(self, z: np.ndarray, cosmo: CosmologyParams) -> np.ndarray:
        """
        Compute distance modulus for flat w0wa cosmology.

        Uses numerical integration for accuracy.
        """
        c_light = 299792.458  # km/s

        def E_inv(zp):
            """1/E(z) for integration."""
            Ode = 1.0 - cosmo.Om
            # w(z) = w0 + wa * z/(1+z)
            w_z = cosmo.w0 + cosmo.wa * zp / (1 + zp)
            de_term = (
                Ode
                * (1 + zp) ** (3 * (1 + cosmo.w0 + cosmo.wa))
                * np.exp(-3 * cosmo.wa * zp / (1 + zp))
            )
            return 1.0 / np.sqrt(cosmo.Om * (1 + zp) ** 3 + de_term)

        # Comoving distance
        dc = np.zeros_like(z)
        for i, zi in enumerate(z):
            if zi > 0:
                integral, _ = quad(E_inv, 0, zi, limit=100)
                dc[i] = (c_light / cosmo.H0) * integral

        # Luminosity distance (flat universe)
        dl = dc * (1 + z)

        # Distance modulus (avoiding log of zero)
        with np.errstate(divide="ignore", invalid="ignore"):
            mu = np.where(dl > 0, 5 * np.log10(dl) + 25, 0)

        return mu

    def standardized_magnitude(
        self, std: StandardizationParams, evo: Optional[EvolutionParams] = None
    ) -> np.ndarray:
        """
        Compute standardized absolute magnitude for each SN.

        M_std = M0 + alpha*x1 - beta*c + gamma*max(log_mass - 10, 0) + evolution
        """
        # Host mass step (only for massive hosts)
        mass_step = std.gamma * np.maximum(self.host_mass - 10.0, 0)

        # Baseline standardization
        M_std = std.M0 + std.alpha * self.x1 - std.beta * self.c + mass_step

        # Population evolution (Spandrel terms)
        if evo is not None:
            M_std += evo.dM_dz * self.z

        return M_std

    def log_likelihood_baseline(self, theta: np.ndarray) -> float:
        """
        Log-likelihood for baseline model (no population evolution).

        theta = [M0, alpha, beta, gamma, sigma_int, H0, Om, w0, wa]
        """
        if len(theta) != 9:
            raise ValueError(f"Expected 9 parameters, got {len(theta)}")

        M0, alpha, beta, gamma, sigma_int, H0, Om, w0, wa = theta

        # Parameter bounds
        if sigma_int <= 0 or H0 <= 0 or Om <= 0 or Om >= 1:
            return -np.inf

        cosmo = CosmologyParams(H0=H0, Om=Om, w0=w0, wa=wa)
        std = StandardizationParams(M0=M0, alpha=alpha, beta=beta, gamma=gamma, sigma_int=sigma_int)

        # Model prediction
        mu_theory = self.distance_modulus(self.z, cosmo)
        M_std = self.standardized_magnitude(std, evo=None)
        m_pred = mu_theory - M_std

        # Residuals
        residuals = self.m_obs - m_pred

        # Covariance (diagonal approximation for now)
        # Full treatment would propagate x1, c uncertainties
        var_total = (
            self.m_err**2
            + (std.alpha * self.x1_err) ** 2
            + (std.beta * self.c_err) ** 2
            + sigma_int**2
        )

        # Gaussian log-likelihood
        log_lik = -0.5 * np.sum(residuals**2 / var_total + np.log(2 * np.pi * var_total))

        return log_lik

    def log_likelihood_spandrel(self, theta: np.ndarray) -> float:
        """
        Log-likelihood for Spandrel model (with population evolution).

        theta = [M0, alpha, beta, gamma, sigma_int, dM_dz, dx1_dz, dc_dz, H0, Om, w0, wa]
        """
        if len(theta) != 12:
            raise ValueError(f"Expected 12 parameters, got {len(theta)}")

        (M0, alpha, beta, gamma, sigma_int, dM_dz, dx1_dz, dc_dz, H0, Om, w0, wa) = theta

        # Parameter bounds
        if sigma_int <= 0 or H0 <= 0 or Om <= 0 or Om >= 1:
            return -np.inf

        cosmo = CosmologyParams(H0=H0, Om=Om, w0=w0, wa=wa)
        std = StandardizationParams(M0=M0, alpha=alpha, beta=beta, gamma=gamma, sigma_int=sigma_int)
        evo = EvolutionParams(dM_dz=dM_dz, dx1_dz=dx1_dz, dc_dz=dc_dz)

        # Model prediction with evolution
        mu_theory = self.distance_modulus(self.z, cosmo)
        M_std = self.standardized_magnitude(std, evo=evo)
        m_pred = mu_theory - M_std

        # Residuals
        residuals = self.m_obs - m_pred

        # Covariance
        var_total = (
            self.m_err**2
            + (std.alpha * self.x1_err) ** 2
            + (std.beta * self.c_err) ** 2
            + sigma_int**2
        )

        # Gaussian log-likelihood
        log_lik = -0.5 * np.sum(residuals**2 / var_total + np.log(2 * np.pi * var_total))

        return log_lik

    def log_prior_baseline(self, theta: np.ndarray) -> float:
        """Priors for baseline model parameters."""
        M0, alpha, beta, gamma, sigma_int, H0, Om, w0, wa = theta

        log_p = 0.0

        # M0: Gaussian prior
        log_p += norm.logpdf(M0, loc=-19.3, scale=0.5)

        # alpha, beta: Gaussian priors centered on typical values
        log_p += norm.logpdf(alpha, loc=0.14, scale=0.05)
        log_p += norm.logpdf(beta, loc=3.1, scale=0.5)

        # gamma: Gaussian centered on zero
        log_p += norm.logpdf(gamma, loc=0.05, scale=0.05)

        # sigma_int: Half-normal
        if sigma_int <= 0:
            return -np.inf
        log_p += norm.logpdf(sigma_int, loc=0, scale=0.2)

        # H0: Gaussian
        log_p += norm.logpdf(H0, loc=70, scale=10)

        # Om: Uniform on [0.1, 0.5]
        if Om < 0.1 or Om > 0.5:
            return -np.inf

        # w0: Uniform on [-2, 0]
        if w0 < -2 or w0 > 0:
            return -np.inf

        # wa: Gaussian centered on zero
        log_p += norm.logpdf(wa, loc=0, scale=1.0)

        return log_p

    def log_prior_spandrel(self, theta: np.ndarray) -> float:
        """Priors for Spandrel model parameters."""
        (M0, alpha, beta, gamma, sigma_int, dM_dz, dx1_dz, dc_dz, H0, Om, w0, wa) = theta

        # Start with baseline priors
        baseline_theta = [M0, alpha, beta, gamma, sigma_int, H0, Om, w0, wa]
        log_p = self.log_prior_baseline(baseline_theta)

        if not np.isfinite(log_p):
            return log_p

        # Evolution parameters: Gaussian priors centered on zero
        log_p += norm.logpdf(dM_dz, loc=0, scale=0.2)
        log_p += norm.logpdf(dx1_dz, loc=0, scale=0.5)
        log_p += norm.logpdf(dc_dz, loc=0, scale=0.1)

        return log_p

    def log_posterior_baseline(self, theta: np.ndarray) -> float:
        """Log posterior for baseline model."""
        log_prior = self.log_prior_baseline(theta)
        if not np.isfinite(log_prior):
            return -np.inf
        return log_prior + self.log_likelihood_baseline(theta)

    def log_posterior_spandrel(self, theta: np.ndarray) -> float:
        """Log posterior for Spandrel model."""
        log_prior = self.log_prior_spandrel(theta)
        if not np.isfinite(log_prior):
            return -np.inf
        return log_prior + self.log_likelihood_spandrel(theta)


def compute_model_comparison(
    likelihood: SpandrelLikelihood, baseline_samples: np.ndarray, spandrel_samples: np.ndarray
) -> Dict:
    """
    Compute model comparison statistics between baseline and Spandrel.

    Args:
        likelihood: SpandrelLikelihood instance
        baseline_samples: MCMC samples from baseline model (n_samples, 9)
        spandrel_samples: MCMC samples from Spandrel model (n_samples, 12)

    Returns:
        Dictionary with comparison statistics
    """
    results = {}

    # Maximum likelihood estimates
    baseline_logliks = np.array([likelihood.log_likelihood_baseline(s) for s in baseline_samples])
    spandrel_logliks = np.array([likelihood.log_likelihood_spandrel(s) for s in spandrel_samples])

    results["baseline_max_loglik"] = np.max(baseline_logliks)
    results["spandrel_max_loglik"] = np.max(spandrel_logliks)

    # BIC approximation
    n_data = likelihood.n_sn
    k_baseline = 9
    k_spandrel = 12

    results["BIC_baseline"] = k_baseline * np.log(n_data) - 2 * results["baseline_max_loglik"]
    results["BIC_spandrel"] = k_spandrel * np.log(n_data) - 2 * results["spandrel_max_loglik"]
    results["delta_BIC"] = results["BIC_spandrel"] - results["BIC_baseline"]

    # Parameter constraints
    results["w0_baseline"] = {
        "mean": np.mean(baseline_samples[:, 7]),
        "std": np.std(baseline_samples[:, 7]),
    }
    results["wa_baseline"] = {
        "mean": np.mean(baseline_samples[:, 8]),
        "std": np.std(baseline_samples[:, 8]),
    }
    results["w0_spandrel"] = {
        "mean": np.mean(spandrel_samples[:, 10]),
        "std": np.std(spandrel_samples[:, 10]),
    }
    results["wa_spandrel"] = {
        "mean": np.mean(spandrel_samples[:, 11]),
        "std": np.std(spandrel_samples[:, 11]),
    }

    # Shift in w0, wa
    delta_w0 = results["w0_spandrel"]["mean"] - results["w0_baseline"]["mean"]
    delta_wa = results["wa_spandrel"]["mean"] - results["wa_baseline"]["mean"]

    # Combined uncertainty
    sigma_w0 = np.sqrt(results["w0_baseline"]["std"] ** 2 + results["w0_spandrel"]["std"] ** 2)
    sigma_wa = np.sqrt(results["wa_baseline"]["std"] ** 2 + results["wa_spandrel"]["std"] ** 2)

    results["delta_w0"] = {"value": delta_w0, "significance": np.abs(delta_w0) / sigma_w0}
    results["delta_wa"] = {"value": delta_wa, "significance": np.abs(delta_wa) / sigma_wa}

    # Evolution parameter constraints
    results["dM_dz"] = {
        "mean": np.mean(spandrel_samples[:, 5]),
        "std": np.std(spandrel_samples[:, 5]),
    }

    return results


if __name__ == "__main__":
    # Test with simulated data
    np.random.seed(42)
    n_sn = 500

    # Simulate SN sample
    z = np.random.uniform(0.01, 1.2, n_sn)
    x1 = np.random.normal(0, 1, n_sn)
    c = np.random.normal(0, 0.1, n_sn)
    host_mass = np.random.normal(10.5, 0.5, n_sn)

    # True parameters
    M0_true = -19.3
    alpha_true = 0.14
    beta_true = 3.1
    gamma_true = 0.05

    # Distance modulus (simplified)
    mu_true = 5 * np.log10((1 + z) * 4285.7 * z) + 25

    # Simulated magnitudes
    M_std = M0_true + alpha_true * x1 - beta_true * c + gamma_true * np.maximum(host_mass - 10, 0)
    m_true = mu_true - M_std
    m_obs = m_true + np.random.normal(0, 0.15, n_sn)

    # Create likelihood
    lik = SpandrelLikelihood(
        z=z,
        m_obs=m_obs,
        m_err=np.full(n_sn, 0.15),
        x1=x1,
        x1_err=np.full(n_sn, 0.1),
        c=c,
        c_err=np.full(n_sn, 0.02),
        host_mass=host_mass,
        host_mass_err=np.full(n_sn, 0.1),
    )

    # Test likelihood evaluation
    theta_baseline = [-19.3, 0.14, 3.1, 0.05, 0.1, 70.0, 0.3, -1.0, 0.0]
    theta_spandrel = [-19.3, 0.14, 3.1, 0.05, 0.1, 0.0, 0.0, 0.0, 70.0, 0.3, -1.0, 0.0]

    loglik_base = lik.log_likelihood_baseline(theta_baseline)
    loglik_span = lik.log_likelihood_spandrel(theta_spandrel)

    print("Spandrel Likelihood Test")
    print("=" * 50)
    print(f"N_SN: {n_sn}")
    print(f"Baseline log-likelihood: {loglik_base:.2f}")
    print(f"Spandrel log-likelihood: {loglik_span:.2f}")
    print(f"Log-posterior baseline: {lik.log_posterior_baseline(theta_baseline):.2f}")
    print(f"Log-posterior spandrel: {lik.log_posterior_spandrel(theta_spandrel):.2f}")
