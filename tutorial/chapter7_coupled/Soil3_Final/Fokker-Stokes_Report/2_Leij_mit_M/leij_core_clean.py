
"""
Lean, paper-faithful implementation (Leij et al., 2002).
- Kosugi PSD
- VG→Kosugi via Eqs. (6–8) in head space + Jurin
- Analytical r-independent kernel (constant coefficients) with log-space quadrature
- Analytical r-dependent kernel (Eq. 24), reflecting boundary, exact term structure
Units: r in µm; head h in cm.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Tuple
import numpy as np
import math

# --------------------------
# Core parameterization
# --------------------------

@dataclass
class SoilParams:
    phi_0: float   # total pore volume
    sigma: float   # ln-space stddev in r
    r_m: float     # median pore radius (µm)

    @staticmethod
    def from_van_genuchten(theta_s: float, theta_r: float, alpha: float, n: float,
                           A_cm2: float = 0.149) -> "SoilParams":
        """
        Leij et al. (2002) Eqs. (6)–(8) in head space, then Jurin r=A/h.
        m = 1 - 1/n, α in cm⁻¹.
        h0 = m**(1-m)/α
        σ_h^2 = (1-m) ln( (2^(1/m) - 1)/m )
        h_m = h0 * exp(σ_h^2)
        r_m(µm) = (A_cm2 / h_m) * 1e4,   σ_r = σ_h,   φ0 = θs - θr
        """
        if n <= 1.0 or alpha <= 0.0:
            raise ValueError("Require n>1 and alpha>0 for VG -> Kosugi.")
        phi_0 = float(theta_s - theta_r)
        m = 1.0 - 1.0/n
        h0 = (m**(1.0 - m)) / alpha
        two_pow_1_over_m_minus_1 = math.expm1(math.log(2.0) / m)
        sigma2 = (1.0 - m) * math.log(max(two_pow_1_over_m_minus_1, 1e-12) / m)
        sigma = math.sqrt(max(0.0, sigma2))
        hm = h0 * math.exp(sigma2)
        r_m_um = (A_cm2 / hm) * 1.0e4
        return SoilParams(phi_0=phi_0, sigma=float(sigma), r_m=float(r_m_um))

def kosugi_psd(r: np.ndarray | float, p: SoilParams) -> np.ndarray | float:
    """Lognormal PSD in r (µm)."""
    r = np.asarray(r, dtype=float)
    mu = math.log(p.r_m)
    inv = 1.0 / (math.sqrt(2.0*math.pi) * p.sigma)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        out = (p.phi_0 * inv / r) * np.exp(-(np.log(r) - mu)**2 / (2.0*p.sigma**2))
    out = np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)
    return float(out) if out.ndim == 0 else out

# ------------------------------------
# Analytical r-independent (constant-coeff)
# ------------------------------------

# =========================  Eq. 13 (linear-r ADE with reflecting boundary)  =========================
# Paper: Leij, Ghezzehei & Or (2002), Soil & Tillage Research
# Governing ADE in transformed time T with V = -1, D = λ, M = 0, and reflecting BC at r = 0.

def _erfc_fast(x: np.ndarray) -> np.ndarray:
    """Vectorized Abramowitz–Stegun approximation; good to ~1e-7."""
    t = 1.0 / (1.0 + 0.5 * np.abs(x))
    tau = t*np.exp(-x*x - 1.26551223 + t*(1.00002368 + t*(0.37409196 +
          t*(0.09678418 + t*(-0.18628806 + t*(0.27886807 +
          t*(-1.13520398 + t*(1.48851587 + t*(-0.82215223 + t*0.17087277)))))))))
    return np.where(x >= 0, tau, 2.0 - tau)

def _greens_kernel_eq13(r: float, xi: np.ndarray, T: float, lam: float) -> np.ndarray:
    """
    Green's function G(r,ξ,T) for Eq. 13 (linear-r ADE) with reflecting boundary at r=0:
      G = [(e^{-(r-ξ+T)^2/(4λT)} + e^{-r/λ} e^{-(r+ξ-T)^2/(4λT)}) / sqrt(4πλT)]
          + e^{-r/λ}/(2λ) * erfc((r+ξ-T)/sqrt(4λT))

    Notes:
      - V = -1 (left shift), D = λ (constant), M = 0
      - r ≥ 0, ξ ≥ 0, T ≥ 0, λ > 0
    """
    xi = np.asarray(xi, dtype=float)
    if lam <= 0.0 or T <= 0.0:
        return np.zeros_like(xi, dtype=float)

    sqrt_4lamT = math.sqrt(4.0*lam*T)
    pref = 1.0 / (math.sqrt(4.0*math.pi*lam*T))        # 1/sqrt(4πλT)
    e_r_over_lam = math.exp(-float(r)/lam)

    a1 = -((r - xi + T)**2) / (4.0*lam*T)              # direct image
    a2 = -((r + xi - T)**2) / (4.0*lam*T)              # reflected image

    with np.errstate(over="ignore"):
        term1 = np.exp(a1)
        term2 = e_r_over_lam * np.exp(a2)
        term3 = (0.5/lam) * e_r_over_lam * _erfc_fast((r + xi - T)/sqrt_4lamT)

    return pref*(term1 + term2) + term3

def analytical_solution(r: np.ndarray,
                        p0: SoilParams,
                        T: float,
                        lam: float,
                        *,
                        xi_bounds: tuple[float, float] | None = None,
                        n_quad: int = 1201) -> np.ndarray:
    """
    Paper-faithful Eq. 13 solution (linear r):
      f(r,T) = ∫_0^∞ f0(ξ) * G(r,ξ,T; λ) dξ,
    with G given by _greens_kernel_eq13 (reflecting boundary at r=0).

    Numerical evaluation:
      - We integrate in log-space (η = ln ξ) for robustness:
            dξ = ξ dη  →  ∫ f0(ξ) G dξ = ∫ f0(e^η) G(r,e^η,T) e^η dη
      - This does NOT change the PDE; it’s only a quadrature trick.

    Parameters
    ----------
    r : array of radii [µm]
    p0 : SoilParams (initial Kosugi/lognormal PSD)
    T : transformed time (≥ 0)
    lam : dispersivity λ (µm)
    xi_bounds : optional (lo, hi) bounds for ξ; if None we use ±6σ about r_m.
    n_quad : number of log-ξ points for trapezoidal integration (default 1201)
    """
    r = np.asarray(r, dtype=float)
    if T <= 0.0 or lam <= 0.0:
        return kosugi_psd(r, p0)

    # log-ξ grid covering the initial lognormal (≈ ±6σ around r_m)
    mu_ln = math.log(p0.r_m)
    if xi_bounds is None:
        L = 6.0 * p0.sigma
        lo = max(1e-6, math.exp(mu_ln - L))
        hi = max(lo*1.001, math.exp(mu_ln + L))
    else:
        lo, hi = xi_bounds
        lo = max(1e-6, float(lo))
        hi = max(lo*1.001, float(hi))

    xis = np.logspace(math.log10(lo), math.log10(hi), int(n_quad))
    dlogxi = math.log(hi/lo) / (int(n_quad) - 1)

    f0 = kosugi_psd(xis, p0)

    out = np.empty_like(r, dtype=float)
    for i, rv in enumerate(r):
        G = _greens_kernel_eq13(float(rv), xis, float(T), float(lam))
        # integrate over η = ln ξ  (Jacobian ξ)
        integrand = f0 * G * xis
        out[i] = np.trapz(integrand, dx=dlogxi)
    return out
# ===================================================================================================


# --------------------------------------
# Analytical r-dependent (Eq. 24) with xi^A placement toggle
# --------------------------------------
def analytical_solution_rdep(r: np.ndarray, p_init: SoilParams, T: float, *,
                             eps0: float = 0.10,
                             M_over_u0w: float = 0.0,
                             xi_bounds: Tuple[float,float] | None = None,
                             n_quad: int = 1201,
                             r_min_um: float = 0.05,
                             xi_power_on_term2: bool = False) -> np.ndarray:
    """
    Leij et al. (2002) Eq. 24 (paper-faithful structure):
      - u(r)=u0 r,  D=eps0*u, left drift implicit
      - reflecting boundary at r=a=r_min_um
      - mass factor exp( + (M/u0w) * T )
      - terms: two Gaussians + erfc with ξ^A factor
    Toggle:
      xi_power_on_term2=False -> use corrected form (ξ^A only on erfc term)
      xi_power_on_term2=True  -> paper-literal (ξ^A on reflected Gaussian and erfc)
    """
    if T <= 0 or eps0 <= 0:
        return kosugi_psd(r, p_init)

    r = np.asarray(r, dtype=float)
    a = float(max(r_min_um, 1e-8))

    # log-ξ grid
    if xi_bounds is None:
        lo = max(a, p_init.r_m*math.exp(-6.0*p_init.sigma))
        hi = p_init.r_m*math.exp(6.0*p_init.sigma)
    else:
        lo, hi = xi_bounds
        lo = max(lo, a)
        hi = max(hi, lo*1.001)

    xis = np.logspace(np.log10(lo), np.log10(hi), int(n_quad))
    dlogxi = math.log(hi/lo) / (int(n_quad)-1)

    f0 = kosugi_psd(xis, p_init)

    den = math.sqrt(4.0 * eps0 * T)
    gauss_pref = 1.0 / math.sqrt(4.0 * math.pi * eps0 * T)
    mu = (1.0 - eps0) * T     # Eq. 24 center shift
    Aexp  = (1.0 - eps0) / eps0
    C_erfc = (1.0 - eps0) / (2.0 * eps0)
    log_a2 = 2.0 * math.log(a)

    def erfc_fast(x):
        t = 1.0 / (1.0 + 0.5 * np.abs(x))
        tau = t*np.exp(-x*x - 1.26551223 + t*(1.00002368 + t*(0.37409196 +
              t*(0.09678418 + t*(-0.18628806 + t*(0.27886807 +
              t*(-1.13520398 + t*(1.48851587 + t*(-0.82215223 + t*0.17087277)))))))))
        return np.where(x >= 0, tau, 2.0 - tau)

    out = np.empty_like(r, dtype=float)
    mass_fac = math.exp(M_over_u0w * T)  # Eq. 24 sign

    xi_pow = np.exp(Aexp * np.log(xis))  # ξ^A

    for i, rr0 in enumerate(r):
        rr = float(max(rr0, a))
        z1 = (math.log(rr) - np.log(xis) - mu) / den
        z2 = (math.log(rr) + np.log(xis) - log_a2 - mu) / den

        term1 = gauss_pref * np.exp(-z1*z1)
        if xi_power_on_term2:
            term2 = gauss_pref * xi_pow * np.exp(-z2*z2)   # paper-literal
        else:
            term2 = gauss_pref * np.exp(-z2*z2)            # corrected
        term3 = C_erfc * xi_pow * erfc_fast(z2)            # ξ^A always here

        kernel = (term1 + term2 + term3) / rr
        integrand = kernel * f0 * xis
        out[i] = mass_fac * np.trapz(integrand, dx=dlogxi)

    return out
