# =============================================================================
#
#       Implementation of Leij et al. (2002) Pore Evolution Model - Phase 1
#       Author: Gemini AI
#       Date: 2025-10-17
#       Description: This script models the evolution of soil pore-size 
#                    distribution (PSD) after tillage, assuming no pore
#                    degradation (M=0). It uses the Fokker-Planck analytical
#                    solution and fits for the dispersivity parameter (lambda).
#
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import curve_fit

# =============================================================================
# ## PART 1: HELPER FUNCTIONS
# =============================================================================

def vg_to_kosugi(alpha, n, theta_s, theta_r):
    """
    Converts Van Genuchten parameters to Kosugi parameters.
    Based on Eqs. 6-103 from Chandrasekhar et al. (2019).
    """
    m = 1 - (1 / n)
    # Proportionality constant from Young-Laplace equation [cm^2]
    A = 0.149
    
    sigma_sq = (1 - m) * np.log((2**(1/m) - 1) / m)
    sigma = np.sqrt(sigma_sq)
    
    h0 = (m**(1 - m)) / alpha
    hm = h0 * np.exp(sigma_sq)
    
    # Convert head (hm) to median pore radius (rm) in micrometers (um)
    r_m_cm = A / hm
    r_m_um = r_m_cm * 10000
    
    return {'sigma': sigma, 'r_m': r_m_um, 'theta_s': theta_s, 'theta_r': theta_r}

def kosugi_psd(r, params):
    """
    Calculates the pore-size distribution f(r) using Kosugi's lognormal model.
    Based on Eq. 5 from Chandrasekhar et al. (2019).
    """
    phi_0 = params['theta_s'] - params['theta_r']
    sigma = params['sigma']
    r_m = params['r_m']
    
    r = np.maximum(r, 1e-9) # Avoid log(0) or division by zero
    
    term1 = phi_0 / (r * sigma * np.sqrt(2 * np.pi))
    exponent = -((np.log(r / r_m))**2) / (2 * sigma**2)
    
    return term1 * np.exp(exponent)

def calculate_moments(params):
    """
    Calculates the zero-order (m0) and first-order (<r>) moments of the PSD.
    Based on Eqs. 10 & 11 from Chandrasekhar et al. (2019).
    """
    phi_0 = params['theta_s'] - params['theta_r']
    sigma = params['sigma']
    r_m = params['r_m']
    
    m0 = phi_0
    r_mean = r_m * np.exp(0.5 * sigma**2)
    
    return m0, r_mean

# =============================================================================
# ## PART 2: THE CORE MODEL (ANALYTICAL SOLUTION)
# =============================================================================

from scipy.special import erfc # <-- Make sure this import is at the top of your script

def analytical_solution(r, initial_params, T, lambda_val):
    """
    The DEFINITIVE analytical solution of the Fokker-Planck equation.
    This corrected version faithfully transcribes the complete Green's function.
    """
    
    def integrand_func(xi, r_val):
        """ The function inside the integral, multiplied by f0(xi). """
        f0_xi = kosugi_psd(xi, initial_params)
        
        # Guard against division by zero or invalid math operations
        if lambda_val <= 1e-9 or T <= 1e-9:
            return 0
        
        # --- Green's Function G(r, xi, T) ---
        
        # Part 1 & 2: The exponential terms
        sqrt_4_pi_lambda_T = np.sqrt(4 * np.pi * lambda_val * T)
        
        exp1_num = -((r_val - xi + T)**2)
        exp1_den = 4 * lambda_val * T
        term1 = np.exp(exp1_num / exp1_den)
        
        exp2_num = -((r_val + xi - T)**2)
        exp2_den = 4 * lambda_val * T
        term2 = np.exp(-(r_val / lambda_val) + (exp2_num / exp2_den))

        # Part 3: The erfc (complementary error function) term
        sqrt_4_lambda_T = np.sqrt(4 * lambda_val * T)
        erfc_arg = (r_val + xi - T) / sqrt_4_lambda_T
        term3_factor = 1 / (2 * lambda_val)
        term3 = term3_factor * np.exp(-r_val / lambda_val) * erfc(erfc_arg)

        # Combine all parts of the Green's function
        G = ((term1 + term2) / sqrt_4_pi_lambda_T) + term3
        
        return f0_xi * G

    final_psd = np.zeros_like(r)
    for i, r_val in enumerate(r):
        # We integrate over a reasonable range of pore sizes
        result, _ = quad(integrand_func, 0.01, 5000, args=(r_val,), limit=100)
        final_psd[i] = result
        
    return final_psd

# =============================================================================
# ## PART 3: FITTING AND EXECUTION SCRIPT
# =============================================================================
def run_simulation():
    """
    Runs simulations for both panels of Leij et al. (2002) Fig. 3,
    prints metrics, and plots a two-panel figure.
    """
    
    # =========================================================================
    # --- SIMULATION FOR PANEL (a) ---
    # =========================================================================
    print("Running Simulation for Figure 3(a)...")
    
    # Data from Leij et al. (2002), Soil & Tillage Research, Table 1
    # Case: Moldboard (1350) -> Zero Tillage (1340)
    initial_a = {'theta_s': 0.469, 'theta_r': 0.104, 'sigma': 2.27, 'r_m': 6.36}
    final_a = {'theta_s': 0.392, 'theta_r': 0.097, 'sigma': 2.33, 'r_m': 2.15}
    T_a = 20.5  # From Table 1

    # Define radius values for calculation
    r_values = np.logspace(-1, 2, 50)
    observed_psd_a = kosugi_psd(r_values, final_a)

    def model_to_fit_a(r, lambda_val):
        return analytical_solution(r, initial_a, T_a, lambda_val)
    
    popt_a, _ = curve_fit(model_to_fit_a, r_values, observed_psd_a, p0=[1.3], bounds=(0, np.inf))
    lambda_a = popt_a[0]
    predicted_psd_a = model_to_fit_a(r_values, lambda_a)

    # --- Print metrics for Panel (a) ---
    print("\n" + "="*50)
    print("      METRICS FOR FIGURE 3(a)")
    print("="*50)
    print(f"  - Target Lambda λ (µm):   1.30")
    print(f"  - Fitted Lambda λ (µm):   {lambda_a:.4f}")
    print("="*50 + "\n")

    # =========================================================================
    # --- SIMULATION FOR PANEL (b) ---
    # =========================================================================
    print("Running Simulation for Figure 3(b)...")
    
    # Data from Leij et al. (2002), Soil & Tillage Research, Table 1
    # Case: Moldboard (2010) -> Reduced Tillage (2000)
    initial_b = {'theta_s': 0.394, 'theta_r': 0.142, 'sigma': 2.60, 'r_m': 3.97}
    final_b = {'theta_s': 0.367, 'theta_r': 0.128, 'sigma': 2.44, 'r_m': 0.997}
    T_b = 24.9  # From Table 1
    
    observed_psd_b = kosugi_psd(r_values, final_b)
    
    def model_to_fit_b(r, lambda_val):
        return analytical_solution(r, initial_b, T_b, lambda_val)

    popt_b, _ = curve_fit(model_to_fit_b, r_values, observed_psd_b, p0=[0.96], bounds=(0, np.inf))
    lambda_b = popt_b[0]
    predicted_psd_b = model_to_fit_b(r_values, lambda_b)

    # --- Print metrics for Panel (b) ---
    print("\n" + "="*50)
    print("      METRICS FOR FIGURE 3(b)")
    print("="*50)
    print(f"  - Target Lambda λ (µm):   0.96")
    print(f"  - Fitted Lambda λ (µm):   {lambda_b:.4f}")
    print("="*50 + "\n")

    # =========================================================================
    # --- CREATE THE TWO-PANEL PLOT ---
    # =========================================================================
    import matplotlib.ticker as mticker

    plt.style.use('seaborn-v0_8-whitegrid')
    # Create a figure with 2 rows, 1 column of subplots
    fig, ax = plt.subplots(2, 1, figsize=(8, 12))

    # --- Plot Panel (a) on the top axis (ax[0]) ---
    ax[0].plot(r_values, observed_psd_a, 'o', color='black', markerfacecolor='black', markersize=6, label='PSD 1340 observed')
    ax[0].plot(r_values, kosugi_psd(r_values, initial_a), color='black', linestyle='--', label='1350')
    ax[0].plot(r_values, predicted_psd_a, color='black', linestyle='-', linewidth=1.5, label='1340 predicted')
    
    # --- Plot Panel (b) on the bottom axis (ax[1]) ---
    ax[1].plot(r_values, observed_psd_b, 'o', color='black', markerfacecolor='black', markersize=6, label='PSD 2000 observed')
    ax[1].plot(r_values, kosugi_psd(r_values, initial_b), color='black', linestyle='--', label='2010')
    ax[1].plot(r_values, predicted_psd_b, color='black', linestyle='-', linewidth=1.5, label='2000 predicted')

    # --- Format both plots to match the paper ---
    for i, panel_label in enumerate(['(a)', '(b)']):
        ax[i].set_xscale('log')
        ax[i].set_xlim(0.1, 100)
        ax[i].xaxis.set_major_locator(mticker.LogLocator(base=10.0, numticks=4))
        ax[i].xaxis.set_major_formatter(mticker.ScalarFormatter())
        ax[i].set_xlabel('r (μm)', fontsize=12)
        ax[i].set_ylabel('f (μm⁻¹)', fontsize=12)
        ax[i].legend(fontsize=11)
        # Add panel label (a) or (b)
        ax[i].text(0.05, 0.9, panel_label, transform=ax[i].transAxes, fontsize=16, fontweight='bold')

    # Set specific y-limits for each panel
    ax[0].set_ylim(0, 0.25)
    ax[1].set_ylim(0, 0.3)
    ax[0].yaxis.set_major_locator(mticker.MultipleLocator(0.05))
    ax[1].yaxis.set_major_locator(mticker.MultipleLocator(0.05))
    
    plt.tight_layout(pad=2.0)
    
    # Save the figure to a file
    plt.savefig('leij_figure_3_recreation.png', dpi=300)
    print("Plot saved as 'leij_figure_3_recreation.png'")
    
    plt.show()
# --- RUN THE SCRIPT ---
if __name__ == "__main__":
    run_simulation()