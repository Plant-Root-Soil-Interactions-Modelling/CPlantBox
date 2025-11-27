import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.integrate import quad
from scipy.special import erfc

# =========================================================================
# ## PART 1: HELPER FUNCTIONS (Defined only once)
# =========================================================================

def kosugi_psd(r, params):
    """Calculates the pore-size distribution f(r) using Kosugi's lognormal model."""
    phi_0 = params['phi_0']
    sigma = params['sigma']
    r_m = params['r_m']
    r = np.maximum(r, 1e-9)
    term1 = phi_0 / (r * sigma * np.sqrt(2 * np.pi))
    exponent = -((np.log(r / r_m))**2) / (2 * sigma**2)
    return term1 * np.exp(exponent)

def analytical_solution(r, initial_params, T, lambda_val):
    """
    A more robust version of the analytical solution with dynamic integration limits
    to prevent numerical artifacts.
    """
    def integrand_func(xi, r_val):
        f0_xi = kosugi_psd(xi, initial_params)
        if lambda_val <= 1e-9 or T <= 1e-9 or xi <= 1e-9: return 0
        
        sqrt_4_pi_lambda_T = np.sqrt(4 * np.pi * lambda_val * T)
        exp1_num = -((r_val - xi + T)**2); exp1_den = 4 * lambda_val * T
        term1 = np.exp(exp1_num / exp1_den)
        exp2_num = -((r_val + xi - T)**2); exp2_den = 4 * lambda_val * T
        term2 = np.exp(-(r_val / lambda_val) + (exp2_num / exp2_den))
        sqrt_4_lambda_T = np.sqrt(4 * lambda_val * T)
        erfc_arg = (r_val + xi - T) / sqrt_4_lambda_T
        term3 = (1 / (2 * lambda_val)) * np.exp(-r_val / lambda_val) * erfc(erfc_arg)
        G = ((term1 + term2) / sqrt_4_pi_lambda_T) + term3
        return f0_xi * G

    final_psd = np.zeros_like(r)
    r_m = initial_params['r_m']; sigma = initial_params['sigma']
    low_bound = r_m * np.exp(-5 * sigma)
    high_bound = r_m * np.exp(5 * sigma)
    
    for i, r_val in enumerate(r):
        result, _ = quad(integrand_func, low_bound, high_bound, args=(r_val,), limit=200, epsabs=1e-9)
        final_psd[i] = result
    return final_psd

# =========================================================================
# ## MAIN SCRIPT TO RECREATE FIGURE 1
# =========================================================================

def run_figure_1_simulation():
    """Runs a forward prediction to recreate Leij et al. (2002) Fig. 1."""
    
    # 1. DERIVE INITIAL STATE from Fig 1(a) table data
    r_mean_t0 = 7.50; variance_t0 = 3.71; m0_t0 = 0.278
    sigma_t0 = np.sqrt(np.log((variance_t0 / r_mean_t0**2) + 1))
    r_m_t0 = r_mean_t0 / np.exp(0.5 * sigma_t0**2)
    initial_params = {'phi_0': m0_t0, 'sigma': sigma_t0, 'r_m': r_m_t0}
    
    print("--- Derived Initial State Parameters ---")
    print(f"  Initial sigma (σ): {sigma_t0:.4f}")
    print(f"  Initial r_m (µm): {r_m_t0:.4f}")
    
    # 2. DEFINE EVOLUTION PARAMETERS from Fig 1 caption
    a = 0.01; b = 5; times = [0, 30, 120]
    lambda_a = 1.0; lambda_b = 0.1

    # 3. CALCULATE PSD at each time step
    r_values = np.linspace(0.01, 20, 200) # Start from 0.01 to avoid issues at r=0
    results_a = {}; results_b = {}
    
    def thornley_r_mean(t):
        return b / (1 + (b/r_mean_t0 - 1) * np.exp(-a*t))

    for t in times:
        r_mean_t = thornley_r_mean(t); T_t = r_mean_t - r_mean_t0
        print(f"\nCalculating for t={t} days... <r(t)> = {r_mean_t:.2f} µm, T(t) = {T_t:.2f} µm")

        if t == 0:
            psd = kosugi_psd(r_values, initial_params)
            results_a[t] = psd; results_b[t] = psd
        else:
            results_a[t] = analytical_solution(r_values, initial_params, abs(T_t), lambda_a)
            results_b[t] = analytical_solution(r_values, initial_params, abs(T_t), lambda_b)

    # 4. PLOT THE RESULTS
    fig, ax = plt.subplots(2, 1, figsize=(8, 11))

    # --- Plot Curves for Both Panels ---
    for t in times:
        ax[0].plot(r_values, results_a[t], 'k-')
    for t in times:
        ax[1].plot(r_values, results_b[t], 'k-')

    # --- NEW: Dynamic Annotation Method ---
    # This method finds the peak of each curve and places the label above it.
    
    # Add dynamic labels to Panel (a)
    for t in times:
        curve_data = results_a[t]
        peak_index = np.argmax(curve_data)
        peak_r = r_values[peak_index]
        peak_f = curve_data[peak_index]
        ax[0].text(peak_r, peak_f + 0.002, f't={t} d', ha='center', va='bottom', fontsize=10)

    # Add dynamic labels to Panel (b)
    for t in times:
        curve_data = results_b[t]
        peak_index = np.argmax(curve_data)
        peak_r = r_values[peak_index]
        peak_f = curve_data[peak_index]
        ax[1].text(peak_r, peak_f + 0.002, f't={t} d', ha='center', va='bottom', fontsize=10)

    # --- General Formatting ---
    for i, panel_label in enumerate(['a)', '(b)']):
        ax[i].set_xlim(0, 20)
        ax[i].set_ylim(0, 0.08)
        ax[i].set_xlabel('r (μm)')
        ax[i].set_ylabel('f (r,t) (μm⁻¹)')
        ax[i].text(0.05, 0.9, panel_label, transform=ax[i].transAxes, fontsize=14, fontweight='bold')
        ax[i].xaxis.set_major_locator(mticker.MultipleLocator(4))
        ax[i].yaxis.set_major_locator(mticker.MultipleLocator(0.02))
    
    # Add Descriptive Caption
    caption = (
        "Figure Difference: The panels show the effect of the dispersivity parameter (λ), which controls\n"
        "how much the pore sizes spread out over time.\n\n"
        "• Panel (a) uses a larger λ = 1.0 µm, causing the distribution to broaden and flatten.\n"
        "• Panel (b) uses a smaller λ = 0.1 µm, causing the distribution to shift with minimal change in shape."
    )
    fig.text(0.5, 0.01, caption, ha='center', va='bottom', fontsize=11, wrap=True)
    fig.subplots_adjust(bottom=0.15)
    
    plt.savefig('leij_figure_1_recreation.png', dpi=300)
    print("\nPlot saved as 'leij_figure_1_recreation.png'")
    plt.show()

if __name__ == "__main__":
    run_figure_1_simulation()