import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.integrate import quad
from scipy.optimize import curve_fit
from scipy.special import erfc

# =========================================================================
# ## PART 1: A ROBUST PARAMETER CLASS AND HELPER FUNCTIONS
# =========================================================================

class SoilParams:
    """A standardized container for soil hydraulic parameters."""
    def __init__(self, phi_0=None, sigma=None, r_m=None, theta_s=None, theta_r=None):
        self.sigma = sigma
        self.r_m = r_m
        if phi_0 is not None:
            self.phi_0 = phi_0
        else:
            self.phi_0 = theta_s - theta_r

def kosugi_psd(r, params: SoilParams):
    """Calculates the pore-size distribution f(r)."""
    r = np.maximum(r, 1e-9)
    term1 = params.phi_0 / (r * params.sigma * np.sqrt(2 * np.pi))
    exponent = -((np.log(r / params.r_m))**2) / (2 * params.sigma**2)
    return term1 * np.exp(exponent)
def analytical_solution(r, initial_params: SoilParams, T, lambda_val):
    """
    The definitive, robust analytical solution with a constrained dynamic
    integration range that works for both narrow and broad distributions.
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
    
    # --- BEST OF BOTH WORLDS: Constrained Dynamic Integration ---
    r_m = initial_params.r_m
    sigma = initial_params.sigma
    # Calculate a dynamic range based on the distribution's properties...
    dynamic_low = r_m * np.exp(-5 * sigma)
    dynamic_high = r_m * np.exp(5 * sigma)
    # ...but constrain it within reasonable, fixed limits.
    low_bound = max(0.01, dynamic_low)
    high_bound = min(5000.0, dynamic_high)
    
    for i, r_val in enumerate(r):
        result, _ = quad(integrand_func, low_bound, high_bound, args=(r_val,), limit=100)
        final_psd[i] = result
    return final_psd
    
    
def vg_to_kosugi(alpha, n, theta_s, theta_r):
    """
    Convert van Genuchten (θs, θr, α, n) -> Kosugi via Chandrasekhar et al. (2019), Eqs. (6)–(8).
    α in 1/cm, n > 1. Returns Kosugi σ on ln h, inflection head h0, median head hm, and rm (µm).
    A_cm2: Jurin proportionality constant (≈0.149 cm² @ ~20°C; ≈0.151–0.152 cm² @ ~10°C).
    """
    m = 1 - (1 / n)
    A = 0.149 # Proportionality constant from Young-Laplace [cm^2]
    sigma_sq = (1 - m) * np.log((2**(1/m) - 1) / m)
    sigma = np.sqrt(sigma_sq)
    h0 = (m**(1 - m)) / alpha
    hm = h0 * np.exp(sigma_sq)
    r_m_cm = A / hm
    r_m_um = r_m_cm * 10000
    return SoilParams(theta_s=theta_s, theta_r=theta_r, sigma=sigma, r_m=r_m_um)
# =========================================================================
# ## PART 2: FIGURE RECREATION FUNCTIONS
# =========================================================================

def recreate_figure_1():
    """Runs a forward prediction to recreate Leij et al. (2002) Fig. 1."""
    print("--- Recreating Figure 1 (Forward Prediction) ---")
    
    r_mean_t0 = 7.50; variance_t0 = 3.71; m0_t0 = 0.278
    sigma_t0 = np.sqrt(np.log((variance_t0 / r_mean_t0**2) + 1))
    r_m_t0 = r_mean_t0 / np.exp(0.5 * sigma_t0**2)
    initial_params = SoilParams(phi_0=m0_t0, sigma=sigma_t0, r_m=r_m_t0)
    
    a = 0.01; b = 5; times = [0, 30, 120]
    lambda_a = 1.0; lambda_b = 0.1
    r_values = np.linspace(0.01, 20, 200)
    results_a = {}; results_b = {}
    
    thornley = lambda t: b / (1 + (b/r_mean_t0 - 1) * np.exp(-a*t))
    for t in times:
        T_t = thornley(t) - r_mean_t0
        if t == 0:
            psd = kosugi_psd(r_values, initial_params)
            results_a[t] = psd; results_b[t] = psd
        else:
            results_a[t] = analytical_solution(r_values, initial_params, abs(T_t), lambda_a)
            results_b[t] = analytical_solution(r_values, initial_params, abs(T_t), lambda_b)

    fig, ax = plt.subplots(2, 1, figsize=(8, 11))
    for t in times:
        ax[0].plot(r_values, results_a[t], 'k-')
        ax[1].plot(r_values, results_b[t], 'k-')
    
    for i, panel_label in enumerate(['a)', '(b)']):
        ax[i].set_xlim(0, 20); ax[i].set_ylim(0, 0.08)
        ax[i].set_xlabel('r (μm)'); ax[i].set_ylabel('f (r,t) (μm⁻¹)')
        ax[i].text(0.05, 0.9, panel_label, transform=ax[i].transAxes, fontsize=14, fontweight='bold')
        ax[i].xaxis.set_major_locator(mticker.MultipleLocator(4))
        ax[i].yaxis.set_major_locator(mticker.MultipleLocator(0.02))

    for t in [0, 30, 120]:
        curve_a = results_a[t]; peak_idx_a = np.argmax(curve_a)
        ax[0].text(r_values[peak_idx_a], curve_a[peak_idx_a] + 0.002, f't={t} d', ha='center', va='bottom', fontsize=10)
        curve_b = results_b[t]; peak_idx_b = np.argmax(curve_b)
        ax[1].text(r_values[peak_idx_b], curve_b[peak_idx_b] + 0.002, f't={t} d', ha='center', va='bottom', fontsize=10)
        
    caption = ( "Figure Difference: The panels show the effect of the dispersivity parameter (λ).\n"
                "Panel (a) uses λ = 1.0 µm, causing the distribution to broaden and flatten.\n"
                "Panel (b) uses λ = 0.1 µm, causing a shift with minimal change in shape." )
    fig.text(0.5, 0.01, caption, ha='center', va='bottom', fontsize=11, wrap=True)
    fig.subplots_adjust(bottom=0.15)
    
    plt.savefig('leij_figure_1_recreation.png', dpi=300); print("Plot saved as 'leij_figure_1_recreation.png'")
    plt.show()

def recreate_figure_3():
    """Runs an inverse fitting to recreate Leij et al. (2002) Fig. 3."""
    print("--- Recreating Figure 3 (Inverse Fitting) ---")
    
    initial_a = SoilParams(theta_s=0.469, theta_r=0.104, sigma=2.27, r_m=6.36)
    final_a = SoilParams(theta_s=0.392, theta_r=0.097, sigma=2.33, r_m=2.15)
    T_a = 20.5
    r_values = np.logspace(-1, 2, 50)
    observed_psd_a = kosugi_psd(r_values, final_a)
    popt_a, _ = curve_fit(lambda r, l: analytical_solution(r, initial_a, T_a, l), r_values, observed_psd_a, p0=[1.3], bounds=(0, np.inf))
    predicted_psd_a = analytical_solution(r_values, initial_a, T_a, popt_a[0])

    initial_b = SoilParams(theta_s=0.394, theta_r=0.142, sigma=2.60, r_m=3.97)
    final_b = SoilParams(theta_s=0.367, theta_r=0.128, sigma=2.44, r_m=0.997)
    T_b = 24.9
    observed_psd_b = kosugi_psd(r_values, final_b)
    popt_b, _ = curve_fit(lambda r, l: analytical_solution(r, initial_b, T_b, l), r_values, observed_psd_b, p0=[0.96], bounds=(0, np.inf))
    predicted_psd_b = analytical_solution(r_values, initial_b, T_b, popt_b[0])

    fig, ax = plt.subplots(2, 1, figsize=(8, 12))
    ax[0].plot(r_values, observed_psd_a, 'o', c='k', mfc='k', ms=6, label='PSD 1340 observed')
    ax[0].plot(r_values, kosugi_psd(r_values, initial_a), 'k--', label='1350')
    ax[0].plot(r_values, predicted_psd_a, 'k-', lw=1.5, label='1340 predicted')
    ax[1].plot(r_values, observed_psd_b, 'o', c='k', mfc='k', ms=6, label='PSD 2000 observed')
    ax[1].plot(r_values, kosugi_psd(r_values, initial_b), 'k--', label='2010')
    ax[1].plot(r_values, predicted_psd_b, 'k-', lw=1.5, label='2000 predicted')
    
    for i, panel_label in enumerate(['(a)', '(b)']):
        ax[i].set_xscale('log')
        ax[i].set_xlim(0.1, 100)
        ax[i].xaxis.set_major_locator(mticker.LogLocator(base=10.0, numticks=4))
        ax[i].xaxis.set_major_formatter(mticker.ScalarFormatter())
        ax[i].set_xlabel('r (μm)'); ax[i].set_ylabel('f (μm⁻¹)')
        ax[i].legend(); ax[i].text(0.05, 0.9, panel_label, transform=ax[i].transAxes, fontsize=16, fontweight='bold')

    ax[0].set_ylim(0, 0.25); ax[1].set_ylim(0, 0.3)
    ax[0].yaxis.set_major_locator(mticker.MultipleLocator(0.05))
    ax[1].yaxis.set_major_locator(mticker.MultipleLocator(0.05))
    
    plt.tight_layout(pad=2.0)
    plt.savefig('leij_figure_3_recreation.png', dpi=300); print("Plot saved as 'leij_figure_3_recreation.png'")
    plt.show()

def recreate_figure_from_vg_data():
    """
    Performs inverse fitting using VG parameter sets as input and plots the result.
    """
    print("--- Running Inverse Fitting from van Genuchten Data ---")
    
    # --- 1. USER INPUT: Provide your initial and final VG sets here ---
    # Using data for the "loamy sand" from Leij's SSSA J. paper as an example
    initial_vg_params = {'alpha': 0.075, 'n': 1.89, 'theta_s': 0.41, 'theta_r': 0.09}
    final_vg_params   = {'alpha': 0.124, 'n': 2.26, 'theta_s': 0.39, 'theta_r': 0.09}

    # --- 2. CONVERT VG sets to our internal SoilParams format ---
    initial_params = vg_to_kosugi(**initial_vg_params)
    final_params = vg_to_kosugi(**final_vg_params)

    # --- 3. RUN THE FITTING ---
    r_mean_initial = initial_params.r_m * np.exp(0.5 * initial_params.sigma**2)
    r_mean_final = final_params.r_m * np.exp(0.5 * final_params.sigma**2)
    T = abs(r_mean_final - r_mean_initial)

    r_values = np.logspace(-1, 3, 50)
    observed_psd = kosugi_psd(r_values, final_params)
    
    popt, _ = curve_fit(lambda r, l: analytical_solution(r, initial_params, T, l), 
                        r_values, observed_psd, p0=[1.0], bounds=(0, np.inf))
    fitted_lambda = popt[0]
    predicted_psd = analytical_solution(r_values, initial_params, T, fitted_lambda)
    
    print(f"✅ Fitting complete. Fitted Dispersivity (λ): {fitted_lambda:.4f} µm")

    # --- 4. PLOT THE RESULTS ---
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(r_values, observed_psd, 'o', c='k', mfc='k', ms=6, label='Final PSD (from VG data)')
    ax.plot(r_values, kosugi_psd(r_values, initial_params), 'k--', label='Initial PSD (from VG data)')
    ax.plot(r_values, predicted_psd, 'r-', lw=2, label=f'Predicted Final PSD (λ={fitted_lambda:.2f})')

    ax.set_xscale('log')
    ax.set_xlim(0.1, 1000)
    ax.xaxis.set_major_locator(mticker.LogLocator(base=10.0, numticks=5))
    ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
    ax.set_xlabel('r (μm)'); ax.set_ylabel('f (μm⁻¹)')
    ax.legend(); ax.set_title('Inverse Fitting Result from VG Parameter Sets')
    ax.grid(True, which="both", ls="--")
    
    plt.tight_layout()
    plt.savefig('vg_data_fit_result.png', dpi=300)
    print("Plot saved as 'vg_data_fit_result.png'")
    plt.show()

# =========================================================================
# ## PART 3: MAIN EXECUTION BLOCK
# =========================================================================

if __name__ == "__main__":
    recreate_figure_1()
    recreate_figure_3()
    recreate_figure_from_vg_data()