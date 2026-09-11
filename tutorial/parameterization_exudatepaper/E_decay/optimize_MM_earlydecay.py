import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib as mpl
from matplotlib.pyplot import figure
from matplotlib.ticker import FormatStrFormatter

fontsize = 20
mpl.rcParams.update({
    'font.size': fontsize,
    'axes.titlesize': fontsize,
    'axes.labelsize': fontsize,
    'xtick.labelsize': fontsize,
    'ytick.labelsize': fontsize,
    'legend.fontsize': fontsize,
    'figure.titlesize': fontsize,
    'mathtext.default': 'regular'
})


# ============================================================
# READ CSV FILE
# ============================================================

df = pd.read_csv("data/early_exudates.csv")

# concentration [mol C / cm^3 water]
Cl = df["Concentration"].values

# decay rate [mol C / cm^3 water / day]
V_obs = df["Decay"].values


# ============================================================
# MICHAELIS-MENTEN MODEL
# ============================================================
#
# Standard Michaelis-Menten:
#
# V = Vmax * Cl / (Km + Cl)
#
# Vmax : maximum decay rate
# Km   : half-saturation concentration
#
# ============================================================

def michaelis_menten(Cl, Vmax, Km):
    return Vmax * Cl / (Km + Cl)

# ============================================================
# INITIAL GUESSES
# ============================================================

Vmax0 = np.max(V_obs)
Km0 = np.median(Cl)

p0 = [Vmax0, Km0]

# ============================================================
# FIT MODEL
# ============================================================

params, covariance = curve_fit(
    michaelis_menten,
    Cl,
    V_obs,
    p0=p0,
    bounds=(0, np.inf)
)

Vmax_fit, Km_fit = params

# ============================================================
# MODEL PREDICTIONS
# ============================================================

Cl_plot = np.linspace(
    np.min(Cl),
    np.max(Cl),
    500
)

V_fit = michaelis_menten(Cl_plot, Vmax_fit, Km_fit)

# ============================================================
# MEAN SQUARED ERROR
# ============================================================

V_pred_data = michaelis_menten(Cl, Vmax_fit, Km_fit)

mse = np.mean((V_obs - V_pred_data)**2)

# ============================================================
# PRINT RESULTS
# ============================================================

print("\n===== FITTED PARAMETERS =====")
print(f"Vmax = {Vmax_fit:.4e} mol C/cm³ water/day")
print(f"Km   = {Km_fit:.4e} mol C/cm³ water")

print(f"\nMSE  = {mse:.4e}")

# ============================================================
# PLOT
# ============================================================

fig, ax = plt.subplots(figsize=(10,8))

# original data
ax.scatter(
    Cl,
    V_obs,
    s=60,
    label='Observed data'
)

# fitted curve
ax.plot(
    Cl_plot,
    V_fit,
    linewidth=2,
    label='Michaelis-Menten fit'
)

# ------------------------------------------------------------
# SCIENTIFIC NOTATION ON X-AXIS
# ------------------------------------------------------------

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1e'))

# ------------------------------------------------------------
# LABELS
# ------------------------------------------------------------

ax.set_xlabel("Concentration Cl (mol C/cm³ water)")
ax.set_ylabel("Decay rate V (mol C/cm³ water/day)")



# ------------------------------------------------------------
# ADD FITTED EQUATION TO PLOT
# ------------------------------------------------------------

equation_text = (
    r"$V = \frac{V_{max} \cdot C_l}{K_m + C_l}$" "\n"
    + rf"$V_{{max}} = {Vmax_fit:.2e}$ mol C cm$^{{-3}}$ d$^{{-1}}$" "\n"
    + rf"$K_m = {Km_fit:.2e}$ mol C cm$^{{-3}}$"
)

ax.text(
    0.05,
    0.95,
    equation_text,
    transform=ax.transAxes,
    fontsize=20,
    verticalalignment='top',
    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
)

ax.legend()

plt.tight_layout()
plt.show()