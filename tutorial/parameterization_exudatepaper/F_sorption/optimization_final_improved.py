import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
import pandas as pd
import sys
import matplotlib as mpl

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

def get_axis_limits(ax, scalex, scaley):
    ymin = ax.get_ylim()[0]
    ymax = ax.get_ylim()[1]
    return ax.get_xlim()[1]*scalex, 10**(np.log10(ymin) + scaley * (np.log10(ymax) - np.log10(ymin)))
    

# =========================================================
# CONSTANTS
# =========================================================
V = np.array([2.25, 2.25, 1])   # solution in mL (cm3)
m = np.array([1.0, 1, 0.2])    # soil in g
mV_ = m / V  # factor to convert Cl in Cs
molarMassC = 12
k_claysilt = 0.67
Qmax = k_claysilt * 0.079
Cs_max = Qmax / molarMassC #mol C/g soil 
min_to_day = 1/(60*24)
#Cl (mol C/cm³ water)
#Cs (mol C/g soil)
 
# =========================================================
# LANGMUIR MODEL
# =========================================================
def model(t, y, ka, kd):
    Cs, Cl = y
 
    dCs = ka * Cl * (Cs_max - Cs) - kd * Cs
    dCl = -mV * dCs
 
    return [dCs, dCl]
 
# =========================================================
# SOLVER (one time series)
# =========================================================
def solve(t_eval, Cl0, params):
 
    t_eval = np.asarray(t_eval, float)
 
    sol = solve_ivp(
        model,
        (t_eval.min(), t_eval.max()),
        [0.0, Cl0],     # initial conditions for Cs and Cl
        method="Radau",  # for stiff odes
        args=tuple(params),    # passes values for K and kd
        dense_output=True,
        rtol=1e-6,
        atol=1e-10
    )
 
    y = sol.sol(t_eval)
 
    return y[0], y[1]
 
# =========================================================
# RESIDUALS (3 time series)
# =========================================================
def residuals(params, datasets):

    Cl_data = []
    Cl_pred = []
    for d in datasets: 
        t = d["t"]
        Cl0 = d["Cl0"]
        Cs_pred_, Cl_pred_ = solve(t, Cl0, params)
        Cl_pred.extend(Cl_pred_)
        Cl_data.extend(d["Cl"])
    
    Cl_pred = np.array(Cl_pred) 
    Cl_data = np.array(Cl_data)
    eps = 1e-12

    return np.log(Cl_pred + eps) - np.log(Cl_data + eps)
 
# =========================================================
# FIT
# =========================================================
def fit_model(datasets):
 
    best = None
    best_cost = np.inf

    for ka0 in [1e-4, 1e-3, 1e-2, 1e-1]:
        for kd0 in [1e-6, 1e-5, 1e-4, 1e-3]:
            
            result = least_squares(
                residuals,
                [ka0, kd0],
                args=(datasets,),
                bounds=([1e-12, 1e-12], [1e3, 1e3]),
                method="trf"
            )

            if result.cost < best_cost:
                best = result
                best_cost = result.cost

    return best
 
# =========================================================
# PLOT (ALL 3 CURVES)
# =========================================================
def plot_all(datasets, params, i):
 
    ax = plt.subplot(1,3,i+1)
    color = 'r'
    for d in datasets: 
        
        t = d["t"]
        Cl = d["Cl"]
        Cl0 = d["Cl0"]
        
        ts = np.linspace(t.min(), t.max(), 300)
     
        _, Cl_fit = solve(ts, Cl0, params)
     
        plt.scatter(t, Cl, color=color)
        plt.plot(ts, Cl_fit, linestyle = '--', color=color)
        color = 'k'
    
    plt.yscale("log")
    plt.xlabel("Time (min)")
    if i == 0: 
        plt.ylabel("$C_l$ (mol $cm^{-3} water$)")
        
    res = residuals(result.x,datasets) 
    rmse_log = np.sqrt(np.mean(res**2))
    print(f"RMSE(log) = {rmse_log:.2e}")
    ax.annotate(annots[i], xy=get_axis_limits(ax, 0.8, 0.85))
    # ax.annotate(f"RMSE(log) = {rmse_log:.2e}", xy=get_axis_limits(ax, 0.3, .7))
 
# =========================================================
# Main
# =========================================================
if __name__ == "__main__":
 
    plt.figure(figsize=(22, 8))
    data = ['alanine Fischer', 'glucose Fischer', 'alanine Dippold']
    annots = ['(a)', '(b)', '(c)']
    
    
    # data
    kin = pd.read_csv("data/removal_data.csv")  
    
    for i in range(0, len(data)): 
        datasets = []
        mV = mV_[i]
        
        t = kin["Time "+data[i]].values             # minutes
        Cl = kin["Cl "+data[i]].values              # mol/cm³ water
        Cs = kin["Cs "+data[i]].values              # mol/g soil
        m = ~np.isnan(t)
        t, Cl, Cs = t[m], Cl[m], Cs[m]
        sort_idx = np.argsort(t)
        t = t[sort_idx]
        Cl = Cl[sort_idx]
        Cs = Cs[sort_idx]
        Cl0 = Cl[0]
        
        datasets.append({
        "t": t,
        "Cl": Cl,
        "Cs": Cs,
        "Cl0": Cl0,
        })

        eq = pd.read_csv("data/eq_"+data[i]+".csv")  
        
        for j in range(0, int((np.shape(eq)[1]-1)/2)): 
            
            t = eq["Time"].values
            Cl = eq["Cl "+data[i]+" "+str(int(j+1))].values              # mol C/cm³ water
            Cs = eq["Cs "+data[i]+" "+str(int(j+1))].values              # mol C/g soil
            sort_idx = np.argsort(t)
            t = t[sort_idx]
            Cl = Cl[sort_idx]
            Cs = Cs[sort_idx]
            Cl0 = Cl[0]

            datasets.append({
            "t": t,
            "Cl": Cl,
            "Cs": Cs,
            "Cl0": Cl0,
            })
     
     
        # =====================================================
        # FIT
        # =====================================================
        result = fit_model(datasets)
     
        print("\nFitted parameters:")
        print(f"kads = {result.x[0]/min_to_day:.2e} ($cm^{3}$ water $mol^{-1}$ C $d^{-1}$)")
        print(f"kdes = {result.x[1]/min_to_day:.2e} ($d^{-1}$)")
     
        # =====================================================
        # PLOT
        # =====================================================
        plot_all(datasets,
                 result.x, i)

handles = [
    plt.Line2D([], [], color='r', linestyle='None', marker='o',
               label='kinetic sorption data'),
    plt.Line2D([], [], color='k', linestyle='None', marker='o',
               label='equilibrium sorption data'),
    plt.Line2D([], [], color='k', linestyle='--',
               label='Langmuir model fit')
]

plt.gcf().legend(
    handles=handles,
    loc='lower center',
    ncol=3,                    # horizontal layout
    bbox_to_anchor=(0.5, -0.02),  # BELOW plots
    frameon=False
)

plt.subplots_adjust(
    left=0.08,
    right=0.98,
    top=0.9,
    bottom=0.2,
    wspace=0.3,
    hspace=0.5
)
plt.tight_layout(rect=[0, 0.1, 1, 1])
plt.show()