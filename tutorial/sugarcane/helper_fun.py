import numpy as np
import pandas as pd

# custom function for stress reduction factor
def Stress_Reduction_Function(A: float, B: float, C: float, x0: float, y0: float, S: float, Qp: float) -> float:
    '''
    Calculate stress reduction factor based on soil water potential (S, in hPa) and soil penetration resistance (Qp, in MPa)

    @param A: Maximum reduction factor (dimensionless)
    @param B: Parameter related to the sensitivity to soil water potential (hPa)
    @param C: Parameter related to the sensitivity to soil penetration resistance (MPa)
    @param x0: Soil water potential at which the reduction factor is maximum (hPa)
    @param y0: Soil penetration resistance at which the reduction factor is maximum (MPa)
    @param S: Current soil water potential (hPa)
    @param Qp: Current soil penetration resistance (MPa)
    '''
    SRF = A*np.e**(-0.5*(((S-x0)/B)**2+((Qp-y0)/C)**2))  # Stress reduction function
    return SRF


def asDataFrame(rs, simdate):
    ''' 
    Function to get root outputs from PlantBox 
    @author: Murilo Vianna <mvianna@uni-bonn.de>
    '''
    
    nTypes = 5

    # get pb-related outputs
    t        = np.array(rs.getParameter("type"))
    rs_len   = np.array(rs.getParameter("length"))
    rs_ang   = np.array(rs.getParameter("theta"))
    rs_age   = np.array(rs.getParameter("age"))
    
    len_tot, ang_tot = 0., 0.
    len_r, ang_r = np.zeros(nTypes), np.zeros(nTypes)
    
    if len(rs_len) > 0: 
      len_tot = np.sum(rs_len)
      ang_tot = np.mean(rs_ang)
    
    for i in range(0,nTypes):
      if len(rs_len[t == i+1]) > 0:
        len_r[i] = np.mean(rs_len[t == i+1])
        ang_r[i] = np.mean(rs_ang[t == i+1])
    
    roots = rs.getPolylines()
    tip_depth = 0.
    max_depth = 0.
    for r in roots:
      tip_depth += r[-1].z
      max_depth = min(max_depth, r[-1].z)
    if len(roots) == 0:
      tip_depth = 0.
    else:
      tip_depth = tip_depth / len(roots)
    
    if len(t) == 0:
      n_tot, n_t1, n_t2, n_t3, n_t4, n_t5 = 0., 0., 0., 0., 0., 0.
      age_t1  = 0.
    else:
      n_tot, n_t1, n_t2, n_t3, n_t4, n_t5 = len(t), len(t[t == 1]), len(t[t == 2]), len(t[t == 3]), len(t[t == 4]), len(t[t == 5])
      age_t1  = np.max(rs_age[t == 1]) # oldest type
    
    # write pb outputs
    res_pb = pd.DataFrame({'CURRENT.DATE':np.array(simdate), 
                           'RootSystem_Age_days':np.array(age_t1),
                           'TotalRootLen_cm':np.array(len_tot), 
                           'SubType1_Len_cm':np.array(len_r[0]), 
                           'SubType2_Len_cm':np.array(len_r[1]), 
                           'SubType3_Len_cm':np.array(len_r[2]), 
                           'SubType4_Len_cm':np.array(len_r[3]),
                           'SubType5_Len_cm':np.array(len_r[4]),
                           'TotalRootAng_rad':np.array(ang_tot),
                           'SubType1_Ang_rad':np.array(ang_r[0]), 
                           'SubType2_Ang_rad':np.array(ang_r[1]), 
                           'SubType3_Ang_rad':np.array(ang_r[2]), 
                           'SubType4_Ang_rad':np.array(ang_r[3]),                               
                           'SubType5_Ang_rad':np.array(ang_r[4]),                               
                           'TotalRootNumber':np.array(n_tot), 
                           'SubType1_Number':np.array(n_t1), 
                           'SubType2_Number':np.array(n_t2), 
                           'SubType3_Number':np.array(n_t3), 
                           'SubType4_Number':np.array(n_t4),                                                              
                           'SubType5_Number':np.array(n_t5),                                                              
                           'RootTipDepthMean_cm':np.array(tip_depth),
                           'RootTipDepthMax_cm':np.array(max_depth)}, index=[0])
    return res_pb

def is_running_in_notebook():
    try:
        shell = get_ipython().__class__.__name__
        return shell in ('ZMQInteractiveShell',)  # Jupyter Notebook or JupyterLab
    except NameError:
        return False  # Probably standard Python script