import numpy as np
import pandas as pd

def asDataFrame(rs, simdate):
    ''' 
    Function to get root outputs from PlantBox 
    @author: Murilo Vianna <mvianna@uni-bonn.de>
    '''
    
    nTypes = 4

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
        len_r[i] = np.sum(rs_len[t == i+1])
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
      n_tot, n_t1, n_t2, n_t3, n_t4 = 0., 0., 0., 0., 0.
      age_t1  = 0.
    else:
      n_tot, n_t1, n_t2, n_t3, n_t4 = len(t), len(t[t == 1]), len(t[t == 2]), len(t[t == 3]), len(t[t == 4])
      age_t1  = np.max(rs_age[t == 1]) # oldest type
    
    # write pb outputs
    res_pb = pd.DataFrame({'CURRENT.DATE':np.array(simdate), 
                           'RootSystem_Age_days':np.array(age_t1),
                           'TotalRootLen_cm':np.array(len_tot), 
                           'SubType1_Len_cm':np.array(len_r[0]), 
                           'SubType2_Len_cm':np.array(len_r[1]), 
                           'SubType3_Len_cm':np.array(len_r[2]), 
                           'SubType4_Len_cm':np.array(len_r[3]),
                           'TotalRootAng_rad':np.array(ang_tot),
                           'SubType1_Ang_rad':np.array(ang_r[0]), 
                           'SubType2_Ang_rad':np.array(ang_r[1]), 
                           'SubType3_Ang_rad':np.array(ang_r[2]), 
                           'SubType4_Ang_rad':np.array(ang_r[3]),                               
                           'TotalRootNumber':np.array(n_tot), 
                           'SubType1_Number':np.array(n_t1), 
                           'SubType2_Number':np.array(n_t2), 
                           'SubType3_Number':np.array(n_t3), 
                           'SubType4_Number':np.array(n_t4),                                                              
                           'RootTipDepthMean_cm':np.array(tip_depth),
                           'RootTipDepthMax_cm':np.array(max_depth)}, index=[0])
    return res_pb    