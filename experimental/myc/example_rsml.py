"""root system length over time"""
# %%
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
# from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
import numpy as np
import pandas as pd
import visualisation.vtk_plot as vp
name = "RSHA_WT_30d_export"
dfs = pd.read_csv(name+".csv", low_memory=False)  
RSA_ids = list(set(dfs['RSA_id'].values))
print('RSA_ids',RSA_ids)
for RSA_id in range(1,len(RSA_ids)+1):
    df = dfs[(dfs['RSA_id'].values == 'Tomato_60_'+str(RSA_id))]
    
    assert df['node1ID'].values[0] == 0
    nodesID = np.concatenate(([df['node1ID'].values[0]] , 
                              df['node2ID'].values), dtype =int)
    nodesX = np.concatenate(([df['x1'].values[0]] , df['x2'].values))
    nodesY = np.concatenate(([df['y1'].values[0]] , df['y2'].values))
    nodesZ = np.concatenate(([df['z1'].values[0]] , df['z2'].values))
    nodes =np.array( [nodesX,nodesY,nodesZ]).T.tolist()
    nodes = [pb.Vector3d(node) for node in nodes]
    segments = np.array([df['node1ID'].values, df['node2ID'].values]).T.tolist()
    segments_ = np.array([df['node1ID'].values, df['node2ID'].values]).T


    assert len( df['node2ID'].values.tolist()) == len(set(df['node2ID'].values.tolist())), "Duplicate values found in the list!"

    segments = [pb.Vector2i(seg) for seg in segments]
    segCTs = df['time'].values.tolist()
    radii = df['radius'].values.tolist()
    
    ### troubleshooting, todo: delete
    length_th = df['length'].values.tolist()
    length_ob = [np.sqrt(sum((np.array(nodes[seg[0]]) - np.array(nodes[seg[1]]))**2)) for seg in segments_]
    lenDiffs = np.array(abs(np.array(length_ob) - np.array(length_th)))
    maxDiff= max(lenDiffs)
    maxDiffIndex = np.where(lenDiffs == maxDiff)[0][0]
    print( maxDiff,maxDiffIndex)
    print('len makes sense?',maxDiff,length_th[maxDiffIndex], 
          length_ob[maxDiffIndex])
    print(segments_[maxDiffIndex], nodes[segments_[maxDiffIndex][0]], 
          nodes[segments_[maxDiffIndex][1]])
    ###
    
    ana = pb.SegmentAnalyser(nodes,segments,segCTs, radii)
    ana.addData("kr",df['kr'].values.tolist())
    ana.addData("Kx",df['Kx'].values.tolist())
    ana.addData("SUF",df['SUF'].values.tolist())
    ana.addData("subType",df['type'].values.tolist())
    ana.addData("branchID",df['branchID'].values.tolist())
    ana.addData("nodesID",nodesID)
   
    ana.write(name+str(RSA_id)+".vtp", ['kr', 'Kx', 'SUF', 'subType', 'branchID','nodesID'])
    vp.plot_roots(ana, p_name= "subType")
    raise Exception
# %%
