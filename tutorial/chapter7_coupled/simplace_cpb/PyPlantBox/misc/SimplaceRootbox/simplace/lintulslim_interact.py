'''
Functions to interact with the LintulSlim simulation
running within Simplace

Created on 25.07.2018
Adapted on 09.08.2026 [Murilo Vianna]

@author: Gunther Krauss <guntherkrauss@uni-bonn.de>
'''
import numpy as np
import math

import plantbox as pb

from ..util import rootelongationreduction as rer


def setSimplaceRoots(sim, RLD, unused_lenght, gram_per_cm, layerthickness = 0.01):
    '''Set root length density'''
    unused_bm = max(0., unused_lenght*gram_per_cm)     
    i=0
    layers = len(RLD)
    FRR = np.zeros(layers)
    depth_nr = 1
    for i in range(0, layers):
        FRR[i] = min(1, 1- math.exp(-.3*RLD[i]))
        if(RLD[i]>0):
            depth_nr = i + 1
    par = {'vFRR': FRR.tolist(), 
           'vMD95': depth_nr, 
           'vUnusedRootBiomass': unused_bm, 
           'vRootLengthDensityPerLayer': RLD.tolist()}
    sim.setSimulationValues(par) 


def calculateRLD(rs, soildepth, layers, area = 10000):
    '''Calculate root length density'''
    analysis = pb.SegmentAnalyser(rs)
    vRLD=np.zeros(layers)
    #vRLD[:] = analysis.distribution(rb.ScalarType.length,0,soildepth,layers,True) 
    vRLD[:] = analysis.distribution('length',0.,-soildepth,layers,True) #/mv: plantbox now uses negative values for vertical discretization of soil
    vRLD[:] /= (soildepth/layers)*area
    return vRLD


def getSimplaceValues(sim, gram_per_cm, area = 10000):
    '''Run a simulation step in Simplace and retrieve the values'''
    flt = ['LintulBiomass.rRWRT','CURRENT.DATE','DefaultManagement.DoHarvest',
           'LintulWaterStress.TRANRF','LintulBiomass.sWSO',
           'SlimRoots.RootLengthDensityPerLayer']
    result = sim.stepSimulation(varFilter = flt)
    res = result.toList() 
    rwrt = res['LintulBiomass.rRWRT']
    date = res['CURRENT.DATE']
    tranrf = res['LintulWaterStress.TRANRF']
    yld = res['LintulBiomass.sWSO']
    doharvest = res['DefaultManagement.DoHarvest']
    rld_s = res['SlimRoots.RootLengthDensityPerLayer']
    return (date, (area/10000) * rwrt / gram_per_cm, doharvest, tranrf, yld, rld_s)


def getSimplaceValuesExtended(sim, gram_per_cm, h1, h2, h3, h4, cont_mp = False, area=10000):
    '''
    Run a simulation step in Simplace and retrieve the values
    Calculate root elongation reduction
    '''
    flt = ['LintulBiomass.rRWRT','CURRENT.DATE','DefaultManagement.DoHarvest',
           'LintulWaterStress.TRANRF','LintulBiomass.sWSO',
           'SlimRoots.RootLengthDensityPerLayer',
           'soil.SoilPenetrationResistance_a',
           'soil.SoilPenetrationResistance_b',
           'soil.SoilPenetrationResistance_c',
           'soil_transform.SoilPenetrationResistance',
           'soil_transform.n',
           'soil_transform.alpha',     
           'soil_transform.soilwater_res',     
           'soil_transform.soilwater_sat',
           'soil_transform.bulkdensity',
           'SlimWater.TotalVolumetricWaterContentPerLayer',
           'vCalculateSoilStrength'       
           ]
    result = sim.stepSimulation(varFilter = flt)
    res = result.toList() 
    rwrt = res['LintulBiomass.rRWRT']
    date = res['CURRENT.DATE']
    tranrf = res['LintulWaterStress.TRANRF']
    yld = res['LintulBiomass.sWSO']
    doharvest = res['DefaultManagement.DoHarvest']
    rld_s = res['SlimRoots.RootLengthDensityPerLayer']

    tr = res['soil_transform.soilwater_res']
    ts = res['soil_transform.soilwater_sat']
    ta = res['SlimWater.TotalVolumetricWaterContentPerLayer']
    alpha = res['soil_transform.alpha']
    n = res['soil_transform.n']
    
    calc_soil_strength = (res['vCalculateSoilStrength'] == True)
    if calc_soil_strength:
        a = res['soil.SoilPenetrationResistance_a']
        b = res['soil.SoilPenetrationResistance_b']
        c = res['soil.SoilPenetrationResistance_c']
        bd = res['soil_transform.bulkdensity']  
        q = rer.SoilStrengthFromBDandWC(ta, bd, a, b, c);
    else:
        q = res['soil_transform.SoilPenetrationResistance']
   
    re_reduction = rer.ReductionRE(q, ta, ts, tr, h1, h2, h3, h4, n, alpha, cont_mp = cont_mp)
    return (date, (area/10000) * rwrt / gram_per_cm, doharvest, tranrf, yld, rld_s, re_reduction)
