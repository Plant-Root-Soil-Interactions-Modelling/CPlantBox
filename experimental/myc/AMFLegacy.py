import sys; sys.path.append("../.."); sys.path.append("../../src/")
import numpy as np
import plantbox as pb

import AMFAnalysis as amf
import AMFDensities as amfd

# Old functions that are currently not called by any script, kept for reference.
# Replaced by AMFDensities.getParasPerRing (one SegmentAnalyser for all rings and parameters, much faster).

def getLengthPerSubtype(plant):
    hyphae = pb.SegmentAnalyser(plant)
    hyphae.filter("organType",5.0)
    hyphae.pack()
    lenSubtype = []
    for i in range(1,4):
        hyphae.filter("subType",i)
        hyphae.pack()
        lenSubtype.append(hyphae.getSummed("length"))
        hyphae = pb.SegmentAnalyser(plant)
        hyphae.filter("organType",5.0)
        hyphae.pack()
    return lenSubtype

def getParaSumperRing(parameter, times, plant, rings):
    return amfd.getParasPerRing([parameter], times, plant, rings)[parameter]

def getParaSumperRing_old(parameter, times, plant, rings):
    # print("Parameter sum for "+ parameter)
    paradenmat = np.zeros((len(rings),len(times[1:])))
    flipped = np.flip(np.asarray(times))
    for k, ring in enumerate(rings):
        ana = amf.getMycSegmentAnalyser(plant)
        for j in range(len(times[1:])-1):
            ana.filter("creationTime", 0, flipped[j])
            ana.pack()
            distrib = ana.getSummed(parameter, ring)
            ana.filter("creationTime",0,flipped[j+1])
            ana.pack()
            summed = 0
            if parameter == "nodeTips": # we want to look only at actively growing tips. (Oyarte Galvez 2025 page 3 when discussing Fig 2b)
                summed = ana.getSummed(parameter, ring)
            paradenmat[k, len(times[1:])-1 -j] = np.array(distrib-summed).sum()
        ana.filter("creationTime",0,flipped[len(times[1:])-1])
        ana.pack()
        summed = ana.getSummed(parameter, ring)
        ana.filter("creationTime",0,flipped[len(times[1:])])
        ana.pack()
        summed = summed - ana.getSummed(parameter, ring)
        paradenmat[k, 0] = np.array(summed).sum()
    return paradenmat

def getParameterOverTime(parameter, times, plant, subType):
    paradenmat = np.zeros((len(subType), len(times[1:])))
    flipped = np.flip(np.asarray(times))
    for k in subType:
        ana = pb.SegmentAnalyser(plant)
        ana.filter("organType", 5, 5)
        for j in range(len(times[1:])-1):
            ana.filter("creationTime", 0, flipped[j])
            ana.pack()
            distrib = ana.getSummed(parameter)
            ana.filter("creationTime",0,flipped[j+1])
            ana.pack()
            summed = ana.getSummed(parameter)
            paradenmat[k-1,len(times[1:])-1 -j] = np.array(distrib-summed).sum() #/np.array(allTips-overTips).sum() if np.array(allTips-overTips).sum() != 0 else 0
        ana.filter("creationTime",0,flipped[len(times[1:])-1])
        ana.pack()
        summed = ana.getSummed(parameter)
        ana.filter("creationTime",0,flipped[len(times[1:])])
        ana.pack()
        summed = summed - ana.getSummed(parameter)
        paradenmat[k-1,-1] = np.array(summed).sum()
    return paradenmat
