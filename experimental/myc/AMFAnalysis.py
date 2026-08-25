import sys; sys.path.append("../.."); sys.path.append("../../src/")
import numpy as np
import plantbox as pb

def getMycSegmentAnalyser(plant):
        ana = pb.SegmentAnalyser(plant)
        ana.addData("colonization", plant.getNodeColonizations(2))
        ana.addData("colonizationTime", plant.getNodeColonizationTime(2))
        ana.addData("anastomosis", plant.getAnastomosisPoints(5))
        ana.addData("nodeTips", plant.getNodeTips(5))
        return ana


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

def getParaDistperRing(parameter, times, plant, rings):
        paradenmat = np.zeros((len(rings),len(times[1:])))
        flipped = np.flip(np.asarray(times))
        for k, ring in enumerate(rings):
            ringana = pb.SegmentAnalyser(plant) # need to copy the whole plant for segment analyzer since cropping to one ring removes all information outside
            ringana.crop(ring)
            for j in range(len(times[1:])-1):
                ringana.filter("creationTime", 0, flipped[j])
                ringana.pack()
                distrib = ringana.getSummed(parameter)
                ringana.filter("creationTime",0,flipped[j+1])
                ringana.pack()
                summed = ringana.getSummed(parameter)
                paradenmat[k, len(times[1:])-1 -j] = np.array(distrib-summed).sum() 
            ringana.filter("creationTime",0,flipped[len(times[1:])-1])
            ringana.pack()
            summed = ringana.getSummed(parameter)
            ringana.filter("creationTime",0,flipped[len(times[1:])])
            ringana.pack()
            summed = summed - ringana.getSummed(parameter)
            paradenmat[k, -1] = np.array(summed).sum() #/(np.pi*(4.7**2)/2)
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

def makedishes(diameter, height, barrier_thickness, barrier_height, opening_length, opening_height, rootradius,nRings):
    radius = diameter / 2
    # petri dish has a radius of 9.4 cm and a height of 1.6 cm
    petri_dish = pb.SDF_PlantContainer(radius,radius,height,False)
    # the helper dish is used to cut the petri dish in half, it has the same radius and height as the petri dish but is rotated and translated to cut the petri dish in half
    helper_dish = pb.SDF_PlantContainer(radius,radius,height,True)
    # moving the helper dish such that it halves the petri dish and removes a bit more to restrict roots to the correct side of barrier
    moved_helper_dish = pb.SDF_RotateTranslate(helper_dish, 0, 0, pb.Vector3d(-(radius+barrier_thickness+rootradius), 0, 0))
    half_dish = pb.SDF_Intersection(petri_dish, moved_helper_dish)
    # helper container for barrier
    helper_staff = pb.SDF_PlantBox(barrier_thickness,barrier_height,diameter)
    # helper container for opening in barrier
    helper_staff2 = pb.SDF_PlantBox(barrier_thickness,opening_height,opening_length)
    # have to  move the helper staff for the right position of the opening and barrier, the midpoint is the distance from the center of the helper staff to the center of the petri dish, the bottompoint is the distance from the center of the helper staff to the center of the petri dish in the y direction, and the helper staff is moved to the position of the opening and barrier
    midpoint = opening_length / 2 - diameter / 2
    bottompoint = opening_height / 2 - barrier_height / 2
    # container for opening moved to the right position
    helper_staff2 = pb.SDF_RotateTranslate(helper_staff2, 0,0,pb.Vector3d(0, bottompoint, midpoint))
    # opening made in the barrier
    helper_dish2 = pb.SDF_Difference(helper_staff, helper_staff2)
    # barrier moved to the right position
    moved_helper_dish2 = pb.SDF_RotateTranslate(helper_dish2, 90, pb.SDF_Axis.xaxis , pb.Vector3d(0, -radius, -height/2))
    petri_dish = pb.SDF_Difference(petri_dish, moved_helper_dish2)

    moved_helper_dish_hyphae = pb.SDF_RotateTranslate(helper_dish, 0, 0, pb.Vector3d(-(radius+barrier_thickness/2), 0, 0))
    small_dish = pb.SDF_PlantContainer(radius,radius,height,False)
    small_hyphae_dish = pb.SDF_Difference(small_dish, moved_helper_dish_hyphae)

    # nRings = 25
    centrepoint = [0, 1.5, 0] ## Set a different centre for the rings for analysis, so that growth is more centred
    small_dish = pb.SDF_PlantContainer(radius*np.sqrt(1/nRings),radius*np.sqrt(1/nRings),height,False)
    ringone = pb.SDF_Difference(small_dish, moved_helper_dish_hyphae)
    moved_ringone = pb.SDF_RotateTranslate(ringone, 0, 0, pb.Vector3d(centrepoint[0], centrepoint[1], centrepoint[2]))
    rings = []
    rings.append(moved_ringone)
    for i in range(2, nRings+1):
        small_dish = pb.SDF_PlantContainer(radius*np.sqrt(i/nRings),radius*np.sqrt(i/nRings),height,False)
        small_hyphae_dish2 = pb.SDF_Difference(small_dish, moved_helper_dish_hyphae)
        old_dish = pb.SDF_Difference(pb.SDF_PlantContainer(radius*np.sqrt((i-1)/nRings),radius*np.sqrt((i-1) /nRings),height,False),moved_helper_dish_hyphae)
        small_hyphae_dish2 = pb.SDF_Difference(small_hyphae_dish2,old_dish)
        moved_small_hyphae_dish = pb.SDF_RotateTranslate(small_hyphae_dish2, 0, 0, pb.Vector3d(centrepoint[0], centrepoint[1], 0))
        rings.append(moved_small_hyphae_dish)
    
    return petri_dish, small_hyphae_dish, half_dish, rings