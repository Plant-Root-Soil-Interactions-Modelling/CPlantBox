"""" rhizotubes as obstacles for multirootsystems"""
import sys;

sys.path.append("../../../")
sys.path.append("../../../src")
sys.path.append("../../../modelparameter/structural/rootsystem/")
import plantbox as pb
import math
import numpy as np       
import matplotlib.pyplot as plt


def run_benchmark(num, name,M, N, distr, distp, distl, simtime, w_field, l_field, hlayer, depth, planes, tube_diam, fieldbox_wo_rhizotubes, rhizotubes, rhizotubeso, sc1, sc2, sc3,sc4,sc_volume, imgw, imgl, size, x_standard, x_lowr, x_highr, x_layer, img_cont, y_, z_, tropismN, tropismS): 

    # Initialize and simulate N*M root systems
    allRS = []
    allRS_cores = []
    axroots_ = [1,4,5]
    for i in range(0, M):
        for j in range(0, N):
            rs = pb.Plant()
            rs_cores = pb.Plant()
            rs.readParameters('../../../modelparameter/structural/rootsystem/'+name + ".xml")
            rs_cores.readParameters('../../../modelparameter/structural/rootsystem/'+name + ".xml")
            
            params = rs.getOrganRandomParameter(pb.root)
            params_cores = rs_cores.getOrganRandomParameter(pb.root)
           
            if len(params)>4: 
                axroots = axroots_
                for k in range(0, len(axroots)): 
                    params[axroots[k]].tropismN=tropismN
                    params[axroots[k]].tropismS=tropismS
                    params_cores[axroots[k]].tropismN=tropismN
                    params_cores[axroots[k]].tropismS=tropismS
            else:                   
                params[1].tropismN=tropismN
                params[1].tropismS=tropismS
                params_cores[1].tropismN=tropismN                                 
                params_cores[1].tropismS=tropismS

            seed = rs.getOrganRandomParameter(pb.seed)[0]
            seed_cores = rs_cores.getOrganRandomParameter(pb.seed)[0]
            seed.seedPos = pb.Vector3d(distr * j-((N-1)*distr/2), distp * i-(distp*M/2), -3.) # ields perpendicular to tubes
            seed_cores.seedPos = pb.Vector3d(distr * j-((N-1)*distr/2), distp * i-(distp*M/2), -3.) # fields perpendicular to tubes
            rs.setGeometry(fieldbox_wo_rhizotubes)
            rs.initialize(False)
            rs_cores.initialize(False)
            rs.simulate(simtime, False)
            rs_cores.simulate(simtime, False)
            allRS.append(rs) 
            allRS_cores.append(rs_cores) 
            if i+ j == 0:
                allAna = pb.SegmentAnalyser(rs) 
                allAna_cores = pb.SegmentAnalyser(rs_cores) 
            else:
                allAna.addSegments(rs)  
                allAna_cores.addSegments(rs_cores)  # collect all in a segAna object
    allAna.mapPeriodic(w_field, l_field) #periodic outer boundaries 
    allAna_cores.mapPeriodic(w_field, l_field) #periodic outer boundaries 
    
    #fieldboxes
    fieldbox = pb.SDF_PlantBox(w_field, l_field, hlayer) #fieldbox for vRLD and soil cores
    fieldbox_tot = pb.SDF_PlantBox(w_field, l_field, depth) #fieldbox for RLD of the complete field 
    
    # Computation of coefficient of variation 
    cubelen = 20 #cm
    rld_CV = np.zeros((len(z_), int(w_field/cubelen), int(l_field/cubelen)))
    fieldbox_CV = pb.SDF_PlantBox(cubelen, cubelen, hlayer) #fieldbox for vRLD and soil cores
    CV = np.zeros((len(z_)))
    for i in range(0,len(z_)):
        for j in range(int(w_field/cubelen)):  
            for k in range(int(l_field/cubelen)): 
                ana_CV= pb.SegmentAnalyser(allAna)
                fieldbox_CV_ = pb.SDF_RotateTranslate(fieldbox_CV, pb.Vector3d(-w_field/2+cubelen/2+j*cubelen, -l_field/2+cubelen/2+k*cubelen, z_[i]+hlayer /2 )) 
                rs.setGeometry(fieldbox_CV_)
                ana_CV.crop(fieldbox_CV_)
                rld_CV[i,j,k]=ana_CV.getSummed("length")/(cubelen* cubelen* hlayer) #root length density cm/cm3
        CV[i] = np.std(rld_CV[i,:,:])/np.mean(rld_CV[i,:,:])     
        
    #Computation of anisotropy factor 
    An = np.zeros(len(z_))
    for i in range(len(z_)):
        ana_aniso = pb.SegmentAnalyser(allAna)
        fieldbox_ = pb.SDF_RotateTranslate(fieldbox, 0, pb.SDF_Axis.yaxis, pb.Vector3d(0, 0, z_[i]+hlayer /2 )) 

        Na = 0
        for j in range(0, int(np.floor((w_field/distl))+1)): 
            soil_layer_ = pb.SDF_HalfPlane(planes[0][0], planes[0][1], planes[0][2]) 
            soil_layer = pb.SDF_RotateTranslate(soil_layer_, pb.Vector3d(-w_field/2+j*distl, -l_field/2, -depth)) 
            ana_x = ana_aniso.cut(soil_layer)
            Na = Na + len(ana_x.getParameter("length"))
            
        Nb = 0
        for j in range(0, int(np.floor((l_field/distl))+1)): 
            soil_layer_ = pb.SDF_HalfPlane(planes[1][0], planes[1][1], planes[1][2])  
            soil_layer = pb.SDF_RotateTranslate(soil_layer_, pb.Vector3d(-w_field/2, -l_field/2+j*distl, -depth)) 
            ana_y = ana_aniso.cut(soil_layer)
            Nb = Nb + len(ana_y.getParameter("length"))
            
        Nc = 0
        for j in range(0, int(np.floor((hlayer/distl))+1)): 
            soil_layer_ = pb.SDF_HalfPlane(planes[2][0], planes[2][1], planes[2][2]) 
            soil_layer = pb.SDF_RotateTranslate(soil_layer_, pb.Vector3d(-w_field/2, -l_field/2, z_[i]+hlayer /2-distl*j)) 
            ana_z = ana_aniso.cut(soil_layer)
            Nc = Nc + len(ana_z.getParameter("length"))
        
        Nm = (Na+Nb+Nc)/3
        if Nm>0: 
            An[i] = (((Na-Nm)**2+(Nb-Nm)**2+(Nc-Nm)**2)/Nm**2)**0.5*1/6**0.5
        else: 
            An[i] = np.nan
            
    #Computation of relative RLD (to describe a systematic trend) and coefficient of variation (CV, to describe clustering) 
    rrld = np.zeros(len(z_))
    ana_tot = pb.SegmentAnalyser(allAna)
    ana_tot.crop(fieldbox_tot)
    RLD_tot = ana_tot.getSummed("length")/(w_field* l_field* depth) #root length density cm/cm3
    
    for i in range(len(z_)):
        ana_vrld = pb.SegmentAnalyser(allAna)
        fieldbox_ = pb.SDF_RotateTranslate(fieldbox, pb.Vector3d(0, 0, z_[i]+hlayer /2 )) 
        ana_vrld.crop(fieldbox_)
        
        #compute rRLD
        rrld[i]=(ana_vrld.getSummed("length")/(w_field* l_field* hlayer))/RLD_tot #relative root length density (-)
    # Computation of vRLD, whole field 
    vrld = np.zeros(len(z_))
    vrvd = np.zeros(len(z_))
    for i in range(len(z_)):
        ana_vrld = pb.SegmentAnalyser(allAna)
        fieldbox_ = pb.SDF_RotateTranslate(fieldbox, pb.Vector3d(0, 0, z_[i]+hlayer /2 )) 
        ana_vrld.crop(fieldbox_)
        vrld[i]=ana_vrld.getSummed("length")/(w_field* l_field* hlayer) #root length density cm/cm3
        vrvd[i]=ana_vrld.getSummed("volume")/(w_field* l_field* hlayer) #root volume density cm3/cm3
        
    #Computation of vRLD, 2 and 4  soil cores 
    vrld_sc2 = np.zeros(len(z_))
    vrld_sc4 = np.zeros(len(z_))
    vrld_sc_ir = np.zeros(len(z_))
    vrld_sc_ip = np.zeros(len(z_))
    vrvd_sc2 = np.zeros(len(z_))
    vrvd_sc4 = np.zeros(len(z_))
    vrvd_sc_ir = np.zeros(len(z_))
    vrvd_sc_ip = np.zeros(len(z_))
    for i in range(len(z_)):
        ana_vrld_sc1 = pb.SegmentAnalyser(allAna_cores)
        sc1_ = pb.SDF_RotateTranslate(sc1, pb.Vector3d(0, 0, z_[i]+hlayer /2 )) 
        ana_vrld_sc1.crop(sc1_)
        vrld_sc2[i] = ana_vrld_sc1.getSummed("length")/(sc_volume*2)
        vrvd_sc2[i] = ana_vrld_sc1.getSummed("volume")/(sc_volume*2)
        
        ana_vrld_sc2 = pb.SegmentAnalyser(allAna_cores)
        sc2_ = pb.SDF_RotateTranslate(sc2, pb.Vector3d(0, 0, z_[i]+hlayer /2 )) 
        ana_vrld_sc2.crop(sc2_)
        vrld_sc4[i] = ana_vrld_sc2.getSummed("length")/(sc_volume*4)
        vrvd_sc4[i] = ana_vrld_sc2.getSummed("volume")/(sc_volume*4)
        
        ana_vrld_sc3 = pb.SegmentAnalyser(allAna_cores)
        sc3_ = pb.SDF_RotateTranslate(sc3, pb.Vector3d(0, 0, z_[i]+hlayer /2 )) 
        ana_vrld_sc3.crop(sc3_)
        vrld_sc_ir[i] = ana_vrld_sc3.getSummed("length")/(sc_volume*2)
        vrvd_sc_ir[i] = ana_vrld_sc3.getSummed("volume")/(sc_volume*2)        
        
        ana_vrld_sc4 = pb.SegmentAnalyser(allAna_cores)
        sc4_ = pb.SDF_RotateTranslate(sc4, pb.Vector3d(0, 0, z_[i]+hlayer /2 )) 
        ana_vrld_sc4.crop(sc4_)
        vrld_sc_ip[i] = ana_vrld_sc4.getSummed("length")/(sc_volume*2)
        vrvd_sc_ip[i] = ana_vrld_sc4.getSummed("volume")/(sc_volume*2)                  

    #Computation of pRLD 
    imgsa = imgw*imgl #image surface area 
    imgbox = pb.SDF_PlantBox(imgw, imgl,10)  # box inhouse camera system
    if size == 'complete':
        imgbox = pb.SDF_PlantBox(imgw, imgl,30)
    
    #rhizotron images (taken 80° clockwise and counterclockwise) 
    imgs_standard = []; imgs_lowr = []; imgs_highr = []; imgs_cont = []
    for i in range(0, len(z_)): 
        dummy1 = []; dummy2 = []; dummy3 = []; dummy4 = []
        for j in range(0, len(x_standard)): 
            if size == 'complete': 
                dummy1.append(pb.SDF_RotateTranslate(imgbox, 90, pb.SDF_Axis.xaxis, pb.Vector3d(x_standard[j], y_[i]-15, z_[i])))
            else: 
                dummy1.append(pb.SDF_RotateTranslate(imgbox, 180-80, pb.SDF_Axis.xaxis, pb.Vector3d(x_standard[j], y_[i], z_[i])))
                dummy1.append(pb.SDF_RotateTranslate(imgbox, 180+80, pb.SDF_Axis.xaxis, pb.Vector3d(x_standard[j], y_[i], z_[i])))
        for j in range(0, len(x_lowr)):   
            if size == 'complete': 
                dummy2.append(pb.SDF_RotateTranslate(imgbox, 90, pb.SDF_Axis.xaxis, pb.Vector3d(x_highr[j], y_[i]-15, z_[i])))
            else: 
                dummy2.append(pb.SDF_RotateTranslate(imgbox, 180-80, pb.SDF_Axis.xaxis, pb.Vector3d(x_highr[j], y_[i], z_[i])))
                dummy2.append(pb.SDF_RotateTranslate(imgbox, 180+80, pb.SDF_Axis.xaxis, pb.Vector3d(x_highr[j], y_[i], z_[i])))
        for j in range(0, len(x_highr)):   
            if size == 'complete': 
                dummy3.append(pb.SDF_RotateTranslate(imgbox, 90, pb.SDF_Axis.xaxis, pb.Vector3d(x_highr[j], y_[i]-15, z_[i])))
            else: 
                dummy3.append(pb.SDF_RotateTranslate(imgbox, 180-80, pb.SDF_Axis.xaxis, pb.Vector3d(x_highr[j], y_[i], z_[i])))
                dummy3.append(pb.SDF_RotateTranslate(imgbox, 180+80, pb.SDF_Axis.xaxis, pb.Vector3d(x_highr[j], y_[i], z_[i])))
        if size == 'complete': 
            dummy4.append(pb.SDF_RotateTranslate(pb.SDF_PlantBox(img_cont, imgl,30), 90, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[i]-15, z_[i])))            
        else: 
            dummy4.append(pb.SDF_RotateTranslate(pb.SDF_PlantBox(img_cont, imgl,10), 180-80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[i], z_[i])))
            dummy4.append(pb.SDF_RotateTranslate(pb.SDF_PlantBox(img_cont, imgl,10), 180+80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[i], z_[i])))
        imgs_standard.append(dummy1); imgs_lowr.append(dummy2); imgs_highr.append(dummy3); imgs_cont.append(dummy4)
    imgs_standard = np.array(imgs_standard); imgs_lowr = np.array(imgs_lowr); imgs_highr = np.array(imgs_highr);  imgs_cont = np.array(imgs_cont)
    
    ana_prld = pb.SegmentAnalyser(allAna)
    ana_prld.crop(pb.SDF_Difference(rhizotubeso, rhizotubes))
    prld_stand = np.zeros(len(z_)); prld_lowr = np.zeros(len(z_)); prld_highr = np.zeros(len(z_)); prld_cont = np.zeros(len(z_))
    prvd_stand = np.zeros(len(z_)); prvd_lowr = np.zeros(len(z_)); prvd_highr = np.zeros(len(z_)); prvd_cont = np.zeros(len(z_))
    prcd_stand = np.zeros(len(z_)); prcd_lowr = np.zeros(len(z_)); prcd_highr = np.zeros(len(z_)); prcd_cont = np.zeros(len(z_))
    for i in range(np.shape(imgs_standard)[0]):
        prld_layer_stand = []; prld_layer_lowr = []; prld_layer_highr = []; prld_layer_cont = []
        prvd_layer_stand = []; prvd_layer_lowr = []; prvd_layer_highr = []; prvd_layer_cont = []
        prcd_layer_stand = []; prcd_layer_lowr = []; prcd_layer_highr = []; prcd_layer_cont = []
        #prld / prvd, standard resolution imgs
        for j in range(np.shape(imgs_standard)[1]):
            ana_prld1= pb.SegmentAnalyser(ana_prld)
            ana_prld1.crop(imgs_standard[i,j])
            ana_prld1.pack()
            #ana_prld1.write('test_results/prld_stand_'+str(j)+'.vtp')
            segs = ana_prld1.segments
            nodes = ana_prld1.nodes
            radii = ana_prld1.getParameter("radius")
            length = ana_prld1.getParameter("length")
            volume = np.array(radii)**2*np.pi*np.array(length)
            leng_sum = 0; vol_sum = 0
            dist = np.zeros((len(segs)))
            count = 0
            for m, s in enumerate(segs):
                s1 = segs[m]
                n1, n2 = nodes[s1.x], nodes[s1.y]
                D1 = rhizotubeso.getDist(pb.Vector3d(n1.x, n1.y, n1.z))
                D2 = rhizotubeso.getDist(pb.Vector3d(n2.x, n2.y, n2.z))
                if D1<radii[m] and D2<radii[m]: 
                    leng_sum = leng_sum + length[m]
                    vol_sum = vol_sum + volume[m]
                    count = count +1
                elif D1<radii[m] or D2<radii[m]: 
                    leng_sum = leng_sum + length[m]/2
                    vol_sum = vol_sum +volume[m]/2
                    count = count +0.5
            prld_layer_stand.append(leng_sum/imgsa) #root length density
            prvd_layer_stand.append(vol_sum/imgsa) #root volume density
            prcd_layer_stand.append(vol_sum/imgsa) #root count density
            
        #prld / prvd, low resolution imgs    
        for j in range(np.shape(imgs_lowr)[1]):
            ana_prld1= pb.SegmentAnalyser(ana_prld)
            ana_prld1.crop(imgs_lowr[i,j])
            ana_prld1.pack()
            segs = ana_prld1.segments
            nodes = ana_prld1.nodes
            radii = ana_prld1.getParameter("radius")
            length = ana_prld1.getParameter("length")
            volume = np.array(radii)**2*np.pi*np.array(length)
            leng_sum = 0; vol_sum = 0
            dist = np.zeros((len(segs)))
            count = 0
            for m, s in enumerate(segs):
                s1 = segs[m]
                n1, n2 = nodes[s1.x], nodes[s1.y]
                D1 = rhizotubeso.getDist(pb.Vector3d(n1.x, n1.y, n1.z))
                D2 = rhizotubeso.getDist(pb.Vector3d(n2.x, n2.y, n2.z))
                if D1<radii[m] and D2<radii[m]: 
                    leng_sum = leng_sum + length[m]
                    vol_sum = vol_sum + volume[m]
                    count = count +1
                elif D1<radii[m] or D2<radii[m]: 
                    leng_sum = leng_sum + length[m]/2
                    vol_sum = vol_sum +volume[m]/2
                    count = count +0.5
            prld_layer_lowr.append(leng_sum/imgsa) #root length density
            prvd_layer_lowr.append(vol_sum/imgsa) #root volume density
            prcd_layer_lowr.append(count/(img_cont*imgl)) #root count density 

        #prld / prvd, high resolution imgs    
        for j in range(np.shape(imgs_highr)[1]):
            ana_prld1= pb.SegmentAnalyser(ana_prld)
            ana_prld1.crop(imgs_highr[i,j])
            ana_prld1.pack()
            segs = ana_prld1.segments
            nodes = ana_prld1.nodes
            radii = ana_prld1.getParameter("radius")
            length = ana_prld1.getParameter("length")
            volume = np.array(radii)**2*np.pi*np.array(length)
            leng_sum = 0; vol_sum = 0
            dist = np.zeros((len(segs)))
            count = 0
            for m, s in enumerate(segs):
                s1 = segs[m]
                n1, n2 = nodes[s1.x], nodes[s1.y]
                D1 = rhizotubeso.getDist(pb.Vector3d(n1.x, n1.y, n1.z))
                D2 = rhizotubeso.getDist(pb.Vector3d(n2.x, n2.y, n2.z))
                if D1<radii[m] and D2<radii[m]: 
                    leng_sum = leng_sum + length[m]
                    vol_sum = vol_sum + volume[m]
                    count = count +1
                elif D1<radii[m] or D2<radii[m]: 
                    leng_sum = leng_sum + length[m]/2
                    vol_sum = vol_sum +volume[m]/2
                    count = count +0.5
            prld_layer_highr.append(leng_sum/imgsa) #root length density
            prvd_layer_highr.append(vol_sum/imgsa) #root volume density
            prcd_layer_highr.append(count/(img_cont*imgl)) #root count density 
            
        #prld / prvd, continuous imgs
        for j in range(np.shape(imgs_cont)[1]):        
            ana_prld1= pb.SegmentAnalyser(ana_prld)
            ana_prld1.crop(imgs_cont[i,j])
            ana_prld1.pack()
            segs = ana_prld1.segments
            nodes = ana_prld1.nodes
            radii = ana_prld1.getParameter("radius")
            length = ana_prld1.getParameter("length")
            volume = np.array(radii)**2*np.pi*np.array(length)
            leng_sum = 0; vol_sum = 0
            dist = np.zeros((len(segs)))
            count = 0
            for m, s in enumerate(segs):
                s1 = segs[m]
                n1, n2 = nodes[s1.x], nodes[s1.y]
                D1 = rhizotubeso.getDist(pb.Vector3d(n1.x, n1.y, n1.z))
                D2 = rhizotubeso.getDist(pb.Vector3d(n2.x, n2.y, n2.z))
                if D1<radii[m] and D2<radii[m]: 
                    leng_sum = leng_sum + length[m]
                    vol_sum = vol_sum + volume[m]
                    count = count +1
                elif D1<radii[m] or D2<radii[m]: 
                    leng_sum = leng_sum + length[m]/2
                    vol_sum = vol_sum +volume[m]/2
                    count = count +0.5
            prld_layer_cont.append(leng_sum/(img_cont*imgl)) #root length density
            prvd_layer_cont.append(vol_sum/(img_cont*imgl)) #root volume density
            prcd_layer_cont.append(count/(img_cont*imgl)) #root count density 

        prld_stand[i] = np.mean(prld_layer_stand); prld_lowr[i] = np.mean(prld_layer_lowr); prld_highr[i] = np.mean(prld_layer_highr); prld_cont[i] = np.mean(prld_layer_cont)    
        prvd_stand[i] = np.mean(prvd_layer_stand); prvd_lowr[i] = np.mean(prvd_layer_lowr); prvd_highr[i] = np.mean(prvd_layer_highr); prvd_cont[i] = np.mean(prvd_layer_cont) 
        prcd_stand[i] = np.mean(prcd_layer_stand); prcd_lowr[i] = np.mean(prcd_layer_lowr); prcd_highr[i] = np.mean(prcd_layer_highr); prcd_cont[i] = np.mean(prcd_layer_cont); 
        
    #Results 
    result=np.zeros((len(z_),26))
    
    #general params 
    result[:,0] = z_
    result[:,1] = tropismN
    result[:,2] = tropismS
    result[:,3] = rrld #relative root elngth density 
    result[:,4] = An #anisotropy factor 
    result[:,5] = CV #anisotropy factor 
    
    #root length density 
    result[:,6] = prld_stand #pRLD standard images
    result[:,7] = prld_lowr #pRLD lowr images
    result[:,8] = prld_highr #pRLD highr images
    result[:,9] = prld_cont #pRLD continuous image
   
   #volumetric root length density 
    result[:,10] = vrld #vRLD
    result[:,11] = vrld_sc2 #vRLD
    result[:,12] = vrld_sc4 #vRLD
    result[:,13] = vrld_sc_ir #vRLD
    result[:,14] = vrld_sc_ip #vRLD
       
    #planar root volume density
    result[:,15] = prvd_stand 
    result[:,16] = prvd_lowr 
    result[:,17] = prvd_highr 
    result[:,18] = prvd_cont 
   
    #volumetric root volume density
    result[:,19] = vrvd 
    result[:,20] = vrvd_sc2 
    result[:,21] = vrvd_sc4          
    
    #planar root count density 
    result[:,22] = prcd_stand #pRCD standard images
    result[:,23] = prcd_lowr #pRCD lowr images
    result[:,24] = prcd_highr #pRCD highr images
    result[:,25] = prcd_cont #pRCD continuous image
    
    return result

