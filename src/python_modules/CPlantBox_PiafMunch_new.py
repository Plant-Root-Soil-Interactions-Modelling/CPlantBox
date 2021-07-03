################### Convert CPlantBox -> PiafMunch and post processing ###################
#\   \ /\  ____ CPlantBox: creates the plant structure. Link: cplantbox.com              #
# \___\|| /   / PiafMunch: carbon and water flow simulation in structure. Link:          #
#      ||/___/  https://www6.ara.inrae.fr/piaf_eng/Methods-and-Models/PiafMunch          #
#  ___ ||                                                                                #  
#  \__\|| ___   Environment:                                                             #  
#      ||/__/               Python==3.8                                                  #
#      ||                                                                                #  
#	  //|                                                                                #   
#   _//|\\                                                                               #
# __// ||\\___                                                                           #  
#/ // //| \\_ \___                                                                       #
# /|\//|\\ \ \ \_ \                                                                      #
#/ |// ||\\ \ \__\                                                                       #
#  /\ \| /|\| \\ \__ Author: Xiao-Ran Zhou  zxrzxr@gmail.com  (adapted by m giraud)      #
##########################################################################################

############################## Import system libraries ###################################
import os
import sys
import math
import timeit
import numpy as np
from scipy.linalg import norm
import datetime
import pandas as pd
import matplotlib.pylab as plt
import plotly
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import xml.etree.ElementTree as ET
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d
from scipy import stats
############################## system libraries imported ###################################

############################## Load precompiled CPlantBox ###################################
sys.path.append("../..") # adding the path of the python binding
import plantbox as pb # CPlantBox Python Binding

################################# Defined functions ########################################



# python_nodes
#returns : node num, org type of node 1,num connections, node coord 



def getLengthsurfaces( plant ):
	nodes = plant.get_nodes()
	segs = plant.get_segments()
	lengths = [norm(nodes[seg[0]]-nodes[seg[1]]) for seg in segs]
	surfacesPhloem = (np.pi * (np.full(len(lengths), 0.0013))**2) #0.0013  or other phliem radius?
	surfacesOrg = (np.pi * (np.array(plant.rs.radii))**2) #1e-2 to go from stem to phloem/xylem surface
	#print(np.concatenate(([lengths[0]],lengths)), np.concatenate(([surfaces[0]],surfaces)))
	return np.concatenate(([lengths[0]],lengths)), np.concatenate(([surfacesPhloem[0]],surfacesPhloem)), np.concatenate(([surfacesOrg[0]],surfacesOrg)) #add surface and length for seed node

def write_PiafMunch_parameter(plant,initphi=[], name = "mogi_test", time=36, nopsi=False): #function of creating parameters of PiafMunch
	nodes = plant.get_nodes()
	Zed = [xi[2] for xi in nodes]
	f = open("coords.txt",'w')
	f.write(str(Zed)[1:-1]+"\n")
	f.close()
	segments = plant.get_segments()
	organTypes = plant.get_organ_types()
	subTypes = plant.get_subtypes()
	nodes_organtype = np.concatenate(([1], organTypes))#seg_Indx = node_y_Indx - 1. just need to add organ type for seed (set to type 'root')     #np.zeros(len(nods))#-1)
	nodes_subtype = np.concatenate(([1], subTypes))#
	f = open("ots.txt",'w')
	f.write(','.join([num for num in map(str,nodes_organtype)])  +'\n')
	f.close()
	lengths, surfaces,surfacesOrg = getLengthsurfaces( plant )
	Nt = len(nodes) #total number of nodes, here the 0 and 1st node are the same, so minus 1 => not true anymore
	Nc = len(segments) #total number of connections = segments
    
	tiproots, tipstems, tipleaves = plant.get_organ_nodes_tips()
	N1L_node = np.concatenate((tipstems, tipleaves)) 
	N1R_node = tiproots 
	N3 = plant.get_nodes_child_base() 
	N2 = np.setdiff1d(np.array([range(len(nodes))]), np.concatenate((N1L_node, N1R_node, N3)))
	print(N1L_node, N1R_node,N2,N3 )

	#'******** CARBON Lateral FLUX - RELATED PARAMETERS *********\n'
	#initialization of the parameters
	kML = np.zeros(len(nodes_organtype))
	vML = np.zeros(len(nodes_organtype))
	kMU = np.zeros(len(nodes_organtype))
	vMU = np.zeros(len(nodes_organtype))
	kMParMb = np.zeros(len(nodes_organtype))
	vMParMb = np.zeros(len(nodes_organtype))
	kM = np.zeros(len(nodes_organtype)) #kinetic parameter / Michaelis - starch Synthesis
	Vmax = np.zeros(len(nodes_organtype)) # kinetic parameter / starch Synthesis
	C_targ = np.zeros(len(nodes_organtype)) #kinetic parameter / starch/sugar equilibrium. (regul. par. sugar conc.)			 (mmol / ml)
	kHyd = np.zeros(len(nodes_organtype))
	k1 = np.zeros(len(nodes_organtype))
	k2 = np.zeros(len(nodes_organtype))
	k3 = np.zeros(len(nodes_organtype))
	#k3 = Soil_water
	StructC = np.zeros(len(nodes_organtype))
	vol_ST = np.zeros(len(nodes_organtype))
	volPhlApo = np.zeros(len(nodes_organtype))
	volParApo = np.zeros(len(nodes_organtype))
	k_Lockhart = np.zeros(len(nodes_organtype))
	P_thr = np.zeros(len(nodes_organtype))
	vol_Sympl_max = np.zeros(len(nodes_organtype))

    # resistance = resistivity * L / A
	r_Trsv = 100#np.full(len(nodes_organtype), r_trsv)
	r_PhlMb = 135.785#np.full(len(nodes_organtype), r_phl_mb ) #135.785
	r_ParMb = 1e+025#np.full(len(nodes_organtype), 1e+025)
	r_Apo = 1e+025#np.full(len(nodes_organtype), 1e+025)
	r_Sympl =1e+025# np.full(len(nodes_organtype), 1e+025)
	if(nopsi):
		r_Trsv = 1e-025#np.full(len(nodes_organtype), r_trsv)
		r_PhlMb =1e-025#np.full(len(nodes_organtype), r_phl_mb ) #135.785
        
	radConductivity = 0.3

	#print(len(nodes_organtype), len(lengths), len(vol_ST))
	for i, ot in enumerate(nodes_organtype): #given different value based on whether it is source, sink or connection
		StructC[i]	= plant.rhoC * lengths[i] * surfacesOrg[i]*1000 #mmol C.
		vol_ST[i]	 = lengths[i] * surfaces[i]
		k1[i]		 = 0.
		k2[i]		 = 0.	 #manually set it to 0.4
		kMParMb[i]	= 1
		vMParMb[i]	= 0
		kM[i]		 = 1e-100
		Vmax[i]	   = 0
		C_targ[i]	 = 0.1
		kHyd[i]	   = 0
		k3[i]		 = 0
		volPhlApo[i]  = 2.6e-05
		volParApo[i]  = 2.6e-05
		k_Lockhart[i] = 0
		P_thr[i]	  = 1
		vol_Sympl_max[i] = 0.00018		
		kML[i]	 =  0
		kMU[i]		= 0	  #different in source, sink or connection of piafmunch2
		vML[i]		= 0	  #different in source, sink or connection of piafmunch2
		vMU[i]		= 0	  #different in source, sink or connection of piafmunch2
		if (ot == 4 ): #leafs	   
			vML[i]		= 6.40314e-005#vml*1e3*60*60	  #mol/s to mmol/hr	  #kinetic parameter / phloem loading (mmol /h) different in source, sink or connection of piafmunch2 oringinal value is 6.40314e-006 
			kMU[i] = plant.rhoC *surfacesOrg[i]*1/24*1000*(1/0.75)/2 #growth sink: cm3_dot * rho/M
			#print(ot, i,kMU[i])#,lengths[i] * surfacesOrg[i])
		elif ot == 2:   #roots
			st = nodes_subtype[i]
			vMU[i]		= 0.12*0.2*lengths[i] *2*np.pi * np.array(plant.rs.radii[i-1])#vmu#radConductivity * ((np.array(plant.rs.radii)[i])**2)* lengths[i]# # default is 2.82627e+95
			print(i, st, vMU[i]*1e-3, )
			if i in N1R_node:
				r = 4.4*(st == 1)+ 1 *(st==2)+0.2*(st==3)
				kMU[i] = plant.rhoC *surfacesOrg[i]*r/24*1000 *(1/0.75)
				#print(ot, i,kMU[i])#,lengths[i] * surfacesOrg[i], r, st)
 				#print(ot, i,surfacesOrg[i],plant.rhoC,((np.pi*surfacesOrg[i]**2)*1)/(60*60)*100*(1/0.75), kMU[i])
		elif ot == 3: #stem
			#kMU[i] = plant.rhoC *surfacesOrg[i]*1/24*1000*(1/0.75)/sum(nodes_organtype==3) #growth sink: (cm3 hr-1) * rho * 1000/M
			#print(sum(nodes_organtype==3),nodes_organtype==3)
			if i in N1L_node :
				kMU[i] = plant.rhoC *surfacesOrg[i]*4.4/24*1000*(1/0.75) #growth sink: (cm3 hr-1) * rho * 1000
				print(ot, i,kMU[i])

	#'******** INITIAL VALUES *********\n'
	#initialization of the parameters
	if(len(initphi)==0):
		Q_ST = np.full(len(nodes_organtype), 0)
		Q_Sympl = np.full(len(nodes_organtype), 0)
	else:
		Q_ST = initphi*vol_ST*1000 #mmol to mol
		Q_Sympl = initphi*vol_ST*1000 #mmol to mol
	#Q_Sympl =0#4.4e-006# np.full(len(nodes_organtype), 4.4e-006)
	Starch =0#1# np.full(len(nodes_organtype), 1)
	Q_PhlApo = 0#4.4e-006#np.full(len(nodes_organtype), 4.4e-006)
	Q_ParApo = 0#4.4e-006#np.full(len(nodes_organtype), 4.4e-006)
	Tr_Q_ST = 0#np.full(len(nodes_organtype), 0)
	Tr_Q_Sympl =0# 4.4e-006#np.full(len(nodes_organtype), 4.4e-006)
	Tr_Starch = 0#1#np.full(len(nodes_organtype), 1)
	Tr_Q_PhlApo = 0#np.full(len(nodes_organtype), 0)
	Tr_Q_ParApo = 0#np.full(len(nodes_organtype), 0)
	vol_Sympl =2.6e-005# np.full(len(nodes_organtype), 2.6e-005)

	#******** SIMULATION SOLVING PARAMETERS *********


	#'******** CARBON Lateral FLUX - RELATED PARAMETERS *********\n'
	#initialization of the parameters
	#k3 = Soil_water ?? pas k3 = kinetic parameter / starch/sugar equilibrium. (regul. par. sugar conc.) 			(h-1)

	Q_ST_Abs =  1e-015
	Q_Sympl_Abs =  1e-015
	Starch_Abs =  1e-012
	Q_PhlApo_Abs =  1e-015
	Q_ParApo_Abs =  1e-015
	Tr_Q_ST_Abs =  1e-012
	Tr_Q_Sympl_Abs =  1e-012
	Tr_Starch_Abs =  1e-012
	Tr_Q_PhlApo_Abs =  1e-015
	Tr_Q_ParApo_Abs =  1e-015
	vol_Sympl_Abs = 1e-012



	sys.path.append(".")
	if(nopsi):
		name = name + "nopsi" 
	f = open(name+".ini",'w')
	f.write('******** DESCRIPTION OF ARCHITECTURE *********\n\n')

	f.write("Total number of Nodes : {0} = {1}\n".format('Nt', Nt))
	f.write("number of Internode Connections : {0} = {1}\n\n".format('Nc', Nc))

	f.write("Nodes Of Connectivity Order 1, Transpiring Leaf Ends : N1L = {0}\n".format(len(N1L_node)))
	f.write("node#\tc.node\tconn.#\n")
	for n in N1L_node: #numerotation starts at 1 in PiafMunch
		f.write(f'{n + 1}\t{segments[n - 1][0] +1 }\t{n - 1  +1}\n') #node num before, after and seg num
	f.write('\n')

	f.write("Nodes Of Connectivity Order 1, Absorbing Root Ends : N1R = {0}\n".format(len(N1R_node)))
	f.write("node#\tc.node\tconn.\tr_abs#\n")
	for n in N1R_node:
		f.write(f'{n +1}\t{segments[n -1][0] +1}\t{-(n - 1 +1)}\t{1e-025}\n') #no flux BC?
	f.write('\n')

	f.write('Nodes Of Connectivity Order 2 :  N2 = {0}\n' .format(len(N2) ))
	f.write("node#\tc.nd.1\tconn.1\tc.nd.2\tconn.2\n")
	for n in N2:    
		supSegIdx = np.where([seg[0] for seg in segments] == n)[0]
		if(len(supSegIdx) == 2): #seed node
			f.write(f'{1}\t{2}\t{1}\t{supSegIdx[1] +2}\t{-(supSegIdx[1] +1)}\n')
		else:    
			ot = organTypes[n - 1]
			sign = -int(ot == 2) + int(ot != 2) #is root?
			f.write(f'{n +1}\t{segments[n -1][0]+1}\t{n  * sign}\t{supSegIdx[0] +2}\t{-(supSegIdx[0] +1) * sign}\n')
	f.write('\n')

	f.write("Nodes Of Connectivity Order 3 :  N3 = {0}\n".format(len(N3)))
	f.write("node#\tc.nd.1\tconn.1\tc.nd.2\tconn.2\tc.nd.3\tconn.3\n")
	for n in N3:
		supSegIdx = np.where([seg[0] for seg in segments] == n)[0]
		ot = organTypes[n - 1]
		sign = -int(ot == 2) + int(ot != 2) #is root?
		f.write(f'{n +1}\t{segments[n-1 ][0]+1}\t{n  * sign}\t{supSegIdx[0] + 2}\t{-(supSegIdx[0] +1) * sign}\t{supSegIdx[1] + 2}\t{-(supSegIdx[1] +1) * sign}\n')	
	f.write('\n')

	string = "node#\tc.nd.1\tconn.1\tc.nd.2\tconn.2\tc.nd.3\tconn.3"
	for i in range(4,9): #max connection at node in CPlantBox = 3
		f.write(f'Nodes Of Connectivity Order {i} :  N{i} = 0\n')
		string = string + f'\tc.nd.{i}\tconn.{i}'
		f.write(string)
		f.write('\n\n')
	

	f.write('******** WATER FLUX - RELATED PARAMETERS *********\n\n')
	f.write("T\xb0 K\tvisco=f(C)\tNonLin.Psi+NonZeroSugarVol.\n")
	f.write(f'{plant.TairK}\ttrue\tfalse\t\n')

	f.write(f'InterNode Connections -- Axial Resistances (MPa h / ml) : Nc= {Nc}\n')
	f.write("conn.#\tupfl.node\tdnfl.node\tr_Xyl\tr_ST\n")
    
	for i, ot in enumerate(organTypes):
		upNode = segments[i][0]*(ot == 2) + segments[i][1]*(ot != 2)
		downNode = segments[i][1]*(ot == 2) + segments[i][0]*(ot != 2)
		st = subTypes[i]
		kx = 0.0005356511 # cm2
		muw =  1.005e-9 / 3600 # viscosity pure water (MPa.h)
		r_Xyl = muw/kx * lengths[i+1] / surfaces[i+1] #+1 because of added length for seed node
		konz = 0.5e-3 # mol ml-1
		waterViscosity = 2.414e-5 * (10**(247.8/(293 - 140)))
		sucroseDensity = 1.59 #g/cm³, https://pubchem.ncbi.nlm.nih.gov/compound/Sucrose#section=Density
		sucroseMolarMass = 342.3 #g/mol https://pubchem.ncbi.nlm.nih.gov/compound/Sucrose
		sucroseMolarVolume = sucroseMolarMass/sucroseDensity #cm³/mol like the one used by lacointe 2019
		VolFractSucrose = konz* sucroseMolarVolume
		phloemViscosity = waterViscosity * np.exp((4.68 * 0.956 * VolFractSucrose)/(1 - 0.956 * VolFractSucrose)) #in [Pa.s = N-s/m^2]
		mu = phloemViscosity * 1e-6 /( 60*60 )#Pa s-1 to MPa hr-1
		k = plant.phloemConductivity ## conductivity (cm2)
		r_ST = mu/k * lengths[i+1] / surfaces[i+1] #MPa hr ml-1, +1=> take of the segments
		print(i, ot, r_ST)
		f.write(f'{i + 1}\t{upNode + 1}\t{downNode + 1}\t{r_Xyl}\t{r_ST}\n')#no radial water flux
	f.write('\n')

	f.write('Individual Node : Lateral Resistances (MPa h / ml)\n')
	f.write("{:s}\t{:s}\t{:s}\t{:s}\t{:s}\t{:s}\n".format('node#','r_Trsv','r_PhlMb','r_ParMb','r_Apo', 'r_Sympl' ))
	for i, j in enumerate(nodes_organtype):#check for  trsv and phlmb
		f.write(f'{ 1+i}\t{r_Trsv}\t{r_PhlMb}\t{1e25}\t{1e25}\t{1e25}\n')#.format(nodes_organtype[i][0],r_Trsv[i],r_PhlMb[i],r_ParMb[i],r_Apo[i],r_Sympl[i]))
	f.write('\n')

	f.write('******** CARBON Lateral FLUX - RELATED PARAMETERS *********\n')	
	f.write("{:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s} \n".format('node#','kML(M)','vML(mmol/h)','kMU(M)','vMU(mmol/h)', 'kMParMb(M)','vMParMb(mmol/h)','kM(M)','Vmax(M/h)','C_targ(M)','kHyd(h-1)','k1(h-1)','k2','k3(h-1)','StructC','vol_ST(ml)','volPhlApo,ml','volParApo,ml','k_Lockhart','P_thr(MPa)','vol_Sympl_max,ml' ))	
	for i in range(len(nodes_organtype)):
		f.write(f'{i+1}\t{kML[i]}\t{vML[i]}\t{ kMU[i]}\t{vMU[i]}\t{kMParMb[i]}\t{vMParMb[i]}\t{kM[i]}\
        \t{Vmax[i]}\t{C_targ[i]}\t{kHyd[i]}\t{k1[i]}\t{k2[i]}\t{ k3[i]}\t{StructC[i]}\t{vol_ST[i]}\
        \t{volPhlApo[i]}\t{volParApo[i]}\t{k_Lockhart[i]}\t{P_thr[i]}\t{vol_Sympl_max[i]}\n')
	f.write('\n')



	f.write('******** INITIAL VALUES  *********\n')
	f.write("{:s}  {:s}  {:s}  {:s}  {:s}  {:s} {:s}  {:s}  {:s}  {:s}  {:s}  {:s}\n".format('node#','Q.ST(mmol)','Q.Sympl(mmol)','Starch','Q.PhlApo(mmol)', 'Q.ParApo(mmol)',
												'Tr.Q.ST(mmol)','Tr.Q.Sympl(mmol)','Tr.Starch','Tr.Q.PhlApo(mmol)', 'Tr.Q.ParApo(mmol)','vol_Sympl(ml)' ))
	for i in range(len(nodes_organtype)):
		f.write("{:.0f}  {:e}  {:e} {:e} {:e} {:e}  {:e}  {:e} {:e} {:e} {:e} {:e}\n"
				.format(i+1, Q_ST[i],Q_Sympl[i],Starch,Q_PhlApo,Q_ParApo, 
						Tr_Q_ST,Tr_Q_Sympl,Tr_Starch,Tr_Q_PhlApo,Tr_Q_ParApo,vol_Sympl))
	f.write('\n')	

	f.write('******** SIMULATION SOLVING PARAMETERS *********\n')
	f.write('{:s}  {:s}  {:s}   {:s}  {:s}\n'.format('StartTime','EndTime','OutputStep', 'TracerHalfLife','Rel_Tol'))
	f.write('{:s}  {:s}  {:s}   {:s}  {:s}\n'.format('0', str(time), '0.166667', '0.33967', '1e-007'))
	f.write('\n')	 

	f.write('***Abs_Tols for individual  nodes ***\n')
	f.write("{:s}  {:s}  {:s}  {:s}  {:s}  {:s} {:s}  {:s}  {:s}  {:s}  {:s}  {:s}\n".format('node#','Q.ST(mmol)','Q.Sympl(mmol)','Starch','Q.PhlApo(mmol)', 'Q.ParApo(mmol)',
												'Tr.Q.ST(mmol)','Tr.Q.Sympl(mmol)','Tr.Starch','Tr.Q.PhlApo(mmol)', 'Tr.Q.ParApo(mmol)','vol_Sympl(ml)' ))
	for i in range(len(nodes_organtype)):
		f.write("{:.0f}  {:e}  {:e} {:e} {:e} {:e}  {:e}  {:e} {:e} {:e} {:e} {:e}\n"
				.format(i + 1, Q_ST_Abs,Q_Sympl_Abs,Starch_Abs,Q_PhlApo_Abs,Q_ParApo_Abs, 
						Tr_Q_ST_Abs,Tr_Q_Sympl_Abs,Tr_Starch_Abs,Tr_Q_PhlApo_Abs,Tr_Q_ParApo_Abs,vol_Sympl_Abs))
	f.write('\n')	 

	out_nodes = np.concatenate((N1L_node, N3,N1R_node))
    
	f.write('******** OUTPUT SETTINGS : INDIVIDUAL NODE - LATERAL FLUXES-RELATED VARIABLES *********\n')	
	f.write(f'Nodes selected for plotting  : nsp = {len(out_nodes)}\n' )	
	for i in out_nodes:
		f.write(f'{i +1}\n')
	f.write('\n')	
    
	f.write('individual-Node-related variables selected for plotting : nvp = 4\n')	
	f.write('C_ST (mmol / ml)\n')		
	f.write('JS_PhlMb (mmol / h)\n')		
	f.write('Input (mmol / h)\n')		
	f.write('Resp_Maint (mmol / h)\n\n')		


	f.write(f'Nodes selected for saving  : nss = {len(nodes_organtype)}\n' )	
	for i in range(len(nodes_organtype)):
		f.write(f'{i + 1}\n')
	f.write('\n')	


	f.write('individual-Node-related variables selected for saving : nvs = 18\n')	
	f.write('C_ST (mmol / ml)\n')		
	f.write('JS_PhlMb (mmol / h)\n')			
	f.write('JS_Sympl (mmol / h)\n')			
	f.write('TracerInput (MBq / h)\n')			
	f.write('TracerJS_Apo (MBq / h)\n')			
	f.write('TracerJS_ParMb (MBq / h)\n')			
	f.write('TracerJS_Sympl (MBq / h)\n')			
	f.write('TracerJS_PhlMb (MBq / h)\n')			
	f.write('StarchSyn (mmol eq. Glu / h)\n')			
	f.write('TracerStarch (MBq)\n')			
	f.write('Input (mmol / h)\n')	
	f.write('TracerStarch (MBq)\n')			
	f.write('TracerQ_ST (MBq)\n')			
	f.write('P_ST (MPa)\n')			
	f.write('Psi_ST (MPa)\n')			
	f.write('Starch (mmol eq. Glu)\n')			
	f.write('C_Sympl (mmol / ml)\n')	
	f.write('Resp_Maint (mmol / h)\n\n')		

	f.write('******** OUTPUT SETTINGS : INTERNODE CONNECTION - AXIAL FLUXES-RELATED VARIABLES *********\n')  
	f.write('node-to-node Fluxes selected for plotting  : fsp = {:.0f} \n'.format(len(segments)))
	for i in range(len(segments)):
		f.write("{:.0f}\n".format((i+1)))
	f.write('\n')

	f.write('node-to-node-Fluxes-related variables selected for plotting : fvp = 3\n')	
	f.write('JS_ST (mmol / h)\n')		
	f.write('JW_ST (ml / h)\n')			
	f.write('C_Upflow (mmol / ml)\n')	


	f.write('node-to-node Fluxes selected for saving  : fss = {} \n' .format(len(segments)))
	for i in range(len(segments)):
		f.write("{:.0f}\n".format((i+1)))
	f.write('\n')

	f.write('node-to-node-Fluxes-related variables selected for saving : fvs = 3\n')	
	f.write('JS_ST (mmol / h)\n')		
	f.write('JW_ST (ml / h)\n')  
	f.write('C_Upflow (mmol / ml)\n')

	print('output successful')
	