""" coupling with DuMux as solver for the soil part, dumux-rosi must be installed & compiled """
import sys; sys.path.append("../.."); sys.path.append("../../src/")
sys.path.append("../../../dumux-rosi/build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../../dumux-rosi/python/modules/")  # python wrappers
import argparse
import plantbox as pb
import numpy as np
import timeit
from scenario_static_setup_mpi import *
import multiprocessing as mp
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()

def seg2cell(j): 
    return plant.seg2cell[j]

def simulate_sra(name, sim_time, out_time, rs_age, trans, wilting_point, soil, s, hm, plant, soil_, target_x, target_y, res, cell_number, cellvol, X, Y, Z, wc_root,soil_depth, save_npz:bool, save_vtr:bool):
    
    
    """ Numerical solution """
    dt = 360. / (24 * 3600)  # [days]  
    x_, y_, z_ = [], [], []
    
    hs = s.getSolutionHead()  # [cm] matric potential 
    if rank == 0: 
        
        peri = PerirhizalPython(hm.ms)  
        peri.open_lookup("lookup/"+soil_) 
        outer_r = peri.get_outer_radii("length") 
        inner_r = peri.ms.radii
        rho_ = np.divide(outer_r, np.array(inner_r))  
        
        hs_ = hm.ms.getHs(hs)  # [cm] matric potential per segment 
        hsr = hs_.copy()  # initial values for fix point iteration

    source_water = None 
    hx = None
    N = round(sim_time / dt)
    t = 0.
    dummy = 0
    for i in range(0, N):  
        
        if rank == 0: 
            hx = hm.solve(rs_age + t, -trans * sinusoidal(t), hsr, cells = False)  
            hx_old = hx.copy()

            kr_ = hm.params.getKr(rs_age + t)  
            inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up; here const

            err = 1.e6
            c = 0
            while err > 100. and c < 100:  

                """ interpolation """
                hsr = peri.soil_root_interface_potentials(hx[1:], hs_, inner_kr_, rho_)  

                """ xylem matric potential """
                hx = hm.solve_again(rs_age + t, -trans * sinusoidal(t), hsr, cells = False)  
                err = np.linalg.norm(hx - hx_old)
                hx_old = hx.copy()
                c += 1  
            
         
            fluxes = hm.radial_fluxes(rs_age + t, hx, hsr, cells = False) 
            source_water = hm.sumSegFluxes(fluxes)
        
        water = s.getWaterVolume() 
        source_water = comm.bcast(source_water, root = 0)
        s.setSource(source_water)  
        s.solve(dt)
        soil_water = (s.getWaterVolume() - water) / dt 

        hs = s.getSolutionHead()  # per cell
        wc = s.getWaterContent()
        if rank == 0:
            x_.append(t)  
            y_.append(hm.get_transpiration(rs_age + t, hx.copy(), hsr.copy()))  # cm3/day
            z_.append(soil_water)  # cm3/day 

            n = round(float(i) / float(N) * 100.) 
            print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], {:g} iterations, soil hs [{:g}, {:g}], interface [{:g}, {:g}] cm, root [{:g}, {:g}] cm, {:g} days"
                  .format(c, np.min(hs), np.max(hs), np.min(hsr), np.max(hsr), np.min(hx), np.max(hx), s.simTime))

            if i == 0 or (np.around(s.simTime,2) in out_time and np.around(s.simTime,2) != dummy): #should be in outtime, but the same time should not be saved twice (thus dummy)
                
                dummy = np.around(s.simTime,2)

                #reshape  water content
                wc = np.reshape(wc, (cell_number[2],cell_number[1],cell_number[0]))
                wc = np.swapaxes(wc,0,2)
                
                #reshape SWP
                hs = np.reshape(hs, (cell_number[2],cell_number[1],cell_number[0]))
                hs = np.swapaxes(hs,0,2)

                #get root params
                ana = pb.SegmentAnalyser(plant.mappedSegments())
                segs = ana.segments
                radius = np.array(ana.getParameter("radius"))
                seglen = np.array(ana.getParameter("length"))

                #compute root volume and root volume fraction
                start = timeit.default_timer()
                x = np.array([plant.seg2cell[j] for j in range(0, len(segs))])      #slow!!!! 
                # pool_obj = mp.Pool()
                # a = (pool_obj.map(seg2cell,range(0,len(segs))))
                # x = np.array(a)
                # pool_obj.close()

                sc = np.unique(x).astype(int) #all voxels with a segment 
                rootvol = np.zeros((cell_number[0]*cell_number[1]*cell_number[2]))
                for k in range(0,len(sc)):
                    sg = plant.cell2seg[sc[k]] #numbers of all segments within that voxel
                    rootvol_ind = np.array([radius[sg[j]] ** 2 * np.pi * seglen[sg[j]] for j in range(0,len(sg))])
                    rootvol[sc[k]]= np.sum(rootvol_ind)
                    
                frac = rootvol / cellvol
                frac[frac>1] = 1
                print('Computation of root volume and root volume fraction took ' +str((timeit.default_timer()-start)/60)+ ' minutes')
                
               
                #reshaping is needed
                frac = np.reshape(frac, (cell_number[2],cell_number[1],cell_number[0]))
                frac = np.swapaxes(frac,0,2)
                rootvol = np.reshape(rootvol, (cell_number[2],cell_number[1],cell_number[0]))
                rootvol = np.swapaxes(rootvol,0,2)

                #stitch the single plant domain to the target domain
                x_stitch = int(np.ceil((target_x/res)/np.shape(wc)[0]))
                y_stitch = int(np.ceil((target_y/res)/np.shape(wc)[1]))
                wc_stitch = np.tile(wc, (x_stitch, y_stitch, 1))
                hs_stitch = np.tile(hs, (x_stitch, y_stitch, 1))
                frac_stitch = np.tile(frac, (x_stitch, y_stitch, 1))
                rootvol_stitch = np.tile(rootvol, (x_stitch, y_stitch, 1))
                
                #check if stitched domain has the correct size and crop it if not
                if np.shape(wc_stitch)[0]>(target_x/res) or np.shape(wc_stitch)[1]>(target_y/res) or np.shape(wc_stitch)[2]>(-soil_depth/res): 
                    wc_stitch = wc_stitch[:int(target_x/res),:int(target_y/res), :int(-soil_depth/res)]
                    hs_stitch = hs_stitch[:int(target_x/res),:int(target_y/res), :int(-soil_depth/res)]
                    frac_stitch = frac_stitch[:int(target_x/res),:int(target_y/res), :int(-soil_depth/res)]
                    rootvol_stitch = rootvol_stitch[:int(target_x/res),:int(target_y/res), :int(-soil_depth/res)]
                
                if save_npz: 
                    write_npz(name, t, wc, hs, frac, rootvol, wc_stitch, hs_stitch, frac_stitch, rootvol_stitch, target_x,target_y, soil_depth)
                if save_vtr: 
                    write_vtr(name, t, target_x,target_y, X, Y, Z, wc_root, wc, hs, rootvol, wc_stitch, hs_stitch, rootvol_stitch, plant, res, soil_depth) 
                    
                
        t += dt


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = 'Simulation options')
    parser.add_argument('plant', type = str, help = 'maize or wheat')
    parser.add_argument('res', type = str, help = 'high or low')
    parser.add_argument('soil', type = str, help = 'soil type (hydrus_loam, hydrus_clay or hydrus_sandyloam)')

    args = parser.parse_args()

    initial = -200 #cm #initial can be 'initial water content (-)' or initial soil water potential (cm)'- the model recognizes if it is the one or the other                          
    trans = 0.5 #cm/day
    infiltration = False #scenario with infiltration only
    evaporation = True #scenario with evaporation only                                 
    rs_age = 70 #d
    sim_time = 14.5  # day
    out_time = np.arange(0.5, sim_time,1) #additionally, the first time step is always saved (--> corresponds to root system with static swc)
    save_npz = True
    save_vtr = False
    
    name = args.plant + "_" + args.res + "_resolution_" + args.soil+'_age'+str(rs_age)
    if infiltration: 
        name = name+'_inf'
    elif evaporation: 
        name = name+'_evap'        
    print()
    print(name, "\n")
    

    wilting_point, soil, s, hm, plant, target_x,target_y, cell_number, cellvol, X, Y, Z, wc_root, res , soil_depth, soil_ = set_scenario(args.plant, args.res, args.soil, initial, trans, rs_age, infiltration, evaporation)

    simulate_sra(name, sim_time, out_time, rs_age, trans, wilting_point, soil, s, hm, plant, soil_, target_x,target_y, res, cell_number, cellvol, X, Y, Z, wc_root, soil_depth, save_npz, save_vtr)
    
    
    
    
    # mpiexec -n 4 python3 scenario_static_mpi.py maize high hydrus_loam