import scipy.sparse.linalg as LA
from scipy import sparse

import py_plantbox as pb
from rb_tools import *

import xylem_flux

# Simulate a root system
name = "Sorghum_bicolor_NA_NA"
dt = 30 # days
plant = pb.Plant()
plant.openXML(name)
plant.initialize()
plant.simulate(dt)

# Create graph
nodes = vv2a(plant.getNodes())/100 # convert from cm to m
rseg = seg2a(plant.getSegments(2)) # root system segments
sseg = seg2a(plant.getSegments(4)) # root system segments

#sseg = seg2a(plant.getShootSegments()) # additional shoot segments
seg = np.vstack((sseg,rseg))

print(sseg)
print(nodes[0,:])
print(rseg)

# Adjacency matrix
A = sparse.coo_matrix((np.ones(seg.shape[0]),(seg[:,0],seg[:,1])))

# Parameters for flux model
rs_Kr = np.array([ 2.e-10, 2.e-10, 2.e-10, 2.e-10, 2.e-10, 2.e-11, 2.e-11 ]) # s/m; root hydraulic radial conductivity per root type
rs_Kz = np.array([ 5.e-14, 5.e-14, 5.e-14, 5.e-14, 5e-14, 5e-14, 5e-14 ]) # m2*s; root hydraulic axial conductivity per root type

soil_psi = -700 # static soil pressure J kg^-1

rho = 1e3 # kg / m^3
g = 1.e-3*9.8065 # m / s^2

pot_trans = np.array([-1.15741e-10]) # # m^3 s^-1 potential transpiration

# Conversions
plant_ana = pb.SegmentAnalyser(plant)
radius = v2a(plant_ana.getScalar("radius"))/100. # convert from cm to m
type = v2a(plant_ana.getScalar("subtype"))
kr = np.array(list(map(lambda t: rs_Kr[int(t)-1], type))) # convert from 'per type' to 'per segment'
kr.resize((kr.shape[0],1))
kz = np.array(list(map(lambda t: rs_Kz[int(t)-1], type)))
kz.resize((kz.shape[0],1))

print(radius)

shoot1 = np.ones((sseg.shape[0],1))
shoot0 = np.ones((sseg.shape[0],1))

# Call back function for soil potential
soil = lambda x,y,z : soil_psi

# glue together shoot and root segments
radius = np.vstack((shoot1,radius))
kr =  np.vstack((shoot0,kr))
kz =  np.vstack((shoot1,kz))

# Calculate fluxes within the root system
Q, b = xylem_flux.linear_system(seg, nodes, radius, kr, kz, rho, g, soil)
Q, b = xylem_flux.bc_neumann(Q, b, np.array([0]), np.array([pot_trans]))
x = LA.spsolve(Q, b) # direct

# Save results into vtp
segP = nodes2seg(nodes,seg,x)# save vtp
axial_flux = xylem_flux.axial_flux(x, seg, nodes, kz, rho, g)
radial_flux = xylem_flux.radial_flux(x, seg, nodes, radius, kr, soil)
net_flux = axial_flux+radial_flux

print(segP)
print(len(segP))


# plant_ana.addUserData(a2v(segP[sseg.shape[0]:-1]),"pressure")
plant_ana.addUserData(a2v(axial_flux[sseg.shape[0]:]),"axial_flux")
plant_ana.addUserData(a2v(radial_flux[sseg.shape[0]:]),"radial_flux")
plant_ana.addUserData(a2v(net_flux[sseg.shape[0]:]),"net_flux")
plant_ana.write("results/example_5a.vtp")
