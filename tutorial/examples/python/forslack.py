from fipy import *
from fipy.meshes.nonUniformGrid1D import NonUniformGrid1D as Grid1D
from fipy.meshes.nonUniformGrid2D import NonUniformGrid2D as Grid2D
from fipy.tools import numerix as np

ny = 6 
nx = 1 
dy = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2] 
dx = 0.2 
translate = [[0.0], [-1.2]]
meshes = Grid2D(ny = ny,nx = nx, dy = dy, dx = dx) + translate


nx = 2 
ny = 1 
dx = [0.05, 0.05]
dy = 0.2 
translate = [[0.2], [-0.2]] 
mesh = Grid2D(ny = ny,nx = nx, dy = dy, dx = dx) + translate

meshes = meshes + mesh

nx = 2 
ny = 1 
dx = [0.05 ,0.05]
dy = 0.2 
translate = [[0.2], [-0.6]] 
mesh = Grid2D(ny = ny,nx = nx, dy = dy, dx = dx) + translate

meshes = meshes + mesh

nx = 2 
ny = 1 
dx = [0.05, 0.05]
dy = 0.2 
translate = [[0.2], [-1.0]] 
mesh = Grid2D(ny = ny,nx = nx, dy = dy, dx = dx) + translate

meshes = meshes + mesh

ny = 3 
nx = 1 
dy = [0.2 ,0.2, 0.2]
dx = 0.2 
mesh = Grid2D(ny = ny,nx = nx, dy = dy, dx = dx) 

meshes = meshes + mesh


ny = 1 
nx = 1 
dx = [1.]
dy = 0.2 
translate = [[0.2], [0.2]]
mesh = Grid2D(ny = ny,nx = nx, dy = dy, dx = dx) + translate

meshes = meshes + mesh
vol = meshes.cellVolumes
phi = CellVariable(mesh= meshes, value = 0.)
Source = CellVariable(mesh= meshes, value = 0.)
Source.constrain(10., where = abs(vol - 0.2) < 0.001) 
#Source.setValue(10., where = abs(vol - 0.2) < 0.001) 
#addition of Source * dt to meshes at each step on the selected cell
dt = 1e-6#0.9*( min(meshes.cellVolumes) ** 2) / (2 *max(Source.faceValue ))
eq = (TransientTerm() ==  DiffusionTerm(coeff= phi.faceValue) + Source )

cumulatedIn = 0. # cumulated input
steps = 10
for _ in range(steps):
    print('steps ', _, dt)
    res= 1e5
    loop =0 
    eq.solve(var = phi, dt= dt)
    
    while loop< 1:#(res > 1e-3* max(abs(Source.value)) or res > 1e-3* max(abs(phi.value)) )and loop < 1000:
            eq.solve(var = phi, dt= dt)
            loop += 1
            #print('res: ',res)
            print('maxRes: ',1e-3* max(abs(Source.value)) , 1e-3* max(abs(phi.value)) )
    
    cumulatedIn += sum(Source * dt * meshes.cellVolumes)
    print('mass balance ', sum(phi * meshes.cellVolumes) - cumulatedIn,'\nphi content in mesh ', sum(phi * meshes.cellVolumes))
    print('phi matrix\n', phi)






