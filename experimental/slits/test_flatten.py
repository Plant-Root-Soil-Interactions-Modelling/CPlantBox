import numpy as np

nx = 2
ny = 3
nz = 4

A = np.zeros((nx, ny, nz))

for i in range(0, nx):
    for j in range(0, ny):
        for k in range(0, nz):
            lind = nx * ny * k + j * nx + i  # // k*(nx*ny)+j*nx+i, same ordering as RectilinearGrid3D
            A[i, j, k] = lind

print(A)

print()

print(A.flatten('C')) 
print(A.flatten('F'))  # we do column major order (Fortran style) 
