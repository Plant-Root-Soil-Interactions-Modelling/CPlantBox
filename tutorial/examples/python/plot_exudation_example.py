# from pyevtk.hl import gridToVTK
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors
import py_rootbox as rb
from rb_tools import *


def v2a(vd):  # rb.std_vector_double_ to numpy array
    l = np.zeros(len(vd))
    for i in range(0, len(vd)):
        l[i] = vd[i]
    return l


width = 6; depth = 26; xres = 0.1; yres = 0.1; zres = 0.1;
nx = int(width / xres);
ny = int(width / yres);
nz = int(depth / zres);
X = np.linspace(-width / 2, width / 2, nx)
Y = np.linspace(-width / 2, width / 2, ny)
Z = np.linspace(-depth, 0, nz)
#X_, Y_, Z_ = np.meshgrid(X, Y, Z, indexing = "ij")


axialc_mps_straight = np.loadtxt('axialc_mps_straight.csv', delimiter='\t')   # reading in from CSV file
axialc_mps_straight = a2v(axialc_mps_straight) 

axialc_mps = np.loadtxt('axialc_mps.csv', delimiter='\t') 
axialc_mps = a2v(axialc_mps)

data_mps_scipy = np.loadtxt('axialc_mps_scipy.csv', delimiter='\t')
z_scipy = a2v(data_mps_scipy[0,:])
axialc_mps_scipy = a2v(data_mps_scipy[1,:])
data_mls_scipy = np.loadtxt('axialc_mls_scipy.csv', delimiter='\t') 
axialc_mls_scipy = a2v(data_mls_scipy[1,:])

axialc_mls = np.loadtxt('axialc_mls.csv', delimiter='\t') 
axialc_mls = a2v(axialc_mls)

axialc_mls_slow100 = np.loadtxt('axialc_mls_slow100.csv', delimiter='\t') 
axialc_mls_slow100 = a2v(axialc_mls_slow100)

axialc_mls_slow1000 = np.loadtxt('axialc_mls_slow1000.csv', delimiter='\t') 
axialc_mls_slow1000 = a2v(axialc_mls_slow1000)

radialc_mps_straight = np.loadtxt('radialc_mps_straight.csv', delimiter='\t')   
radialc_mps_straight = a2v(radialc_mps_straight) 

radialc_mps = np.loadtxt('radialc_mps.csv', delimiter='\t')   
radialc_mps = a2v(radialc_mps_straight)

radialc_mls = np.loadtxt('radialc_mls.csv', delimiter='\t') 
radialc_mls = a2v(radialc_mls)

radialc_mls_slow100 = np.loadtxt('radialc_mls_slow100.csv', delimiter='\t') 
radialc_mls_slow100 = a2v(radialc_mls_slow100)

radialc_mls_slow1000 = np.loadtxt('radialc_mls_slow1000.csv', delimiter='\t') 
radialc_mls_slow1000 = a2v(radialc_mls_slow1000)

fig1 = plt.figure()
ax = plt.axes()
plt.plot(-Z, axialc_mps_straight, 'k-',label="MPS, velocity-based")
plt.plot(-Z, axialc_mps, 'r--',label="MPS, segment-based")
plt.plot(z_scipy, axialc_mps_scipy, 'g.',label="MPS, scipy")
plt.xlabel('z')
plt.ylabel('c')
plt.legend()

fig2 = plt.figure()
ax = plt.axes()
plt.plot(-Z, axialc_mls, 'y-',label="MLS, 1 cm/d")
plt.plot(z_scipy, axialc_mls_scipy, 'g.',label="MLS, scipy")
plt.plot(-Z, axialc_mls_slow100,'g-',label="MLS, 0.1 cm/d")
plt.plot(-Z, axialc_mls_slow1000,'b-',label="MLS, 0.01 cm/d")
plt.xlabel('z')
plt.ylabel('c')
plt.legend()

fig3 = plt.figure()
ax = plt.axes()
plt.plot(X, radialc_mps_straight, 'k-',label="MPS, velocity-based")
plt.plot(X, radialc_mps, 'r--',label="MPS, segment-based")
plt.plot(X, radialc_mls, 'y-',label="MLS, 1 cm/d")
plt.plot(X, radialc_mls_slow100,'g-',label="MLS, 0.1 cm/d")
plt.plot(X, radialc_mls_slow1000,'b-',label="MLS, 0.01 cm/d")
plt.xlabel('z')
plt.ylabel('c')
plt.legend()

plt.show()
