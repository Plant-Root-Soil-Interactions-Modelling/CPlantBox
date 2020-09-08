#
# Example Exudation
#
# 1) Opens plant and root parameters from a file
# 2) Simulates root growth
# 3) Outputs a VTP (for vizualisation in ParaView)
#
#  Computes analytical solution of moving point/line sources based on Carslaw and Jaeger
#
import sys; sys.path.append("../..")

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

import plantbox as rb

rootsystem = rb.RootSystem()
name = "anagallis_straight"

#
# Open plant and root parameter from a file
#
rootsystem.openFile(name)
# rootsystem.writeParameters() # not exposed to python yet

#
# Initialize
#
rootsystem.initialize()

#
# Simulate
#
simtime = 1  # or 20, 40, 60 days
dt = 1  # try other values here
N = round(simtime / dt)  # steps

for i in range(0, int(N)):
    rootsystem.simulate(dt);

tips = rootsystem.getRootTips()  # 3D vector containing the coordinates of all the nodes
bases = rootsystem.getRootBases()
ctimes = [0, 0, 0, 0, 0, 0]  # TODO: replace with real calculation
#
# Export final result (as vtp)
#
rootsystem.write(name + ".vtp", rb.OutputType.segments)  # use ot_polylines for nicer visualization, ot_segments for animations

#
# Export segments for Matlab analysis
#
# analysis = rb.AnalysisSDF(rootsystem)
# analysis.write(name+".txt")

#
# Total length and surface
#
# l = analysis.getSummed(rb.ScalarType.length)
# print("Total root system length is "+str(l)+" cm")

# print("Finished with a total of "+str(rootsystem.getNumberOfNodes())+ " nodes")

# Parameter definition for exudation
simtime = simtime * 3600 * 24
ctimes = ctimes * 3600 * 24
M = 1e-5;
Dt = 1e-5;  # cm2/s
Dl = Dt;
theta = 0.3;
R = 1;
lambda_ = 1e-6;
L = 1;

x = np.linspace(-2, 2, 10)
y = np.linspace(-2, 2, 10)
z = np.linspace(-4, 0, 10);
X, Y, Z = np.meshgrid(x, y, z)
output_times = [1.0]

tend = 1.0;
Csum = np.zeros((len(x), len(y), len(z)))

# loop over all the tips
for i in range(len(tips)):
    base = bases[i]
    tip = tips[i]  # plume position
    ctime = ctimes[i]
    age_r = tend - ctime
    C = np.zeros((len(x), len(y), len(z)))

    if age_r > 0:

        vx = -(tip.x - base.x) / age_r
        vy = -(tip.y - base.y) / age_r
        vz = -(tip.z - base.z) / age_r
        vl = np.linalg.norm([vx, vy, vz])
        vxn = vx / vl; vyn = vy / vl; vzn = vz / vl
        L = min(L, vl)  # if root is not yet as long as L
        tipx = tip.x; tipy = tip.y; tipz = tip.z;

        def funC(t, l, theta, Dt, Dl, age_r, tipx, tipy, tipz, vx, vxn, vy, vyn, vz, vzn, R, xa, yb, zd):
            exp_x = -((xa - tipx - vxn * l) * R - vx * (age_r - t)) ** 2 / (4 * R * Dl * (age_r - t))
            exp_y = -((yb - tipy - vyn * l) * R - vy * (age_r - t)) ** 2 / (4 * R * Dl * (age_r - t))
            exp_z = -((zd - tipz - vzn * l) * R - vz * (age_r - t)) ** 2 / (4 * R * Dl * (age_r - t))
            C_ = M / (8 * theta * np.sqrt(math.pi ** 3 * Dt ** 2 * Dl * (age_r - t) ** 3)) * np.exp(exp_x + exp_y + exp_z - lambda_ * (age_r - t) / R);
            return C_

        # loop over all points of the domain
        for a in range(len(x)):
            for b in range(len(y)):
                for d in range(len(z)):
                    xa = x[a]; yb = y[b]; zd = z[d];
                    params = (theta, Dt, Dl, age_r, tipx, tipy, tipz, vx, vxn, vy, vyn, vz, vzn, R, xa, yb, zd)
                    I = integrate.dblquad(funC, 0.0, L, lambda l: 0.0, lambda l: age_r, args = params)
                    C[a, b, d] = I[0]
                    print(I[0])

    # sum up concentrations due to all tips : principle of superposition
    Csum = Csum + C;

    gridToVTK("./Exudates", X, Y, Z, pointData = {"Exudates":Csum})

