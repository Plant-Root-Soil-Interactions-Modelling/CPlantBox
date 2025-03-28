import sys; sys.path.append(".."); sys.path.append("../src")

import unittestimport plantbox as pb

from functional.PlantHydraulicModel import HydraulicModel_Doussan
from functional.PlantHydraulicModel import HydraulicModel_Meunier
from functional.PlantHydraulicParameters import PlantHydraulicParameters
import visualisation.vtk_plot as vp

import numpy as np
import matplotlib.pyplot as plt

"""
    test are based on the Benchmark example M31 (single root, static soil, comparison to analytic solution)
    
    * tests finished: neumann, dirichlet (for single root)
    * todo: root system: M32 compare to Meunier(?)
    * todo: suf, krs, cached solution (including timings?)
    
"""


class TestPlantHydraulicModel(unittest.TestCase):

    def setup_analytic(self):
        """ Parameters """
        g = 9.8065 * 100.*24.*3600.*24.*3600.  # gravitational acceleration [cm day-2]
        rho = 1.  # density of water, [g/cm^3]
        L = 50  # length of single straight root [cm]
        self.L = L
        a = 0.2  # radius [cm] <---- ??? rather large
        self.a = a
        kz0 = 4.32e-2  # [cm^3/day]
        self.kx0 = kz0
        kz = kz0 / (rho * g)  # axial conductivity [cm^5 s / g]
        kr0 = 1.728e-4  # [1/day]
        self.kr0 = kr0
        kr = kr0 / (rho * g)  # radial conductivity per root type [cm^2 s / g]
        self.soil_matricpotential = -200  # static soil pressure [cm]
        self.collar_potential = -1000  # dircichlet bc at top
        """ Analytical solution """
        c = 2 * a * np.pi * kr / kz
        p_r = lambda z: self.soil_matricpotential + d[0] * np.exp(np.sqrt(c) * z) + d[1] * np.exp(-np.sqrt(c) * z)  #
        AA = np.array([[1, 1], [np.sqrt(c) * np.exp(-np.sqrt(c) * L), -np.sqrt(c) * np.exp(np.sqrt(c) * L)] ])  # # Boundary conditions dirichlet top, neumann bot
        bb = np.array([self.collar_potential - self.soil_matricpotential, -1])  # -rho * g
        d = np.linalg.solve(AA, bb)  # compute constants d_1 and d_2 from bc
        N = 100
        self.z_a = np.linspace(0, -L, N)  # Evaluate function
        self.rx_a = list(map(p_r, self.z_a))

    def setup_numeric(self):
        """ Numeric solution (run setup_analytic before)"""
        N = self.z_a.shape[0]
        nodes, segs, radii = [], [], []
        for z in self.z_a:
            nodes.append(pb.Vector3d(0, 0, z))
        for s in range(0, N - 1):
            segs.append(pb.Vector2i(s, s + 1))
            radii.append(self.a)
        self.rs = pb.MappedSegments(nodes, segs, radii)
        soil_index = lambda x, y, z: 0
        self.rs.setSoilGrid(soil_index)
        self.params = PlantHydraulicParameters(self.rs)
        self.params.setKrConst(self.kr0,0)
        self.params.setKxConst(self.kx0,0)

    def equilibrium_(self, solver):
        t = 0.
        N = self.z_a.shape[0]
        total_potential = np.ones((N - 1,)) * -500
        matric_potential = self.rs.total2matric(total_potential)
        # matric_potential = total_potential - 0.5 * (self.z_a[1:] + self.z_a[:-1])
        # Dirichlet
        rx_d = solver.solve_dirichlet(t, -500, matric_potential, cells = False)
        print("test_equilibrium(), dirichlet, rmse", np.linalg.norm(rx_d[1:] - matric_potential) / np.sqrt(N))
        trans = solver.get_transpiration(t, rx_d, matric_potential, cells = False)
        print("test_equilibrium(), dirichlet, transpiration", trans, "cm3/day")
        # Neumann
        rx_d = solver.solve_neumann(t, 0., matric_potential, cells = False)
        print("test_equilibrium(), neumann, rmse", np.linalg.norm(rx_d[1:] - matric_potential) / np.sqrt(N))
        trans = solver.get_transpiration(t, rx_d, matric_potential, cells = False)
        print("test_equilibrium(), neumann, transpiration", trans, "cm3/day")
        return None

    def test_equilibrium(self):
        """ For constant total potential and zero uptake, all fluxes should vanish. 
            Depending on the method there might be some numerical noise
        """
        self.setup_analytic()
        self.setup_numeric()
        print("\nDoussan model")
        self.equilibrium_(HydraulicModel_Doussan(self.rs, self.params, cached = False))

        print("\nMeunier model")
        self.equilibrium_(HydraulicModel_Meunier(self.rs, self.params, cached = False))

        print("")
        return None

    def test_dirichlet(self):
        """ Compare to analytic solution for a single root, and plausibility for static root system"""
        self.setup_analytic()
        self.setup_numeric()
        t = 0.  # simulation time (can be neglected since kr and kx are constant)
        N = self.z_a.shape[0]

        solver = HydraulicModel_Meunier(self.rs, self. params, cached = False)
        rx = solver.solve_dirichlet(t, self.collar_potential, [self.soil_matricpotential], cells = True)
        error = np.linalg.norm(rx - self.rx_a) / np.sqrt(N)
        print("test_dirichlet(), Meunier, rmse ", error)
        trans = solver.get_transpiration(t, rx, [self.soil_matricpotential], cells = True)
        print("test_dirichlet(), Meunier, transpiration", trans, "cm3/day")
        self.assertAlmostEqual(rx[0], self.collar_potential)

        solver = HydraulicModel_Doussan(self.rs, self.params, cached = False)  # or HydraulicModel_Doussan, HydraulicModel_Meunier
        rx = solver.solve_dirichlet(t, self.collar_potential, [self.soil_matricpotential], cells = True)
        error = np.linalg.norm(rx - self.rx_a) / np.sqrt(N)
        print("test_dirichlet(), Doussan, rmse ", error)
        trans = solver.get_transpiration(t, rx, [self.soil_matricpotential], cells = True)
        print("test_dirichlet(), Doussan, transpiration", trans, "cm3/day")
        self.assertAlmostEqual(rx[0], self.collar_potential)

        print("")
        return None

    def test_neumann(self):
        """ Compare to analytic solution for a single root, and plausibility for static root system"""
        self.setup_analytic()
        self.setup_numeric()
        t = 0.  # simulation time (can be neglected since kr and kx are constant)
        N = self.z_a.shape[0]
        trans_ = -2.4054512  # should result again in a collar potential of -1000

        solver = HydraulicModel_Meunier(self.rs, self. params, cached = False)
        rx = solver.solve_neumann(t, trans_, [self.soil_matricpotential], cells = True)
        trans = solver.get_transpiration(t, rx, [self.soil_matricpotential], cells = True)
        print("test_neumann(), Meunier, transpiration", trans, "cm3/day")
        print("test_neumann(), Meunier, collar potential", rx[0], "cm3/day")
        self.assertAlmostEqual(trans_, trans)

        solver = HydraulicModel_Doussan(self.rs, self.params, cached = False)  # or HydraulicModel_Doussan, HydraulicModel_Meunier
        rx = solver.solve_neumann(t, trans_, [self.soil_matricpotential], cells = True)
        trans = solver.get_transpiration(t, rx, [self.soil_matricpotential], cells = True)
        print("test_neumann(), Doussan, transpiration", trans, "cm3/day")
        print("test_neumann(), Doussan, collar potential", rx[0], "cm3/day")
        self.assertAlmostEqual(trans_, trans)

        print("")
        return None

    def test_consistancy(self):
        """ Test if switching between Dirichlet and Neumann yield same results """
        self.setup_analytic()
        self.setup_numeric()
        t = 0.  # simulation time (can be neglected since kr and kx are constant)

        solver = HydraulicModel_Meunier(self.rs, self.params, cached = False)  # or HydraulicModel_Doussan, HydraulicModel_Meunier
        collar_pot = -1000
        rx = solver.solve_dirichlet(t, collar_pot, [self.soil_matricpotential], cells = True)
        trans = solver.get_transpiration(t, rx, [self.soil_matricpotential], cells = True)
        print("test_consistancy(), Meunier, transpiration", trans, "cm3/day", collar_pot, "cm collar potential")
        rx = solver.solve_neumann(t, trans, [self.soil_matricpotential], cells = True)
        trans = solver.get_transpiration(t, rx, [self.soil_matricpotential], cells = True)
        print("test_consistancy(), Meunier, transpiration", trans, "cm3/day", rx[0], "cm collar potential")
        self.assertAlmostEqual(rx[0], collar_pot)

        solver = HydraulicModel_Doussan(self.rs, self.params, cached = False)  # or HydraulicModel_Doussan, HydraulicModel_Meunier
        collar_pot = -1000
        rx = solver.solve_dirichlet(t, collar_pot, [self.soil_matricpotential], cells = True)
        trans = solver.get_transpiration(t, rx, [self.soil_matricpotential], cells = True)
        print("test_consistancy(), Doussan, transpiration", trans, "cm3/day", collar_pot, "cm collar potential")
        rx = solver.solve_neumann(t, trans, [self.soil_matricpotential], cells = True)
        trans = solver.get_transpiration(t, rx, [self.soil_matricpotential], cells = True)
        print("test_consistancy(), Doussan, transpiration", trans, "cm3/day", rx[0], "cm collar potential")
        self.assertAlmostEqual(rx[0], collar_pot)

        print("")
        return None

    def test_krs(self):
        """ """
        return None

    def test_suf(self):
        """ """
        return None


if __name__ == '__main__':
    unittest.main()
    print("done.")
