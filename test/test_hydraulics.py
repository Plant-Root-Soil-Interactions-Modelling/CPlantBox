import sys

sys.path.append("..")
sys.path.append("../src")

import os
import unittest

import matplotlib.pyplot as plt
import numpy as np

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from plantbox.functional.PlantHydraulicModel import HydraulicModel_Doussan, HydraulicModel_Meunier
from plantbox.functional.PlantHydraulicParameters import PlantHydraulicParameters

"""
    test are based on the Benchmark example M31 (single root, static soil, comparison to analytic solution)
    and M32 (root system, static soil, comparison between Doussan and Meunier solvers)
    
    * tests finished: neumann, dirichlet (for single root), m32 (root system, Doussan vs. Meunier)    
"""

_HERE = os.path.dirname(os.path.abspath(__file__))
M32_RSML = os.path.join(_HERE, "..", "tutorial", "examples", "root_hydraulics", "RootSystem.rsml")
M32_REFERENCE = os.path.join(_HERE, "m32_reference.npz")  # baseline results, regenerate if a model is intentionally changed


class TestPlantHydraulicModel(unittest.TestCase):

    def setup_analytic(self):
        """Parameters"""
        g = 9.8065 * 100.0 * 24.0 * 3600.0 * 24.0 * 3600.0  # gravitational acceleration [cm day-2]
        rho = 1.0  # density of water, [g/cm^3]
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
        AA = np.array([[1, 1], [np.sqrt(c) * np.exp(-np.sqrt(c) * L), -np.sqrt(c) * np.exp(np.sqrt(c) * L)]])  # # Boundary conditions dirichlet top, neumann bot
        bb = np.array([self.collar_potential - self.soil_matricpotential, -1])  # -rho * g
        d = np.linalg.solve(AA, bb)  # compute constants d_1 and d_2 from bc
        N = 100
        self.z_a = np.linspace(0, -L, N)  # Evaluate function
        self.rx_a = list(map(p_r, self.z_a))

    def setup_numeric(self):
        """Numeric solution (run setup_analytic before)"""
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
        self.params.set_kr_const(self.kr0)
        self.params.set_kx_const(self.kx0)

    def equilibrium_(self, solver):
        t = 0.0
        N = self.z_a.shape[0]
        total_potential = np.ones((N - 1,)) * -500
        matric_potential = self.rs.total2matric(total_potential)
        # matric_potential = total_potential - 0.5 * (self.z_a[1:] + self.z_a[:-1])
        # Dirichlet
        rx_d = solver.solve_dirichlet(t, -500, matric_potential, cells=False)
        print("test_equilibrium(), dirichlet, rmse", np.linalg.norm(rx_d[1:] - matric_potential) / np.sqrt(N))
        trans = solver.get_transpiration(t, rx_d, matric_potential, cells=False)
        print("test_equilibrium(), dirichlet, transpiration", trans, "cm3/day")
        # Neumann
        rx_d = solver.solve_neumann(t, 0.0, matric_potential, cells=False)
        print("test_equilibrium(), neumann, rmse", np.linalg.norm(rx_d[1:] - matric_potential) / np.sqrt(N))
        trans = solver.get_transpiration(t, rx_d, matric_potential, cells=False)
        print("test_equilibrium(), neumann, transpiration", trans, "cm3/day")
        return None

    def test_equilibrium(self):
        """For constant total potential and zero uptake, all fluxes should vanish.
        Depending on the method there might be some numerical noise
        """
        self.setup_analytic()
        self.setup_numeric()
        print("\nDoussan model")
        self.equilibrium_(HydraulicModel_Doussan(self.rs, self.params, cached=False))

        print("\nMeunier model")
        self.equilibrium_(HydraulicModel_Meunier(self.rs, self.params, cached=False))

        print("")
        return None

    def test_dirichlet(self):
        """Compare to analytic solution for a single root, and plausibility for static root system"""
        self.setup_analytic()
        self.setup_numeric()
        t = 0.0  # simulation time (can be neglected since kr and kx are constant)
        N = self.z_a.shape[0]

        solver = HydraulicModel_Meunier(self.rs, self.params, cached=False)
        rx = solver.solve_dirichlet(t, self.collar_potential, [self.soil_matricpotential], cells=True)
        error = np.linalg.norm(rx - self.rx_a) / np.sqrt(N)
        print("test_dirichlet(), Meunier, rmse ", error)
        trans = solver.get_transpiration(t, rx, [self.soil_matricpotential], cells=True)
        print("test_dirichlet(), Meunier, transpiration", trans, "cm3/day")
        self.assertAlmostEqual(rx[0], self.collar_potential)

        solver = HydraulicModel_Doussan(self.rs, self.params, cached=False)  # or HydraulicModel_Doussan, HydraulicModel_Meunier
        rx = solver.solve_dirichlet(t, self.collar_potential, [self.soil_matricpotential], cells=True)
        error = np.linalg.norm(rx - self.rx_a) / np.sqrt(N)
        print("test_dirichlet(), Doussan, rmse ", error)
        trans = solver.get_transpiration(t, rx, [self.soil_matricpotential], cells=True)
        print("test_dirichlet(), Doussan, transpiration", trans, "cm3/day")
        self.assertAlmostEqual(rx[0], self.collar_potential)

        print("")
        return None

    def test_neumann(self):
        """Compare to analytic solution for a single root, and plausibility for static root system"""
        self.setup_analytic()
        self.setup_numeric()
        t = 0.0  # simulation time (can be neglected since kr and kx are constant)
        N = self.z_a.shape[0]
        trans_ = -2.4054512  # should result again in a collar potential of -1000

        solver = HydraulicModel_Meunier(self.rs, self.params, cached=False)
        rx = solver.solve_neumann(t, trans_, [self.soil_matricpotential], cells=True)
        trans = solver.get_transpiration(t, rx, [self.soil_matricpotential], cells=True)
        print("test_neumann(), Meunier, transpiration", trans, "cm3/day")
        print("test_neumann(), Meunier, collar potential", rx[0], "cm3/day")
        self.assertAlmostEqual(trans_, trans)

        solver = HydraulicModel_Doussan(self.rs, self.params, cached=False)  # or HydraulicModel_Doussan, HydraulicModel_Meunier
        rx = solver.solve_neumann(t, trans_, [self.soil_matricpotential], cells=True)
        trans = solver.get_transpiration(t, rx, [self.soil_matricpotential], cells=True)
        print("test_neumann(), Doussan, transpiration", trans, "cm3/day")
        print("test_neumann(), Doussan, collar potential", rx[0], "cm3/day")
        self.assertAlmostEqual(trans_, trans)

        print("")
        return None

    def test_consistancy(self):
        """Test if switching between Dirichlet and Neumann yield same results"""
        self.setup_analytic()
        self.setup_numeric()
        t = 0.0  # simulation time (can be neglected since kr and kx are constant)

        solver = HydraulicModel_Meunier(self.rs, self.params, cached=False)  # or HydraulicModel_Doussan, HydraulicModel_Meunier
        collar_pot = -1000
        rx = solver.solve_dirichlet(t, collar_pot, [self.soil_matricpotential], cells=True)
        trans = solver.get_transpiration(t, rx, [self.soil_matricpotential], cells=True)
        print("test_consistancy(), Meunier, transpiration", trans, "cm3/day", collar_pot, "cm collar potential")
        rx = solver.solve_neumann(t, trans, [self.soil_matricpotential], cells=True)
        trans = solver.get_transpiration(t, rx, [self.soil_matricpotential], cells=True)
        print("test_consistancy(), Meunier, transpiration", trans, "cm3/day", rx[0], "cm collar potential")
        self.assertAlmostEqual(rx[0], collar_pot)

        solver = HydraulicModel_Doussan(self.rs, self.params, cached=False)  # or HydraulicModel_Doussan, HydraulicModel_Meunier
        collar_pot = -1000
        rx = solver.solve_dirichlet(t, collar_pot, [self.soil_matricpotential], cells=True)
        trans = solver.get_transpiration(t, rx, [self.soil_matricpotential], cells=True)
        print("test_consistancy(), Doussan, transpiration", trans, "cm3/day", collar_pot, "cm collar potential")
        rx = solver.solve_neumann(t, trans, [self.soil_matricpotential], cells=True)
        trans = solver.get_transpiration(t, rx, [self.soil_matricpotential], cells=True)
        print("test_consistancy(), Doussan, transpiration", trans, "cm3/day", rx[0], "cm collar potential")
        self.assertAlmostEqual(rx[0], collar_pot)

        print("")
        return None

    def m32_(self, param, sim_time, rmse_tol, trans_tol, krs_tol, suf_tol, name):
        """Benchmark M3.2: root system, static soil, Dirichlet bc
        compares Doussan and Meunier solvers, and both against a saved reference solution
        """
        p_s = -200  # static soil pressure [cm]
        p0 = -500  # dirichlet bc at top
        soil_index = lambda x, y, z: 0

        solver_m = HydraulicModel_Meunier(M32_RSML, param, cached=False)
        solver_m.ms.setSoilGrid(soil_index)
        rx_m = solver_m.solve_dirichlet(sim_time, p0, [p_s], cells=True)
        trans_m = solver_m.get_transpiration(sim_time, rx_m, [p_s], cells=True)
        krs_m, _ = solver_m.get_krs(sim_time)
        suf_m = solver_m.get_suf(sim_time)

        solver_d = HydraulicModel_Doussan(M32_RSML, param, cached=False)
        solver_d.ms.setSoilGrid(soil_index)
        rx_d = solver_d.solve_dirichlet(sim_time, p0, [p_s], cells=True)
        trans_d = solver_d.get_transpiration(sim_time, rx_d, [p_s], cells=True)
        krs_d, _ = solver_d.get_krs(sim_time)
        suf_d = solver_d.get_suf(sim_time)

        rmse = np.linalg.norm(rx_d - rx_m) / np.sqrt(len(rx_d))
        suf_rmse = np.linalg.norm(suf_d - suf_m) / np.sqrt(len(suf_d))
        print("test_m32(), Doussan vs. Meunier, rmse", rmse)
        print("test_m32(), Doussan transpiration", trans_d, "cm3/day, Meunier transpiration", trans_m, "cm3/day")
        print("test_m32(), Doussan krs", krs_d, "cm2/day, Meunier krs", krs_m, "cm2/day")
        print("test_m32(), Doussan vs. Meunier, suf rmse", suf_rmse)

        self.assertAlmostEqual(rx_m[0], p0)
        self.assertAlmostEqual(rx_d[0], p0)
        self.assertLess(rmse, rmse_tol)
        self.assertAlmostEqual(trans_d, trans_m, delta=trans_tol)
        self.assertAlmostEqual(suf_m.sum(), 1.0, places=6)
        self.assertAlmostEqual(suf_d.sum(), 1.0, places=6)
        self.assertAlmostEqual(krs_d, krs_m, delta=krs_tol)
        self.assertLess(suf_rmse, suf_tol)

        # regression check against saved baseline (regenerate m32_reference.npz if a model is intentionally changed)
        ref = np.load(M32_REFERENCE)
        rmse_ref_m = np.linalg.norm(rx_m - ref[f"rx_m_{name}"]) / np.sqrt(len(rx_m))
        rmse_ref_d = np.linalg.norm(rx_d - ref[f"rx_d_{name}"]) / np.sqrt(len(rx_d))
        suf_rmse_ref_m = np.linalg.norm(suf_m - ref[f"suf_m_{name}"]) / np.sqrt(len(suf_m))
        suf_rmse_ref_d = np.linalg.norm(suf_d - ref[f"suf_d_{name}"]) / np.sqrt(len(suf_d))
        print("test_m32(), Meunier vs. reference, rmse", rmse_ref_m)
        print("test_m32(), Doussan vs. reference, rmse", rmse_ref_d)
        self.assertLess(rmse_ref_m, 1e-6)
        self.assertLess(rmse_ref_d, 1e-6)
        self.assertAlmostEqual(trans_m, float(ref[f"trans_m_{name}"]), places=6)
        self.assertAlmostEqual(trans_d, float(ref[f"trans_d_{name}"]), places=6)
        self.assertAlmostEqual(krs_m, float(ref[f"krs_m_{name}"]), places=6)
        self.assertAlmostEqual(krs_d, float(ref[f"krs_d_{name}"]), places=6)
        self.assertLess(suf_rmse_ref_m, 1e-6)
        self.assertLess(suf_rmse_ref_d, 1e-6)
        return None

    def test_m32a(self):
        """Benchmark M3.2a: root system with constant conductivities"""
        kz = 4.32e-2  # [cm^3/day]
        kr = 1.728e-4  # [1/day]
        param = PlantHydraulicParameters()
        param.set_kr_const(kr)
        param.set_kx_const(kz)
        self.m32_(param, sim_time=0.0, rmse_tol=0.5, trans_tol=0.05, krs_tol=1e-4, suf_tol=1e-4, name="a")
        print("")
        return None

    def test_m32b(self):
        """Benchmark M3.2b: root system with age dependent conductivities"""
        kr0 = np.array(
            [
                [0, 1.14e-03],
                [2, 1.09e-03],
                [4, 1.03e-03],
                [6, 9.83e-04],
                [8, 9.35e-04],
                [10, 8.90e-04],
                [12, 8.47e-04],
                [14, 8.06e-04],
                [16, 7.67e-04],
                [18, 7.30e-04],
                [20, 6.95e-04],
                [22, 6.62e-04],
                [24, 6.30e-04],
                [26, 5.99e-04],
                [28, 5.70e-04],
                [30, 5.43e-04],
                [32, 5.17e-04],
            ]
        )
        kr1 = np.array(
            [
                [0, 4.11e-03],
                [1, 3.89e-03],
                [2, 3.67e-03],
                [3, 3.47e-03],
                [4, 3.28e-03],
                [5, 3.10e-03],
                [6, 2.93e-03],
                [7, 2.77e-03],
                [8, 2.62e-03],
                [9, 2.48e-03],
                [10, 2.34e-03],
                [11, 2.21e-03],
                [12, 2.09e-03],
                [13, 1.98e-03],
                [14, 1.87e-03],
                [15, 1.77e-03],
                [16, 1.67e-03],
                [17, 1.58e-03],
            ]
        )
        kx0 = np.array(
            [
                [0, 6.74e-02],
                [2, 7.48e-02],
                [4, 8.30e-02],
                [6, 9.21e-02],
                [8, 1.02e-01],
                [10, 1.13e-01],
                [12, 1.26e-01],
                [14, 1.40e-01],
                [16, 1.55e-01],
                [18, 1.72e-01],
                [20, 1.91e-01],
                [22, 2.12e-01],
                [24, 2.35e-01],
                [26, 2.61e-01],
                [28, 2.90e-01],
                [30, 3.21e-01],
                [32, 3.57e-01],
            ]
        )
        kx1 = np.array(
            [
                [0, 4.07e-04],
                [1, 5.00e-04],
                [2, 6.15e-04],
                [3, 7.56e-04],
                [4, 9.30e-04],
                [5, 1.14e-03],
                [6, 1.41e-03],
                [7, 1.73e-03],
                [8, 2.12e-03],
                [9, 2.61e-03],
                [10, 3.21e-03],
                [11, 3.95e-03],
                [12, 4.86e-03],
                [13, 5.97e-03],
                [14, 7.34e-03],
                [15, 9.03e-03],
                [16, 1.11e-02],
                [17, 1.36e-02],
            ]
        )
        param = PlantHydraulicParameters()
        param.set_kr_age_dependent(kr0[:, 0], kr0[:, 1], subType=1)
        param.set_kr_age_dependent(kr1[:, 0], kr1[:, 1], subType=[2, 3])
        param.set_kx_age_dependent(kx0[:, 0], kx0[:, 1], subType=1)
        param.set_kx_age_dependent(kx1[:, 0], kx1[:, 1], subType=[2, 3])
        self.m32_(param, sim_time=14, rmse_tol=8.0, trans_tol=0.5, krs_tol=5e-4, suf_tol=5e-4, name="b")
        print("")
        return None


if __name__ == "__main__":
    unittest.main()
    print("done.")
