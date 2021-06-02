"""small example"""
import sys
sys.path.append("../../.."); sys.path.append("../../../src/python_modules"); sys.path.append("../../../test")
import plantbox as pb
import vtk_plot as vp
import test_root as tr
rs = tr.TestRoot()
rs.root_dxMin_test([1., 2., 10.])