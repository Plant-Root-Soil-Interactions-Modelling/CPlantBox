"""
# @Testing RootSystem::simulate(double dt, double maxinc_, ProportionalElongation* se, bool verbose)
"""
import sys
import numpy as np
import unittest

sys.path.append('..')

import plantbox as pb

rootparfiles = ['../modelparameter/structural/rootsystem/Zea_mays_3_Postma_2011.xml']

class TestMaxInc(unittest.TestCase):
    def _short_sim(self, dt: int, soildepth: float, layers: float, rootparfile: str, maxinc: float, re_reduction: np.array):
        """ A short root growth simulation (1 step) using RootSystem::simulate(double dt, double maxinc_, ProportionalElongation* se, bool verbose) """
        #--- read some root pars
        rootparname = rootparfile

        #--- Set scale elongaton according to tutorial example "example5b_scaleelongation.py"
        scale_elongation = pb.EquidistantGrid1D(0, -soildepth, layers)
        scale_elongation.data = re_reduction
        se = pb.ProportionalElongation()
        se.setBaseLookUp(scale_elongation)

        #--- Update root parameters
        plant = pb.Plant()
        plant.setSeed(0)
        plant.readParameters(rootparname)
        for p in plant.getOrganRandomParameter(pb.OrganTypes.root):
            p.f_se = se  # set scale elongation function
            p.r = 5. # set a high initial growth
        
        #--- intialize root system
        plant.initialize()

        #--- set a low maxinc and simulate        
        ol = 0.
        plant.simulateLimited(dt, maxinc, "lengthTh", [1., 1., 1., 1., 1.], se, True)
        
        #--- get actual root increment
        l = np.sum(plant.getParameter('length'))
        self.inc =  l - ol

    def test_WithoutReduction(self):
        """ Tests RootSystem::simulate(double dt, double maxinc_, ProportionalElongation* se, bool verbose) without re_reduction"""
        tol = 0.1 # set tolerance as maxinc runs binary search [cm]
        maxinc = 4.
        layers=60
        for rootparfile in rootparfiles:
            self._short_sim(dt=1, soildepth=180, layers=layers, rootparfile=rootparfile, maxinc=maxinc, re_reduction=np.ones((layers-1,)))
            self.assertGreater(maxinc+tol, self.inc, 'maxinc not taking an effect')

    def test_WithReduction(self):
        """ Tests RootSystem::simulate(double dt, double maxinc_, ProportionalElongation* se, bool verbose) with re_reduction"""
        tol = 0.1 # set tolerance as maxinc runs binary search [cm]
        maxinc = 4.
        layers=60
        for rootparfile in rootparfiles:
            self._short_sim(dt=1, soildepth=180, layers=layers, rootparfile=rootparfile, maxinc=maxinc, re_reduction=np.ones((layers-1,))*0.9) # -10% reduction
            self.assertGreater(maxinc+tol, self.inc, 'maxinc not taking an effect')
 
if __name__ == '__main__':
    unittest.main(verbosity=2)