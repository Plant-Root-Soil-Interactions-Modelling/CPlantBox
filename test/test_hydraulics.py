import sys; sys.path.append("..")

import unittestimport plantbox as pb


class TestPlantHydraulicModel(unittest.TestCase):

    def set_up_Meunier_(self):
        return None

    def set_up_Doussan_(self):
        return None

    def equilibrium_(self):
        """ for constant total potential and zero uptake fluxes should vanish. 
            Depending on the method there might be some numerical noise
        """
        return None

    def dirichlet_(self):
        """ Compare to analytic solution for a single root, and plausibility for static root system"""

        return None

    def neumann_(self):
        """ Compare to analytic solution for a single root, and plausibility for static root system"""
        return None

    def consistancy_(self):
        """ Test if switching between Dirichlet and Neumann yield same results """
        return None

    def switch_(self):
        """ Test if switching between Dirichlet and Neumann yield same results """
        return None

    def test_equilibrium_Meunier(self):
        return None

    def test_equilibrium_Doussan(self):
        return None

    def test_switch(self):
        return None


if __name__ == '__main__':
    unittest.main()
    print("done.")
