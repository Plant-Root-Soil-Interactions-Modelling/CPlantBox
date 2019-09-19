import unittest
import py_rootbox as rb
from rsml import *
from rb_tools import *

class TestRootSystem(unittest.TestCase):
  
    def test_rsml(self):
        """ checks rsml functionality with Python rsml reader """
        name = "Anagallis_femina_Leitner_2010"
        rs = rb.RootSystem()
        rs.readParameters("modelparameter2/" + name + ".xml")
        rs.initialize()
        simtime = 60
        rs.simulate(simtime)
        rs.writeRSML(name + ".rsml")


if __name__ == '__main__':
#     test = TestRootSystem()
#     test.test_copy()
#     print("fin.")
    unittest.main()
