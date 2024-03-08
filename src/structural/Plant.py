import sys; sys.path.append(".."); sys.path.append("../..");

import plantbox as pb
from plantbox import Plant

import numpy as np


class PlantPython(Plant):

    def initialize_static(self, rsml_file, initialSubType):

        # i guess
        # 1. open the file ...
        # 2. parse topology (create mean radius)
        # 3. add StaticRoot to self

        pass
