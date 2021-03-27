from fipy.meshes import mesh1D
from fipy.tools.numerix import MA
import numpy as np

class Mesh1Dmod(mesh1D.Mesh1D):
    """ adapted mesh1D to axcept 3D network
    """
    
    def _calcScaleArea(self):
        return 1.

    def _calcScaleVolume(self):
        return self.scale['length']

    def _calcFaceAreas(self):
        return 1. * np.ones(self.numberOfFaces, 'd')  

    def _calcFaceNormals(self):
        # faceNormals = np.array((np.ones(self.numberOfFaces, 'd'),))        
        n = self.numberOfFaces         
        faceNormals = -np.vstack((np.zeros((n,)), np.zeros((n,)), np.ones((n,)))) 
        # The left-most face has neighboring cells None and the left-most cell.
        # We must reverse the normal to make fluxes work correctly.
#         if self.numberOfFaces > 0:
#             faceNormals[:, 0] = -faceNormals[:, 0]
        return faceNormals

    def _calcCellVolumes(self):
        return np.array((np.ones(self.numberOfCells, 'd'),))  

    def _calcFaceTangents(self):
        n = self.numberOfFaces 
        faceTangents2 = np.vstack((np.ones((n,)), np.zeros((n,)), np.zeros((n,)))) #vstack: add a row to matrix
        faceTangents1 = np.vstack((np.zeros((n,)), np.ones((n,)), np.zeros((n,))))
        return faceTangents1, faceTangents2
        