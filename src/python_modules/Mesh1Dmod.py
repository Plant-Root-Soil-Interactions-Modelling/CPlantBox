from fipy.meshes import mesh1D
from fipy.tools.numerix import MA
import numpy as np
from fipy.tools import numerix

class Mesh1Dmod(mesh1D.Mesh1D):
    """ adapted mesh1D to axcept 3D network
    """
    def __init__(self, facesNorm, vertexCoords, faceVertexIDs, cellFaceIDs):
        self.facesNorm = facesNorm
        super(Mesh1Dmod, self).__init__(vertexCoords=vertexCoords, faceVertexIDs=faceVertexIDs, cellFaceIDs=cellFaceIDs)
    
    def _calcScaleArea(self): #adapt? what is this for?
        return 1.

    def _calcScaleVolume(self): #adapt? what is this for?
        return self.scale['length']

    def _calcFaceAreas(self): #change to make dependent on segment length and radius
        return 1. * np.ones(self.numberOfFaces, 'd')  

    def _calcFaceNormals(self):
        norm_x_coord = np.array([self.facesNorm[xi][0] for xi in range(len(self.faceVertexIDs[0]))], dtype=np.float64) 
        norm_y_coord = np.array([self.facesNorm[xi][1] for xi in range(len(self.faceVertexIDs[0]))], dtype=np.float64) 
        norm_z_coord = np.array([self.facesNorm[xi][2] for xi in range(len(self.faceVertexIDs[0]))], dtype=np.float64) 
        norm_coord =  np.vstack((norm_x_coord,norm_y_coord,norm_z_coord)) #change how coords are stored 
        #print(norm_coord)
        return norm_coord

    def _calcCellVolumes(self):#changed to avoid error when computing source term. ATT! => change to make dependent on segment length and radius
        return numerix.ones(self.numberOfCells, 'd')#np.array((np.ones(self.numberOfCells, 'd'),))  

    def _calcFaceTangents(self): #adapt? what is this for?
        n = self.numberOfFaces 
        faceTangents2 = np.vstack((np.ones((n,)), np.zeros((n,)), np.zeros((n,)))) #vstack: add a row to matrix
        faceTangents1 = np.vstack((np.zeros((n,)), np.ones((n,)), np.zeros((n,))))
        return faceTangents1, faceTangents2
        
    @property
    def faces(self):
        return self.interiorFaces|self.exteriorFaces
        