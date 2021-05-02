from fipy.meshes import mesh1D
from fipy.tools.numerix import MA
import numpy as np
from fipy.tools import numerix

class Mesh1Dmod(mesh1D.Mesh1D):
    """ adapted mesh1D to axcept 3D network
    """
    def __init__(self, radiiVertices, radiiCells, length, vertexCoords, faceVertexIDs, cellFaceIDs, facesNorm = []):
        self.facesNorm = facesNorm
        self.radiiVertices = radiiVertices
        self.radiiCells = radiiCells
        self.length = length
        #self.tryval = 2
        self.valArea = 1
        self.volChoice = 1
        super(Mesh1Dmod, self).__init__(vertexCoords=vertexCoords, faceVertexIDs=faceVertexIDs, cellFaceIDs=cellFaceIDs)
        
    
    def _calcScaleArea(self): #adapt? what is this for?
        return 1.

    def _calcScaleVolume(self): #adapt? what is this for?
        return self.scale['length']

    def _calcFaceAreas(self): #able to handle when segment face vary?
        areas = np.array([(xi**2) * np.pi for xi in self.radiiVertices])
        return areas * (self.valArea == 1) + (1. * np.ones(self.numberOfFaces, 'd')  ) * (self.valArea == 2)

    def _calcFaceNormals(self): #used when?
        #norm_x_coord = np.array([self.facesNorm[xi][0] for xi in range(len(self.faceVertexIDs[0]))], dtype=np.float64) 
        #norm_y_coord = np.array([self.facesNorm[xi][1] for xi in range(len(self.faceVertexIDs[0]))], dtype=np.float64) 
        #norm_z_coord = np.array([self.facesNorm[xi][2] for xi in range(len(self.faceVertexIDs[0]))], dtype=np.float64) 
        #norm_coord =  -np.vstack((norm_x_coord,norm_y_coord,norm_z_coord)) #change how coords are stored
        n = self.numberOfFaces         
        faceNormals = -np.vstack((np.zeros((n,)), np.zeros((n,)), np.ones((n,)))) 
        #The left-most face has neighboring cells None and the left-most cell.
        # We must reverse the normal to make fluxes work correctly.
#         if self.numberOfFaces > 0:
#             faceNormals[:, 0] = -faceNormals[:, 0]
            
        return (faceNormals)# *(self.tryval == 2) #+ (norm_coord)*(self.tryval == 1) 



    def _calcCellVolumes(self):#leads to func not cunverging. no effect on mesh.cellVolumes
        volumes = np.array([(np.pi * self.radiiCells[xi]**2) *self.length[xi] for xi in range(self.numberOfCells) ])
        return self.length* (self.volChoice == 2) +volumes * (self.volChoice == 1) + (self.volChoice == 0) * numerix.ones(self.numberOfCells, 'd')#np.array((np.ones(self.numberOfCells, 'd'),))  

    def _calcFaceTangents(self): #adapt? what is this for?
        n = self.numberOfFaces 
        faceTangents2 = np.vstack((np.ones((n,)), np.zeros((n,)), np.zeros((n,)))) #vstack: add a row to matrix
        faceTangents1 = np.vstack((np.zeros((n,)), np.ones((n,)), np.zeros((n,))))
        return faceTangents1, faceTangents2
        
    def _calcCellCenters(self): #otherwise, wrong cellCenter for cell with 3 faces
        tmp = numerix.take(self._faceCenters, self.cellFaceIDs[:-1], axis=1)
        return MA.filled(MA.average(tmp, 1))  
        
    @property
    def faces(self):
        return self.interiorFaces|self.exteriorFaces
        
    @property
    def _cellDistances(self): #default function gave wrong results for some reason
        distance = np.take(self.length/2,self._adjacentCellIDs[0] ) \
            + self.interiorFaces * np.take(self.length/2,self._adjacentCellIDs[1] )
        return distance
        
        
        
        