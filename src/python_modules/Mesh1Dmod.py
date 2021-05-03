from fipy.meshes import mesh1D
from fipy.tools.numerix import MA
import numpy as np
from fipy.tools import numerix

class Mesh1Dmod(mesh1D.Mesh1D):
    """ adapted mesh1D to axcept 3D network
        already ok/apparently no need to adapt:
        cellFaceIDs
        cellToCellIDs
        _adjacentCellIDs
        faceCenters
        _cellToFaceOrientations
        _interiorCellIDs
        _exteriorCellIDs
        _interiorFaces
        _exteriorFaces
        _cellToFaceDistanceVectors : shape =( 3[x,y,z coord] , 2 [each side of face] , Num of faces )
        _internalCellDistances : ok
        _internalFaceToCellDistances
        cellDistanceVectors : ok
        cellToFaceDistanceVectors: ok
        _cellToFaceOrientations: faces not facing same way as would for fipy mesh but still gives ok results for faceGrad.divergence
        cellvar.grad: bcause tangents = 0, not used to compute faceVAr.grad. so ok if equation is wrong
    """
    def __init__(self, radiiVertices, radiiCells, length, vertexCoords, faceVertexIDs, cellFaceIDs, facesNorm = []):
        self.facesNorm = facesNorm
        self.radiiVertices = radiiVertices
        self.radiiCells = radiiCells
        self.length = length
        #self.tryval = 2
        self.valArea = 1
        self.volChoice = 1
        #self.faceContribution = self._calcFaceContribution()
        super(Mesh1Dmod, self).__init__(vertexCoords=vertexCoords, faceVertexIDs=faceVertexIDs, cellFaceIDs=cellFaceIDs)

    def _calcScaleArea(self): #adapt? what is this for?
        return 1.

    def _calcScaleVolume(self): #adapt? what is this for?
        return self.scale['length']

    def _calcFaceAreas(self): #able to handle when segment face vary?
        areas = np.array([(xi**2) * np.pi for xi in self.radiiVertices])
        return areas * (self.valArea == 1) + (1. * np.ones(self.numberOfFaces, 'd')  ) * (self.valArea == 2)

    def _calcFaceNormals(self): #used for computation gradient/divergence
    #need to set last or first normal to 1 instead of -1?
        n = self.numberOfFaces         
        faceNormals = -np.vstack((np.zeros((n,)), np.zeros((n,)), np.ones((n,)))) 
        #The left-most face has neighboring cells None and the left-most cell.
        # We must reverse the normal to make fluxes work correctly.
#         if self.numberOfFaces > 0:
#             faceNormals[:, 0] = -faceNormals[:, 0]
            
        return (faceNormals)# *(self.tryval == 2) #

    def _calcCellVolumes(self):
        volumes = np.array([(np.pi * self.radiiCells[xi]**2) *self.length[xi] for xi in range(self.numberOfCells) ])
        return self.length* (self.volChoice == 2) +volumes * (self.volChoice == 1) + (self.volChoice == 0) * numerix.ones(self.numberOfCells, 'd')#np.array((np.ones(self.numberOfCells, 'd'),))  

        
        
    def _calcFaceTangents(self): #adapt? what is this for?
        n = self.numberOfFaces 
        #faceTangents2 = np.vstack((np.ones((n,)), np.zeros((n,)), np.zeros((n,)))) #vstack: add a row to matrix
        #faceTangents1 = np.vstack((np.zeros((n,)), np.ones((n,)), np.zeros((n,))))
        faceTangents2 = np.vstack((np.zeros((n,)), np.zeros((n,)), np.zeros((n,)))) #vstack: add a row to matrix
        faceTangents1 = np.vstack((np.zeros((n,)), np.zeros((n,)), np.zeros((n,))))
        return faceTangents1, faceTangents2
        #set every thing to 0 instead, like for 1D grid of fipy
        
    def _calcCellCenters(self):
        tmp = numerix.take(self._faceCenters, self.cellFaceIDs[:2,:], axis=1)
        #print('\n\ncss ',self._faceCenters, self.cellFaceIDs, self.cellFaceIDs[:-1])
        return MA.filled(MA.average(tmp, 1))  
        
    def _calcCellToCellDist(self): #round otherwise adds a slight error for some reason
        return np.around(numerix.take(self._cellDistances, self.cellFaceIDs), 5)
        
    
    def _calcCellDistAndVec(self):
        tmp = numerix.take(self._cellCenters, self.faceCellIDs, axis=1)
        tmp = tmp[..., 1,:] - tmp[..., 0,:]
        tmp = MA.filled(MA.where(MA.getmaskarray(tmp), self._cellToFaceDistanceVectors[:, 0], tmp))
        cellDistanceVectors = tmp
        #cellDistances = MA.filled(MA.sqrt(MA.sum(tmp * tmp, 0)))
        return self._cellDistances, cellDistanceVectors #check celldistance vectors
        
    def _calcFaceToCellDistAndVec(self):
        tmp = MA.repeat(self._faceCenters[..., numerix.NewAxis,:], 2, 1)
        # array -= masked_array screws up masking for on numpy 1.1

        tmp = tmp - numerix.take(self._cellCenters, self.faceCellIDs, axis=1)
        cellToFaceDistanceVectors = tmp
        faceToCellDistances = MA.sqrt(MA.sum(tmp * tmp, 0))
        return np.around(faceToCellDistances, 5), cellToFaceDistanceVectors
        
    @property
    def faces(self):
        return self.interiorFaces|self.exteriorFaces
        
    @property
    def _cellDistances(self): #default function gave wrong results for some reason
        distance = np.take(self.length/2,self._adjacentCellIDs[0] ) \
            + self.interiorFaces * np.take(self.length/2,self._adjacentCellIDs[1] )
        return distance
        
        
        
        