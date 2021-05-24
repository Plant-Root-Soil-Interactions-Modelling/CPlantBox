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
        super(Mesh1Dmod, self).__init__(vertexCoords=vertexCoords, faceVertexIDs=faceVertexIDs, cellFaceIDs=cellFaceIDs)

    def _calcScaleArea(self): #adapt? what is this for?
        return 1.

    def _calcScaleVolume(self): #adapt? what is this for?
        return self.scale['length']

    def _calcFaceAreas(self): #able to handle when segment face vary?
        areas = np.array([(xi**2) * np.pi for xi in self.radiiVertices])
        return areas 

    def _calcFaceNormals(self): #used for computation gradient/divergence
        n = self.numberOfFaces         
        faceNormals = -np.vstack((np.zeros((n,)), np.zeros((n,)), np.ones((n,)))) 
            
        return (faceNormals)

    def _calcCellVolumes(self):#fipy just return length
        volumes = np.array([(np.pi * self.radiiCells[xi]**2) *self.length[xi] for xi in range(self.numberOfCells) ])
        return volumes 
     
        
    def _calcFaceTangents(self): 
        n = self.numberOfFaces 
        faceTangents2 = np.vstack((np.zeros((n,)), np.zeros((n,)), np.zeros((n,)))) #vstack: add a row to matrix
        faceTangents1 = np.vstack((np.zeros((n,)), np.zeros((n,)), np.zeros((n,))))
        return faceTangents1, faceTangents2
        #set every thing to 0, like for 1D grid of fipy
        
    def _calcCellCenters(self):#cell center from bot1 and top1 (ignore coordinates of eventual bot2 and top2)
        tmp = numerix.take(self._faceCenters, self.cellFaceIDs[:2,:], axis=1)
        return MA.filled(MA.average(tmp, 1))  
        
    def _calcCellToCellDist(self): #round otherwise adds a slight error for some reason
        return np.around(numerix.take(self._cellDistances, self.cellFaceIDs), 5)
    
    def _calcCellDistAndVec(self):
        tmp = numerix.take(self._cellCenters, self.faceCellIDs, axis=1)
        tmp = tmp[..., 1,:] - tmp[..., 0,:]
        tmp = MA.filled(MA.where(MA.getmaskarray(tmp), self._cellToFaceDistanceVectors[:, 0], tmp))
        cellDistanceVectors = tmp
        return self._cellDistances, cellDistanceVectors 
        
    def _calcFaceToCellDistAndVec(self):
        tmp = MA.repeat(self._faceCenters[..., numerix.NewAxis,:], 2, 1)
        tmp = tmp - numerix.take(self._cellCenters, self.faceCellIDs, axis=1)
        cellToFaceDistanceVectors = tmp
        faceToCellDistances = MA.sqrt(MA.sum(tmp * tmp, 0))
        return np.around(faceToCellDistances, 5), cellToFaceDistanceVectors
        
    @property
    def _cellDistances(self): #default function gave wrong results for some reason
        #for internal faces: distance between cell center on each side
        #exernal faces : distance between face and internal cell center
        distance = np.take(self.length/2,self._adjacentCellIDs[0] ) \
            + self.interiorFaces * np.take(self.length/2,self._adjacentCellIDs[1] )
        return distance
        
        
        
        