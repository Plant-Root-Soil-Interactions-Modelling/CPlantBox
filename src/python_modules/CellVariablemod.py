
from fipy import CellVariable, FaceVariable, Variable
from fipy.variables.addOverFacesVariable import _AddOverFacesVariable
from fipy.variables.arithmeticCellToFaceVariable import _ArithmeticCellToFaceVariable
from fipy.variables.faceGradVariable import _FaceGradVariable
from fipy.tools.numerix import MA
import numpy as np
from fipy.tools import numerix
from fipy.tools import inline
from fipy.variables.meshVariable import _MeshVariable
from functools import reduce

class CellVariablemod(CellVariable):
    """ adapted cellVariable
    """        
    @property
    def arithmeticFaceValue(self):
        if not hasattr(self, '_arithmeticFaceValue'):
            
            self._arithmeticFaceValue = _ArithmeticCellToFaceVariablemod(self)

        return self._arithmeticFaceValue

    faceValue = arithmeticFaceValue
        
    @property
    def faceGrad(self):
        r"""
        Return :math:`\nabla \phi` as a rank-1 `FaceVariable` using differencing
        for the normal direction(second-order gradient).
        """
        if not hasattr(self, '_faceGrad'):
            self._faceGrad = _FaceGradVariablemod(self)

        return self._faceGrad
        
    def _getArithmeticBaseClass(self, other = None):
        """
        Given `self` and `other`, return the desired base
        class for an operation result.
        """
        return CellVariablemod


class _FaceGradVariablemod(_FaceGradVariable):   
    
    @property
    def _variableClass(self):
        return FaceVariablemod

    def _getArithmeticBaseClass(self, other = None):
        """
        Given `self` and `other`, return the desired base class for an operation
        result.
        """
        if other is None:
            return FaceVariablemod
        return _MeshVariable._getArithmeticBaseClass(self, other)
        
    @property
    def divergence(self): 
        raise Exception('FaceGradVariablemod::divergence')
        s = (slice(0, None, None),) + (numerix.newaxis,) * (len(self.shape) - 2) + (slice(0, None, None),)
        self._divergence = _AddOverFacesVariablemod((self * self.mesh._orientedAreaProjections[s]).sum(0))

        return self._divergence


        
class _ArithmeticCellToFaceVariablemod(_ArithmeticCellToFaceVariable):
    if inline.doInline:
        def _calcValue_(self, alpha, id1, id2):
            raise Exception("cellVariablemod::_ArithmeticCellToFaceVariablemod unexpected not inline calculation")
    else:
        def _calcValue_(self, alpha, id1, id2):
            cell1 = numerix.take(self.var, id1, axis=-1)
            cell2 = numerix.take(self.var, id2, axis=-1)
            return (cell2 - cell1) * alpha + cell1       
    
    @property
    def _variableClass(self):
        return FaceVariablemod
        
    def _getArithmeticBaseClass(self, other = None):
        """
        Given `self` and `other`, return the desired base class for an operation
        result.
        """
        return FaceVariablemod

    def decimal(self,num):
        num_str = str(num)
        if(num < 1):
            parts = num_str.split('e-', 2)
        else:
            parts = num_str.split('e', 2)
        decimal = parts[1] if len(parts) > 1 else '0'
        decimal = int(decimal)
        return decimal
            
    @property
    def divergence(self): 
        raise Exception('_ArithmeticCellToFaceVariablemod::divergence')
        s = (slice(0, None, None),) + (numerix.newaxis,) * (len(self.shape) - 2) + (slice(0, None, None),)
        self._divergence = _AddOverFacesVariablemod((self * self.mesh._orientedAreaProjections[s]).sum(0))

        return self._divergence

class FaceVariablemod(FaceVariable):
    @property
    def _variableClass(self):
        return FaceVariablemod


    @property
    def divergence(self): #this is the one used for computation of carbon flow
        s = (slice(0, None, None),) + (numerix.newaxis,) * (len(self.shape) - 2) + (slice(0, None, None),)
        self._divergence = _AddOverFacesVariablemod((self * self.mesh._orientedAreaProjections[s] ).sum(0))#=> change here the orinted area projection sent
        #print('oriented Area ', self.mesh._orientedAreaProjections[s],  self.mesh.weightedNodeFaces)
        return self._divergence

    def _getArithmeticBaseClass(self, other = None):
        """
        Given `self` and `other`, return the desired base class for an operation
        result.
        """
        if other is None:
            return FaceVariablemod

        return _MeshVariable._getArithmeticBaseClass(self, other)        


        
class _AddOverFacesVariablemod(_AddOverFacesVariable):
    r"""surface integral of `self.faceVariable`, :math:`\phi_f`

    .. math:: \int_S \phi_f\,dS \approx \frac{\sum_f \phi_f A_f}{V_P}

    Returns
    -------
    integral : CellVariable
        volume-weighted sum
    """
    def _calcValueNoInlineSTOP(self):
        ids = self.mesh.cellFaceIDs
        contributions = numerix.take(self.faceVariable, ids, axis=-1)
        
        s = (numerix.newaxis,) * (len(contributions.shape) - 2) + (slice(0, None, None),) + (slice(0, None, None),)
        
        faceContributions = contributions * self.mesh._cellToFaceOrientations[s]*self.mesh.weightedNodeFaces #changes sign + does weighing by area surface
        
        if(len(faceContributions)>1):
            faceContributions = numerix.vstack((numerix.array(faceContributions[0] + faceContributions[2]), numerix.array(faceContributions[1]+faceContributions[3])))

        return numerix.tensordot(numerix.ones(faceContributions.shape[-2], 'd'),
                                faceContributions, (0, -2)) / self.mesh.cellVolumes
                                                               

    def _calcValueInline(self):
        print('_calcValueInline')
        NCells = self.mesh.numberOfCells
        ids = self.mesh.cellFaceIDs

        val = self._array.copy()

        inline._runInline("""
        int i;

        for(i = 0; i < numberOfCells; i++)
          {
          int j;
          value[i] = 0.;
          for(j = 0; j < numberOfCellFaces; j++)
            {
              // cellFaceIDs can be masked, which caused subtle and
              // unreproducible problems on OS X (who knows why not elsewhere)
              long id = ids[i + j * numberOfCells];
              if (id >= 0) {
                  value[i] += orientations[i + j * numberOfCells] * faceVariable[id] + orientations[i +2 + j * numberOfCells] * faceVariable[id] ;
              }
            }
            value[i] = value[i] / cellVolume[i];
          }
          """,
                          numberOfCellFaces = self.mesh._maxFacesPerCell -2,
                          numberOfCells = NCells,
                          faceVariable = self.faceVariable.numericValue,
                          ids = numerix.array(ids),
                          value = val,
                          orientations = numerix.array(self.mesh._cellToFaceOrientations),
                          cellVolume = numerix.array(self.mesh.cellVolumes))

        return self._makeValue(value = val)
        
    def _calcValueNoInline(self):
        ids = self.mesh.cellFaceIDs
        contributions = numerix.take(self.faceVariable, ids, axis=-1)
        
        s = (numerix.newaxis,) * (len(contributions.shape) - 2) + (slice(0, None, None),) + (slice(0, None, None),)
        #print('contrib ',contributions , conductivity2conductance)
        faceContributions = contributions * self.mesh._cellToFaceOrientations[s] #just changes sign
        #print('contrib ',contributions * self.mesh._cellToFaceOrientations[s],faceContributions)
                                 
        temp = faceContributions.copy()
        maskedValues = np.vstack(( [not isinstance(xi, float) for xi in temp[0]],[not isinstance(xi, float) for xi in temp[1]],
            [not isinstance(xi, float) for xi in temp[2]],[not isinstance(xi, float) for xi in temp[3]]))
        faceAreas = np.take(self.mesh._faceAreas,ids)
        faceAreas = MA.masked_where(maskedValues,faceAreas)
        #faceAreas = numerix.MA.filled(faceAreas, 0)
        
        if(len(temp)>1):
            temptemp = temp.copy()
            temp = temp * faceAreas
            temp[2] = numerix.MA.filled(temp[2], temp[0])
            temp[3] = numerix.MA.filled(temp[3], temp[1])
            
            faceAreas[2] = numerix.MA.filled(faceAreas[2], faceAreas[0])
            faceAreas[3] = numerix.MA.filled(faceAreas[3], faceAreas[1])
            faceAreas[0] = (faceAreas[0] + faceAreas[2])
            faceAreas[1] = (faceAreas[1] + faceAreas[3])
            
            temp[0] = (temp[0] + temp[2])/faceAreas[0]
            temp[1] = (temp[1] + temp[3])/faceAreas[1]
            temp = numerix.vstack((numerix.array(temp[0]), numerix.array(temp[1])))
            
            #temptemp = temptemp * self.mesh.weightedNodeFaces
            #temptemp = numerix.vstack((numerix.array(temptemp[0] + temptemp[2]), numerix.array(temptemp[1]+temptemp[3])))
            #print(self.mesh.cellFaceIDs)
            #print("temptemps, ", temp[0]-temptemp[0],  temp[0],temptemp[0], self.mesh.weightedNodeFaces[0], self.mesh.weightedNodeFaces[2], contributions[0])

        return numerix.tensordot(numerix.ones(temp.shape[-2], 'd'),
                                temp, (0, -2)) / self.mesh.cellVolumes