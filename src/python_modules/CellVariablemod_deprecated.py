
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

            
    @property
    def divergence(self): 
        s = (slice(0, None, None),) + (numerix.newaxis,) * (len(self.shape) - 2) + (slice(0, None, None),)
        self._divergence = _AddOverFacesVariablemod((self * self.mesh._orientedAreaProjections[s]).sum(0))

        return self._divergence

class FaceVariablemod(FaceVariable):
    @property
    def _variableClass(self):
        return FaceVariablemod

    @property
    def divergence(self): 
        s = (slice(0, None, None),) + (numerix.newaxis,) * (len(self.shape) - 2) + (slice(0, None, None),)
        self._divergence = _AddOverFacesVariablemod((self * self.mesh._orientedAreaProjections[s]).sum(0))

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
    def _calcValue(self):
        ids = self.mesh.cellFaceIDs
        contributions = numerix.take(self.faceVariable, ids, axis=-1)
        s = (numerix.newaxis,) * (len(contributions.shape) - 2) + (slice(0, None, None),) + (slice(0, None, None),)

        faceContributions = contributions * self.mesh._cellToFaceOrientations[s] #just changes sign
        
                                 
        temp = faceContributions.copy()
        maskedValues = np.vstack(( [not isinstance(xi, float) for xi in temp[0]],[not isinstance(xi, float) for xi in temp[1]],
            [not isinstance(xi, float) for xi in temp[2]],[not isinstance(xi, float) for xi in temp[3]]))
        faceAreas = np.take(self.mesh._faceAreas,ids)
        faceAreas = MA.masked_where(maskedValues,faceAreas)
        if(len(temp)>1):
        
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
            
        return numerix.tensordot(numerix.ones(temp.shape[-2], 'd'),
                                temp, (0, -2)) / self.mesh.cellVolumes
                                
