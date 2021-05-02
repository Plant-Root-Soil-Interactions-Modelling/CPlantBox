
from fipy import CellVariable, FaceVariable, Variable
from fipy.variables.addOverFacesVariable import _AddOverFacesVariable
from fipy.variables.arithmeticCellToFaceVariable import _ArithmeticCellToFaceVariable
from fipy.variables.faceGradVariable import _FaceGradVariable
from fipy.tools.numerix import MA
import numpy as np
from fipy.tools import numerix
from fipy.tools import inline

class CellVariablemod(CellVariable):
    """ adapted cellVariable
    """

    def __init__(self, mesh, name='', value=0., rank=None, elementshape=None, unit=None, hasOld=0):
        CellVariable.__init__(self, mesh=mesh, name=name, value=value,
                               rank=rank, elementshape=elementshape, unit=unit, hasOld=hasOld)
        
    def __mulmod__(self, other):
        print('\n\ntypes__mul__mod', type(self) , type(other))
        if(type(self) == type(other)):
            return self._BinaryOperatorVariablemod(lambda a, b: a*b, other)
        else:
            return CellVariable.__mul__
    __rmul__ = __mulmod__

    @property
    def _variableClass(self):
        return CellVariablemod
        
    def _shapeClassAndOther(self, opShape, operatorClass, other):
        """
        Determine the shape of the result, the base class of the result, and (if
        necessary) a modified form of `other` that is suitable for the
        operation.
        """
        # If the caller has not specified a base class for the binop,
        # check if the member Variables know what type of Variable should
        # result from the operation.
        print('\nbaseclassBB ',operatorClass, other)
        baseClass = operatorClass or self._getArithmeticBaseClass(other)
        print('\nbaseclassAA ',baseClass)

        # If the caller has not specified a shape for the result, determine the
        # shape from the base class or from the inputs
        if opShape is None:
            opShape = self._broadcastShape(other)

        return (opShape, baseClass, other)


    def _getArithmeticBaseClass(self, other = None):
        """
        Given `self` and `other`, return the desired base
        class for an operation result.
        """
        if other is None:
            print('_getArithmeticBaseClass::if other is None')
            return CellVariablemod

        return _MeshVariable._getArithmeticBaseClass(self, other)
    
    def _BinaryOperatorVariablemod(self, op, other, operatorClass = None, opShape=None, canInline=True, unit=None,
                                    value0mattersForUnit=False, value1mattersForUnit=False):
            opShape, baseClass, other = self._shapeClassAndOther(opShape, operatorClass, other)
            baseClass = type(self)
            print('type beginning ',type(self), type(other))
            if opShape is None or baseClass is None:
                return NotImplemented

            canInline = False

            # obtain a general operator class with the desired base class
            print('\n\noperatorClass, ',operatorClass, self._OperatorVariableClass(baseClass),opShape,  other)
            operatorClass = operatorClass or self._OperatorVariableClass(baseClass)
            print('\n\noperatorClass, ',operatorClass)
            #binOp = _BinaryOperatorVariablemod(CellVariablemod)

            return op(self.var[0].value, val1)   
        
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

    @property
    def divergence(self): 
        r"""the divergence of `self`, :math:`\vec{u}`,

        .. math:: \nabla\cdot\vec{u} \approx \frac{\sum_f (\vec{u}\cdot\hat{n})_f A_f}{V_P}

        Returns
        -------
        divergence : fipy.variables.cellVariable.CellVariable
            one rank lower than `self`

        """
        if not hasattr(self, '_divergence'):
            
            print('\n\nfacegradVariablemod::dvergence2')
            #slice(0, None, None): takes everything
            s = (slice(0, None, None),) + (numerix.newaxis,) * (len(self.shape) - 2) + (slice(0, None, None),)
            #print()
            self._divergence = _AddOverFacesVariablemod((self * self.mesh._orientedAreaProjections[s]).sum(0))
            #creates a new cellVar
            ##value * surface of faces
            #sum(iterable, start) : sum the dimentions (?): for us changes nothing as everything along z axis
            #print('s\n',self,'\n',(self * self.mesh._orientedAreaProjections[s]).sum(0))

        return self._divergence

class _FaceGradVariablemod(FaceVariable):   
    def __init__(self, var):
        FaceVariable.__init__(self, mesh=var.mesh, elementshape=(var.mesh.dim,) + var.shape[:-1])
        self.var = self._requires(var)
    
    
    def __mul__(self, other):
        print('\n\ntypes__mulfacevar__', type(self) , type(other))
        if(type(self) == type(other)):
            return self._BinaryOperatorVariablemod(lambda a, b: a*b, other)
        else:
            return CellVariable.__mul__
    __rmul__ = __mul__
    

        
    def _calcValue(self):
        #print(' _calcValueNoInline')
        dAP = self.mesh._cellDistances
        id1, id2 = self.mesh._adjacentCellIDs

        N2 = numerix.take(self.var.value, id2, axis=-1)

        faceMask = numerix.array(self.mesh.exteriorFaces)
        #print('dAP\n',dAP,'\n',id1,'\n',id2, '\n',N2, '\n',faceMask)

        ## The following conditional is required because empty
        ## indexing is not altogether functional.  This
        ## numpy.empty((0,))[[]] and this numpy.empty((0,))[...,[]]
        ## both work, but this numpy.empty((3, 0))[...,[]] is
        ## broken.

        if self.var.faceValue.shape[-1] != 0:
            s = (Ellipsis, faceMask)
        else:
            s = (faceMask,)

        N2[s] = self.var.faceValue[s]

        N = (N2 - numerix.take(self.var, id1, axis=-1)) / dAP

        normals = self.mesh._orientedFaceNormals
        #print('s\n',s,N2[s],N2,self.var.faceValue[s], N,numerix.take(self.var, id1, axis=-1))  
        tangents1 = self.mesh._faceTangents1
        tangents2 = self.mesh._faceTangents2
        cellGrad = self.var.grad.numericValue

        grad1 = numerix.take(cellGrad, id1, axis=-1)
        grad2 = numerix.take(cellGrad, id2, axis=-1)
        #print('grad4facegrad\n',cellGrad,'\n',grad1,'\n',grad2)
        s = (slice(0, None, None),) + (numerix.newaxis,) * (len(grad1.shape) - 2) + (slice(0, None, None),)
        t1grad1 = numerix.sum(tangents1[s] * grad1, 0)
        t1grad2 = numerix.sum(tangents1[s] * grad2, 0)
        t2grad1 = numerix.sum(tangents2[s] * grad1, 0)
        t2grad2 = numerix.sum(tangents2[s] * grad2, 0)

        T1 = (t1grad1 + t1grad2) / 2.
        T2 = (t2grad1 + t2grad2) / 2.
        
        return normals[s] * N[numerix.newaxis] + tangents1[s] * T1[numerix.newaxis] + tangents2[s] * T2[numerix.newaxis]
        
    @property
    def divergence(self): 
        r"""the divergence of `self`, :math:`\vec{u}`,

        .. math:: \nabla\cdot\vec{u} \approx \frac{\sum_f (\vec{u}\cdot\hat{n})_f A_f}{V_P}

        Returns
        -------
        divergence : fipy.variables.cellVariable.CellVariable
            one rank lower than `self`

        """
        if not hasattr(self, '_divergence'):
            
            #print('facegradVariablemod::dvergence')
            #slice(0, None, None): takes everything
            s = (slice(0, None, None),) + (numerix.newaxis,) * (len(self.shape) - 2) + (slice(0, None, None),)
            #print()
            self._divergence = _AddOverFacesVariablemod((self * self.mesh._orientedAreaProjections[s]).sum(0))
            #creates a new cellVar
            ##value * surface of faces
            #sum(iterable, start) : sum the dimentions (?): for us changes nothing as everything along z axis
            #print('s\n',self,'\n',(self * self.mesh._orientedAreaProjections[s]).sum(0))

        return self._divergence
        
class _ArithmeticCellToFaceVariablemod(_ArithmeticCellToFaceVariable):
    
    def divergence(self, other):
        print('\n\ntypes__mul_arithm_mod', type(self))
        return self._BinaryOperatorVariablemod(lambda a, b: a*b, other)
    
    
    @property
    def divergencemod(self): 
        r"""the divergence of `self`, :math:`\vec{u}`,

        .. math:: \nabla\cdot\vec{u} \approx \frac{\sum_f (\vec{u}\cdot\hat{n})_f A_f}{V_P}

        Returns
        -------
        divergence : fipy.variables.cellVariable.CellVariable
            one rank lower than `self`

        """
        if not hasattr(self, '_divergencemod arithm'):
            
            #print('ArithmeticCellToFaceVariablemod::divergence')
            #slice(0, None, None): takes everything
            s = (slice(0, None, None),) + (numerix.newaxis,) * (len(self.shape) - 2) + (slice(0, None, None),)
            self._divergence = _AddOverFacesVariablemod((self * self.mesh._orientedAreaProjections[s]).sum(0))
            #creates a new cellVar
            ##value * surface of faces
            #sum(iterable, start) : sum the dimentions (?): for us changes nothing as everything along z axis
            #print('s\n',self,'\n',(self * self.mesh._orientedAreaProjections[s]).sum(0))

        return self._divergence

    def _BinaryOperatorVariablemod(self, op, other):


        binOp = _BinaryOperatorVariablemod()

        return binOp(op=op, var=[self, other], mesh = self.mesh)       

        
class _AddOverFacesVariablemod(_AddOverFacesVariable):
    r"""surface integral of `self.faceVariable`, :math:`\phi_f`

    .. math:: \int_S \phi_f\,dS \approx \frac{\sum_f \phi_f A_f}{V_P}

    Returns
    -------
    integral : CellVariable
        volume-weighted sum
    """
    def _calcValue_(self):
        print('\n\n_calcValueNoInline\n\n')
        ids = self.mesh.cellFaceIDs
        contributions = numerix.take(self.faceVariable, ids, axis=-1)
        # FIXME: numerix.MA.filled casts away dimensions
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
                                

def _BinaryOperatorVariablemod():
    
    # declare a binary operator class with the desired base class
    class binOp(_FaceGradVariablemod):

        def _calcValue_(self):
            val1 = self.var[1].value
            s = (slice(0, None, None),) + (numerix.newaxis,) * (len(self.shape) - 2) + (slice(0, None, None),)
            print('\n\n_BinaryOperatorVariablemod ',type(self))
            a = (self.op(self.var[0].value, val1)* self.mesh._orientedAreaProjections[s]).sum(0)
            print('\n\n_calcValueNoInline\n\n')
            ids = self.mesh.cellFaceIDs
            contributions = numerix.take(a, ids, axis=-1)
            # FIXME: numerix.MA.filled casts away dimensions
            s = (numerix.newaxis,) * (len(contributions.shape) - 2) + (slice(0, None, None),) + (slice(0, None, None),)

            faceContributions = contributions * self.mesh._cellToFaceOrientations[s] #just changes sign
            
                                     
            temp = faceContributions#.copy()
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
            #)
            
        


    return binOp
    
    
    