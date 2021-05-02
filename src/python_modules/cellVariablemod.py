
#from fipy.variables.cellVariable import cellVariable
from fipy.tools.numerix import MA
import numpy as np
from fipy.tools import numerix

class cellVariablemod(cellVariable):
    """ adapted cellVariable
    """

    def __init__(self, mesh, name='', value=0., rank=None, elementshape=None, unit=None, hasOld=0):
        cellVariable.__init__(self, mesh=mesh, name=name, value=value,
                               rank=rank, elementshape=elementshape, unit=unit, hasOld=hasOld)
        
    
    @property
    def arithmeticFaceValue(self):
        if not hasattr(self, '_arithmeticFaceValue'):
            from fipy.variables.arithmeticCellToFaceVariable import _ArithmeticCellToFaceVariable
            self._arithmeticFaceValue = _ArithmeticCellToFaceVariablemod(self)

        return self._arithmeticFaceValue

    faceValue = arithmeticFaceValue
        
        
        
class _ArithmeticCellToFaceVariablemod(_ArithmeticCellToFaceVariable):
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
            from fipy.variables.addOverFacesVariable import _AddOverFacesVariable
            #print('faceVariable::dvergence')
            #slice(0, None, None): takes everything
            s = (slice(0, None, None),) + (numerix.newaxis,) * (len(self.shape) - 2) + (slice(0, None, None),)
            self._divergence = _AddOverFacesVariablemod((self * self.mesh._orientedAreaProjections[s]).sum(0))
            #creates a new cellVar
            ##value * surface of faces
            #sum(iterable, start) : sum the dimentions (?): for us changes nothing as everything along z axis
            #print('s\n',self,'\n',(self * self.mesh._orientedAreaProjections[s]).sum(0))

        return self._divergence
        
class _AddOverFacesVariablemod(_AddOverFacesVariable):
    r"""surface integral of `self.faceVariable`, :math:`\phi_f`

    .. math:: \int_S \phi_f\,dS \approx \frac{\sum_f \phi_f A_f}{V_P}

    Returns
    -------
    integral : CellVariable
        volume-weighted sum
    """
    def _calcValueNoInline(self):
        print('_calcValueNoInline')
        ids = self.mesh.cellFaceIDs
        contributions = numerix.take(self.faceVariable, ids, axis=-1)
        # FIXME: numerix.MA.filled casts away dimensions
        s = (numerix.newaxis,) * (len(contributions.shape) - 2) + (slice(0, None, None),) + (slice(0, None, None),)

        faceContributions = contributions * self.mesh._cellToFaceOrientations[s] #just changes sign
        print('ids ', ids, '\n\n',contributions,'\n\n', faceContributions,'\n\n',
        numerix.tensordot(numerix.ones(faceContributions.shape[-2], 'd'),
                                 numerix.MA.filled(faceContributions, 0.), (0, -2)))
        print('cell to face orientation ',self.mesh._cellToFaceOrientations)
        
        print('_calcValueNoInline end')
        print('_calcValueNoInline try')
        temp = faceContributions.copy()
        if(len(temp)>1):
            temp[2] = numerix.MA.filled(temp[2], temp[0])
            temp[3] = numerix.MA.filled(temp[3], temp[1])
            print(temp[2],temp[3],'\n', temp)
            temp[0] = (temp[0] + temp[2])/2
            temp[1] = (temp[1] + temp[3])/2
            #temp = numerix.array([numerix.array(temp[0])[0], numerix.array(temp[1])[0]])
            print(numerix.vstack((numerix.array(temp[0]), numerix.array(temp[1]))))#, numerix.MA.filled(temp[2], temp[0]))
            temp = numerix.vstack((numerix.array(temp[0]), numerix.array(temp[1])))
            print(numerix.tensordot(numerix.ones(temp.shape[-2], 'd'),
                                    temp, (0, -2)))
        print('_calcValueNoInline try end')
        #TODO:
        #mean of bot 1, bot2 on one side and mean of top 1 top2 2 on the other side.
        #
        #
        
        
        return numerix.tensordot(numerix.ones(faceContributions.shape[-2], 'd'),
                                 numerix.MA.filled(faceContributions, 0.), (0, -2)) / self.mesh.cellVolumes
        #return numerix.tensordot(numerix.ones(temp.shape[-2], 'd'),
         #                       temp, (0, -2)) / self.mesh.cellVolumes