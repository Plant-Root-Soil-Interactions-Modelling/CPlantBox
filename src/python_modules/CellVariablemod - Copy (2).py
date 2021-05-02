
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

    def __init__(self, mesh, name='', value=0., rank=None, elementshape=None, unit=None, hasOld=0):
        CellVariable.__init__(self, mesh=mesh, name=name, value=value,
                               rank=rank, elementshape=elementshape, unit=unit, hasOld=hasOld)
        
   
        
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



class _FaceGradVariablemod(_FaceGradVariable):   
    
    @property
    def _variableClass(self):
        return FaceGradVariablemod    
    @property
    def divergence(self): 
        if not hasattr(self, '_divergence'):
            s = (slice(0, None, None),) + (numerix.newaxis,) * (len(self.shape) - 2) + (slice(0, None, None),)
            self._divergence = _AddOverFacesVariablemod((self * self.mesh._orientedAreaProjections[s]).sum(0))

        return self._divergence

'''
def _OperatorVariableClass(_OperatorVariableClass):
    class _OperatorVariable(baseClass):
        def __init__(self, op, var, opShape=(), canInline=True, unit=None, inlineComment=None, valueMattersForUnit=None, *args, **kwargs):
            self.op = op
            self.var = var
            self.opShape = opShape
            self._unit = unit
            if valueMattersForUnit is None:
                self.valueMattersForUnit = [False for v in var]
            else:
                self.valueMattersForUnit = valueMattersForUnit
            self.canInline = canInline  #allows for certain functions to opt out of --inline
            baseClass = 
            baseClass.__init__(self, value=None, *args, **kwargs)
            self.name = ''
            for var in self.var:    #C does not accept units
                if not var.unit.isDimensionless():
                    self.canInline = False
                    break

            for aVar in self.var:
                self._requires(aVar)

            self.dontCacheMe()

            self.comment = inlineComment
'''

        
class _ArithmeticCellToFaceVariablemod(_ArithmeticCellToFaceVariable):
    @property
    def _variableClass(self):
        return FaceVariablemod#_ArithmeticCellToFaceVariablemod
        
    def _getArithmeticBaseClass(self, other=None):
        """
        Given `self` and `other`, return the desired base class for an operation
        result.
        """
        if other is None:
            return _ArithmeticCellToFaceVariablemod

        if self._broadcastShape(other) is not None:
            # If self and other have the same base class, result has that base class.
            # If self derives from other, result has self's base class.
            # If other derives from self, result has other's base class.
            # If self and other don't have a common base, we don't know how to combine them.
            from fipy.variables.constant import _Constant
            if isinstance(self, other._getArithmeticBaseClass()) or isinstance(other, _Constant):
                return self._getArithmeticBaseClass()
            else:
                return None
        else:
            # If self and other have un-broadcastable shapes, we don't know how to combine them.
            return None
    def _OperatorVariableClass(self, baseClass=None):
        baseClass = Variable._OperatorVariableClass(self, baseClass=baseClass)

        class _MeshOperatorVariable(baseClass):
            def __init__(self, op, var, opShape=None, canInline=True,
                         *args, **kwargs):
                mesh = reduce(lambda a, b: a or b,
                              [getattr(v, "mesh", None) for v in var])
                for shape in [opShape] + [getattr(v, "opShape", None) for v in var]:
                    if shape is not None:
                        opShape = shape
                        break
##                 opShape = reduce(lambda a, b: a or b,
##                                  [opShape] + [getattr(v, "opShape", None) for v in var])
                if opShape is not None:
                    elementshape = opShape[:-1]
                else:
                    elementshape = reduce(lambda a, b: a or b,
                                          [getattr(v, "elementshape", None) for v in var])

                baseClass.__init__(self, mesh=mesh, op=op, var=var,
                                   opShape=opShape, canInline=canInline,
                                   elementshape=elementshape,
                                   *args, **kwargs)

            @property
            def rank(self):
                return len(self.opShape) - 1

        return _MeshOperatorVariable 


class FaceVariablemod(_MeshVariable):
    @property
    def _variableClass(self):
        return FaceVariablemod

    @property
    def divergence(self): 
        if not hasattr(self, '_divergence'):
            s = (slice(0, None, None),) + (numerix.newaxis,) * (len(self.shape) - 2) + (slice(0, None, None),)
            self._divergence = _AddOverFacesVariablemod((self * self.mesh._orientedAreaProjections[s]).sum(0))

        return self._divergence
        
    def _OperatorVariableClass(self, baseClass=None):
        baseClass = Variable._OperatorVariableClass(self, baseClass=baseClass)

        class _FaceVariablemodOperatorVariable(baseClass):
            def __init__(self, op, var, opShape=None, canInline=True,
                         *args, **kwargs):
                mesh = reduce(lambda a, b: a or b,
                              [getattr(v, "mesh", None) for v in var])
                for shape in [opShape] + [getattr(v, "opShape", None) for v in var]:
                    if shape is not None:
                        opShape = shape
                        break
##                 opShape = reduce(lambda a, b: a or b,
##                                  [opShape] + [getattr(v, "opShape", None) for v in var])
                if opShape is not None:
                    elementshape = opShape[:-1]
                else:
                    elementshape = reduce(lambda a, b: a or b,
                                          [getattr(v, "elementshape", None) for v in var])

                baseClass.__init__(self, mesh=mesh, op=op, var=var,
                                   opShape=opShape, canInline=canInline,
                                   elementshape=elementshape,
                                   *args, **kwargs)

            @property
            def rank(self):
                return len(self.opShape) - 1

        return _FaceVariablemodOperatorVariable

        
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
                                

def _BinaryOperatorVariable(operatorClass=None):
    class binOp(operatorClass):

        def _calcValue_(self):
            from fipy.variables.variable import Variable
            print('binop normal',type(self))
            if isinstance(self.var[1], Variable):
                val1 = self.var[1].value
            else:
                if isinstance(self.var[1], type('')):
                    self.var[1] = physicalField.PhysicalField(value=self.var[1])
                val1 = self.var[1]

            return self.op(self.var[0].value, val1)

        @property
        def unit(self):
            if self._unit is None:
                try:
                    var = self._varProxy
                    return self._extractUnit(self.op(var[0], var[1]))
                except:
                    return self._extractUnit(self._calcValue_())
            else:
                return self._unit

        def _getRepresentation(self, style="__repr__", argDict={}, id=id, freshen=False):
            self.id = id
            if (style == "__repr__") and hasattr(self, '_name') and len(self._name) > 0:
                return self._name
            else:
                return "(" + operatorClass._getRepresentation(self, style=style, argDict=argDict, id=id, freshen=freshen) + ")"

    return binOp
    
def _OperatorVariableClass(baseClass=object):
    class _OperatorVariable(baseClass):
        def __init__(self, op, var, opShape=(), canInline=True, unit=None, inlineComment=None, valueMattersForUnit=None, *args, **kwargs):
            self.op = op
            self.var = var
            self.opShape = opShape
            self._unit = unit
            if valueMattersForUnit is None:
                self.valueMattersForUnit = [False for v in var]
            else:
                self.valueMattersForUnit = valueMattersForUnit
            self.canInline = canInline  #allows for certain functions to opt out of --inline
            print('opVaropVar 226', baseClass, args, kwargs)
            baseClass.__init__(self, value=None, *args, **kwargs)
            self.name = ''
            for var in self.var:    #C does not accept units
                if not var.unit.isDimensionless():
                    self.canInline = False
                    break

            for aVar in self.var:
                self._requires(aVar)

            self.dontCacheMe()

            self.comment = inlineComment

        def __setitem__(self, index, value):
            raise TypeError("The value of an `_OperatorVariable` cannot be assigned")

        def setValue(self, value, unit=None, where=None):
            raise TypeError("The value of an `_OperatorVariable` cannot be assigned")

        def _calcValue(self):
            if not self.canInline:
                return self._calcValue_()
            else:
                from fipy.tools import inline
                if inline.doInline:
                    return self._execInline(comment=self.comment)
                else:
                    return self._calcValue_()

        def _calcValue_(self):
            pass

        def _isCached(self):
            return (Variable._isCached(self)
                    or (len(self.subscribedVariables) > 1 and not self._cacheNever))

        def _getCstring(self, argDict={}, id="", freshen=False):
            if self.canInline: # and not self._isCached():
                s = self._getRepresentation(style="C", argDict=argDict, id=id, freshen=freshen)
            else:
                s = baseClass._getCstring(self, argDict=argDict, id=id)
            if freshen:
                self._markFresh()

            return s

        def _getRepresentation(self, style="__repr__", argDict={}, id=id, freshen=False):
            """

            Parameters
            ----------
            style : {'__repr__', 'name', 'TeX', 'C'}
               desired formatting for representation
            """
            if isinstance(self.op, numerix.ufunc):
                return "%s(%s)" % (self.op.__name__, ", ".join([self.__var(i, style, argDict, id, freshen)
                                                               for i in range(len(self.var))]))

            try:
                instructions = dis.get_instructions(self.op.__code__)
                parseInstructions = self._py3kInstructions
            except AttributeError:
                instructions = [ord(byte) for byte in self.op.__code__.co_code]
                parseInstructions = self._py2kInstructions

            return parseInstructions(instructions, style=style, argDict=argDict, id=id, freshen=freshen)

        def __var(self, i, style, argDict, id, freshen):
            v = self.var[i]
            if style == "__repr__":
                result = repr(v)
            elif style == "name":
                if isinstance(v, Variable):
                    result = v.name
                    if len(result) == 0:
                        # The string form of a variable
                        # would probably be too long and messy.
                        # Just give shorthand.
                        result = "%s(...)" % v.__class__.__name__
                elif type(v) in (type(1), type(1.)):
                    result = repr(v)
                else:
                    # The string form of anything but a
                    # number would be too long and messy.
                    # Just give shorthand.
                    result = "<...>"

            elif style == "TeX":
                raise Exception("TeX style not yet implemented")
            elif style == "C":
                if not v._isCached():
                    result = v._getCstring(argDict, id=id + str(i), freshen=freshen)
                    if isinstance(v, Variable):
                        v._value = None
                    else:
                        v.value = None
                else:
                    result = v._variableClass._getCstring(v, argDict,
                                                               id=id + str(i),
                                                               freshen=False)
            else:
                raise SyntaxError("Unknown style: %s" % style)

            return result

        _unop = {
            10: "+", 11: "-", 12: "not ", 15: "~"
        }

        _binop = {
            19: "**", 20: "*", 21: "/", 22: "%", 23: "+", 24: "-", 26: "//", 27: "/",
                    62: "<<", 63: ">>", 64: "&", 65: "^", 66: "|", 106: "=="
        }

        def _py2kInstructions(self, bytecodes, style, argDict, id, freshen):
            def _popIndex():
                return bytecodes.pop(0) + bytecodes.pop(0) * 256

            allbytecodes = bytecodes[:]

            stack = []

            while len(bytecodes) > 0:
                bytecode = bytecodes.pop(0)
                if dis.opname[bytecode] == 'UNARY_CONVERT':
                    stack.append("`" + stack.pop() + "`")
                elif dis.opname[bytecode] == 'BINARY_SUBSCR':
                    stack.append(stack.pop(-2) + "[" + stack.pop() + "]")
                elif dis.opname[bytecode] == 'RETURN_VALUE':
                    s = stack.pop()
                    if style == 'C':
                        return s.replace('numerix.', '').replace('arc', 'a')
                    else:
                        return s
                elif dis.opname[bytecode] == 'LOAD_CONST':
                    stack.append(self.op.__code__.co_consts[_popIndex()])
                elif dis.opname[bytecode] == 'LOAD_ATTR':
                    stack.append(stack.pop() + "." + self.op.__code__.co_names[_popIndex()])
                elif dis.opname[bytecode] == 'COMPARE_OP':
                    stack.append(stack.pop(-2) + " " + dis.cmp_op[_popIndex()] + " " + stack.pop())
                elif dis.opname[bytecode] == 'LOAD_GLOBAL':
                    counter = _popIndex()
                    stack.append(self.op.__code__.co_names[counter])
                elif dis.opname[bytecode] == 'LOAD_FAST':
                    stack.append(self.__var(_popIndex(), style=style, argDict=argDict, id=id, freshen=freshen))
                elif dis.opname[bytecode] == 'CALL_FUNCTION':
                    args = []
                    for j in range(bytecodes.pop(1)):
                        # keyword parameters
                        args.insert(0, stack.pop(-2) + " = " + stack.pop())
                    for j in range(bytecodes.pop(0)):
                        # positional parameters
                        args.insert(0, stack.pop())
                    stack.append(stack.pop() + "(" + ", ".join(args) + ")")
                elif dis.opname[bytecode] == 'LOAD_DEREF':
                    free = self.op.__code__.co_cellvars + self.op.__code__.co_freevars
                    stack.append(free[_popIndex()])
                elif bytecode in self._unop:
                    stack.append(self._unop[bytecode] + '(' + stack.pop() + ')')
                elif bytecode in self._binop:
                    stack.append(stack.pop(-2) + " " + self._binop[bytecode] + " " + stack.pop())
                else:
                    raise SyntaxError("Unknown bytecode: %s in %s of %s, got %s" % (
                       repr(bytecode),
                       repr([bytecode] + bytecodes),
                       repr(allbytecodes),
                       repr(stack)))

        def _py3kInstructions(self, instructions, style, argDict, id, freshen):
            stack = []
            
            for ins in instructions:
                if ins.opname == 'UNARY_CONVERT':
                    stack.append("`" + stack.pop() + "`")
                elif ins.opname == 'BINARY_SUBSCR':
                    stack.append(stack.pop(-2) + "[" + stack.pop() + "]")
                elif ins.opname == 'RETURN_VALUE':
                    s = stack.pop()
                    if style == 'C':
                        return s.replace('numerix.', '').replace('arc', 'a')
                    else:
                        return s
                elif ins.opname == 'LOAD_CONST':
                    stack.append(ins.argval)
                elif ins.opname == 'LOAD_ATTR':
                    stack.append(stack.pop() + "." + ins.argval)
                elif ins.opname == 'COMPARE_OP':
                    stack.append(stack.pop(-2) + " " + dis.cmp_op[ins.arg] + " " + stack.pop())
                elif ins.opname == 'LOAD_GLOBAL':
                    stack.append(ins.argval)
                elif ins.opname == 'LOAD_FAST':
                    stack.append(self.__var(ins.arg, style=style, argDict=argDict, id=id, freshen=freshen))
                elif ins.opname == 'CALL_FUNCTION':
                    # args are last ins.arg items on stack
                    args, stack = stack[-ins.arg:], stack[:-ins.arg]
                    stack.append(stack.pop() + "(" + ", ".join(args) + ")")
                elif ins.opname == 'CALL_FUNCTION_KW':
                    kws = list(stack.pop())
                    # args are last ins.arg items on stack
                    args, stack = stack[-ins.arg:], stack[:-ins.arg]
                    kwargs = []
                    while kws:
                        kwargs.append(kws.pop() + "=" + args.pop())
                    stack.append(stack.pop() + "(" + ", ".join(args + kwargs) + ")")
                elif ins.opname == 'LOAD_DEREF':
                    stack.append(ins.argval)
                elif ins.opcode in self._unop:
                    stack.append(self._unop[ins.opcode] + '(' + stack.pop() + ')')
                elif ins.opcode in self._binop:
                    stack.append(stack.pop(-2) + " " + self._binop[ins.opcode] + " " + stack.pop())
                else:
                    raise SyntaxError("Unknown instruction: %s" % repr(ins))

        @property
        def _varProxy(self):
            """list of dimensional scalars that stand in for `self.var`

            Used for determining units of result without doing expensive computation"""

            return [v if valueMatters else v._unitAsOne for v, valueMatters in zip(self.var, self.valueMattersForUnit)]

        def __repr__(self):
            return self._getRepresentation()

        def __reduce__(self):
            """
            Allows `_OperatorVariables` to be pickled
            """
            state =  self.__getstate__()
            if 'mesh' in list(state.keys()):
                args = (state['mesh'],)
            else:
                args = ()

            return (self._variableClass, args, self.__getstate__())

        def _getName(self):
            name = baseClass._getName(self)
            if len(name) == 0:
                name = self._getRepresentation(style="name")
            return name

        name = property(_getName, baseClass._setName)

        @property
        def shape(self):
            if self.opShape is not None:
                return self.opShape
            else:
                return baseClass.getShape(self)
##             return baseClass.getShape(self) or self.opShape
    print('end opvaropvar HERE 481', _OperatorVariable)
    return _OperatorVariable    