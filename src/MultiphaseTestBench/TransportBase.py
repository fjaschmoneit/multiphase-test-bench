import Fields
import numpy as np
import MeshConfig

class transportBase:

    def __init__(self, mesh, field, linSystem, depFieldShape):
        self.mesh = mesh
        self.linSystem = linSystem
        self.depFieldShape = depFieldShape
        self.constSourceField = None
        self.boundary = {}
        self.depField = field

        self.constSourceField = Fields.newDataField(shape=MeshConfig.SHAPE_SCALAR_CV, value=0.0)
        self.sourceField_c = Fields.newDataField(shape=MeshConfig.SHAPE_SCALAR_CV, value=0.0)
        self.a_e = Fields.newDataField(shape=MeshConfig.SHAPE_SCALAR_CV)
        self.a_w = Fields.newDataField(shape=MeshConfig.SHAPE_SCALAR_CV)
        self.a_s = Fields.newDataField(shape=MeshConfig.SHAPE_SCALAR_CV)
        self.a_n = Fields.newDataField(shape=MeshConfig.SHAPE_SCALAR_CV)
        self.a_p = Fields.newDataField(shape=MeshConfig.SHAPE_SCALAR_CV)

    def getField(self):
        return self.depField

    # def updateFluxes(self):
    #     self.reset()
    #
    #     F_u, F_v = self.calcConvFlux()
    #     D_u, D_v = self.calcDiffFlux()
    #
    #     self.a_w[east] = DifferenceSchemes.centralDifference(D_u, F_u, 'west' )[internal_u]
    #     self.a_e[west] = DifferenceSchemes.centralDifference(D_u, F_u, 'east')[internal_u]
    #
    #     self.a_n = DifferenceSchemes.centralDifference(D_v, F_v, 'north')[north]
    #     self.a_s = DifferenceSchemes.centralDifference(D_v, F_v, 'south')[south]
    #
    #     self.correctBCs()
    #
    #     self.a_p += self.a_w + self.a_e + self.a_s + self.a_n
    # # self.a_p[internal_u] += F_u[east] - F_u[west] + F_v[south] - F_v[north]    # do I need to include these terms when u itself is transported?

    def reset(self):
        self.linSystem.reset(shape=self.depFieldShape)
        self.sourceField_c.fill(0.0)
        self.a_e.fill(0.0)
        self.a_w.fill(0.0)
        self.a_s.fill(0.0)
        self.a_n.fill(0.0)
        self.a_p.fill(0.0)

    def setConstSourceField(self, value):
        self.constSourceField.fill(value)

    def setDiffusionCoefficient(self, value):
        self.diffusionCoefficient = value

    def getCentreMatrixCoeffs(self):
        return self.a_p

    def getCoefficientsInDirection(self,direction):
        directionCoeffDict = {
            'west' : self.a_w,
            'east': self.a_e,
            'north': self.a_n,
            'south': self.a_s
        }
        return directionCoeffDict[direction]

    def getCoefficientsInOppositeDirection(self, direction):
        oppositeDirectionCoeffDict = {
            'west': self.a_e,
            'east': self.a_w,
            'north': self.a_s,
            'south': self.a_n
        }
        return oppositeDirectionCoeffDict[direction]
    # def getDepField(self):
    #     return self.depField

    def correctBCs(self):
        for argDict in self.boundary.values():
            bcType = argDict['type']
            self.boundaryModels[bcType](self, **argDict)
#            self.boundaryModels[bcType](self.depField.data, **argDict)

    def updateLinSystem(self):
        self.linSystem.reset(shape=self.depFieldShape)
        self.linSystem.set_e_coeffs(self.a_e)
        self.linSystem.set_w_coeffs(self.a_w)
        self.linSystem.set_n_coeffs(self.a_n)
        self.linSystem.set_s_coeffs(self.a_s)
        self.linSystem.set_p_coeffs(self.a_p)
        self.linSystem.set_b(self.sourceField_c)

    def solve(self):
        return np.copy(self.linSystem.solve())

#------------------ abstract member function ------------------

    def setVonNeumann(self, direction):
        pass

    def setDerichlet(self, direction, value):
        pass

    def calcConvFlux(self):
        pass

    def calcDiffFlux(self):
        pass