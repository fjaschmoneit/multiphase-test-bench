import BoundaryConditions
import Fields
from fieldAccess import *
import MeshConfig
import ObjectRegistry as objReg
import numpy as np
class Pressure():

    def __init__(self, mesh, linSystem):
        self.p = Fields.newField(shape=MeshConfig.SHAPE_SCALAR_CV, governingModel=self, value=0)
        self.depField = self.p

        # self.pcorr = Fields.newDataField(shape=MeshConfig.SHAPE_SCALAR_CV_GHOST, value=0)
        # objReg.FIELDS['pcorr'] = self.pcorr

        self.mesh = mesh
        self.depFieldShape = MeshConfig.SHAPE_SCALAR_CV
        self.sourceField_c = Fields.newDataField(shape=self.depFieldShape, value=0.0)
        self.a_e = Fields.newDataField(shape=self.depFieldShape, value=0.0)
        self.a_w = Fields.newDataField(shape=self.depFieldShape, value=0.0)
        self.a_s = Fields.newDataField(shape=self.depFieldShape, value=0.0)
        self.a_n = Fields.newDataField(shape=self.depFieldShape, value=0.0)
        self.a_p = Fields.newDataField(shape=self.depFieldShape, value=0.0)
        self.linSystem = linSystem

        # pointers to corresp. boundary functions:
        self.boundaryModels = {
            'freeFlow' : BoundaryConditions.pressure.freeFlow,
            'totalPressure' : BoundaryConditions.pressure.totalPressure,
            'constantPressure' : BoundaryConditions.pressure.derichlet
        }

    def getCentreMatrixCoeffs(self):
        return self.a_p

    def updateLinSystem(self):
        self.linSystem.reset(shape=self.depFieldShape)
        self.linSystem.set_e_coeffs(self.a_e)
        self.linSystem.set_w_coeffs(self.a_w)
        self.linSystem.set_n_coeffs(self.a_n)
        self.linSystem.set_s_coeffs(self.a_s)
        self.linSystem.set_p_coeffs(self.a_p)
        self.linSystem.set_b(self.sourceField_c)

    def getCoefficientsInDirection(self,direction):
        directionCoeffDict = {
            'west' : self.a_w,
            'east': self.a_e,
            'north': self.a_n,
            'south': self.a_s
        }
        return directionCoeffDict[direction]

    def correctBCs(self):
        for argDict in self.depField.boundary.values():
            bcType = argDict['type']
            self.boundaryModels[bcType](self, **argDict)

    def getField(self):
        return self.p

    def linkOtherFields(self, closure):
        if len(closure) == 0:
            self.v = objReg.FIELDS['v']
            self.u = objReg.FIELDS['u']
        else:
            self.v = objReg.FIELDS[closure[1]]
            self.u = objReg.FIELDS[closure[0]]

    def solve(self):
        # while np.linalg.det(self.linSystem.A) < 1e-20:
        # self.linSystem.A *= 1e6
        # self.linSystem.b *= 1e6
#        print("condition number: ", np.linalg.cond(self.linSystem.A))
#        print("determinant ", np.linalg.det(self.linSystem.A))

        return self.linSystem.solve()

    def reset(self):
        self.linSystem.reset(shape=self.depFieldShape)
        self.sourceField_c.fill(0.0)
        self.a_e.fill(0.0)
        self.a_w.fill(0.0)
        self.a_s.fill(0.0)
        self.a_n.fill(0.0)
        self.a_p.fill(0.0)

    def updateFluxes(self):
        self.reset()

        transModel_u = self.u.govModel
        transModel_v = self.v.govModel

        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas()

        ap_u = self.u.govModel.a_p
        ap_v = self.v.govModel.a_p

        # SIMPLE ALGORITHM: comment following SIMPLEC block out.
        #
        # # SIMPLEC ALGORITHM
        # centreCoeffs_u -= transModel_u.getCoefficientsInDirection('east')
        # centreCoeffs_u -= transModel_u.getCoefficientsInDirection('west')
        # centreCoeffs_u -= transModel_u.getCoefficientsInDirection('north')
        # centreCoeffs_u -= transModel_u.getCoefficientsInDirection('south')
        #
        # centreCoeffs_v -= transModel_v.getCoefficientsInDirection('east')
        # centreCoeffs_v -= transModel_v.getCoefficientsInDirection('west')
        # centreCoeffs_v -= transModel_v.getCoefficientsInDirection('north')
        # centreCoeffs_v -= transModel_v.getCoefficientsInDirection('south')

        self.d_u = faceAreas_u/ap_u
        self.d_v = faceAreas_v/ap_v

        self.a_e = self.d_u[east] * faceAreas_u[east]
        self.a_w = self.d_u[west] * faceAreas_u[west]
        self.a_n = self.d_v[north] * faceAreas_v[north]
        self.a_s = self.d_v[south] * faceAreas_v[south]

        # I don't think i have to correct these:
#        self.correctBCs()

        self.a_p += self.a_w + self.a_e + self.a_s + self.a_n

    def updateSourceField(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas()

        u = self.u.data
        v = self.v.data
        self.sourceField_c = u[west] *faceAreas_u[west]\
                             -u[east] *faceAreas_u[east]\
                             + v[south] *faceAreas_v[south] \
                             - v[north] *faceAreas_v[north]

