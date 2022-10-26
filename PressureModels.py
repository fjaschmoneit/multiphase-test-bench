import BoundaryConditions
import TransportBase
from fieldAccess import *

class Pressure(TransportBase.transportBase):

    def __init__(self, mesh, fieldCreator, fieldReg, linSystem):
        self.depFieldShape = fieldCreator.typeShapeDict['scalarCV']
        super().__init__(mesh, fieldCreator, fieldReg, linSystem, self.depFieldShape)
        self.p = self.depField

        self.boundaryModels = {
            'freeFlow' : BoundaryConditions.pressure.freeFlow,
            'totalPressure' : BoundaryConditions.pressure.totalPressure,
            'constantPressure' : BoundaryConditions.pressure.derichlet
        }

    def linkOtherFields(self, fieldReg):
        self.u = fieldReg['u']
        self.v = fieldReg['v']

    def solve(self):
        return self.linSystem.solve()

    def updateFluxes(self):
        self.reset()

        transModel_u = self.u.govModel
        transModel_v = self.v.govModel

        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas(self.fc)

        centreCoeffs_u = transModel_u.getCentreMatrixCoeffs()
        centreCoeffs_v = transModel_v.getCentreMatrixCoeffs()

        self.d_u = faceAreas_u/centreCoeffs_u
        self.d_v = faceAreas_v/centreCoeffs_v

        self.a_e = self.d_u[east] * faceAreas_u[east]
        self.a_w = self.d_u[west] * faceAreas_u[west]
        self.a_n = self.d_v[north] * faceAreas_v[north]
        self.a_s = self.d_v[south] * faceAreas_v[south]

        self.correctBCs()

        self.a_p += self.a_w + self.a_e + self.a_s + self.a_n

    def updateSourceField(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas(self.fc)

        u = self.u.data
        v = self.v.data
        self.sourceField_c = u[west] *faceAreas_u[west]\
                             -u[east] *faceAreas_u[east]\
                             + v[south] *faceAreas_v[south] \
                             - v[north] *faceAreas_v[north]