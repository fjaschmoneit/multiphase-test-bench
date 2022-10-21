import Fields
import BoundaryConditions
import TransportBase
from fieldAccess import *

class Pressure(TransportBase.transportBase):

    def __init__(self, mesh, fieldCreator, fieldReg, linSystem):
        self.depFieldShape = fieldCreator.typeShapeDict['scalarCV']
        super().__init__(mesh, fieldCreator, fieldReg, linSystem, self.depFieldShape)
        self.setSourceField(0.0)
        self.p = self.depField

    def linkOtherFields(self, fieldReg):
        self.u = fieldReg['u']
        self.v = fieldReg['v']

    def updateSourceField(self):
        pass

    def updateFluxes(self):
        return
#        faceAreas = self.mesh.calcFaceAreas(self.fc)
    #     d = self.calcDCoeffs()
    #     self.presFluxes = faceAreas*faceAreas*d

    #def calcDCoeffs(self):
#
#         d = Fields.fieldContainer(
#             u=self.fc.newField(type='faces_u'),
#             v=self.fc.newField(type='faces_v')
#         )
#         #uFlowModel = self.flowmodel_u.getTotalFlux_e
# #not wirking yet
#         uFlowModel = self.u.getGoverningTransportModel()
#         vFlowModel = self.v.getGoverningTransportModel()
#
#         a = uFlowModel.getTotalFlux_e
#         #b = vFlowModel.getTotalFlux().e
#
#         d.u.data.fill(1)
#         d.v.data.fill(0)
#         return d


    # def correctBCs(self):
    #     for loc, typeValueTupel in fieldReg[self.depVariableName].boundary.items():
    #         (type, value) = typeValueTupel
    #         if type == 'zeroGradient':
    #             self.setVonNeumann(loc)
    #         else:
    #             self.setDerichlet(loc, value)

    def solve(self):
        return self.linSystem.solve()

    def updateFluxes(self):
        self.reset()

        transModel_u = self.u.govModel
        transModel_v = self.v.govModel

        faceAreas = self.mesh.calcFaceAreas(self.fc)

        A_u = faceAreas.u.data
        centreCoeffs_u = transModel_u.getCentreMatrixCoeffs()
        d_u = A_u/centreCoeffs_u

        self.a_e = d_u[east]*A_u[east]
        # self.linSystem.set_e_coeffs( a_i )
        # self.a_p += a_i

        self.a_w = d_u[west] * A_u[west]
        # self.linSystem.set_w_coeffs(a_i)
        # a_p+=a_i

        A_v = faceAreas.v.data
        centreCoeffs_v = transModel_v.getCentreMatrixCoeffs()
        d_v = faceAreas.v / centreCoeffs_v

        self.a_n = d_v[north] * A_v[north]
        # self.linSystem.set_n_coeffs(a_i)
        # a_p+=a_i

        self.a_s = d_v[south] * A_v[south]
        # self.linSystem.set_s_coeffs(a_i)
        # a_p += a_i

        self.a_p += self.a_w + self.a_e + self.a_s + self.a_n

        self.linSystem.set_p_coeffs( a_p.data )

    def updateSourceField(self):
        faceAreas = self.mesh.calcFaceAreas(self.fc)
        self.sourceField_c = self.u.west*faceAreas.u.west-self.u.east*faceAreas.u.east \
                              + self.v.south*faceAreas.v.south - self.v.north*faceAreas.v.north

        linLength = self.linSystem.shape[0] * self.linSystem.shape[1]
        # make this also a funtion

        self.linSystem.b = np.reshape(self.sourceField_c.data, linLength)
        #self.linSystem.b

    def setVonNeumann(self, location):
        return
    #     BoundaryConditions.scalarBC.vonNeumann(location, D=self.diffFluxes)   Why do I have to pass these parameters?

    def setDerichlet(self, location, value):
        return
        self.sourceField_c.bw.fill(0)
        #BoundaryConditions.scalarBC.derichlet(location, value, F=self.convFluxes, D=self.diffFluxes,
