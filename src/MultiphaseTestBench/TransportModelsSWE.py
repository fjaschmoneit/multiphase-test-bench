import Differentiation
import Interpolation
import TransportBase
import BoundaryConditions
from fieldAccess import *
import MeshConfig
import ObjectRegistry as objReg
import Fields
import numpy as np

LARGE = 1e20

class staggeredTransport_u(TransportBase.transportBase):

    def __init__(self, mesh, linSystem):
        self.u = Fields.newField(shape=MeshConfig.SHAPE_FACES_U)
        self.grav = 9.81
        self.seabed = None
        super().__init__(mesh=mesh, field=self.u, linSystem=linSystem, depFieldShape=MeshConfig.SHAPE_FACES_U)

        self.boundaryModels = {
            'fixedValue' : BoundaryConditions.staggered_u.derichlet,
            'zeroGradient' : BoundaryConditions.staggered_u.vonNeumann
        }

    def setGravitation(self, value):
        self.grav = value

    def setSeaBed(self, field):
        self.seabed = field

    def linkOtherFields(self, closure):
        if len(closure)==0:
            self.v = objReg.FIELDS['v']
            self.h = objReg.FIELDS['h']
        else:
            self.v = objReg.FIELDS[closure[0]]
            self.h = objReg.FIELDS[closure[1]]

    def updateConvectiveFluxes(self, diffMethod='CDS'):
        faceAreas_ew, faceAreas_sn = self.mesh.calcFaceAreas()
        u = self.u.data
        v = self.v.data
        hu = Interpolation.toStaggered(self.h.data, 'u')
        hv = Interpolation.toStaggered(self.h.data, 'v')

        Mc_ew = Interpolation.toStaggered(hu* u * faceAreas_ew,'u')
        Mc_sn = Interpolation.toStaggered(hv* v * faceAreas_sn,'u')

        if diffMethod=='CDS':
            lmbd = 0.5
            self.a_w += lmbd*Mc_ew[west]
            self.a_e += -lmbd*Mc_ew[east]
            self.a_n += -lmbd*Mc_sn[north]
            self.a_s += lmbd*Mc_sn[south]
        elif diffMethod=='UDS':
            self.a_w += np.maximum(Mc_ew[west],0.0)
            self.a_e += -np.maximum(Mc_ew[east],0.0)
            self.a_n += -np.maximum(Mc_sn[north],0.0)
            self.a_s += np.maximum(Mc_sn[south],0.0)

    def updateFluxes(self):
        self.reset()

        self.updateConvectiveFluxes(diffMethod='CDS')
        self.correctBCs()

        self.a_p += (self.a_w + self.a_e + self.a_s + self.a_n )

    def calcTimeStep(self):
        q_ew = self.u.data*Interpolation.toStaggered(self.h.data, 'u')
        velMax = np.max(np.abs(0.5*(q_ew[east] + q_ew[west])/self.h.data) +np.sqrt(9.81*self.h.data))
        return 0.5*self.mesh.uniformSpacing/velMax

    def updateSourceField(self):
        gradh =  Differentiation.grad_u(self.h.data, self.mesh)
        g = 9.81
        self.sourceField_c[internal_u] -= g*gradh*self.mesh.getCellVolumes()

    def updateTransientCoeffs(self,dt, method):
        h_u = Interpolation.toStaggered(self.h.data, 'u')
        a_p0 = h_u * self.mesh.getCellVolumes() / dt
        if method == 'implicit':
            self.a_p += a_p0
            self.sourceField_c += a_p0 * self.u.data
        elif method == 'explicit':
            self.a_p = a_p0
            self.sourceField_c += (a_p0 + self.a_w + self.a_e + self.a_s + self.a_n) * self.u.data
        elif method == 'semi-implicit':
            hnew = self.h.data
            hold = self.h.govModel.phiOldData
            uold = self.u.data.copy()

            q_ew = uold.copy()
            q_ew[west] = np.maximum(hold*uold[west], np.zeros(hold.shape))
            q_ew[east] += np.minimum(hold*uold[east], np.zeros(hold.shape))
            q_c = 0.5*(q_ew[east] + q_ew[west])
            u_c = np.maximum( uold[west], np.zeros(hold.shape)) + np.minimum( uold[east], np.zeros(hold.shape))

            self.u.data[internal_u] = Interpolation.toStaggered(hold,'u')[internal_u]*uold[internal_u] \
                          -dt/self.mesh.getCellVolumes()*(
                                              q_c[east]*u_c[east] - q_c[west]*u_c[west]
                                              + 0.5*self.grav*(np.power(hnew[east],2) -np.power(hnew[west],2))
                                              +self.grav*Interpolation.toStaggered(hnew,'u')[internal_u]*(self.seabed[east] - self.seabed[west])
                                      )

            self.u.data /= Interpolation.toStaggered(hnew, 'u')

            self.a_p.fill(1)
            self.a_e.fill(0)
            self.a_w.fill(0)
            self.a_n.fill(0)
            self.a_s.fill(0)
            self.sourceField_c = self.u.data

#            print(0.1/np.max( 0.5*np.abs( q_ew[east]+q_ew[west]/hold ) ))


    def solve(self):
        return self.linSystem.solve()