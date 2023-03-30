import numpy as np

import Differentiation
#import DifferenceSchemes
import Interpolation
import TransportBase
import BoundaryConditions
from fieldAccess import *
import MeshConfig
import ObjectRegistry as objReg
import Fields


class scalarTransport(TransportBase.transportBase):

    def __init__(self, mesh, linSystem ):

        self.phi = Fields.newField(shape=MeshConfig.SHAPE_SCALAR_CV)
        self.phiOldData = self.phi.data.copy()
        self.heatCapacity = None

        super().__init__(mesh=mesh,
                         field = self.phi,
                         linSystem=linSystem,
                         depFieldShape=MeshConfig.SHAPE_SCALAR_CV)

        # pointer to boundary functions:
        self.boundaryModels = {
            'fixedValue' : BoundaryConditions.scalarBC.derichlet,
            'zeroGradient' : BoundaryConditions.scalarBC.vonNeumann
        }

    def setHeatCapacity(self, value):
        self.heatCapacity = value

    def listAvailableBoundaryModels(self):
        for name, type in self.boundaryModels.items():
            print(name, type)

    def linkOtherFields(self, closure):
        if len(closure)==0:
            self.v = objReg.FIELDS['v']
            self.u = objReg.FIELDS['u']
        else:
            self.v = objReg.FIELDS[closure[1]]
            self.u = objReg.FIELDS[closure[0]]

    def updateDiffusiveFluxes(self):
        rCellDist_ew, rCellDist_sn = self.mesh.calcInvCellDistances()
        faceAreas_ew, faceAreas_sn = self.mesh.calcFaceAreas()

        D = self.diffusionCoefficient

        # mass flow
        Md_ew = D*rCellDist_ew * faceAreas_ew
        Md_sn = D*rCellDist_sn * faceAreas_sn

        self.a_e += Md_ew[east]
        self.a_w += Md_ew[west]
        self.a_n += Md_sn[north]
        self.a_s += Md_sn[south]

    def updateConvectiveFluxes(self, diffMethod='UDS'):
        faceAreas_ew, faceAreas_sn = self.mesh.calcFaceAreas()
        u = self.u.data
        v = self.v.data

        Mc_ew = u * faceAreas_ew
        Mc_sn = v * faceAreas_sn

        if diffMethod=='CDS':
            lmbd = 0.5      # mesh dependent - for structured, regular cells, lambda = 0.5
            self.a_w += lmbd*Mc_ew[west]
            self.a_e += -lmbd*Mc_ew[east]
            self.a_n += -lmbd*Mc_sn[north]
            self.a_s += lmbd*Mc_sn[south]
        elif diffMethod=='UDS':
            self.a_w += np.maximum(Mc_ew[west],0.0)
            self.a_e += -np.minimum(Mc_ew[east],0.0)
            self.a_n += -np.minimum(Mc_sn[north],0.0)
            self.a_s += np.maximum(Mc_sn[south],0.0)

    def updateFluxes(self):
        self.reset()

        self.updateDiffusiveFluxes()
        self.updateConvectiveFluxes('CDS')
        self.correctBCs()

        self.a_p += (self.a_w + self.a_e + self.a_s + self.a_n)

    def updateSourceField(self):
        self.sourceField_c += self.constSourceField

    def updateTransientCoeffs(self, dt, method):
        if method == 'implicit':
            a_p0 = self.heatCapacity*self.mesh.getCellVolumes()/dt
            self.a_p += a_p0
            self.sourceField_c += a_p0*self.phi.data
        elif method == 'explicit':
            # old height data
            self.phiOldData = self.phi.data.copy()
            u = self.u.data
#            h_ew = np.maximum( phi[west]*u[internal_u], 0.0 ) + np.maximum( -phi[east]*u[internal_u], 0.0)
            q_ew = u.copy()
            q_ew[west] = np.maximum(self.phiOldData*u[west], np.zeros(self.phiOldData.shape))
            q_ew[east] += np.minimum(self.phiOldData*u[east], np.zeros(self.phiOldData.shape))

            correction = dt/self.mesh.getCellVolumes()*(q_ew[east] - q_ew[west])

            self.phi.data = self.phiOldData - correction

            # i'm not solving a linear eq system,
            # so I just adjust my coefficients such that they
            # return my already evaluated phi.

            self.a_p.fill(1)
            self.a_e.fill(0)
            self.a_w.fill(0)
            self.a_n.fill(0)
            self.a_s.fill(0)
            self.sourceField_c = self.phi.data


    def solve(self):
        return self.linSystem.solve()


class staggeredTransport_u(TransportBase.transportBase):

    def __init__(self, mesh, linSystem):
        self.u = Fields.newField(shape=MeshConfig.SHAPE_FACES_U)

        super().__init__(mesh=mesh, field=self.u, linSystem=linSystem, depFieldShape=MeshConfig.SHAPE_FACES_U)

        self.boundaryModels = {
            'fixedValue' : BoundaryConditions.staggered_u.derichlet,
            'zeroGradient' : BoundaryConditions.staggered_u.vonNeumann
        }

    def linkOtherFields(self, closure):
        if len(closure)==0:
            self.v = objReg.FIELDS['v']
            self.p = objReg.FIELDS['p']
        else:
            self.v = objReg.FIELDS[closure[0]]
            self.p = objReg.FIELDS[closure[1]]

    def updateDiffusiveFluxes(self):
        rCellDist_ew, rCellDist_sn = self.mesh.calcInvCellDistances()
        faceAreas_ew, faceAreas_sn = self.mesh.calcFaceAreas()

        D = self.diffusionCoefficient

        # mass flow
        Md_ew =  Interpolation.toStaggered(D*rCellDist_ew * faceAreas_ew,'u')
        Md_sn =  Interpolation.toStaggered(D*rCellDist_sn * faceAreas_sn,'u')

        self.a_e += Md_ew[east]
        self.a_w += Md_ew[west]
        self.a_n += Md_sn[north]
        self.a_s += Md_sn[south]

    def updateConvectiveFluxes(self, diffMethod='CDS'):
        faceAreas_ew, faceAreas_sn = self.mesh.calcFaceAreas()
        u = self.u.data
        v = self.v.data

        Mc_ew = Interpolation.toStaggered(u * faceAreas_ew,'u')
        Mc_sn = Interpolation.toStaggered(v * faceAreas_sn,'u')

        if diffMethod=='CDS':
            lmbd = 0.5      # mesh dependent - for structured, regular cells, lambda = 0.5
            self.a_w += lmbd*Mc_ew[west]
            self.a_e += -lmbd*Mc_ew[east]
            self.a_n += -lmbd*Mc_sn[north]
            self.a_s += lmbd*Mc_sn[south]
        elif diffMethod=='UDS':
            self.a_w += np.maximum(Mc_ew[west],0.0)
            self.a_e += -np.minimum(Mc_ew[east],0.0)
            self.a_n += -np.minimum(Mc_sn[north],0.0)
            self.a_s += np.maximum(Mc_sn[south],0.0)

    def updateFluxes(self):
        self.reset()

        self.updateDiffusiveFluxes()
        self.updateConvectiveFluxes()
        self.correctBCs()

        self.a_p += self.a_w + self.a_e + self.a_s + self.a_n

    def updateSourceField(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas()
        gradP_u = Differentiation.grad_u(self.p.data, self.mesh)
        self.sourceField_c[internal_u] += -gradP_u * faceAreas_u[internal_u]

    def solve(self):
        return self.linSystem.solve()

class staggeredTransport_v(TransportBase.transportBase):

    def __init__(self, mesh, linSystem):
        self.v = Fields.newField(shape=MeshConfig.SHAPE_FACES_V)

        super().__init__(mesh=mesh, field=self.v, linSystem=linSystem, depFieldShape=MeshConfig.SHAPE_FACES_V)

        self.boundaryModels = {
            'fixedValue' : BoundaryConditions.staggered_v.derichlet,
            'zeroGradient' : BoundaryConditions.staggered_v.vonNeumann
        }

    def linkOtherFields(self, closure):
        if len(closure)==0:
            self.u = objReg.FIELDS['u']
            self.p = objReg.FIELDS['p']
        else:
            self.u = objReg.FIELDS[closure[0]]
            self.p = objReg.FIELDS[closure[1]]

    def updateSourceField(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas()
        gradP_v = Differentiation.grad_v(self.p.data, self.mesh)
        self.sourceField_c[internal_v] += -gradP_v * faceAreas_v[internal_v]

    def updateDiffusiveFluxes(self):
        rCellDist_ew, rCellDist_sn = self.mesh.calcInvCellDistances()
        faceAreas_ew, faceAreas_sn = self.mesh.calcFaceAreas()

        D = self.diffusionCoefficient

        # mass flow
        Md_ew =  Interpolation.toStaggered(D*rCellDist_ew * faceAreas_ew,'v')
        Md_sn =  Interpolation.toStaggered(D*rCellDist_sn * faceAreas_sn,'v')

        self.a_e += Md_ew[east]
        self.a_w += Md_ew[west]
        self.a_n += Md_sn[north]
        self.a_s += Md_sn[south]

    def updateConvectiveFluxes(self, diffMethod='UDS'):
        faceAreas_ew, faceAreas_sn = self.mesh.calcFaceAreas()
        u = self.u.data
        v = self.v.data

        Mc_ew = Interpolation.toStaggered(u * faceAreas_ew,'v')
        Mc_sn = Interpolation.toStaggered(v * faceAreas_sn,'v')

        if diffMethod=='CDS':
            lmbd = 0.5
            self.a_w += lmbd*Mc_ew[west]
            self.a_e += -lmbd*Mc_ew[east]
            self.a_n += -lmbd*Mc_sn[north]
            self.a_s += lmbd*Mc_sn[south]
        elif diffMethod=='UDS':
            self.a_w += np.maximum(Mc_ew[west],0.0)
            self.a_e += -np.minimum(Mc_ew[east],0.0)
            self.a_n += -np.minimum(Mc_sn[north],0.0)
            self.a_s += np.maximum(Mc_sn[south],0.0)


    def updateFluxes(self):
        self.reset()

        self.updateDiffusiveFluxes()
        self.updateConvectiveFluxes()
        self.correctBCs()

        self.a_p += (self.a_w + self.a_e + self.a_s + self.a_n )

    def solve(self):
        return self.linSystem.solve()

