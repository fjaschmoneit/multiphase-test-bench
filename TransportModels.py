import Differentiation
import DifferenceSchemes
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

        super().__init__(mesh=mesh,
                         field = self.phi,
                         linSystem=linSystem,
                         depFieldShape=MeshConfig.SHAPE_SCALAR_CV)

        # pointer to boundary functions:
        self.boundaryModels = {
            'fixedValue' : BoundaryConditions.scalarBC.derichlet,
            'zeroGradient' : BoundaryConditions.scalarBC.vonNeumann
        }

    def listAvailableBoundaryModels(self):
        for name, type in self.boundaryModels.items():
            print(name, type)

    def linkOtherFields(self):
        self.u = objReg.FIELDS['u']
        self.v = objReg.FIELDS['v']

    def calcConvFlux(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas()

        u = self.u.data
        v = self.v.data

        return(
            u*faceAreas_u,
            v*faceAreas_v
        )

    def calcDiffFlux(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas()
        rCellDist_u, rCellDist_v = self.mesh.calcInvCellDistances()
        D = self.diffusionCoefficient

        return (
            D * faceAreas_u * rCellDist_u,
            D * faceAreas_v * rCellDist_v
        )

    # this includes the advective fluxes when calculating the central coeffs
    def updateFluxes(self):
        self.reset()

        F_u, F_v = self.calcConvFlux()
        D_u, D_v = self.calcDiffFlux()

        self.a_w = DifferenceSchemes.centralDifference( D_u[west], F_u[west], 'west' )
        self.a_e = DifferenceSchemes.centralDifference(D_u[east], F_u[east], 'east')
        self.a_n = DifferenceSchemes.centralDifference(D_v[north], F_v[north], 'north')
        self.a_s = DifferenceSchemes.centralDifference(D_v[south], F_v[south], 'south')

        self.correctBCs()

        self.a_p += ( self.a_w + self.a_e + self.a_s + self.a_n + F_u[east] - F_u[west] + F_v[south] - F_v[north] ) # why do I substract these??


    def updateSourceField(self):
        self.sourceField_c += self.constSourceField



class staggeredTransport_u(TransportBase.transportBase):

    def __init__(self, mesh, linSystem):
        self.u = Fields.newField(shape=MeshConfig.SHAPE_FACES_U)

        super().__init__(mesh=mesh, field=self.u, linSystem=linSystem, depFieldShape=MeshConfig.SHAPE_FACES_U)

        self.boundaryModels = {
            'fixedValue' : BoundaryConditions.staggered_u.derichlet,
            'zeroGradient' : BoundaryConditions.staggered_u.vonNeumann
        }

    def linkOtherFields(self):
        self.v = objReg.FIELDS['v']
        self.p = objReg.FIELDS['p']

    def calcConvFlux(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas()

        u = self.u.data
        v = self.v.data

        return (
            Interpolation.toStaggered(u * faceAreas_u, 'u'),
            Interpolation.toStaggered(v * faceAreas_v, 'u')
        )

    def calcDiffFlux(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas()
        rCellDist_u, rCellDist_v = self.mesh.calcInvCellDistances()
        D = self.diffusionCoefficient

        return(
            Interpolation.toStaggered(D * rCellDist_u * faceAreas_u, 'u'),
            Interpolation.toStaggered(D * rCellDist_v * faceAreas_v, 'u')
        )

    def updateFluxes(self):
        self.reset()

        F_u, F_v = self.calcConvFlux()
        D_u, D_v = self.calcDiffFlux()

        self.a_w = DifferenceSchemes.centralDifference(D_u[west], F_u[west], 'west')
        self.a_e = DifferenceSchemes.centralDifference(D_u[east], F_u[east], 'east')
        self.a_n = DifferenceSchemes.centralDifference(D_v[north], F_v[north], 'north')
        self.a_s = DifferenceSchemes.centralDifference(D_v[south], F_v[south], 'south')

        self.correctBCs()

        self.a_p += (self.a_w + self.a_e + self.a_s + self.a_n )

    def updateSourceField(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas()
        gradP_u = Differentiation.grad_u(self.p.data, self.mesh)
        self.sourceField_c[internal_u] += -gradP_u * faceAreas_u[internal_u]

class staggeredTransport_v(TransportBase.transportBase):

    def __init__(self, mesh, linSystem):
        self.v = Fields.newField(shape=MeshConfig.SHAPE_FACES_V)

        super().__init__(mesh=mesh, field=self.v, linSystem=linSystem, depFieldShape=MeshConfig.SHAPE_FACES_V)

        self.boundaryModels = {
            'fixedValue' : BoundaryConditions.staggered_v.derichlet,
            'zeroGradient' : BoundaryConditions.staggered_v.vonNeumann
        }

    def linkOtherFields(self):
        self.u = objReg.FIELDS['u']
        self.p = objReg.FIELDS['p']

    def calcConvFlux(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas()

        u = self.u.data
        v = self.v.data

        return (
            Interpolation.toStaggered(u * faceAreas_u, 'v'),
            Interpolation.toStaggered(v * faceAreas_v, 'v')
        )

    def calcDiffFlux(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas()
        rCellDist_u, rCellDist_v = self.mesh.calcInvCellDistances()
        D = self.diffusionCoefficient

        return (
            Interpolation.toStaggered(D * rCellDist_u * faceAreas_u, 'v'),
            Interpolation.toStaggered(D * rCellDist_v * faceAreas_v, 'v')
        )

    def updateSourceField(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas()
        gradP_v = Differentiation.grad_v(self.p.data, self.mesh)
        self.sourceField_c[internal_v] += -gradP_v * faceAreas_v[internal_v]


