import Differentiation
import DifferenceSchemes
import Interpolation
import TransportBase
import BoundaryConditions
from fieldAccess import *

class scalarTransport(TransportBase.transportBase):

    def __init__(self, mesh, fieldCreator, fieldReg, linSystem ):
        self._depFieldShape = fieldCreator.typeShapeDict['scalarCV']
        super().__init__(mesh, fieldCreator, fieldReg, linSystem, self._depFieldShape)

        self.boundaryModels = {
            'fixedValue' : BoundaryConditions.scalarBC.derichlet,
            'zeroGradient' : BoundaryConditions.scalarBC.vonNeumann
        }

    def listAvailableBoundaryModels(self):
        for name, type in self.boundaryModels.items():
            print(name, type)

    def linkOtherFields(self, fieldReg):
        self.u = fieldReg['u']
        self.v = fieldReg['v']

    def calcConvFlux(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas(self.fc)

        u = self.u.data
        v = self.v.data

        return(
            u*faceAreas_u,
            v*faceAreas_v
        )

    def calcDiffFlux(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas(self.fc)
        rCellDist_u, rCellDist_v = self.fieldReg['invCellDist']
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

        self.a_p += ( self.a_w + self.a_e + self.a_s + self.a_n + F_u[east] - F_u[west] + F_v[south] - F_v[north] )

    def setVonNeumann(self, direction):
        BoundaryConditions.scalarBC.vonNeumann(direction=direction, transportInstance = self)

    def setDerichlet(self, direction, value):
        BoundaryConditions.scalarBC.derichlet( direction=direction, value=value, transportInstance=self )

    def updateSourceField(self):
        self.sourceField_c += self.constSourceField

class staggeredTransport_u(TransportBase.transportBase):

    def __init__(self, mesh, fieldCreator, fieldReg, linSystem):
        self.depFieldShape = fieldCreator.typeShapeDict['faces_u']
        super().__init__(mesh, fieldCreator, fieldReg, linSystem, self.depFieldShape)
        self.u = self.depField

        self.boundaryModels = {
            'fixedValue' : BoundaryConditions.staggered_u.derichlet,
            'zeroGradient' : BoundaryConditions.staggered_u.vonNeumann
        }

    def linkOtherFields(self, fieldReg):
        self.v = fieldReg['v']
        self.p = fieldReg['p']

    def calcConvFlux(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas(self.fc)

        u = self.u.data
        v = self.v.data

        return (
            Interpolation.toStaggered(u * faceAreas_u, 'u'),
            Interpolation.toStaggered(v * faceAreas_v, 'u')
        )

    def calcDiffFlux(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas(self.fc)
        rCellDist_u, rCellDist_v = self.fieldReg['invCellDist']
        D = self.diffusionCoefficient

        return(
            Interpolation.toStaggered(D * rCellDist_u * faceAreas_u, 'u'),
            Interpolation.toStaggered(D * rCellDist_v * faceAreas_v, 'u')
        )

    def updateSourceField(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas(self.fc)
        gradP_u = Differentiation.grad_u(self.p.data, self.fieldReg)
        self.sourceField_c[internal_u] += -gradP_u * faceAreas_u[internal_u]

class staggeredTransport_v(TransportBase.transportBase):

    def __init__(self, mesh, fieldCreator, fieldReg, linSystem):
        self.depFieldShape = fieldCreator.typeShapeDict['faces_v']
        super().__init__(mesh, fieldCreator, fieldReg, linSystem, self.depFieldShape)
        self.v = self.depField

        self.boundaryModels = {
            'fixedValue' : BoundaryConditions.staggered_v.derichlet,
            'zeroGradient' : BoundaryConditions.staggered_v.vonNeumann
        }

    def linkOtherFields(self, fieldReg):
        self.u = fieldReg['u']
        self.p = fieldReg['p']

    def calcConvFlux(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas(self.fc)

        u = self.u.data
        v = self.v.data

        return (
            Interpolation.toStaggered(u * faceAreas_u, 'v'),
            Interpolation.toStaggered(v * faceAreas_v, 'v')
        )

    def calcDiffFlux(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas(self.fc)
        rCellDist_u, rCellDist_v = self.fieldReg['invCellDist']
        D = self.diffusionCoefficient

        return (
            Interpolation.toStaggered(D * rCellDist_u * faceAreas_u, 'v'),
            Interpolation.toStaggered(D * rCellDist_v * faceAreas_v, 'v')
        )

    def updateSourceField(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas(self.fc)

        # faceAreas_v = self.mesh.calcFaceAreas(self.fc).v.data
        gradP_v = Differentiation.grad_v(self.p.data, self.fieldReg)
        self.sourceField_c[internal_v] += -gradP_v * faceAreas_v[internal_v]


