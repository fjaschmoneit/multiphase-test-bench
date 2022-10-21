import Differentiation
import DifferenceSchemes
import Fields
import Interpolation
import TransportBase
import BoundaryConditions
from fieldAccess import *

class scalarTransport(TransportBase.transportBase):

    def __init__(self, mesh, fieldCreator, fieldReg, linSystem ):
        self._depFieldShape = fieldCreator._typeShapeDict['scalarCV']
        super().__init__(mesh, fieldCreator, fieldReg, linSystem, self._depFieldShape)
        self.setSourceField(0.0)

    def linkOtherFields(self, fieldReg):
        self.u = fieldReg['u']
        self.v = fieldReg['v']

    def calcConvFlux(self, vel, faceAreas):
        return vel*faceAreas

    def calcDiffFlux(self,invCellDist,diffCoeff,faceAreas):
        return invCellDist*diffCoeff*faceAreas

    def updateFluxes(self):
        self.reset()

        faceAreas = self.mesh.calcFaceAreas(self.fc)
        faceAreas_u = faceAreas.u.data
        faceAreas_v = faceAreas.v.data

        rCellDist_u = self.fieldReg['invCellDist'].u.data
        rCellDist_v = self.fieldReg['invCellDist'].v.data

        F_u = self.calcConvFlux(self.u.data, faceAreas_u)
        F_v = self.calcConvFlux(self.v.data, faceAreas_v)

        D_u = self.calcDiffFlux(rCellDist_u, self.diffusionCoefficient, faceAreas_u)
        D_v = self.calcDiffFlux(rCellDist_v, self.diffusionCoefficient, faceAreas_v)

        self.a_w = DifferenceSchemes.centralDifference( D_u[west], F_u[west], 'west' )
        self.a_e = DifferenceSchemes.centralDifference(D_u[east], F_u[east], 'east')
        self.a_n = DifferenceSchemes.centralDifference(D_v[north], F_v[north], 'north')
        self.a_s = DifferenceSchemes.centralDifference(D_v[south], F_v[south], 'south')

        self.correctBCs()

        self.a_p += ( self.a_w + self.a_e + self.a_s + self.a_n + F_u[east] - F_u[west] + F_v[south] - F_v[north] )

    def setVonNeumann(self, loc):
        BoundaryConditions.scalarBC.vonNeumann(loc=loc, transportInstance = self)

    def setDerichlet(self, loc, value):
        BoundaryConditions.scalarBC.derichlet( loc=loc, value=value, transportInstance=self )

    def updateSourceField(self):
        pass

class staggeredTransport_u(TransportBase.transportBase):

    def __init__(self, mesh, fieldCreator, fieldReg, linSystem):
        self.depFieldShape = fieldCreator.typeShapeDict['faces_u']
        super().__init__(mesh, fieldCreator, fieldReg, linSystem, self.depFieldShape)
        self.u = self.depField
        self.setSourceField(0.0)

    def linkOtherFields(self, fieldReg):
        self.v = fieldReg['v']
        self.p = fieldReg['p']

    def calcConvFlux(self, vel, faceAreas, dir):
        return Interpolation.toStaggered(vel*faceAreas, dir)

    def calcDiffFlux(self, invCellDist, diffCoeff, faceAreas, dir):
        return Interpolation.toStaggered(diffCoeff*invCellDist*faceAreas, dir)

    def updateFluxes(self):
        self.reset()

        faceAreas = self.mesh.calcFaceAreas(self.fc)
        faceAreas_u = faceAreas.u.data
        faceAreas_v = faceAreas.v.data

        rCellDist_u = self.fieldReg['invCellDist'].u.data
        rCellDist_v = self.fieldReg['invCellDist'].v.data

        F_u = self.calcConvFlux(self.u.data, faceAreas_u, 'u')
        F_v = self.calcConvFlux(self.v.data, faceAreas_v, 'u')

        D_u = self.calcDiffFlux(rCellDist_u, self.diffusionCoefficient, faceAreas_u, 'u')
        D_v = self.calcDiffFlux(rCellDist_v, self.diffusionCoefficient, faceAreas_v, 'u')

        self.a_w[east] = DifferenceSchemes.centralDifference(D_u, F_u, 'west' )[internal_u]
        self.a_e[west] = DifferenceSchemes.centralDifference(D_u, F_u, 'east')[internal_u]

        self.a_n = DifferenceSchemes.centralDifference(D_v, F_v, 'north')[north]
        self.a_s = DifferenceSchemes.centralDifference(D_v, F_v, 'south')[south]

        self.correctBCs()

        self.a_p += self.a_w + self.a_e + self.a_s + self.a_n
        #self.a_p[internal_u] += F_u[east] - F_u[west] + F_v[south] - F_v[north]    # do I need to include these terms when u itself is transported?

    def setVonNeumann(self, location):
        if location == 'left':
            self.a_w[boundary_west] = 0
        elif location == 'right':
            self.a_e[boundary_east] = 0
        elif location == 'top':
            self.a_n[boundary_north] = 0
        elif location == 'bottom':
            self.a_s[boundary_south] = 0

    def setDerichlet(self, location, value):
        # BoundaryConditions.staggeredBC.derichlet(location, value, self)

        # in flow direction:
        if location == 'left':
            self.a_p[boundary_west] += 1e15
            self.sourceField_c[boundary_west] += 1e15 * value
            self.a_w[boundary_west] = 0.0
        elif location == 'right':
            self.a_p[boundary_east] += 1e15
            self.sourceField_c[boundary_east] += 1e15 * value
            self.a_e[boundary_east] = 0.0

        # across flow direction
        elif location == 'top':
            faceAreas = self.mesh.calcFaceAreas(self.fc)
            faceAreas_v = faceAreas.v.data

            invCellWallDist_A = 2.0 * self.fieldReg['invCellDist'].v.data[boundary_north]*faceAreas_v[boundary_north]
            nodeFluxes = Interpolation.toStaggered(invCellWallDist_A, 'u')

            self.a_p[boundary_north] += self.diffusionCoefficient*nodeFluxes
            self.sourceField_c[boundary_north] += self.diffusionCoefficient*nodeFluxes*value
            self.a_n[boundary_north] = 0

        elif location == 'bottom':
            faceAreas = self.mesh.calcFaceAreas(self.fc)
            faceAreas_v = faceAreas.v.data

            invCellWallDist_A = 2.0 * self.fieldReg['invCellDist'].v.data[boundary_south]*faceAreas_v[boundary_south]
            nodeFluxes = Interpolation.toStaggered(invCellWallDist_A, 'u')

            self.a_p[boundary_south] += self.diffusionCoefficient*nodeFluxes
            self.sourceField_c[boundary_south] += self.diffusionCoefficient*nodeFluxes*value
            self.a_s[boundary_south] = 0

    def updateSourceField(self):
        faceAreas_u = self.mesh.calcFaceAreas(self.fc).u.data
        gradP_u = Differentiation.grad_u(self.p.data, self.fieldReg)
        self.sourceField_c[internal_u] = -gradP_u * faceAreas_u[internal_u]

class staggeredTransport_v(TransportBase.transportBase):

    def __init__(self, mesh, fieldCreator, fieldReg, linSystem):
        self.depFieldShape = fieldCreator.typeShapeDict['faces_v']
        super().__init__(mesh, fieldCreator, fieldReg, linSystem, self.depFieldShape)
        self.v = self.depField
        self.setSourceField(0.0)

    def linkOtherFields(self, fieldReg):
        self.u = fieldReg['u']
        self.p = fieldReg['p']

    def calcConvFlux(self, vel, faceAreas, dir):
        return Interpolation.toStaggered(vel*faceAreas, dir)

    def calcDiffFlux(self, invCellDist, diffCoeff, faceAreas, dir):
        return Interpolation.toStaggered(diffCoeff*invCellDist*faceAreas, dir)

    def updateFluxes(self):
        self.reset()

        faceAreas = self.mesh.calcFaceAreas(self.fc)
        faceAreas_u = faceAreas.u.data
        faceAreas_v = faceAreas.v.data

        rCellDist_u = self.fieldReg['invCellDist'].u.data
        rCellDist_v = self.fieldReg['invCellDist'].v.data

        F_u = self.calcConvFlux(self.u.data, faceAreas_u, 'v')
        F_v = self.calcConvFlux(self.v.data, faceAreas_v, 'v')

        D_u = self.calcDiffFlux(rCellDist_u, self.diffusionCoefficient, faceAreas_u, 'v')
        D_v = self.calcDiffFlux(rCellDist_v, self.diffusionCoefficient, faceAreas_v, 'v')

        self.a_n[south] = DifferenceSchemes.centralDifference(D_v, F_v, 'north')[internal_v]
        self.a_s[north] = DifferenceSchemes.centralDifference(D_v, F_v, 'south')[internal_v]

        self.a_e = DifferenceSchemes.centralDifference(D_u, F_u, 'east')[east]
        self.a_w = DifferenceSchemes.centralDifference(D_u, F_u, 'west')[west]

        self.correctBCs()

        self.a_p += self.a_w + self.a_e + self.a_s + self.a_n
        #self.a_p[internal_u] += F_u[east] - F_u[west] + F_v[south] - F_v[north]    # do I need to include these terms when u itself is transported?

    def setVonNeumann(self, location):
        if location == 'left':
            self.a_w[boundary_west] = 0
        elif location == 'right':
            self.a_e[boundary_east] = 0
        elif location == 'top':
            self.a_n[boundary_north] = 0
        elif location == 'bottom':
            self.a_s[boundary_south] = 0

    def setDerichlet(self, location, value):
        # BoundaryConditions.staggeredBC.derichlet(location, value, self)

        # in flow direction:
        if location == 'top':
            self.a_p[boundary_north] += 1e15
            self.sourceField_c[boundary_north] += 1e15 * value
            self.a_n[boundary_north] = 0.0
        elif location == 'bottom':
            self.a_p[boundary_south] += 1e15
            self.sourceField_c[boundary_south] += 1e15 * value
            self.a_s[boundary_south] = 0.0

        # across flow direction
        elif location == 'left':
            faceAreas = self.mesh.calcFaceAreas(self.fc)
            faceAreas_u = faceAreas.u.data

            invCellWallDist_A = 2.0 * self.fieldReg['invCellDist'].u.data[boundary_west]*faceAreas_u[boundary_west]
            nodeFluxes = Interpolation.toStaggered(invCellWallDist_A, 'v')

            self.a_p[boundary_west] += self.diffusionCoefficient*nodeFluxes
            self.sourceField_c[boundary_west] += self.diffusionCoefficient*nodeFluxes*value
            self.a_w[boundary_west] = 0

        elif location == 'right':
            faceAreas = self.mesh.calcFaceAreas(self.fc)
            faceAreas_u = faceAreas.u.data

            invCellWallDist_A = 2.0 * self.fieldReg['invCellDist'].u.data[boundary_east]*faceAreas_u[boundary_east]
            nodeFluxes = Interpolation.toStaggered(invCellWallDist_A, 'v')

            self.a_p[boundary_east] += self.diffusionCoefficient*nodeFluxes
            self.sourceField_c[boundary_east] += self.diffusionCoefficient*nodeFluxes*value
            self.a_e[boundary_east] = 0

    def updateSourceField(self):
        faceAreas_v = self.mesh.calcFaceAreas(self.fc).v.data
        gradP_v = Differentiation.grad_v(self.p.data, self.fieldReg)
        self.sourceField_c[internal_v] = -gradP_v * faceAreas_v[internal_v]


