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

        #self._depField.setGoverningTransportModel(self)
        #self._phi = self._depField      # defining my local name
        self.setSourceField(0.0)

    def linkOtherFields(self, fieldReg):
        self._u = fieldReg['u']
        self._v = fieldReg['v']


    def calcConvFlux(self, vel, faceAreas):
        return vel*faceAreas

    def calcDiffFlux(self,invCellDist,diffCoeff,faceAreas):
        return invCellDist*diffCoeff*faceAreas

    def updateFluxes(self):

        faceAreas = self._mesh.calcFaceAreas(self._fc)
        faceAreas_u = faceAreas.u.data
        faceAreas_v = faceAreas.v.data

        rCellDist_u = self._fieldReg['invCellDist'].u.data
        rCellDist_v = self._fieldReg['invCellDist'].v.data

        F_u = self.calcConvFlux(self._u.data, faceAreas_u)
        F_v = self.calcConvFlux(self._v.data, faceAreas_v)

        D_u = self.calcDiffFlux(rCellDist_u, self._diffusionCoefficient, faceAreas_u)
        D_v = self.calcDiffFlux(rCellDist_v, self._diffusionCoefficient, faceAreas_v)

        self.a_w = DifferenceSchemes.centralDifference( D_u[west], F_u[west], 'west' )
        self.a_e = DifferenceSchemes.centralDifference(D_u[east], F_u[east], 'east')
        self.a_n = DifferenceSchemes.centralDifference(D_v[north], F_v[north], 'north')
        self.a_s = DifferenceSchemes.centralDifference(D_v[south], F_v[south], 'south')

        self.correctBCs()

        self.a_p += ( self.a_w + self.a_e + self.a_s + self.a_n + F_u[east] - F_u[west] + F_v[south] - F_v[north] )

        # self.a_p[internal_u] = ( self.a_w+self.a_e+( F_u[east] - F_u[west]) )[internal_u]
        # self.a_p[internal_v]+= ( self.a_s+self.a_n+( F_v[south] - F_v[north]))[internal_v]
        # #self.a_p += self._source

    def setVonNeumann(self, loc):
        if loc == 'left':
            self.a_w[boundary_west].fill(0)
        elif loc == 'right':
            self.a_e[boundary_east].fill(0)
        elif loc == 'top':
            self.a_n[boundary_north].fill(0)
        elif loc == 'bottom':
            self.a_s[boundary_south].fill(0)

        #BoundaryConditions.scalarBC.vonNeumann( location, D=self._diffFluxes )

    def setDerichlet(self, loc, value):
        if loc == 'left':
            westFlux = self.a_w[boundary_west]
            self.a_p[boundary_west] += 2*westFlux
            self._sourceField_c[boundary_west] += 2*westFlux*value
            self.a_w[boundary_west] = 0.0

        elif loc == 'right':
            eastFlux = self.a_e[boundary_east]
            self.a_p[boundary_east] += 2*eastFlux
            self._sourceField_c[boundary_east] += 2*eastFlux*value
            self.a_e[boundary_east] = 0.0

        elif loc == 'top':
            northFlux = self.a_n[boundary_north]
            self.a_p[boundary_north] += 2*northFlux
            self._sourceField_c[boundary_north] += 2*northFlux*value
            self.a_n[boundary_north] = 0.0

        elif loc == 'bottom':
            southFlux = self.a_s[boundary_south]
            self.a_p[boundary_south] += 2*southFlux
            self._sourceField_c[boundary_south] += 2*southFlux*value
            self.a_s[boundary_south] = 0.0

        #BoundaryConditions.scalarBC.derichlet( location, value, F=self._convFluxes, D=self._diffFluxes, Sc=self._sourceField_c, Sp=self._sourceField_p )

    def updateLinSystem(self):
        self._linSystem.reset(shape=self._depFieldShape)
#        self._linSystem.reset()
        self._linSystem.set_e_coeffs(self.a_e)
        self._linSystem.set_w_coeffs(self.a_w)
        self._linSystem.set_n_coeffs(self.a_n)
        self._linSystem.set_s_coeffs(self.a_s)
        self._linSystem.set_p_coeffs(self.a_p)

        self._linSystem.set_b(self._sourceField_c)

        # self._centreMatrixCoeffs = self._linSystem.update(F=self._convFluxes, D=self._diffFluxes, Sc=self._sourceField_c, Sp=self._sourceField_p)

    def updateSourceField(self):
        pass

class staggeredTransport_u(TransportBase.transportBase):

    def __init__(self, mesh, fieldCreator, fieldReg, linSystem):
        self._depFieldShape = fieldCreator._typeShapeDict['faces_u']
        super().__init__(mesh, fieldCreator, fieldReg, linSystem, self._depFieldShape)

        self._depField.setGoverningTransportModel(self)

        self._u = self._depField
        self.setSourceField(0.0)

    def linkOtherFields(self, fieldReg):
        self._v = fieldReg['v']
        self._p = fieldReg['p']

    def updateFluxes(self):
        rCellDist_u = self._fieldReg['invCellDist'].u
        rCellDist_v = self._fieldReg['invCellDist'].v

        self._convFluxes.u.internal_u = 0.5 * (self._u.east + self._u.west)
        self._convFluxes.v.internal_u = 0.5 * (self._v.east + self._v.west)

        self._diffFluxes.u.internal_u = 0.5 * (rCellDist_u.east + rCellDist_u.west) * self._diffusionCoefficient
        self._diffFluxes.v.internal_u = 0.5 * (rCellDist_v.east + rCellDist_v.west) * self._diffusionCoefficient

    def updateSourceField(self):
        faceAreas_u = self._mesh.calcFaceAreas(self._fc).u
        gradP_u = Differentiation.grad_u(self._p, self._fieldReg)
        self._sourceField_c.internal_u = -gradP_u * faceAreas_u.internal_u

    def setVonNeumann(self, location):
        BoundaryConditions.staggeredBC.vonNeumann(location, F=self._convFluxes, D=self._diffFluxes,
                                              Sc=self._sourceField_c, Sp=self._sourceField_p)

    def setDerichlet(self, location, value):
        BoundaryConditions.staggeredBC.derichlet(location, value, F=self._convFluxes, D=self._diffFluxes,
                                              Sc=self._sourceField_c, Sp=self._sourceField_p)
    # implement difference schemes here
    def getTotalFlux_e(self):
        return -(self._diffFluxes.u.east - 0.5*self._convFluxes.u.east)

    def getTotalFlux_w(self):
        return -(self._diffFluxes.u.west + 0.5 * self._convFluxes.u.west)

    # def getTotalFlux_n(self):
    #     return -(self._diffFluxes.v.north - 0.5 * self._convFluxes.v.north)
    #
    # def getTotalFlux_s(self):
    #     return -(self._diffFluxes.v.south - 0.5 * self._convFluxes.v.south)

class staggeredTransport_v(TransportBase.transportBase):

    def __init__(self, mesh, fieldCreator, fieldReg, linSystem):
        self._depFieldShape = fieldCreator._typeShapeDict['faces_v']
        super().__init__(mesh, fieldCreator, fieldReg, linSystem, self._depFieldShape)

        self._depField.setGoverningTransportModel(self)

        self._v = self._depField
        self.setSourceField(0.0)

    def linkOtherFields(self, fieldReg):
        self._u = fieldReg['u']
        self._p = fieldReg['p']

        # this is the most important method. The others could possibly be base class methods
    def updateFluxes(self):
        rCellDist_u = self._fieldReg['invCellDist'].u
        rCellDist_v = self._fieldReg['invCellDist'].v

        self._convFluxes.u.internal_v = 0.5 * (self._u.north + self._u.south)
        self._convFluxes.v.internal_v = 0.5 * (self._v.north + self._v.south)

        self._diffFluxes.u.internal_v = 0.5 * (rCellDist_u.north + rCellDist_u.south) * self._diffusionCoefficient
        self._diffFluxes.v.internal_v = 0.5 * (rCellDist_v.north + rCellDist_v.south) * self._diffusionCoefficient

    def updateSourceField(self):
        faceAreas_v = self._mesh.calcFaceAreas(self._fc).v
        gradP_v = Differentiation.grad_v(self._p, self._fieldReg)
        self._sourceField_c.internal_v = -gradP_v*faceAreas_v.internal_v

    def setVonNeumann(self, location):
        BoundaryConditions.staggeredBC.vonNeumann(location, F=self._convFluxes, D=self._diffFluxes,
                                              Sc=self._sourceField_c, Sp=self._sourceField_p)

    def setDerichlet(self, location, value):
        BoundaryConditions.staggeredBC.derichlet(location, value, F=self._convFluxes, D=self._diffFluxes,
                                              Sc=self._sourceField_c, Sp=self._sourceField_p)

    # implement difference schemes here
    # def getTotalFlux_e(self):
    #     return -(self._diffFluxes.u.east - 0.5*self._convFluxes.u.east)
    #
    # def getTotalFlux_w(self):
    #     return -(self._diffFluxes.u.west + 0.5 * self._convFluxes.u.west)

    def getTotalFlux_n(self):
        return -(self._diffFluxes.v.north - 0.5 * self._convFluxes.v.north)

    def getTotalFlux_s(self):
        return -(self._diffFluxes.v.south + 0.5 * self._convFluxes.v.south)
