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

        faceAreas = self._mesh.calcFaceAreas(self._fc)
        faceAreas_u = faceAreas.u.data
        faceAreas_v = faceAreas.v.data

        rCellDist_u = self._fieldReg['invCellDist'].u.data
        rCellDist_v = self._fieldReg['invCellDist'].v.data

        F_u = self.calcConvFlux(self.u.data, faceAreas_u)
        F_v = self.calcConvFlux(self.v.data, faceAreas_v)

        D_u = self.calcDiffFlux(rCellDist_u, self._diffusionCoefficient, faceAreas_u)
        D_v = self.calcDiffFlux(rCellDist_v, self._diffusionCoefficient, faceAreas_v)

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

    def updateLinSystem(self):
        self._linSystem.reset(shape=self._depFieldShape)
        self._linSystem.set_e_coeffs(self.a_e)
        self._linSystem.set_w_coeffs(self.a_w)
        self._linSystem.set_n_coeffs(self.a_n)
        self._linSystem.set_s_coeffs(self.a_s)
        self._linSystem.set_p_coeffs(self.a_p)
        self._linSystem.set_b(self._sourceField_c)

    def updateSourceField(self):
        pass

class staggeredTransport_u(TransportBase.transportBase):

    def __init__(self, mesh, fieldCreator, fieldReg, linSystem):
        self._depFieldShape = fieldCreator._typeShapeDict['faces_u']
        super().__init__(mesh, fieldCreator, fieldReg, linSystem, self._depFieldShape)
        self.u = self.depField
        self.setSourceField(0.0)

    def linkOtherFields(self, fieldReg):
        self.v = fieldReg['v']
        self.p = fieldReg['p']

    def calcInternalConvFlux(self, vel, faceAreas):
        f = vel*faceAreas
        return 0.5*( f[east] + f[west] )

    def calcInternalDiffFlux(self, invCellDist, diffCoeff, faceAreas):
        f = diffCoeff*invCellDist*faceAreas
        return 0.5*( f[east] + f[west] )


    # # what flux do I define in the ghost cells?
    def calcConvFlux_old(self, vel, faceAreas):
        (ny,nx) = vel.shape
        velA = Fields.newDataField((ny,nx+1), value=0.0)
        alpha = vel * faceAreas
        velA[west] += 0.5*alpha
        velA[east] += 0.5*alpha

        # velA[boundary_east] *=200
        # velA[boundary_west] *=200

        return velA

    def calcDiffFlux_old(self, invCellDist, diffCoeff, faceAreas):
        (ny,nx) = faceAreas.shape
        DA = Fields.newDataField((ny,nx+1), value=0.0)
        alpha = diffCoeff*invCellDist  * faceAreas
        DA[west] += 0.5*alpha
        DA[east] += 0.5*alpha
        # DA[boundary_west] *= 2
        # DA[boundary_east] *= 2
        return DA

    def updateFluxes(self):
        self.reset()

        faceAreas = self._mesh.calcFaceAreas(self._fc)
        faceAreas_u = faceAreas.u.data
        faceAreas_v = faceAreas.v.data

        rCellDist_u = self._fieldReg['invCellDist'].u.data
        rCellDist_v = self._fieldReg['invCellDist'].v.data

        # can I not implement my BC in the flux equations here?
        F_u = self.calcInternalConvFlux(self.u.data, faceAreas_u)
        #F_v = self.calcInternalConvFlux(self.v.data, faceAreas_v)

        D_u = self.calcInternalDiffFlux(rCellDist_u, self._diffusionCoefficient, faceAreas_u)
        #D_v = self.calcInternalDiffFlux(rCellDist_v, self._diffusionCoefficient, faceAreas_v)

        self.a_w[east] = DifferenceSchemes.centralDifference( D_u, F_u, 'west' )
        self.a_e[west] = DifferenceSchemes.centralDifference(D_u, F_u, 'east')

        F_v = self.calcConvFlux_old(self.v.data, faceAreas_v)
        D_v = self.calcDiffFlux_old(rCellDist_v, self._diffusionCoefficient, faceAreas_v)

        self.a_n = DifferenceSchemes.centralDifference(D_v, F_v, 'north')[north]
        self.a_s = DifferenceSchemes.centralDifference(D_v, F_v, 'south')[south]

#------------- BOUNDARIES ____________
        #---------WEST__________
        # fixed boundary
        # self.a_p[boundary_west] += 1e15
        # self._sourceField_c[boundary_west] += 1e15 * 2
        # self.a_w[boundary_west] = 0.0


        #---------EAST__________
        #fixed value
        # self.a_p[boundary_east] += 1e15
        # self._sourceField_c[boundary_east] += 1e15 * 2
        # self.a_e[boundary_east] = 0.0

        # zero Grad at east: nothing tobe done in flow direction
      #  self.a_w[boundary_west] = self.a_w[boundary_nb1_west]
        #self.a_e[boundary_east] = -2*self.a_w[boundary_east]
        #self.a_e[boundary_east] = 0

        #---------NORTH__________

        # fixed value
        #self.a_n[boundary_north] -= self._diffusionCoefficient*faceAreas_v[boundary_north]*2.0*rCellDist_v[boundary_north]
        #self.a_n[boundary_north] = 1 * 1 * 2.0 * 1

        # zero Grad
        #self.a_e[boundary_north] = 0.0
        #self.a_n[boundary_north] = 0


        #---------SOUTH__________

        # zeroGrad at south:
#        self.a_s[boundary_south] = self.a_s[boundary_nb1_south]
        #self.a_s[boundary_south] = 9

        # fixedValue at south:
       # self.a_s[boundary_south] = 1 * 1 * 2.0 * 1

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
            self._sourceField_c[boundary_west] += 1e15 * value
            self.a_w[boundary_west] = 0.0
        elif location == 'right':
            self.a_p[boundary_east] += 1e15
            self._sourceField_c[boundary_east] += 1e15 * value
            self.a_e[boundary_east] = 0.0

        # across flow direction
        elif location == 'top':
            # self.a_n[boundary_north] -= self._diffusionCoefficient*faceAreas_v[boundary_north]*2.0*rCellDist_v[boundary_north]
            self.a_n[boundary_north] = 1 * 1 * 2.0 * 1
        elif location == 'bottom':
            # self.a_n[boundary_north] -= self._diffusionCoefficient*faceAreas_v[boundary_north]*2.0*rCellDist_v[boundary_north]
            self.a_s[boundary_south] = 1 * 1 * 2.0 * 1


    def updateLinSystem(self):
        self._linSystem.reset(shape=self._depFieldShape)
        self._linSystem.set_e_coeffs(self.a_e)
        self._linSystem.set_w_coeffs(self.a_w)
        self._linSystem.set_n_coeffs(self.a_n)
        self._linSystem.set_s_coeffs(self.a_s)
        self._linSystem.set_p_coeffs(self.a_p)
        self._linSystem.set_b(self._sourceField_c)

    def updateSourceField(self):
        faceAreas_u = self._mesh.calcFaceAreas(self._fc).u
        gradP_u = Differentiation.grad_u(self.p, self._fieldReg)
        self._sourceField_c.internal_u = -gradP_u * faceAreas_u.internal_u





#        BoundaryConditions.staggeredBC.derichlet(location, value, F=self._convFluxes, D=self._diffFluxes,
 #                                             Sc=self._sourceField_c, Sp=self._sourceField_p)
    # # implement difference schemes here
    # def getTotalFlux_e(self):
    #     return -(self._diffFluxes.u.east - 0.5*self._convFluxes.u.east)
    #
    # def getTotalFlux_w(self):
    #     return -(self._diffFluxes.u.west + 0.5 * self._convFluxes.u.west)

    # def getTotalFlux_n(self):
    #     return -(self._diffFluxes.v.north - 0.5 * self._convFluxes.v.north)
    #
    # def getTotalFlux_s(self):
    #     return -(self._diffFluxes.v.south - 0.5 * self._convFluxes.v.south)

class staggeredTransport_v(TransportBase.transportBase):

    def __init__(self, mesh, fieldCreator, fieldReg, linSystem):
        self._depFieldShape = fieldCreator._typeShapeDict['faces_v']
        super().__init__(mesh, fieldCreator, fieldReg, linSystem, self._depFieldShape)

        self.v = self.depField
        self.setSourceField(0.0)

    def linkOtherFields(self, fieldReg):
        self.u = fieldReg['u']
        self.p = fieldReg['p']

        # this is the most important method. The others could possibly be base class methods
    def updateFluxes(self):
        rCellDist_u = self._fieldReg['invCellDist'].u
        rCellDist_v = self._fieldReg['invCellDist'].v

        self._convFluxes.u.internal_v = 0.5 * (self.u.north + self.u.south)
        self._convFluxes.v.internal_v = 0.5 * (self.v.north + self.v.south)

        self._diffFluxes.u.internal_v = 0.5 * (rCellDist_u.north + rCellDist_u.south) * self._diffusionCoefficient
        self._diffFluxes.v.internal_v = 0.5 * (rCellDist_v.north + rCellDist_v.south) * self._diffusionCoefficient

    def updateSourceField(self):
        faceAreas_v = self._mesh.calcFaceAreas(self._fc).v
        gradP_v = Differentiation.grad_v(self.p, self._fieldReg)
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
