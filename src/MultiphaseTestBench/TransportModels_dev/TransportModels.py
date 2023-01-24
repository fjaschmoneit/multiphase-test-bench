import Differentiation
import DifferenceSchemes
import Interpolation
import BoundaryConditions
from fieldAccess import *
import MeshConfig
import ObjectRegistry as objReg
import Fields
from LinearSystems import LinearEquationSystems_old
#import LinearSystems

class transportBase():
    a = 1

class scalarTransport():

    def __init__(self, mesh ):

        self.phi = Fields.newField(shape=MeshConfig.SHAPE_SCALAR_CV)
        self.diffusionCoefficient = 0.0
        self.constSourceField = Fields.newField(shape=MeshConfig.SHAPE_SCALAR_CV, value=0.0)
        self.mesh = mesh
        #self.sourceField

        # pointer to boundary functions:
        self.boundaryModels = {
            'fixedValue' : BoundaryConditions.scalarBC.derichlet,
            'zeroGradient' : BoundaryConditions.scalarBC.vonNeumann
        }

        self.linSystem = LinearEquationSystems_old.linearSystem(MeshConfig.SHAPE_SCALAR_CV)



    def setConstSourceField(self, value):
        self.constSourceField.fill(value)

    def setDiffusionCoefficient(self, value):
        self.diffusionCoefficient = value

    def getField(self):
        return self.phi

    def solve(self):
        import numpy as np
        return np.copy(self.linSystem.solve())

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

    def updateFluxes(self):
        self.linSystem.reset()

        # F_u, F_v = self.calcConvFlux()
        # D_u, D_v = self.calcDiffFlux()
        #
        # a_w = DifferenceSchemes.centralDifference(D_u[west], F_u[west], 'west')
        # a_e = DifferenceSchemes.centralDifference(D_u[east], F_u[east], 'east')
        # a_n = DifferenceSchemes.centralDifference(D_v[north], F_v[north], 'north')
        # a_s = DifferenceSchemes.centralDifference(D_v[south], F_v[south], 'south')

        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas()
        rCellDist_u, rCellDist_v = self.mesh.calcInvCellDistances()

        D_u = self.diffusionCoefficient *faceAreas_u *rCellDist_u
        D_v = self.diffusionCoefficient *faceAreas_v *rCellDist_v

        # make this a function:
        D_u[boundary_west] = 0
        D_u[boundary_east] = 0
        D_v[boundary_north] = 0
        D_v[boundary_south] = 0

        self.linSystem.a_w = -D_u[west]
        self.linSystem.a_p += D_u[west]

        self.linSystem.a_e = -D_u[east]
        self.linSystem.a_p += D_u[east]

        self.linSystem.a_n = -D_v[north]
        self.linSystem.a_p += D_v[north]

        self.linSystem.a_s = -D_v[south]
        self.linSystem.a_p += D_v[south]

        self.correctBCs()


    # fixing boundaries:
        D_u = self.diffusionCoefficient *faceAreas_u *rCellDist_u
        D_v = self.diffusionCoefficient *faceAreas_v *rCellDist_v

        # I cannot access this by reference
        ap = self.linSystem.a_p.copy()

        self.linSystem.b[boundary_west] = 2*D_u[boundary_west]*100
        ap[boundary_west] += 2*D_u[boundary_west]

        self.linSystem.b[boundary_east] = 2*D_u[boundary_east]*500
        ap[boundary_east] += 2*D_u[boundary_east]
        self.linSystem.a_p = ap

        # self.linSystem.setEastEntries(D_u[east])
        # self.linSystem.setWestEntries(D_u[west])
        # self.linSystem.setSouthEntries(D_v[internal_v])
        # self.linSystem.setNorthEntries(D_v[internal_v])

        #self.a_p += ( self.a_w + self.a_e + self.a_s + self.a_n + F_u[east] - F_u[west] + F_v[south] - F_v[north] ) # why do I substract these??
        # a_p = a_w+a_w+a_s+a_n
        # self.linSystem.setNorthEntries(a_p)

    def correctBCs(self):
        for argDict in self.boundary.values():
            bcType = argDict['type']
            self.boundaryModels[bcType](self, **argDict)

#            self.boundaryModels[bcType](self.depField.data, **argDict)

    # this includes the advective fluxes when calculating the central coeffs
    def updateFluxes_old(self):
        self.reset()

        F_u, F_v = self.calcConvFlux()
        D_u, D_v = self.calcDiffFlux()

        self.a_w = DifferenceSchemes.centralDifference(D_u[west], F_u[west], 'west')
        self.a_e = DifferenceSchemes.centralDifference(D_u[east], F_u[east], 'east')
        self.a_n = DifferenceSchemes.centralDifference(D_v[north], F_v[north], 'north')
        self.a_s = DifferenceSchemes.centralDifference(D_v[south], F_v[south], 'south')

        self.correctBCs()

        self.a_p += ( self.a_w + self.a_e + self.a_s + self.a_n + F_u[east] - F_u[west] + F_v[south] - F_v[north] ) # why do I substract these??

    def updateSourceField(self):
        pass
        #self.sourceField_c += self.constSourceField
        #self.linSystem.b += self.constSourceField


class staggeredTransport_u(transportBase):

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

class staggeredTransport_v(transportBase):

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

