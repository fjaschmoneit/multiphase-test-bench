import Differentiation
import DifferenceSchemes
import Interpolation
import TransportBase
import BoundaryConditions
from fieldAccess import *
import MeshConfig
import ObjectRegistry as objReg
import Fields

LARGE = 1e20


class height(TransportBase.transportBase):

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

    def linkOtherFields(self, closure):
        if len(closure)==0:
            self.v = objReg.FIELDS['v']
            self.u = objReg.FIELDS['u']
        else:
            self.v = objReg.FIELDS[closure[1]]
            self.u = objReg.FIELDS[closure[0]]


    # this includes the advective fluxes when calculating the central coeffs
    def updateFluxes(self):
        self.reset()

        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas()
        u = self.u.data
        v = self.v.data

        F_u = u * faceAreas_u
        F_v = v * faceAreas_v

        self.a_w = -0.5*F_u[west]
        self.a_e = 0.5*F_u[east]
        self.a_n = 0.5*F_v[north]
        self.a_s = -0.5*F_v[south]

        # self.a_w = DifferenceSchemes.centralDifference(D_u[west], F_u[west], 'west')
        # self.a_e = DifferenceSchemes.centralDifference(D_u[east], F_u[east], 'east')
        # self.a_n = DifferenceSchemes.centralDifference(D_v[north], F_v[north], 'north')
        # self.a_s = DifferenceSchemes.centralDifference(D_v[south], F_v[south], 'south')

        self.correctBCs()

        self.a_p -= ( self.a_w + self.a_e + self.a_s + self.a_n  )


    def updateSourceField(self):
        self.sourceField_c += self.constSourceField



class staggeredTransport_u(TransportBase.transportBase):

    def __init__(self, mesh, linSystem):
        self.u = Fields.newField(shape=MeshConfig.SHAPE_FACES_U)

        super().__init__(mesh=mesh, field=self.u, linSystem=linSystem, depFieldShape=MeshConfig.SHAPE_FACES_U)

        self.boundaryModels = {
            'fixedValue' : self.derichlet,
            'zeroGradient' : self.vonNeumann
        }

    def linkOtherFields(self, closure):
        if len(closure)==0:
            self.v = objReg.FIELDS['v']
            self.p = objReg.FIELDS['p']
        else:
            self.v = objReg.FIELDS[closure[0]]
            self.p = objReg.FIELDS[closure[1]]


    def calcConvFlux(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas()

        u = self.u.data
        v = self.v.data

        return (
            Interpolation.toStaggered(u * faceAreas_u, 'u'),
            Interpolation.toStaggered(v * faceAreas_v, 'u')
        )

    def updateFluxes(self):
        self.reset()

        F_u, F_v = self.calcConvFlux()

        self.a_w = -0.5 * F_u[west]
        self.a_e = 0.5 * F_u[east]
        self.a_n = 0.5 * F_v[north]
        self.a_s = -0.5 * F_v[south]

        #self.correctBCs()

        self.a_p = 1e-15 +  (self.a_w + self.a_e + self.a_s + self.a_n )

    def updateSourceField(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas()
        gradP_u = Differentiation.grad_u(self.p.data, self.mesh)
        self.sourceField_c[internal_u] += -gradP_u * faceAreas_u[internal_u]

    def correctBCs(self):
        for argDict in self.boundary.values():
            bcType = argDict['type']
            self.boundaryModels[bcType](self, **argDict)

#            self.boundaryModels[bcType](self.depField.data, **argDict)
    @staticmethod
    def vonNeumann(transportInstance, **argDict):
        direction = argDict.get('direction')

        (boundary_dir, boundary_nb1_dir) = fieldSlice(direction)
        a_dir = transportInstance.getCoefficientsInDirection(direction)
        a_dir[boundary_dir] = 0.0

    @staticmethod
    def derichlet(transportInstance, **argDict):
        direction = argDict.get('direction')
        value = argDict.get('value')

        (boundary_dir, boundary_nb1_dir) = fieldSlice(direction)
        a_dir = transportInstance.getCoefficientsInDirection(direction)
        a_p = transportInstance.getCentreMatrixCoeffs()
        S_c = transportInstance.sourceField_c

        if direction == 'west' or direction == 'east':
            a_p[boundary_dir] += LARGE
            S_c[boundary_dir] += LARGE * value
            a_dir[boundary_dir] = 0.0

        elif direction == 'north' or direction == 'south':
            a_dir[boundary_dir] = 0.0


class staggeredTransport_v(TransportBase.transportBase):

    def __init__(self, mesh, linSystem):
        self.v = Fields.newField(shape=MeshConfig.SHAPE_FACES_V)

        super().__init__(mesh=mesh, field=self.v, linSystem=linSystem, depFieldShape=MeshConfig.SHAPE_FACES_V)

        self.boundaryModels = {
            'fixedValue' : self.derichlet,
            'zeroGradient' : self.vonNeumann
        }

    def linkOtherFields(self, closure):
        if len(closure)==0:
            self.u = objReg.FIELDS['u']
            self.p = objReg.FIELDS['p']
        else:
            self.u = objReg.FIELDS[closure[0]]
            self.p = objReg.FIELDS[closure[1]]

    def calcConvFlux(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas()

        u = self.u.data
        v = self.v.data

        return (
            Interpolation.toStaggered(u * faceAreas_u, 'v'),
            Interpolation.toStaggered(v * faceAreas_v, 'v')
        )

    def updateSourceField(self):
        faceAreas_u, faceAreas_v = self.mesh.calcFaceAreas()
        gradP_v = Differentiation.grad_v(self.p.data, self.mesh)
        self.sourceField_c[internal_v] += -gradP_v * faceAreas_v[internal_v]

    def updateFluxes(self):
        self.reset()

        F_u, F_v = self.calcConvFlux()

        self.a_w = -0.5 * F_u[west]
        self.a_e = 0.5 * F_u[east]
        self.a_n = 0.5 * F_v[north]
        self.a_s = -0.5 * F_v[south]

        self.correctBCs()

        self.a_p = 1e-15+ (self.a_w + self.a_e + self.a_s + self.a_n )

        def correctBCs(self):
            for argDict in self.boundary.values():
                bcType = argDict['type']
                self.boundaryModels[bcType](self, **argDict)

    #            self.boundaryModels[bcType](self.depField.data, **argDict)
    @staticmethod
    def vonNeumann(transportInstance, **argDict):
        direction = argDict.get('direction')

        (boundary_dir, boundary_nb1_dir) = fieldSlice(direction)
        a_dir = transportInstance.getCoefficientsInDirection(direction)
        a_dir[boundary_dir] = 0.0

    @staticmethod
    def derichlet(transportInstance, **argDict):
        direction = argDict.get('direction')
        value = argDict.get('value')

        (boundary_dir, boundary_nb1_dir) = fieldSlice(direction)
        a_dir = transportInstance.getCoefficientsInDirection(direction)
        a_p = transportInstance.getCentreMatrixCoeffs()
        S_c = transportInstance.sourceField_c

        if direction == 'north' or direction == 'south':
            a_p[boundary_dir] += LARGE
            S_c[boundary_dir] += LARGE * value
            a_dir[boundary_dir] = 0.0

        elif direction == 'east' or direction == 'west':

            a_dir[boundary_dir] = 0.0


