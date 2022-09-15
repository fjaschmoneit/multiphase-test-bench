import Fields
import Interpolation


class ScalarConvectionDiffusion():
    def __init__(self, mesh, geom, fieldRegistry, fieldFLowModelLink):

        self._depScalarName = 'T'
        self._velocityName = 'U'
        self._diffCoeffName = 'D'

        # should become a sparse field
        self._sourceField = Fields.scalarField(mesh=mesh, value=0)

        # dependent field
        if self._depScalarName not in fieldRegistry or not isinstance(fieldRegistry['T'], Fields.varScalarField):
            fieldRegistry[self._depScalarName] = Fields.varScalarField(mesh, geom)
        else:
            print("I don't understand. there seems to already exist a varScalarField 'T' ")

        # link to this flowmodel instance
        fieldFLowModelLink[self._depScalarName] = self

        # passive fields, need to be linked
        if self._velocityName not in fieldRegistry:
            fieldRegistry[self._velocityName] = Fields.vectorField(mesh)

        # coefficients
        fieldRegistry[self._diffCoeffName] = 0.0

    def updateLinearEquationSystem(self, mesh, objReg, linSys):
        # defining fluxes F,D and source S in entire mesh
        # passing to linearEquation System instance

        # definings convective flux, vectorField
        F = objReg[self._velocityName]

        # defining diffusive flux, vectorField
        A = mesh.getFaceAreas()
        D = mesh.getInverseCellDistances() * objReg[self._diffCoeffName] * A

        # defining source, scalarField
        S = self._sourceField * mesh._uniformSpacing

        linSys.update(mesh, F, D, S, objReg[self._depScalarName])

class IncompressibleFlow:
    def __init__(self, mesh, geom, fieldRegistry, fieldFLowModelLink):
    # updates velocity, under provided pressure field

        self._varVelocityXName = 'u'
        self._varVelocityYName = 'v'

        self._pressureName = 'p'
        self._densityName = 'rho'
        self._kinViscosityName = 'nu'

        # should become a sparse field
        self._sourceField = Fields.scalarField(mesh=mesh, value=0)

        # Defining fields, passing to obj reg:

        # dependent fields U
        if self._varVelocityXName not in fieldRegistry or not isinstance(fieldRegistry[self._varVelocityXName], Fields.varVectorField):
            fieldRegistry[self._varVelocityXName] = Fields.varVectorField(mesh, geom)   # should be a staggered field, not a vector field
        else:
            print("I don't understand. there seems to already exist a varScalarField 'U' ")

        if self._varVelocityYName not in fieldRegistry or not isinstance(fieldRegistry[self._varVelocityYName],
                                                                        Fields.varVectorField):
            fieldRegistry[self._varVelocityYName] = Fields.varVectorField(mesh, geom)
        else:
            print("I don't understand. there seems to already exist a varScalarField 'U' ")

        # link to this flowmodel instance
        fieldFLowModelLink[self._varVelocityXName] = self
        fieldFLowModelLink[self._varVelocityYName] = self

    # passive fields, need to be linked
        if self._pressureName not in fieldRegistry:
            fieldRegistry[self._pressureName] = Fields.scalarField(mesh)

        # coefficients
        fieldRegistry[self._kinViscosityName] = 0.0
        fieldRegistry[self._densityName] = 1.0

    def calcConvectiveFluxes(self, fieldname, mesh, vel, rho):
        F = Fields.staggeredFluxField_U(mesh=mesh, value=0)

        if fieldname == self._varVelocityXName:
            rho_f = Interpolation.cellToVector(rho)
            F.entries_EW = 0.5*(rho_f.e * vel.e + rho_f.w * vel.w)     # EW entries are a cellfield

            # interpolation of north-south values along EW direction (resulting at vertex interpolation)
            phi = 0.5*( rho_f.entries_NS[:,1:] * vel.entries_NS[:,1:] + rho_f.entries_NS[:,:-1] * vel.entries_NS[:,:-1])
            F.entries_NS = phi
        else:
            print("???")

        return F

    def calcDiffFluxes(self, gamma, mesh):

        gamma_f = Interpolation.cellToVector(gamma)
        q_f =  mesh.getInverseCellDistances()
        return gamma_f*q_f

    def updateLinearEquationSystem(self, mesh, objReg, linSys):
        # defining fluxes F,D and source S in entire mesh

        u = objReg[self._varVelocityXName]
        rho = Fields.scalarField(mesh=mesh, value=objReg[self._densityName])  # scalar or scalarField
        #gamma = Fields.scalarField(mesh=mesh, value=objReg[self._kinViscosityName])  # scalar or scalarField
        #rho = objReg[self._densityName]
        # defining convection flux field for u-staggered mesh:
        F = self.calcConvectiveFluxes(self._varVelocityXName, mesh, u, rho)
        D=F

#        D = self.calcDiffFluxes(gamma, mesh)
            # defining diffusive flux, vectorField
        #A = mesh.getFaceAreas()

            #D = Fields.staggeredFluxField_U(mesh=mesh, value=0)
        #D = Fields.staggeredFluxField_U(mesh=mesh, value=10)
    #        D = mesh.getInverseCellDistances() * objReg[self._kinViscosityName] * A

            # I don't expect a velocity source, maybe from bcs
            # S = self._sourceField * mesh._uniformSpacing

        S = Fields.scalarField(mesh,primitiveField=F.n)
        linSys.update(mesh, F, D, S, objReg[self._varVelocityXName])
