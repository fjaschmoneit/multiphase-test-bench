import Interpolation

class ScalarConvectionDiffusion():

    def __init__(self, depVariableName, **fieldnames):
        self._depVariableName = depVariableName
        self._velocityFieldName = fieldnames['velocityFieldName']
        self._diffusionCoeffName = fieldnames['diffusionCoeffName']
        self._sourceField = None
        self._convFluxes = None
        self._diffFluxes = None

    # I should have global instances of these fields
    def initializeFluxesAndSource(self, fieldReg):
        fGov = fieldReg['governor']
        self._sourceField = fGov.sourceField(fGov)
        self._convFluxes = fGov.newVectorField(shape_u=fGov.shapeFaces_u, shape_v=fGov.shapeFaces_v)  # do I need the boundary flux here at all?
        self._diffFluxes = fGov.newVectorField(shape_u=fGov.shapeFaces_u, shape_v=fGov.shapeFaces_v)

    def updateSourceField(self, mesh, fieldReg):
        constHeatSource = fieldReg[self._depVariableName]._source
        self._sourceField.setConstantSource( constHeatSource*mesh._uniformSpacing )

    def linkDepFieldToModel(self, fieldFlowModelLink):
        fieldFlowModelLink[self._depVariableName] = self

    def initializeFields(self, fieldReg, mesh):
        # adding neccessary field to registry
        # check if fields are not already defined

        fGov = fieldReg['governor']
        fieldReg[self._depVariableName] = fGov.newScalarField(shape=fGov.shapeVarScalarField)
        #fieldReg[self._depVariableName] = fGov.scalarField(shape=fGov.shapeCVField)
        fieldReg[self._velocityFieldName] = fGov.newVectorField(shape_u=fGov.shapeFaces_u, shape_v=fGov.shapeFaces_v)
        fieldReg[self._diffusionCoeffName] = 0.0

    # def updateBoundaryConditions(self, fieldReg):
    #     # top, bottom zeroFlux, left, right Derichlet
    #
    #     phi = fieldReg[self._depVariableName]
    #
    #     phi.gw = 100
    #     phi.ge = 500
    #     phi.gn = phi.nb
    #     phi.gs = phi.bs

    def correctBCs(self):
        D = self._diffFluxes
        Sc = self._sourceField.Sc
        Sp = self._sourceField.Sp

        # west
        Sc.bw += 2* D.u.bw * 100
        Sp.bw += D.u.bw
        D.u.bw =0

        # east
        Sc.be += 2*D.u.be * 500
        Sp.be += D.u.be
        D.u.be = 0

        # north
        Sp.bn += D.v.bn
        D.v.nb = 0

        # south
        Sp.bs += D.v.bs
        D.v.ns = 0

    def updateFluxes(self, fieldReg):
        self._convFluxes = fieldReg[self._velocityFieldName]
        print(fieldReg['invCellDist'].u.data)
        print(fieldReg['invCellDist'].v.data)

        self._diffFluxes = fieldReg['invCellDist'] * fieldReg[self._diffusionCoeffName]

class IncompressibleMomentumComp():

    def __init__(self, depVariableName, **kwargs):
        self._depVariableName = depVariableName
        self._orientation = kwargs['orientation']
        self._otherVelFieldName = kwargs['otherVelocityFieldName']
        self._pressureFieldName = kwargs['pressureFieldName']
        self._kinViscosity = kwargs['kinViscosityName']

    def defineFields(self, fieldReg, mesh):
        if self._orientation == 0:
            staggeredShape = mesh.shapeStaggered_U
        elif self._orientation == 1:
            staggeredShape = mesh.shapeStaggered_V
        else:
            print("error: chose orientation: 0 refers to u-field, 1 refers to v-field")
            staggeredShape = (0,0)

        fieldReg[self._depVariableName] = ScalarField.scalarField.fromShape(staggeredShape)
        fieldReg[self._otherVelFieldName] = ScalarField.scalarField.fromShape( mesh.shapeCVField )
        fieldReg[self._pressureFieldName] = ScalarField.scalarField.fromShape(mesh.shapeCVField)
        fieldReg[self._kinViscosity] = 0.0

    # should this be a cell field? Do I need it at all, or only at the pressure model?
    def calcSourceField(self, mesh, fieldReg):
        constMomentumSource = fieldReg[self._depVariableName]._source
        sourceField = ScalarField.scalarField.fromShape(mesh.shapeCVField, constMomentumSource)
        volSource = sourceField*mesh._uniformSpacing
        return volSource

    def linkDepFieldToModel(self, fieldFlowModelLink):
        fieldFlowModelLink[self._depVariableName] = self

    def calcFluxes(self, fieldReg, mesh):
        #F = fieldReg[self._velocityFieldName]

        F = Fields.vectorField.fromShapes(mesh.shapeStaggered_U, mesh.shapeStaggered_V)

        u_field = fieldReg[self._depVariableName]
        v_field = fieldReg[self._otherVelFieldName]

        F.u = 0.5*( u_field.e + u_field.w )

        # interpolation of north-south values along EW direction (resulting at vertex interpolation)
        F.v = 0.5*( v_field.e + v_field.w )

        #F.u =
        # this is for incompressible flow:
        # rho_f = Interpolation.cellToVector(rho, mesh)
        # F.entries_EW = 0.5 * (rho_f.e * vel.e + rho_f.w * vel.w)  # EW entries are a cellfield



        # phi = 0.5*( rho_f.entries_NS[:,1:] * vel.entries_NS[:,1:] + rho_f.entries_NS[:,:-1] * vel.entries_NS[:,:-1])
        # phi = 0.5 * (rho_f.v.e * vel.v.e + rho_f.v.w * vel.v.w)
        # F.entries_NS = phi

        D = F
#        D = mesh.getInverseCellDistances() * fieldReg[self._diffusionCoeffName]
        return F,D


class _IncompressibleFlow:
    def __init__(self, mesh, geom, fieldRegistry, fieldFLowModelLink):
    # updates velocity, under provided pressure field

        self._varVelocityXName = 'u'
        self._varVelocityYName = 'v'

        self._pressureName = 'p'
        self._densityName = 'rho'
        self._kinViscosityName = 'nu'

        # should become a sparse field
        self._sourceField = ScalarField.scalarField.fromShape( mesh.shapeCVField )

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
            fieldRegistry[self._pressureName] = ScalarField.scalarField.fromShape( mesh.shapeCVField )

        # coefficients
        fieldRegistry[self._kinViscosityName] = 0.0
        fieldRegistry[self._densityName] = 1.0

    def calcConvectiveFluxes(self, fieldname, mesh, vel, rho):
        F = Fields.staggeredFluxField_U(mesh=mesh, value=0)
        #F =

        if fieldname == self._varVelocityXName:
            rho_f = Interpolation.cellToVector(rho, mesh)
            F.entries_EW = 0.5*(rho_f.e * vel.e + rho_f.w * vel.w)     # EW entries are a cellfield

            # interpolation of north-south values along EW direction (resulting at vertex interpolation)
            #phi = 0.5*( rho_f.entries_NS[:,1:] * vel.entries_NS[:,1:] + rho_f.entries_NS[:,:-1] * vel.entries_NS[:,:-1])
            phi = 0.5*( rho_f.v.e * vel.v.e + rho_f.v.w * vel.v.w )
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
        rho = ScalarField.scalarField.fromShape(mesh.shapeCVField, objReg[self._densityName])  # scalar or scalarField
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

        S = ScalarField.scalarField(F.n)
        linSys.update(mesh, F, D, S, objReg[self._varVelocityXName])
