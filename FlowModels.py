import Differentiation
import Fields
import Interpolation

class ConstPresGrad:
    def __init__(self, depVariableName, **fieldnames):
        self._mesh = None
        self._depVariableName = depVariableName

    def initializeFields(self, fieldReg, mesh):
        # adding neccessary field to registry
        # check if fields are not already defined

        fGov = fieldReg['governor']
        fieldReg[self._depVariableName] = fGov.newField(type='scalarCV', includeGhostNodes=True)
        self._mesh = mesh

    def initializeFlowModel(self,  mesh, fieldReg, fieldFlowModelLink):
        pass

    def setConstPresGrad(self, presGrad, fieldReg):
        # [presGrad] = Pa/m
        cellDist = self._mesh._uniformSpacing
        p = fieldReg[self._depVariableName]
        p.data.fill(0.0)
        offset = 0.0
        p.data[:,:1] = offset + 0.5*cellDist*presGrad
        for i in range(1,p.data.shape[1]):
            p.data[:,i:i+1] = p.data[:,i-1:i] + cellDist*presGrad

class ScalarConvectionDiffusion():

    def __init__(self, depVariableName, **fieldnames):
        self._mesh = None   # is thin still needed?
        self._depVariableName = depVariableName
        self._velocityFieldName = fieldnames['velocityFieldName']
        self._diffusionCoeffName = fieldnames['diffusionCoeffName']
        self._sourceField_p = None
        self._sourceField_c = None
        self._convFluxes = None
        self._diffFluxes = None

    # I should have global instances of these fields
    def initializeFlowModel(self, mesh, fieldReg, fieldFlowModelLink):
        self._mesh = mesh
        fGov = fieldReg['governor']

        self._sourceField_p = fGov.newField(type='scalarCV')
        self._sourceField_c = fGov.newField(type='scalarCV')

        self._convFluxes = Fields.fieldContainer(
            u = fGov.newField(type='faces_u'),
            v = fGov.newField(type='faces_v') )

        self._diffFluxes = Fields.fieldContainer(
            u = fGov.newField(type='faces_u'),
            v = fGov.newField(type='faces_v') )

        self.initializeFields(fieldReg)
        self.linkDepFieldToModel(fieldFlowModelLink)

    def linkDepFieldToModel(self, fieldFlowModelLink):
        fieldFlowModelLink[self._depVariableName] = self

    def initializeFields(self, fieldReg):
        # adding neccessary field to registry
        # check if fields are not already defined

        fGov = fieldReg['governor']
        fieldReg[self._depVariableName] = fGov.newField(type='scalarCV', includeGhostNodes=True)
        fieldReg[self._velocityFieldName] = Fields.fieldContainer(
            u=fGov.newField(type='faces_u'),
            v=fGov.newField(type='faces_v') )
        fieldReg[self._diffusionCoeffName] = 0.0

    def setConstSource(self, value):
        self._sourceField_c.fill(value*self._mesh._uniformSpacing)

    def setDerichlet(self, loc, value):
        D = self._diffFluxes
        F = self._convFluxes
        Sc = self._sourceField_c
        Sp = self._sourceField_p

        if loc == 'left':
            Sc.bw += ( 2.0 * D.u.bw + F.u.bw )* value
            Sp.bw += 2.0 * D.u.bw + F.u.bw
            D.u.bw = 0.0
            F.u.bw = 0.0
        elif loc == 'right':
            Sc.be += ( 2.0 * D.u.be - F.u.be) * value
            Sp.be += 2.0 * D.u.be - F.u.be
            D.u.be = 0.0
            F.u.be = 0.0
        elif loc == 'top':
            Sc.bn += ( 2.0 * D.v.bn - F.v.bn) * value
            Sp.bn +=  2.0 * D.v.bn - F.v.bn
            D.v.bn = 0.0
            F.v.bn = 0.0
        elif loc == 'bottom':
            Sc.bs += ( 2.0 * D.v.bs + F.v.bs) * value
            Sp.bs += 2.0 * D.u.bs + F.v.bs
            D.v.bs = 0.0
            F.v.bs = 0.0

    def setVonNeumann(self, loc):
        D = self._diffFluxes
        # Sc = self._sourceField.Sc
        # Sp = self._sourceField.Sp

        if loc == 'left':
            D.u.bw = 0.0
        elif loc == 'right':
            D.u.be = 0.0
        elif loc == 'top':
            D.v.bn = 0.0
        elif loc == 'bottom':
            D.v.bs = 0.0

    def correctBCs(self, fieldReg):
        for loc, typeValueTupel in fieldReg[self._depVariableName]._boundary.items():
            (type, value) = typeValueTupel
            if type == 'zeroGradient':
                self.setVonNeumann(loc)
            else:
                self.setDerichlet(loc, value)

    # this is the most important FM method. The others could possibly be base class methods
    def updateFluxes(self, fieldReg):
        self._convFluxes = fieldReg[self._velocityFieldName] * self._mesh.calcFaceAreas(fieldReg['governor'])
        self._diffFluxes = fieldReg['invCellDist'] * fieldReg[self._diffusionCoeffName] * self._mesh.calcFaceAreas(fieldReg['governor'])

    def updateSourceField(self, fieldReg, mesh):
        pass


# make inherit from some base flowmodel class
class IncompressibleMomentumComp:

    def __init__(self, depVariableName, **fieldnames):
        self._mesh = None
        self._sourceField_c = None
        self._sourceField_p = None
        self._depVariableName = depVariableName
        self._orientation = fieldnames['orientation']
        self._otherVelFieldName = fieldnames['otherVelocityFieldName']
        self._pressureFieldName = fieldnames['pressureFieldName']
        self._kinViscosity = fieldnames['kinViscosityName']
        self._boundaryDict = {}

    # I should have global instances of these fields
    def initializeFlowModel(self, mesh, fieldReg, fieldFlowModelLink):
        self._mesh = mesh
        fGov = fieldReg['governor']

        # these container dimensions are only valid for momentum U component
        # must make abstraction
        self._convFluxes = Fields.fieldContainer(
            u=fGov.newField(type='scalarCV'),
            v=fGov.newField(type='internalVertices_u')
        )

        self._diffFluxes = Fields.fieldContainer(
            u=fGov.newField(type='scalarCV'),
            v=fGov.newField(type='internalVertices_u')
        )

        self._sourceField_c = fGov.newField(type='faceSource_u')
        self._sourceField_p = fGov.newField(type='faceSource_u')

        self.initializeFields(fieldReg)
        self.linkDepFieldToModel(fieldFlowModelLink)

    # this could be a base class
    def linkDepFieldToModel(self, fieldFlowModelLink):
        fieldFlowModelLink[self._depVariableName] = self

    def initializeFields(self, fieldReg):
        # adding neccessary field to registry
        # check if fields are not already defined

        fieldReg[self._depVariableName] = self.linkOrCreateField(self._depVariableName, fieldReg, type='faces_u')
        fieldReg[self._otherVelFieldName] = self.linkOrCreateField(self._otherVelFieldName, fieldReg, type='faces_v')
        fieldReg[self._pressureFieldName] = self.linkOrCreateField(self._pressureFieldName, fieldReg, type='scalarCV')
        fieldReg[self._kinViscosity] = 0.0

    def correctBCs(self, fieldReg):
        for loc, typeValueTupel in fieldReg[self._depVariableName]._boundary.items():
            (type, value) = typeValueTupel
            if type == 'zeroGradient':
                self.setVonNeumann(loc)
            else:
                self.setDerichlet(loc, value)


    # this is the most important FM method. The others could possibly be base class methods
    def updateFluxes(self, fieldReg):

        fGov = fieldReg['governor']

        u = fieldReg[self._depVariableName]
        v = fieldReg[self._otherVelFieldName]
        rCellDist_u = fieldReg['invCellDist'].u
        rCellDist_v = fieldReg['invCellDist'].v

        self._convFluxes.u = fGov.newField( data=0.5*( u.east + u.west ) )
        self._convFluxes.v = fGov.newField( data=0.5 *(v.east + v.west) )

        self._diffFluxes.u = fGov.newField( data=0.5*( rCellDist_u.east + rCellDist_u.west ) * fieldReg[self._kinViscosity] )
        self._diffFluxes.v = fGov.newField( data=0.5*( rCellDist_v.east + rCellDist_v.west )* fieldReg[self._kinViscosity] )

    def updateSourceField(self, fieldReg, mesh):
        p = fieldReg[self._pressureFieldName]
        faceAreas_u = mesh.calcFaceAreas(fieldReg['governor']).u
        gradP_u = Differentiation.grad_u(p, fieldReg)
        self._sourceField_c.data = -gradP_u*faceAreas_u.internal_u

    # this could be a base method
    def linkOrCreateField(self, name, fieldReg, type):
        if name not in fieldReg:
            fGov = fieldReg['governor']
            return fGov.newField(type=type)
        else:
            return fieldReg[name]

    def setDerichlet(self, loc, value):
        D = self._diffFluxes
        F = self._convFluxes
        Sc = self._sourceField_c
        Sp = self._sourceField_p

        # if loc == 'left':
        #     Sc.bw += ( 2.0 * D.u.bw + F.u.bw )* value
        #     Sp.bw += 2.0 * D.u.bw + F.u.bw
        #     D.u.bw = 0.0
        #     F.u.bw = 0.0
        # elif loc == 'right':
        #     Sc.be += ( 2.0 * D.u.be - F.u.be) * value
        #     Sp.be += 2.0 * D.u.be - F.u.be
        #     D.u.be = 0.0
        #     F.u.be = 0.0
        # elif loc == 'top':
        #     Sc.bn += ( 2.0 * D.v.bn - F.v.bn) * value
        #     Sp.bn +=  2.0 * D.v.bn - F.v.bn
        #     D.v.bn = 0.0
        #     F.v.bn = 0.0
        # elif loc == 'bottom':
        #     Sc.bs += ( 2.0 * D.v.bs + F.v.bs) * value
        #     Sp.bs += 2.0 * D.v.bs + F.v.bs
        #     D.v.bs = 0.0
        #     F.v.bs = 0.0

    def setVonNeumann(self, loc):
        D = self._diffFluxes
        F = self._convFluxes
        Sc = self._sourceField_c
        Sp = self._sourceField_p

        # if loc == 'left':
        #     D.u.bw = 0.0
        #     F.u.bw = 0.0
        # elif loc == 'right':
        #     D.u.be = 0.0
        #     F.u.be = 0.0
        # elif loc == 'top':
        #     D.v.bn = 0.0
        #     F.v.bn = 0.0
        # elif loc == 'bottom':
        #     D.v.bs = 0.0
        #     F.v.bs = 0.0









