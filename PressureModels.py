

class Pressure:

    def __init__(self, depVariableName, velocityField_u_name, velocityField_v_name):
        self._depVariableName = depVariableName
        self._velocityField_u_name = velocityField_u_name
        self._velocityField_v_name = velocityField_v_name

        self._fieldReg = None
        self._sourceField_c = None
        self._presFluxes = None

    # I should have global instances of these fields
    def initializeFlowModel(self, mesh, fieldReg, fieldFlowModelLink):
        self._mesh = mesh
        self._fieldReg = fieldReg
        fGov = fieldReg['governor']     # remove this
        self._sourceField_c = fGov.newField(type='scalarCV')
        self._presFluxes = Fields.fieldContainer(
            u = fGov.newField(type='faces_u'),
            v = fGov.newField(type='faces_v') )

        self.initializeFields(fieldReg)
        self.linkDepFieldToModel(fieldFlowModelLink)

    def linkDepFieldToModel(self, flowModelFieldDict):
        flowModelFieldDict[self._depVariableName] = self
        self._flowmodel_u = flowModelFieldDict[self._velocityField_u_name]
        self._flowmodel_v = flowModelFieldDict[self._velocityField_v_name]


    # Should be done by Simulation
    def initializeFields(self, fieldReg):
        # adding neccessary field to registry
        # check if fields are not already defined

        fGov = fieldReg['governor']
        fieldReg[self._depVariableName] = fGov.newField(type='scalarCV', includeGhostNodes=True)
        fieldReg[self._velocityField_u_name] = fGov.newField(type='faces_u')
    #
    # def calcDCoeffs(self, fieldReg):
    #     fGov = fieldReg['governor']
    #     d = Fields.fieldContainer(
    #         u=fGov.newField(type='faces_u'),
    #         v=fGov.newField(type='faces_v')
    #     )
    #     uFlowModel = self._flowmodel_u.getTotalFlux_e
    #     d.u.data.fill(1)
    #     d.v.data.fill(2)
    #     return d

    def solve(self, linSystem):
        faceAreas = self._mesh.calcFaceAreas(self._fieldReg['governor'])
        a_p = self._fieldReg['governor'].newField(type='scalarCV', value=0.0)

        linSystem.reset()

        d_i = self._flowmodel_u.getTotalFluxes_e()[:,:-1]   # removing the ghost node
        A_i = faceAreas.u.east
        a_i = d_i*A_i*A_i
        linSystem.set_e_coeffs( a_i )
        a_p += a_i

        d_i = self._flowmodel_u.getTotalFluxes_w()[:,1:]
        A_i = faceAreas.u.west
        a_i = d_i*A_i*A_i
        linSystem.set_e_coeffs( a_i )
        a_p += a_i

        d_i = self._flowmodel_v.getTotalFluxes_n()
        A_i = faceAreas.v.north
        a_i = d_i*A_i*A_i
        linSystem.set_e_coeffs( a_i )
        a_p += a_i

        d_i = self._flowmodel_v.getTotalFluxes_s()
        A_i = faceAreas.v.south
        a_i = d_i*A_i*A_i
        linSystem.set_e_coeffs( a_i )
        a_p += a_i

        linSystem.set_e_coeffs( a_p )



    # this is the most important FM method. The others could possibly be base class methods
    def updateFluxes_pressure(self, fieldReg):
        faceAreas = self._mesh.calcFaceAreas(fieldReg['governor'])
        d = self.calcDCoeffs(fieldReg)
        self._presFluxes = faceAreas*faceAreas*d

    def setDerichlet(self, loc, value):
        return
        D = self._diffFluxes
        F = self._convFluxes
        Sc = self._sourceField_c
        #Sp = self._sourceField_p

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

    # why are convective fluxes not affected?
    def setVonNeumann(self, loc):
        return
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

    def updateSourceField(self, fieldReg, mesh):
        pass






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
