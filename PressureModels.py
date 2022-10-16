import Fields
import BoundaryConditions
import numpy as np

class Pressure:

    def __init__(self, mesh, fieldCreator, fieldReg, linSystem):
        # self._depVariableName = depVariableName
        # self._velocityField_u_name = velocityField_u_name
        # self._velocityField_v_name = velocityField_v_name
        #
        # self._fieldReg = None
        # self._sourceField_c = None
        # self._presFluxes = None

        self._mesh = mesh
        self._linSystem = linSystem
        self._fieldReg = fieldReg
        self._fc = fieldCreator
        self._depFieldShape = fieldCreator._typeShapeDict['scalarCV']

        self._p = self._fc.newField(shape=self._depFieldShape, governingModel=self)
        self._sourceField_p = self._fc.newField(shape=self._depFieldShape, value=0.0)
        self._sourceField_c = self._fc.newField(shape=self._depFieldShape, value=0.0)

        #self.setSourceField(0.0)

    def getDepField(self):
        return self._p

    def linkOtherFields(self, fieldReg):
        self._u = fieldReg['u']
        self._v = fieldReg['v']

    def updateSourceField(self):
        pass

    def updateFluxes(self):
        return
#        faceAreas = self._mesh.calcFaceAreas(self._fc)
    #     d = self.calcDCoeffs()
    #     self._presFluxes = faceAreas*faceAreas*d

    #def calcDCoeffs(self):
#
#         d = Fields.fieldContainer(
#             u=self._fc.newField(type='faces_u'),
#             v=self._fc.newField(type='faces_v')
#         )
#         #uFlowModel = self._flowmodel_u.getTotalFlux_e
# #not wirking yet
#         uFlowModel = self._u.getGoverningTransportModel()
#         vFlowModel = self._v.getGoverningTransportModel()
#
#         a = uFlowModel.getTotalFlux_e
#         #b = vFlowModel.getTotalFlux().e
#
#         d.u.data.fill(1)
#         d.v.data.fill(0)
#         return d

    def correctBCs(self):
        for loc, typeValueTupel in self._p._boundary.items():
            (type, value) = typeValueTupel
            if type == 'zeroGradient':
                self.setVonNeumann(loc)
            else:
                self.setDerichlet(loc, value)

    # def correctBCs(self):
    #     for loc, typeValueTupel in fieldReg[self._depVariableName]._boundary.items():
    #         (type, value) = typeValueTupel
    #         if type == 'zeroGradient':
    #             self.setVonNeumann(loc)
    #         else:
    #             self.setDerichlet(loc, value)

    def solve(self):
        return self._linSystem.solve()

    def updateLinSystem(self):
        self._linSystem.reset(shape=self._depFieldShape)

        transModel_u = self._u.getGoverningTransportModel()
        transModel_v = self._v.getGoverningTransportModel()

        faceAreas = self._mesh.calcFaceAreas(self._fc)
        a_p = self._fc.newField(shape=self._depFieldShape, value=0.0)

        A_u = faceAreas.u
        a_u = self._fc.newField(data= transModel_u.getCentreMatrixCoeffs() )
        d_u = faceAreas.u/a_u

        a_i = d_u.east*A_u.east
        self._linSystem.set_e_coeffs( a_i )
        a_p += a_i

        a_i = d_u.west * A_u.west
        self._linSystem.set_w_coeffs(a_i)
        a_p+=a_i

        A_v = faceAreas.v
        a_v = self._fc.newField(data=transModel_v.getCentreMatrixCoeffs())
        d_v = faceAreas.v / a_v

        a_i = d_v.north * A_v.north
        self._linSystem.set_n_coeffs(a_i)
        a_p+=a_i

        a_i = d_v.south * A_v.south
        self._linSystem.set_s_coeffs(a_i)
        a_p += a_i
        #
        # A_i = faceAreas.u.east
        # d_i = A_i/transModel_u.getTotalFlux_e()[:,:-1]   # removing the ghost node
        # a_i = d_i*A_i
        # self._linSystem.set_e_coeffs( a_i )
        # a_p += a_i
        #
        # A_i = faceAreas.u.west
        # d_i = A_i/transModel_u.getTotalFlux_w()[:,1:]
        # a_i = d_i*A_i
        # self._linSystem.set_w_coeffs( a_i )
        # a_p += a_i
        #
        # A_i = faceAreas.v.north
        # n = transModel_v.getTotalFlux_n()[:-1,:]
        # d_i = A_i/transModel_v.getTotalFlux_n()[:-1,:]
        # a_i = d_i*A_i
        # self._linSystem.set_n_coeffs( a_i )
        # a_p += a_i
        #
        # A_i = faceAreas.v.south
        # d_i = A_i/transModel_v.getTotalFlux_s()[1:,:]
        # a_i = d_i*A_i
        # self._linSystem.set_s_coeffs( a_i )
        # a_p += a_i

        self._linSystem.set_p_coeffs( a_p.data )

        # setting source Field:
        self._sourceField_c = self._u.west*faceAreas.u.west-self._u.east*faceAreas.u.east \
                              + self._v.south*faceAreas.v.south - self._v.north*faceAreas.v.north

        linLength = self._linSystem._shape[0] * self._linSystem._shape[1]
        # make this also a funtion

        self._linSystem._b = np.reshape(self._sourceField_c.data, linLength)
        #self._linSystem._b

    def setVonNeumann(self, location):
        return
    #     BoundaryConditions.scalarBC.vonNeumann(location, D=self._diffFluxes)   Why do I have to pass these parameters?

    def setDerichlet(self, location, value):
        self._sourceField_c.bw.fill(0)
        #BoundaryConditions.scalarBC.derichlet(location, value, F=self._convFluxes, D=self._diffFluxes,
                                      #        Sc=self._sourceField_c, Sp=self._sourceField_p)
    #
    # def setDerichlet(self, loc, value):
    #     return
    #     D = self._diffFluxes
    #     F = self._convFluxes
    #     Sc = self._sourceField_c
    #     #Sp = self._sourceField_p
    #
    #     if loc == 'left':
    #         Sc.bw += ( 2.0 * D.u.bw + F.u.bw )* value
    #         Sp.bw += 2.0 * D.u.bw + F.u.bw
    #         D.u.bw = 0.0
    #         F.u.bw = 0.0
    #     elif loc == 'right':
    #         Sc.be += ( 2.0 * D.u.be - F.u.be) * value
    #         Sp.be += 2.0 * D.u.be - F.u.be
    #         D.u.be = 0.0
    #         F.u.be = 0.0
    #     elif loc == 'top':
    #         Sc.bn += ( 2.0 * D.v.bn - F.v.bn) * value
    #         Sp.bn +=  2.0 * D.v.bn - F.v.bn
    #         D.v.bn = 0.0
    #         F.v.bn = 0.0
    #     elif loc == 'bottom':
    #         Sc.bs += ( 2.0 * D.v.bs + F.v.bs) * value
    #         Sp.bs += 2.0 * D.u.bs + F.v.bs
    #         D.v.bs = 0.0
    #         F.v.bs = 0.0
    #
    # # why are convective fluxes not affected?
    # def setVonNeumann(self, loc):
    #     return
    #     D = self._diffFluxes
    #     # Sc = self._sourceField.Sc
    #     # Sp = self._sourceField.Sp
    #
    #     if loc == 'left':
    #         D.u.bw = 0.0
    #     elif loc == 'right':
    #         D.u.be = 0.0
    #     elif loc == 'top':
    #         D.v.bn = 0.0
    #     elif loc == 'bottom':
    #         D.v.bs = 0.0
    #




    # # I should have global instances of these fields
    # def initializeFlowModel(self, mesh, fieldReg, fieldFlowModelLink):
    #     self._mesh = mesh
    #     self._fieldReg = fieldReg
    #     fGov = fieldReg['governor']     # remove this
    #     self._sourceField_c = fGov.newField(type='scalarCV')
    #     self._presFluxes = Fields.fieldContainer(
    #         u = fGov.newField(type='faces_u'),
    #         v = fGov.newField(type='faces_v') )
    #
    #     self.initializeFields(fieldReg)
    #     self.linkDepFieldToModel(fieldFlowModelLink)
    #
    # def linkDepFieldToModel(self, flowModelFieldDict):
    #     flowModelFieldDict[self._depVariableName] = self
    #     self._flowmodel_u = flowModelFieldDict[self._velocityField_u_name]
    #     self._flowmodel_v = flowModelFieldDict[self._velocityField_v_name]
    #
    #
    # # Should be done by Simulation
    # def initializeFields(self, fieldReg):
    #     # adding neccessary field to registry
    #     # check if fields are not already defined
    #
    #     fGov = fieldReg['governor']
    #     fieldReg[self._depVariableName] = fGov.newField(type='scalarCV', includeGhostNodes=True)
    #     fieldReg[self._velocityField_u_name] = fGov.newField(type='faces_u')
    # #



#
# class ConstPresGrad:
#     def __init__(self, depVariableName, **fieldnames):
#         self._mesh = None
#         self._depVariableName = depVariableName
#
#     def initializeFields(self, fieldReg, mesh):
#         # adding neccessary field to registry
#         # check if fields are not already defined
#
#         fGov = fieldReg['governor']
#         fieldReg[self._depVariableName] = fGov.newField(type='scalarCV', includeGhostNodes=True)
#         self._mesh = mesh
#
#     def initializeFlowModel(self,  mesh, fieldReg, fieldFlowModelLink):
#         pass
#
#     def setConstPresGrad(self, presGrad, fieldReg):
#         # [presGrad] = Pa/m
#         cellDist = self._mesh._uniformSpacing
#         p = fieldReg[self._depVariableName]
#         p.data.fill(0.0)
#         offset = 0.0
#         p.data[:,:1] = offset + 0.5*cellDist*presGrad
#         for i in range(1,p.data.shape[1]):
#             p.data[:,i:i+1] = p.data[:,i-1:i] + cellDist*presGrad
