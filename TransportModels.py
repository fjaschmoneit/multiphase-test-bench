import Differentiation
import Fields
import Interpolation
import TransportBase
import BoundaryConditions


class scalarTransport(TransportBase.transportBase):

    def __init__(self, mesh, fieldCreator, fieldReg, linSystem ):

        self._depFieldShape = fieldCreator._typeShapeDict['scalarCV']

        super().__init__(mesh, fieldCreator, fieldReg, linSystem, self._depFieldShape)
        self._depField.setGoverningTransportModel(self)
        self._phi = self._depField      # defining my local name
        self.setSourceField(0.0)

    def linkOtherFields(self, fieldReg):
        self._u = fieldReg['u']
        self._v = fieldReg['v']

    def updateFluxes(self):
        A_f = self._mesh.calcFaceAreas(self._fc)
        self._convFluxes = self._convectionField * A_f
        self._diffFluxes = self._fieldReg['invCellDist'] * self._diffusionCoefficient * A_f

    def setVonNeumann(self, location):
        BoundaryConditions.scalarBC.vonNeumann( location, D=self._diffFluxes )

    def setDerichlet(self, location, value):
        BoundaryConditions.scalarBC.derichlet( location, value, F=self._convFluxes, D=self._diffFluxes, Sc=self._sourceField_c, Sp=self._sourceField_p )
    # def updateSourceField(self):
    #     self._sourceField_c = self._fc.newField(data=self._outerSource.data)

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

        # this is the most important method. The others could possibly be base class methods
    def updateFluxes(self):
        rCellDist_u = self._fieldReg['invCellDist'].u
        rCellDist_v = self._fieldReg['invCellDist'].v

        self._convFluxes.u.internal_u = 0.5 * (self._u.east + self._u.west)
        self._convFluxes.v.internal_u = 0.5 * (self._v.east + self._v.west)

        self._diffFluxes.u.internal_u = 0.5 * (rCellDist_u.east + rCellDist_u.west) * self._diffusionCoefficient
        self._diffFluxes.v.internal_u = 0.5 * (rCellDist_v.east + rCellDist_v.west) * self._diffusionCoefficient

    # def updateSourceField(self, fieldReg, mesh):
    #     p = fieldReg[self._pressureFieldName]
    #     faceAreas_u = mesh.calcFaceAreas(fieldReg['governor']).u
    #     gradP_u = Differentiation.grad_u(p, fieldReg)
    #     self._sourceField_c.internal_u = -gradP_u*faceAreas_u.internal_u

    def setVonNeumann(self, location):
        BoundaryConditions.staggeredBC.vonNeumann(location, F=self._convFluxes, D=self._diffFluxes,
                                              Sc=self._sourceField_c, Sp=self._sourceField_p)

    def setDerichlet(self, location, value):
        BoundaryConditions.staggeredBC.derichlet(location, value, F=self._convFluxes, D=self._diffFluxes,
                                              Sc=self._sourceField_c, Sp=self._sourceField_p)

    # implement difference schemes here
    def getTotalFluxes_e(self):
        return -(self._diffFluxes.u.east - 0.5*self._convFluxes.u.east)

    def getTotalFluxes_w(self):
        return -(self._diffFluxes.u.west + 0.5 * self._convFluxes.u.west)

    def getTotalFluxes_n(self):
        return -(self._diffFluxes.v.north - 0.5 * self._convFluxes.v.north)

    def getTotalFluxes_s(self):
        return -(self._diffFluxes.v.south - 0.5 * self._convFluxes.v.south)




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

    # def updateSourceField(self, fieldReg, mesh):
    #     p = fieldReg[self._pressureFieldName]
    #     faceAreas_u = mesh.calcFaceAreas(fieldReg['governor']).u
    #     gradP_u = Differentiation.grad_u(p, fieldReg)
    #     self._sourceField_c.internal_u = -gradP_u*faceAreas_u.internal_u

    def setVonNeumann(self, location):
        BoundaryConditions.staggeredBC.vonNeumann(location, F=self._convFluxes, D=self._diffFluxes,
                                              Sc=self._sourceField_c, Sp=self._sourceField_p)

    def setDerichlet(self, location, value):
        BoundaryConditions.staggeredBC.derichlet(location, value, F=self._convFluxes, D=self._diffFluxes,
                                              Sc=self._sourceField_c, Sp=self._sourceField_p)

    # implement difference schemes here
    def getTotalFluxes_e(self):
        return -(self._diffFluxes.u.east - 0.5*self._convFluxes.u.east)

    def getTotalFluxes_w(self):
        return -(self._diffFluxes.u.west + 0.5 * self._convFluxes.u.west)

    def getTotalFluxes_n(self):
        return -(self._diffFluxes.v.north - 0.5 * self._convFluxes.v.north)

    def getTotalFluxes_s(self):
        return -(self._diffFluxes.v.south - 0.5 * self._convFluxes.v.south)






