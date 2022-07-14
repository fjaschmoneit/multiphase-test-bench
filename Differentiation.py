import Fields


def grad(cellField, mesh):
    # differentialtion of cellField
    # returns a parameterFaceField

    ff = Fields.parameterFaceField(mesh=mesh, value=0)
    invCellDist = mesh.invCellDist

    d_x = (cellField.raw[:,1:] - cellField.raw[:,:-1])*invCellDist.internalEntries_EW
    d_y = (cellField.raw[1:,:] - cellField.raw[:-1,:])*invCellDist.internalEntries_NS

    ff.internalEntries_EW = d_x
    ff.internalEntries_NS = d_y

#    ff.setInternalValues(field_x=d_x, field_y=d_y)


    # d_x.setBoundaryValues_w( 1.5*cellField.field[:,0:1] - 0.5*cellField.field[:,1:2] )
    # d_x.setBoundaryValues_e(1.5 * cellField.field[:,-1:] - 0.5*cellField.field[:,-2:-1])
#    d_y = Fields.faceField_x(mesh)

    return ff


