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

    return ff


