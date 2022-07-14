import Fields


def centralDifferencing(cellField, mesh):
    # returns the mean value of neighboring cell values as two face fields
    # A_e = 0.5*( A_P + A_E )    for f_e inner face
    # A_w = 1.5*A_P - 0.5*A_E    for f_w boundary face

    ff = Fields.parameterFaceField(mesh)

    differnceField_x = 0.5 * (cellField.internal[:, :-1] + cellField.internal[:, 1:])
    differnceField_y = 0.5 * (cellField.internal[:-1, :] + cellField.internal[1:, :])

    ff.setInternalValues(differnceField_x, differnceField_y)

    boundaryValues_w = 1.5 * cellField.internal[:, 0:1] - 0.5 * cellField.internal[:, 1:2]
    boundaryValues_e = 1.5 * cellField.internal[:, -1:] - 0.5 * cellField.internal[:, -2:-1]

    boundaryValues_n = 1.5 * cellField.internal[0:1,:] - 0.5 * cellField.internal[1:2,:]
    boundaryValues_s = 1.5 * cellField.internal[-1:,:] - 0.5 * cellField.internal[-2:-1,:]

    ff.setBoudnaryValues(field_n=boundaryValues_n, field_e=boundaryValues_e, field_w=boundaryValues_w, field_s=boundaryValues_s)

    return ff
