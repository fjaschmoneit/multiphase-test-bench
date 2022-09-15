import Fields
import PrimitiveFields

# def vectorToCell(vectorField, mesh):
#     sf = Fields.scalarField(mesh)
#     sf._raw[:,:] = vectorField.entries_EW[:,:-1] + vectorField.entries_EW[:,1:]
#     sf._raw[:,:] += vectorField.entries_NS[:-1,:] + vectorField.entries_NS[1:,:]
#     sf._raw[:,:] /= 4.0
#     return sf
# #    return Fields.scalarField(mesh, 10)

def scalarToFaces_EW(scalarField):
    # returns a primitive field
    return 0.5*(scalarField._raw[:,1:] + scalarField._raw[:,:-1])


#returns a primitive cell field as linear interpolation of faces in primitive facefield
def getCellInterpolation(primitiveFaceField, direction):
    f = primitiveFaceField
    if direction == 'x':
        return 0.5*( f[:,:-1] + f[:,1:] )
    elif direction == 'y':
        return 0.5*( f[:-1,:] + f[1:,:] )
    else:
        print("error: no direction chosen")

def cellToVertex(cellField, scheme='centralDifference'):

    mesh = cellField._mesh
    vf = Fields.vertexField(mesh)

    if scheme == 'centralDifference':

        ff = cellToVector(cellField)

        vf._raw[:-1,:] = ff.entries_EW
        vf._raw[1:,:] += ff.entries_EW
        vf._raw[:,1:] += ff.entries_NS
        vf._raw[:,:-1] += ff.entries_NS

    else:
        print("interpolation scheme {scheme} not supported.".format(**locals()))
    return vf

def cellToVector(cellField, scheme='centralDifference'):
    # returns the mean value of neighboring cell values as two face fields
    # A_e = 0.5*( A_P + A_E )    for f_e inner face
    # A_w = 1.5*A_P - 0.5*A_E    for f_w boundary face

    mesh = cellField._mesh
    ff = Fields.vectorField(mesh)

    if scheme == 'centralDifference':
        ff.internalEntries_EW = 0.5 * (cellField._raw[:, :-1] + cellField._raw[:, 1:])
        ff.internalEntries_NS = 0.5 * (cellField._raw[:-1, :] + cellField._raw[1:, :])
    # differnceField_x = 0.5 * (cellField.internal[:, :-1] + cellField.internal[:, 1:])
    # differnceField_y = 0.5 * (cellField.internal[:-1, :] + cellField.internal[1:, :])
    # ff.setInternalValues(differnceField_x, differnceField_y)

        ff.bn = 1.5 * cellField._raw[0:1,:] - 0.5 * cellField._raw[1:2,:]
        ff.bs = 1.5 * cellField._raw[-1:,:] - 0.5 * cellField._raw[-2:-1,:]
        ff.be = 1.5 * cellField._raw[:, -1:] - 0.5 * cellField._raw[:, -2:-1]
        ff.bw = 1.5 * cellField._raw[:, 0:1] - 0.5 * cellField._raw[:, 1:2]
    #
    # boundaryValues_w = 1.5 * cellField.internal[:, 0:1] - 0.5 * cellField.internal[:, 1:2]
    # boundaryValues_e = 1.5 * cellField.internal[:, -1:] - 0.5 * cellField.internal[:, -2:-1]
    #
    # boundaryValues_n = 1.5 * cellField.internal[0:1,:] - 0.5 * cellField.internal[1:2,:]
    # boundaryValues_s = 1.5 * cellField.internal[-1:,:] - 0.5 * cellField.internal[-2:-1,:]
    #
    # ff.setBoudnaryValues(field_n=boundaryValues_n, field_e=boundaryValues_e, field_w=boundaryValues_w, field_s=boundaryValues_s)
    else:
        print("interpolation scheme {scheme} not supported.".format(**locals()))
    return ff
