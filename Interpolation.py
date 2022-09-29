from Fields import fieldGovernor
#import ScalarField



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

def cellToVector(cellField, mesh, scheme='centralDifference'):
    # returns the mean value of neighboring cell values as two face fields
    # A_e = 0.5*( A_P + A_E )    for f_e inner face
    # A_w = 1.5*A_P - 0.5*A_E    for f_w boundary face

    if scheme == 'centralDifference':
        (nbC_U, nbC_V) = cellField.shape

        u = ScalarField.scalarField.fromShape( mesh.shapeStaggered_U )
        v = ScalarField.scalarField.fromShape( mesh.shapeStaggered_V )

        u.internalValues_U = 0.5*( cellField.e + cellField.w )
        v.internalValues_V = 0.5*( cellField.n + cellField.s )


        # u = ScalarField.scalarField(  0.5 * (cellField._raw[:, :-1] + cellField._raw[:, 1:]) )
        # v = ScalarField.scalarField( 0.5 * (cellField._raw[:-1, :] + cellField._raw[1:, :]) )

        ff = Fields.vectorField(u,v)

        ff.bn = 1.5 * cellField._raw[0:1,:] - 0.5 * cellField._raw[1:2,:]
        ff.bs = 1.5 * cellField._raw[-1:,:] - 0.5 * cellField._raw[-2:-1,:]
        ff.be = 1.5 * cellField._raw[:, -1:] - 0.5 * cellField._raw[:, -2:-1]
        ff.bw = 1.5 * cellField._raw[:, 0:1] - 0.5 * cellField._raw[:, 1:2]

    else:
        print("interpolation scheme {scheme} not supported.".format(**locals()))
    return ff


# def vectorToCell(vectorField, mesh):
#     sf = Fields.scalarField(mesh)
#     sf._raw[:,:] = vectorField.entries_EW[:,:-1] + vectorField.entries_EW[:,1:]
#     sf._raw[:,:] += vectorField.entries_NS[:-1,:] + vectorField.entries_NS[1:,:]
#     sf._raw[:,:] /= 4.0
#     return sf
# #    return Fields.scalarField(mesh, 10)
#
# def cellToVertex(cellField, fieldReg, scheme='centralDifference'):
#
#     fGov = fieldReg['governor']
#     vf = fGov.newScalarField(shape=fGov.shapeVertices)
#
#
#     if scheme == 'centralDifference':
#
#         ff = cellToVector(cellField)
#
#         vf._raw[:-1,:] = ff.entries_EW
#         vf._raw[1:,:] += ff.entries_EW
#         vf._raw[:,1:] += ff.entries_NS
#         vf._raw[:,:-1] += ff.entries_NS
#
#     else:
#         print("interpolation scheme {scheme} not supported.".format(**locals()))
#     return vf
