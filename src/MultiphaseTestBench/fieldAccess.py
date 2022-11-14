
west = (slice(None,None), slice(None,-1))
east = (slice(None,None), slice(1,None))
north = (slice(None,-1), slice(None,None))
south = (slice(1,None), slice(None,None))

internal = (slice(1,-1), slice(1,-1))
internal_u = (slice(None,None), slice(1,-1))
internal_v = (slice(1,-1), slice(None,None))

boundary_west = (slice(None,None),slice(0,1))
boundary_east = (slice(None,None),slice(-1,None))
boundary_north = (slice(0,1),slice(None,None))
boundary_south = (slice(-1,None),slice(None,None))

boundary_nb1_west = (slice(None,None),slice(1,2))
boundary_nb1_east = (slice(None,None),slice(-2,-1))
boundary_nb1_north = (slice(1,2),slice(None,None))
boundary_nb1_south = (slice(-2,-1),slice(None,None))

def fieldSlice( direction ):
    if direction == 'west':
        return (boundary_west, boundary_nb1_west)
    elif direction == 'east':
        return (boundary_east, boundary_nb1_east)
    elif direction == 'north':
        return (boundary_north, boundary_nb1_north)
    elif direction == 'south':
        return (boundary_south, boundary_nb1_south)
    else:
        print("direction not recognized: \t", direction)

def opposite(direction):
    if direction == 'west':
        return 'east'
    elif direction == 'east':
        return 'west'
    elif direction == 'north':
        return 'south'
    elif direction == 'south':
        return 'north'
    else:
        print("direction not recognized: \t", direction)
