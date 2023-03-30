

def centralDifference(D,F,orientation):
    if orientation not in ['east', 'west', 'north', 'south']:
        print("error: choose an orientation from: 'east', 'west', 'north', 'south'")

    if orientation in ['west', 'south']:
        alpha = 1.0
    else:
        alpha = -1.0

    #return D +alpha*0.5*F
#    return D +alpha*F