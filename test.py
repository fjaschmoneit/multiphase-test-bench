
def printArgs(*args):

    for a in args:
        print(a)


b = 2
c = 4

printArgs(*locals())