import Mesh
import Fields



class objectRegistry(object):

    def __init__(self):
        self.objects = {}

    def newMesh(self, type, **kwargs):
        #make this a mesh method. new()
        if type == 'cartesian2D':
            self.objects['mesh'] = Mesh.cartesian2D( cells_x=kwargs['cells_x'], cells_y=kwargs['cells_y'], dist=kwargs['dist'] )
        else:
            print(type, " is not supported")

    def newParameterCellField(self, name):
        self.objects[name] = Fields.parameterCellField( mesh=self.objects['mesh'], value=0 )

    def newParameterFaceField(self, name, interpolateFrom=None):
        self.objects[name] = Fields.parameterFaceField( mesh=self.objects['mesh'], value=0 )



    def getObjectNames(self):
        print("object registry contains:")
        for object in self.objects:
            print(object)





