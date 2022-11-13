import numpy as np

#
# def newField(type):
#     return field(type)

class field:
    def __init__(self, type):
        self.type = type
        if self.type == 'scalar':
            self.data = np.ones(4)
        else:
            self.data = np.zeros(4)

class tempField(field):
    def __init__(self, type):
        self.isUsed = True
        super().__init__(type)

    def __del__(self):
        print("moin")
        self.isUsed = False

class temporaryField:
    def __init__(self):
        self.tempFieldList = []

    def new(self, type):
        for field in self.tempFieldList:
            if field.type == type and field.isUsed == False:
                #isUsed = True
                return field

        field = tempField(type)
        self.tempFieldList.append( field )
        return field

tf = temporaryField()

a = tf.new('scalar')
b = tf.new('vector')
c = tf.new('scalar')

del a

print(a)

for field in tf.tempFieldList:
    print(field.type, field.isUsed, field)

print('bye bye')