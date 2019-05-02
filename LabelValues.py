import numpy as np
import copy

class LabelValues(object):
    def __init__(self,label=None,values=None):
        self.label=None
        self.values=None
        if label is not None:
            self.label = copy.deepcopy(np.array(label))
        if values is not None:
            self.values = copy.deepcopy(np.array(values))
        # size 0 must be the same
        if self.label is not None and self.values is not None:
            if (self.label.shape[0] != self.values.shape[0] ): 
                print("size is not the same for label and values")
                raise 
    def set_nparray(self,v):
        self.label = copy.deepcopy(v[:,0])
        self.values = copy.deepcopy(v[:,1:])
        return self
    def get_nparray(self):
        label = self.label.reshape(-1,1)
        return np.hstack( (label, self.values) ) # 横につなげる
    def __add__(self, other):
        if (self.label == other.label ).all():
            return LabelValues( label = self.label, values = self.values + other.values)
        else:
            print("error: LabelValues::__add__, label is different")
            raise 
    def __sub__(self, other):
        if (self.label == other.label ).all():
            return LabelValues( label = self.label, values = self.values - other.values)
        else:
            print("error: LabelValues::__sub__, label is different")
            raise 
    def __mul__(self,f):
        return LabelValues( label = self.label, values = self.values*f ) 
    def __div__(self,f):
        return LabelValues( label = self.label, values = self.values/f ) 
    def __eq__(self,other):
        if  (self.label == other.label ).all():
            if (self.values == other.values).all():
                return True
        return False

if __name__ == "__main__":
    label = np.array( [0,1,2] )
    values =  np.array([ [ 10,11,12], [20,21,22], [30,31,32] ] )

    dos = LabelValues(label,values)   
    dos2 = LabelValues(label,values*0.01)
    dos3 = dos+dos2

    print(dos3.label)
    print(dos3.values)

    dos4 = dos3-dos
    print(dos4.label)
    print(dos4.values)


    print(dos.label)
    print(dos.values)
    print(dos2.label)
    print(dos2.values)
    print("==",dos==dos2)

    v =dos.get_nparray()
    print(v)

    dos5 = LabelValues().set_nparray(v)
    print(dos5.label)
    print(dos5.values)

