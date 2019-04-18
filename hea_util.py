from pymatgen import Element,Composition

def heakey2zlist(key):
    if isinstance(key,str):
        # key -> zlist
        n = int(len(key)/2)
        zlist = []
        for i in range(n):
            zlist.append( key[2*i:2*i+2])
        return zlist
    else:
        print("input error",key)
        raise

def heakey2elements(key):
    if isinstance(key,str):
        # key -> zlist
        #n = int(len(key)/2)
        #zlist = []
        #for i in range(n):
        #    zlist.append( key[2*i:2*i+2])
        zlist = heakey2zlist(key)

    if isinstance(key,list) or isinstance(key,tuple):
        zlist = key
        
    # keylist -> elementlist
    elementlist = []
    for z in zlist:
        s = Element("H").from_Z(int(z))
        elementlist.append(str(s))
    return elementlist

def material2heakey(material):
    if isinstance(material,str):
        # material -> elementlist
        c = Composition(material)
        dic = c.as_dict()
        keylist = dic.keys()
    if isinstance(material,list) or isinstance(material,tuple):
        keylist = material
        
    # keylist -> zlist
    zlist = []
    for x in keylist:
        s = Element(x).Z
        zlist.append("{:02d}".format(int(s)))
    # zlist -> heakey
    return "".join(zlist)


