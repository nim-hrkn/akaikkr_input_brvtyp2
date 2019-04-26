import numpy as np
from pymatgen import Element,Composition
import pandas as pd

def heakey2zlist(key):
    if isinstance(key,int):
        key = str(key)
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
    if isinstance(key,int) or isinstance(key,np.int64):
        key = str(key)
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



def make_nlconfig_list(zlist,zconfig,finalcoreconf):
    df_finalcore_list = []
    for coreconfall in finalcoreconf:
        dftmplist = []
        for ispin,coreconf in enumerate(coreconfall):
            dftmp = pd.DataFrame(coreconf,
                                columns=["nl","energy_{}".format(ispin),\
                                         "choice_{}".format(ispin)]
                                )
            dftmp.set_index("nl",inplace=True,drop=True)
            #display(dftmp)
            dftmplist.append(dftmp)
        dfnew = pd.concat(dftmplist,axis=1)
        #display(dfnew)
        df_finalcore_list.append(dfnew)
        

    df_atomconfig_list = []
    for z,atomconfig in zip(zlist,zconfig):

        df_atomconfig = pd.DataFrame(atomconfig,columns=["nl","occupied_atom","energy_atom"])
        df_atomconfig.set_index("nl",inplace=True,drop=True)
        #display(df_atomconfig)
        df_atomconfig_list.append(df_atomconfig)
    
    df_nlconfig_list = []
    for z,df_atomconfig, df_finalcore in zip(zlist,df_atomconfig_list , df_finalcore_list):
        #print("Z",z)
        dftmp = pd.concat([df_atomconfig,df_finalcore],axis=1)
        #display(dftmp)
        df_nlconfig_list.append(dftmp)


    return df_nlconfig_list 





