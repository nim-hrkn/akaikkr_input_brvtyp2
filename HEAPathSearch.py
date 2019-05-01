
# coding: utf-8

# In[60]:

# add KKR script path
import sys
#p = "/home/kino/kino/pythonscript/akaikkr/akaikkr_input_brvtyp2"
p = "/home/kino/kino/work/fukushima_HEA/akaikkr_input_brvtyp2"
sys.path.append(p)
from kkrinput_brvtyp import KKRinput_brv
from hea_util import *
from akaikkrio2 import OutputDOS, OutputGo, OutputJ

import os
import subprocess

import matplotlib.pyplot as plt
import numpy as np
from itertools import combinations
import pandas as pd
import glob
from IPython.display import display
import collections
import copy

os.environ["MKL_NUM_THREADS"]="1"


# In[61]:

polytyplist = ["bcc","fcc"]


# In[62]:

def headirname2prop(dir_):
    s = dir_.split(",")
    dic = collections.OrderedDict()
    for q in s:
        ss = q.split("_")
        dic[ss[0]] =ss[1]
    return dic

class HEAPathSearch:
    def __init__(self,prefix,dic):
        self.prefix = copy.deepcopy(prefix)
        self.dirkeys_default = collections.OrderedDict()
        self.dirkeys_default.update({"key":"*","ew":"*","polytyp":"*","ed":"*","pm":"*"})
        
        self.dirkeys = copy.deepcopy(self.dirkeys_default)
        self.dirkeys.update(dic)
        plist = []
        for key in self.dirkeys:
            plist.append("{}_{}".format(key,self.dirkeys[key]))
        self.dirname = ",".join(plist)
        
        self.targetkey = None
        for key in self.dirkeys:
            if self.dirkeys[key] == "*":
                self.targetkey = key
                break
        self.find_called = False

    def find(self):
        self.find_called = True
        if self.targetkey is None:
            self.targetvalues = None
            return self
        targetlist = []
        for dir_ in glob.glob(os.path.join(self.prefix,self.dirname)):
            headirname = os.path.split(dir_)[-1]
            prop = headirname2prop(headirname)
            targetlist.append(prop[self.targetkey])
        targetlist = list(set(targetlist))
        self.targetvalues = sorted(targetlist)
        return self
        
    def find_last(self):
        if not self.find_called:
            self.find()
        if self.targetkey is None:
            return None
        return self.targetvalues[-1]            
        
                   
if __name__ == "__main__":

    def make_dfresult(keylist):
        prefix = "RUN"
        result = []
        #keylist = ["13244274"]
        for key in keylist:
            ewlist = []
            heaserach = HEAPathSearch(prefix,{"key":key})
            ew = heaserach.find_last()
            eddic = {}    
            for polytyp in polytyplist:
                edlist = []
                heaPathSearch_ew = HEAPathSearch(prefix,{"key":key,"ew":ew,"polytyp":polytyp})
                eddic[polytyp] = heaPathSearch_ew.find_last()
            #print("ed",eddic)
            pmdic= {}
            for polytyp in polytyplist:
                pmlist = []
                heaPathSearch_pm = HEAPathSearch(prefix,{"key":key,"ew":ew,"polytyp":polytyp,"ed":eddic[polytyp]})
                pmdic[polytyp] = heaPathSearch_pm.find_last()
            #print("pm",pmdic)

            r = [key,ew]
            label = ["key","ew"]
            for polytyp in polytyplist:
                r.append(eddic[polytyp])
                r.append(pmdic[polytyp])        
                label.append("ed_{}".format(polytyp))
                label.append("pm_{}".format(polytyp))

            result.append(r)
        df = pd.DataFrame(result,columns=label)
        return df

    def get_uniqkeylist():
        keylist = []
        for dir_ in glob.glob("RUN/key_*"):
            headirname = os.path.split(dir_)[-1]
            prop = headirname2prop(headirname)
            keylist.append(prop["key"])
        keylist = list(set(keylist))
        return keylist


    def main():
        keylist =  get_uniqkeylist()
        df = make_dfresult(keylist)
        display( df.query("ed_bcc != ed_fcc"))

    main()

