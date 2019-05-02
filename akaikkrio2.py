
# coding: utf-8

# In[6]:

import numpy as np
import pickle

import os
import sys
import copy
import glob
import matplotlib.pyplot as plt
import pandas as pd

from pymatgen import Element

from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler



def linearregressionscore(y,plot=False):
    #print(len(y))
    y0 = np.array(y).reshape(-1,1)
    #x0 = np.array([i for i in range(len(y))]).reshape(-1,1)
    x0 = np.arange( len(y0) ).astype(float).reshape(-1,1)
    #print("x0,y0",x0,y0)
    scalerx = StandardScaler()
    scalerx.fit(x0)
    x = scalerx.transform(x0)
    
    scalery = StandardScaler()
    scalery.fit(y0)
    y = scalery.transform(y0)
    
    reg = LinearRegression().fit(x,y)
    yp = reg.predict(x)
    #print( y.ravel(),yp.ravel())

    if plot:
        plt.figure(figsize=(5,5))
        plt.plot( y.ravel(),yp.ravel())
        plt.xlabel("normalized raw")
        plt.ylabel("normalized predict")
        plt.show()

    r2score = r2_score(y.ravel(),yp.ravel())
    return r2score



class OutputJ: 
    def __init__(self,outputfile):
        self.dic = {}
        self.read(outputfile)
        
    def read(self,outputfile ):
        """
          analyze outjij.out
        """
        debug_ = False
        outputlist = open(outputfile).readlines()
        moment_comp=[]
        for output in outputlist:
                if 'total energy' in output:
                        str1=output.split() ; te=float(str1[2])
                elif 'spin moment' in output:
                        str1=output.split() ; moment_comp.append(float(str1[2]))
                elif 'mean field approximation' in output:
                        str1=output.split() ; Tc=str1[6]
                        Tc=Tc.split('K') ; Tc=float(Tc[0])
                elif 'a=' in output:
                        str1=output.split() ; 
                        for i,x in enumerate(str1): 
                            if x =="a=":
                                alen = float(str1[i+1])
                                break
                elif 'unit cell volume=' in output:
                        str1=output.split() ; vol = str1[3]
                        vol = float(vol.split('(a.u.)')[0])
                elif "itr=" in output and "moment=" in output:
                        s = output.replace("=","= ")
                        x = s.split()[5]
                        moment = float(x)

        self.dic = {"totalenergy":te, "moment":moment,"moment_comp": moment_comp, 
                    "Tc": Tc, "alen":alen,"vol":vol  }

        return self.dic

    def __str__(self,sep=" "):
        moment_comp = self.dic["moment_comp"] 
        key = self.dic["key"] 
        moment = self.dic["moment"]
        Tc = self.dic["Tc"]
        alen = self.dic["alen"]
        vol = self.dic["vol"]
        te = self.dic["totalenergy"]

        nelem  = len(moment_comp) 

        if sep==" ":
            s = "{:>10} ".format(key) 
            s += '{:>10} {:>10} '.format(alen,vol)
            s += "{:>10} {:>10} {:>10} {:>10}".format(moment,moment/vol,Tc,te)
            for i in range(nelem):
                s += "{:>10} ".format( moment_comp[i] )
        elif sep==",":
            v = [key,alen,vol,moment,moment/vol,Tc,te]
            v.extend(moment_comp)
            v = list( map(str,v) )
            s = ",".join(v)
        return s


class OutDisp:

    def __init__(self,filename,key,history=False):
        self.dic = {}
        self.read(filename,key,history)

    def short(self):
        sss = []
        for d in self.dic["property"]:
            s = [ d["alen"] ,d["len_h_err"], d["h_err"][-1], d["h_moment"][-1],d["totalenergy"] ]
            s.extend(d["moment"]) 
            sss.append(" ".join(list(map(str,s))))
        return "\n".join(sss)

    def read_file(self,filename):
        with open(filename) as f:
            lines = f.readlines()

        jobs = []
        onejob = []
        for x in lines:
            if 'data read in' in x:
                 if len(onejob)>0:
                     jobs.append(onejob)
                     onejob = []
                     onejob.append(x)
            else:
                 onejob.append(x)
        if len(onejob)>0:
            jobs.append(onejob)
        return jobs

    def read(self,filename,key,history=True):
        jobs = self.read_file(filename)
        # initialize
        diclist = []
        for onejob in jobs[1:]:
            moment_comp  = []
            valence_charge = []
            anclr = []
            conc = []
            dic = {}
            dic["totalenergy"] = None
            for output in onejob:
                if 'total energy' in output:
                        str1=output.split() ; te=str1[2]; dic["totalenergy"] = float(te)

                elif "anclr=" in output:
                        str1=output.split() ;  anclr.append(str1[3]) ; conc.append(float(str1[5]))
                elif 'spin moment' in output:
                        str1=output.split() ; moment_comp.append(float(str1[2]) )
                elif 'valence charge (up/down)=' in output:
                        str1=output.split() ; valence_charge.extend(list(map(float,str1[-2:-1]))) 
                elif 'mean field approximation' in output:
                        str1=output.split() ; Tc=str1[6]
                        Tc=Tc.split('K') ; Tc=Tc[0]; dic["Tc"] = float(Tc)
                elif 'a=' in output:
                        str1=output.split() ;
                        for i,x in enumerate(str1):
                            if x =="a=":
                                alen = str1[i+1]
                                break
                        dic["alen"] = float(alen)
                elif 'unit cell volume=' in output:
                        str1=output.split() ; vol = str1[3]
                        vol = vol.split('(a.u.)')[0]
                        dic["vol"] = float(vol)
            dic["moment"] = moment_comp
            dic["valence_charge"] = valence_charge 
            dic["anclr"] = anclr
            dic["conc"] = conc
            for i,output in enumerate(onejob):
                if 'self-consistent iteration starts' in output:
                   break  

            errlist = []
            h_momentlist = []
            h_telist = []
            kw = '   itr='
            nkw = len(kw)
            for output in onejob[i:]:
                if output[:nkw] == kw:
                    s = output.split()
                    if 'err=' in s:
                        str1=output.split() ;
                        for i,x in enumerate(str1):
                            if x =="err=":
                                err = str1[i+1]
                                errlist.append(float(err))
                                break
                        for i,x in enumerate(str1):
                            if x =="moment=":
                                h_moment = str1[i+1]
                                h_momentlist.append(float(h_moment))
                                break
                        for i,x in enumerate(str1):
                            if x =="te=":
                                h_te = str1[i+1]
                                h_telist.append(float(h_te))
                                break

            # original I planed to use all of the lists, but they are too big.
            if history==True:
                dic["len_h_err"] = len(errlist)
                dic["h_err"] = errlist
                dic["h_moment"] = h_momentlist
                dic["h_te"] = h_telist
            else:
                dic["len_h_err"] = len(errlist)
                dic["h_err"] = errlist[-1:]
                dic["h_moment"] = h_momentlist[-1:]
                dic["h_te"] = h_telist[-1:]

            dic["totalmoment"] = h_momentlist[-1]
            
            diclist.append(dic)

        converged = True
        for d in diclist:
            try:
                v =  d["h_err"][-1]
            except:
               converged = False
               break 
            if d["h_err"][-1] > -6.0:
               converged = False

        self.dic = {"key":key, "converged": converged, "property": diclist }

        return self



# coding: utf-8

# In[ ]:

# In[ ]:

def read_file(filename = "out.log"):
    with open(filename) as f:
        lines = f.readlines()
    lines2 = []
    for x in lines:
        x = x.rstrip()
        lines2.append(x)
    del lines
    lines = lines2
    return lines

def del_null(s):
    s2 = []
    for x in s:
        if len(x)>0:
            s2.append(x)
    return s2

class OutputDOS:
    
    def __init__(self,filename="out.log",target="lastDOS"):
        self.lines = read_file(filename = filename)    
        
        self.ewidth,_ = self.get_value(key="ewidth=",astype=float)
        self.magtyp,_ = self.get_value(key="magtyp=")
        self.ncmpx,_ = self.get_value(key="ncmpx=",astype=int)
        self.dic = { "ewidth":self.ewidth, "magtyp":self.magtyp, "ncmpx":self.ncmpx }

        self.cut_atomiclevel()
        if target=="firstDOS":
            self.totaldos,_ = self.cut_dos(key = " total DOS",normalize=True)
            self.derdos,_ = self.cut_dos(key = ' derivative of DOS',normalize=True)
            self.der2dos,_ = self.cut_dos(key = ' 2nd derivative of DOS',normalize=True)
            self.chargeneutral,_ = self.cut_dos(key = ' charge neutrality')        


        #print(self.magtyp,self.ncmpx)
        if target == "lastDOS":
            if self.magtyp == "mag":
                nspin = 2
            else:
                nspin = 1
            self.nspin = nspin
            self.dic["nspin"] = nspin

            #print("ncmpx",self.ncmpx)
            istart=0
            for icmp in range(self.ncmpx):
                 z,istart =  self.get_value("anclr=",astype=float,start=istart)
                 istart += 1 
                 self.zlist.append(z)

            istart = 0 
            self.totaldos = []
            self.componentdos = []
            self.integrateddos = []
            for ispin in range(nspin):

                componentdos = []
                for icmp in range(self.ncmpx):
                    dos,istart = self.cut_dos(key = " DOS of component", start=istart)
                    componentdos.append(dos)
                componentdos = np.array(componentdos)
                mse = componentdos.shape[1]
                self.componentdos.append(componentdos)

                dos,istart = self.cut_dos(key = " total DOS", start=istart)
                self.totaldos.append(dos)

                dos,istart = self.cut_dos(key = " integrated DOS", start=istart,mse=mse)
                self.integrateddos.append(dos)

        self.dic["componentdos"] = self.componentdos
        self.dic["totaldos"] = self.totaldos
        self.dic["integrateddos"] = self.integrateddos

    def get(self,key,spinsum=True):
        """
        averageを返す
        """

        if spinsum:
            v0 = self.dic[key][0]
            v1 = self.dic[key][1]
            vsum = (v0+v1)*0.5
            return vsum
        else:
            return self.dic[key][0]

    def plot(self,key,spinsum=True):
        if key=="totaldos":
            dos = self.get(key,spinsum)

            figsize = (10,1.5*2)
            fig,ax = plt.subplots(  2,1, figsize=figsize)
            ax[0].plot(dos[:,0],dos[:,1]*2.0)
            ax[0].set_xlabel("E")
            ax[0].set_ylabel("DOS")
            ax[1].plot(dos[:,0],np.log(dos[:,1]*2.0))
            ax[1].set_xlabel("E")
            ax[1].set_ylabel("log10(DOS)")
            plt.tight_layout()
            plt.show()
        elif key=="componentdos":
            dos = self.get(key,spinsum)
            ncmpx = self.dic["ncmpx"]
            figsize = (10,1.5*ncmpx)
            lmx = dos.shape[2]
            emin = dos[:,:,0].ravel().min()
            emax = dos[:,:,0].ravel().max()
            fig,ax = plt.subplots(  ncmpx,1, figsize=figsize)
            for icmpx,z in enumerate(self.dic["zlist"]):
                elem = Element("H").from_Z(z)
                for ilm in range(1,lmx):
                    ax[icmpx].plot(dos[icmpx,:,0],np.log10(dos[icmpx,:,ilm]*2.0),label=str(ilm-1))
                ax[icmpx].set_title("Z={} {}".format(z,str(elem)))
                ax[icmpx].legend()
                ax[icmpx].set_xlim( (emin,emax) )
                ax[icmpx].set_xlabel( "E" )
                ax[icmpx].set_ylabel( "log10(DOS)" )
            plt.tight_layout()

            plt.show()


    def get_value(self,key="ewidth=",astype=str,start=0):
        lines = self.lines
        #for x in lines:
        for iline in range(start,len(lines)):
            x = lines[iline]
            x = x.replace("=","= ")
            s = x.split(" ")
            if key in s:
                s = del_null(s)
                for i in range(len(s)):
                    if s[i]==key:
                        value = s[i+1]
                        return astype(value),iline
        return None,0


    def cut_atomiclevel(self,key = '   nuclear charge=',key2 = '         nl      cnf         energy' ):
        lines = self.lines
        
        flag = False
        config = None
        zlist = []
        zconfig = []
        i = 0
        n = len(lines)
        while i< n:
            line = lines[i]
            if line.startswith(key):
                x = line.split("=")
                zlist.append(float(x[1]))
            if line.startswith(key2):
                i += 2
                config = []
                for j in range(20):
                    line = lines[i]

                    if len(line)==0:
                        zconfig.append(config)
                        config = []
                        break
                    #print(line)
                    x = line.split()
                    config.append( [ x[0], float(x[1]), float(x[2])])                    
                    i+=1
            i = i+1
        if config is not None:
            if len(config)>0:
                zconfig.append(config)
        self.zlist , self.zconfig =  zlist,zconfig
        self.dic.update({"zlist":zlist,"zconfig":zconfig})


    def print_valencelevels(self,ewidth=-1.0):
        zlist , zconfig = self.zlist, self.zconfig
        for z, atomicconfig in zip(zlist,zconfig):
            for config in atomicconfig:
                if config[2]>ewidth:
                    print(z,config)


    def cut_dos(self,key = " total DOS", start=0, mse=None, normalize=False):
        istart = start
        lines = self.lines
        totaldos= []        
        if mse is not None:

            for ia,line in enumerate(self.lines[start:]):
                if line.startswith(key):
                    start += ia+1
                    break

            for ia,x in enumerate(self.lines[start:]):
                if ia==mse:
                    break

                vm =  list(map(float,x.split()))
                totaldos.append( vm ) 
            iline = start+mse

        else:
            started = False
            done = False

            #for x in lines:
            for iline in range(istart,len(lines)):
                x = lines[iline]
                if started and len(x)==0:
                    done = True
                if started and x.startswith(" ***err"):
                    done = True
                if started and x.startswith(" eb="):
                    done = True
                if started and x.startswith("   itr="):
                    done = True
                if started and x.startswith(" **itr="):
                    done = True
                if started and x.startswith("  *itr="):
                    done = True
                if started and x.startswith(" * itr="):
                    done = True
                if started and x.startswith("***itr="):
                    done = True
                if started and x.startswith("   ***msg"):
                    done = True

                if started and done:
                    break

                if started and not done:
                    try:
                        vm =  list(map(float,x.split()))
                        totaldos.append( vm )
                    except:
                        print("warning: conversion to float failed. drop the line.")
                        print(x)
                        vm = [x.split()[0],0,0,0]
                        vm = list(map(float,vm)) 
                        totaldos.append( vm )

                if x.startswith(key):
                    started = True
                
        totaldos = np.array(totaldos)
        if normalize and totaldos.shape[0]>0:
            for ic in range(1,totaldos.shape[1]):
                vm = np.max(totaldos[:,ic])
                totaldos[:,ic] = totaldos[:,ic]/vm 
                
        return np.array(totaldos), iline


"""
(1x,3(2x,f15.7,' Ry(',a2,')',a1))

1+3*(2+15+4+2+1+1) = 1*3*(25)

"""


def finalcorelevel2list(slist0):
    levels = []
    for ss in slist0:
        n = len(ss)
        m = int((n-1)/24)
        #print(m)
        slist = []
        for i in range(m):
            slist.append(  ss[1+25*i:1+25*i+25])
        #print("slist",slist)
        for s in slist:
            s = s.replace("("," ").replace(")"," ")
            x = s.split()
            #print("x",x)
            if len(x)==4:
                op = "valence"
            else:
                op = "core"
            nl = x[2]
            ene = x[0]
            #levelsdic[nl] = [ene, op]
            levels.append([nl,ene, op])
            
    return levels
        
def type2_name_z(s):
    x = s.replace("("," ").replace(")"," ").replace("=","= ")
    s = x.split()
    name =s[2]
    z = s[-2]
    return name,z
    
class OutputGo:
    
    def __init__(self,filename="out.log",heakey=None,action=["default"]):

        self.lines = read_file(filename = filename)    
        self.action = action
        self.dic = {}
        self.dic["heakey"] = heakey        
        
        self.dic["go"],_ = self.get_value(" go=")
        self.dic["file"],_ = self.get_value("file=")
        self.dic["edelt"],_ = self.get_value("edelt=",astype=float)        
        self.ewidth,_ = self.get_value(key="ewidth=",astype=float)
        self.magtyp,_ = self.get_value(key="magtyp=")
        if self.magtyp == "mag":
            self.nspin = 2
        else:
            self.nspin = 1
        self.dic["record"],_ = self.get_value(key="record=")
        self.dic["outtyp"],_ = self.get_value(key="outtyp=")
        self.dic["bzqlty"],_ = self.get_value(key="bzqlty=",astype=int)
        self.dic["maxitr"],_ = self.get_value(key="maxitr=",astype=int)
        self.dic["pmix"],_ = self.get_value(key="pmix=",astype=float)
            
        self.dic["ewidth"] = self.ewidth
        self.dic["magtyp"] = self.magtyp
        self.dic["nspin"] = self.nspin
        
        if "energylevels" in action:
            self.cut_atomiclevel()

        self.ncmpx,_= self.get_value(key="ncmpx=",astype=int)
        self.dic["ncmpx"] = self.ncmpx
        
        _,start = self.get_value(key="lattice constant")
        self.dic["alen"] ,_= self.get_value(key="a=",astype=float,start=start)
        if True:
            vol,_ = self.get_value(key="unit cell volume=")
            s = vol.replace("("," ").strip().split()
            self.dic["vol"] = float(s[0] )
        
        
        istart = 0
        self.zlist = []
        self.conclist = []
        for i in range(self.ncmpx):
            x,istart = self.get_value(key="anclr=",astype=float,start=istart)
            conc, istart = self.get_value(key="conc=",astype=float,start=istart)
            istart += 1
            self.zlist.append(x)
            self.conclist.append(conc)
        self.dic["zlist"] = self.zlist
        self.dic["anclr"] = self.zlist
        self.dic["conc"] = self.conclist
        
        if "energylevels" in action:
            self.finalcoreconf = []
            self.finalcorez = []
            for i in range(self.ncmpx):
                z, finalcoreconfig, start = self.get_finalcorelevel(start)
                #print(z,finalcoreconfig)
                self.finalcoreconf.append(finalcoreconfig)
                self.finalcorez.append(z)
            self.dic["finalcoreconf"] = self.finalcoreconf
            self.dic["finalcorez"] = self.finalcorez
            
        self.dic["totalmoment"],_ = self.get_value(key="spin moment=",astype=float)
        
        
        self.totalenergy,start = self.get_value(key="total energy=",astype=float)
        self.dic["totalenergy"] = self.totalenergy
        moment_comp = []
        for icmpx in range(self.dic["ncmpx"]):
            x,start = self.get_value(key="spin moment=",astype=float)
            moment_comp.append(x)
        self.dic["moment"] = moment_comp
        
        h_moment = []
        h_te = []
        h_err = []
        _,start = self.get_value("   ***** self-consistent")
        for line in self.lines[start+1:]:
            if line.startswith("   itr="):
                if line.find("chr,spn")>=0:
                    break
                x = self.get_value_fromline(line,"moment=",astype=float)
                h_moment.append(x)
                x = self.get_value_fromline(line,"te=",astype=float)
                h_te.append(x)
                x = self.get_value_fromline(line,"err=",astype=float)
                h_err.append(x)
            if line.startswith("   total energy="):
                break
        if "history" in action:

            self.dic["h_err"] = h_err
            self.dic["h_moment"] = h_moment
            self.dic["h_te"] = h_te
        else:
            self.dic["h_err"] = h_err[-1:]
            self.dic["h_moment"] = h_moment[-1:]
            self.dic["h_te"] = h_te[-1:] 

        if self.dic["h_err"][-1]<=-6.:
            self.dic["converged"] = True
        else:
            self.dic["converged"] = False
        
    def plot(self,plotkind="all"):
        if "history" in self.action:
            figsize=(10,3)
            fig,ax = plt.subplots(  2,1, figsize=figsize)
            h_err = self.dic["h_err"]
            h_moment = self.dic["h_moment"]
            ax[0].set_xlabel("iteration")
            ax[0].set_ylabel("h_err")
            ax[0].plot(range(len(h_err)), h_err)
            ax[1].set_xlabel("iteration")
            ax[1].set_ylabel("h_moment")
            ax[1].plot(range(len(h_moment)), h_moment)
            plt.tight_layout()
            plt.show()
                
    def get_finalcorelevel(self,start=0):
        lines = self.lines
        key = "                             *** type"
        for ispin  in range(self.nspin):
            if True:  # 下と見かけのlevelを合わせているだけ。
                found = False
                for iline in range(start,len(self.lines)):
                    x = lines[iline]
                    if x.startswith(key):
                        found = True
                        break
        name, z = type2_name_z(x)
        
        configlist  =  []
        key = "   core level  "
        for ispin  in range(self.nspin):
            if True:  # 下と見かけのlevelを合わせているだけ。
                found = False
                for iline in range(start,len(self.lines)):
                    x = lines[iline]
                    if x.startswith(key):
                        found = True
                        break
            if found:
                done = False
                corelevellines = []
                start = iline+1
                for iline in range(start,len(self.lines)):
                    x = lines[iline]
                    if x.startswith(key) or len(x)==0:
                        #print( x.startswith(key), len(x)==0)
                        done = True
                        start = iline
                        break
                    if not done:
                        corelevellines.append(x)

            #print(ispin,finalcorelevel2dic("".join(corelevellines)))
            configlist.append(finalcorelevel2list(corelevellines))
        return z, configlist, start
        
    def get_value(self,key="ewidth=",astype=str,start=0):
        lines = self.lines
        #for x in line:
        #for iline in range(start,len(lines)):
        for ia,x in enumerate(lines[start:]) :
            #x = lines[iline]
            ipos = x.find(key)
            if ipos>=0:
                x = x[ipos+len(key):]
                x = x.strip()
                s = x.split(" ")
                value = s[0]
                try:
                    return astype(value),start+ia
                except:
                    return None,start+ia
        return None,0
    
    def get_value_fromline(self,line,key="ewidth=",astype=str,start=0):

        x = line
        ipos = x.find(key)
        if ipos>=0:
            x = x[ipos+len(key):]
            x = x.strip()
            s = x.split(" ")
            value = s[0]

            try:
                return astype(value)
            except:
                return np.nan
        return None
    
    def cut_atomiclevel(self,key = '   nuclear charge=',key2 = '         nl      cnf         energy' ):
        lines = self.lines
        
        flag = False
        config = None
        zlist = []
        zconfig = []
        i = 0
        n = len(lines)
        while i< n:
            line = lines[i]
            if line.startswith(key):
                x = line.split("=")
                #print(x)
                zlist.append(float(x[1]))
            if line.startswith(key2):
                i += 2
                config = []
                for j in range(20):
                    line = lines[i]

                    if len(line)==0:
                        zconfig.append(config)
                        config = []
                        break
                    #print(line)
                    x = line.split()
                    config.append( [ x[0], float(x[1]), float(x[2])])                    
                    i+=1
            i = i+1
        if config is not None:
            if len(config)>0:
                zconfig.append(config)
        #self.zlist , self.zconfig =  zlist,zconfig
        self.zconfig =  zconfig
        self.dic["zconfig"] = zconfig

    def print_valencelevels(self,ewidth=-1.0):
        zlist , zconfig = self.zlist, self.zconfig
        for z, atomicconfig in zip(zlist,zconfig):
            for config in atomicconfig:
                if config[2]>ewidth:
                    print(z,config)


    def converging(self):
        last = 100
        dic = self.dic  
        r2moment = linearregressionscore( dic["h_moment"][-last:] )
        r2err = linearregressionscore(  dic["h_err"][-last:] )
        print( "use the last {}, r2moment={:6.2f} r2err={:6.2f}".format(last,r2moment, r2err) )
        if r2err>0.80:
            return True
        elif r2moment > 0.80 and r2err>0.80:
            return True
        else:
            return False


if __name__ == "__main__":

    #outgo = OutputGo("/home/kino/kino/work/fukushima_HEA/set3/err_sample/fcc_13143183_002/out_go.log",
    #                 "11223344",action=["default"])
    #print(outgo.dic)
    
    outdos = OutputDOS("/home/kino/kino/work/fukushima_HEA/set3/run110/bcc_25497683_3/out_dos.log")
    #outj = OutputJ("/home/kino/kino/work/fukushima_HEA/set3/run110/bcc_25497683_3/out_j.log")

    #print(outgo.dic)
    #print(outdos.dic.keys())
    #print(outj.dic)

    outdos.plot("componentdos")

