
# coding: utf-8

# In[99]:

# add KKR script path
import sys
#p = "/home/kino/kino/pythonscript/akaikkr/akaikkr_input_brvtyp2"
p = "/home/kino/kino/work/fukushima_HEA/akaikkr_input_brvtyp2"
sys.path.append(p)
from kkrinput_brvtyp import KKRinput_brv
from hea_util import *
from akaikkrio2 import OutputDOS, OutputGo, OutputJ,OutputCommon

import os
import subprocess

import matplotlib.pyplot as plt
import numpy as np
from itertools import combinations
import pandas as pd
import glob
from IPython.display import display
import shutil

os.environ["MKL_NUM_THREADS"]="1"

# In[104]:


class analyzeDOS:
    def __init__(self,dossumlist,dosth):
        """
        dossum list
        dosth: DOS threshold , dossum[:,1] < dosth
        """
        self.dossumlist = dossumlist
        self.dosth = dosth

        # check tne energy range of the DOS is the same
        eminlist = []
        emaxlist = []
        for dossum in dossumlist:
            #print("dossum",dossum[:5,:])
            eminlist.append(dossum[0,0])
            emaxlist.append(dossum[-1,0])
        eminlist = np.array(eminlist)
        emaxlist = np.array(emaxlist)
        if (eminlist != eminlist[0]).any():
            print("energy range different, for min")
            raise
        if (emaxlist != emaxlist[0]).any():
            print("energy range different, for max")
            raise
        print("emin,emax",eminlist[0],emaxlist[0])
        
        smalldosregion = []
        for dossum in dossumlist:
            #print("dossum",dossum.shape)
            smalldos0 =  dossum[:,1] < dosth
            smalldos0 = list( map(int,smalldos0) )
            smalldosregion.append( smalldos0 )                 
            #for e0,d0, sm in zip(dossum[:,0],dossum[:,1],smalldos0):
            #    if sm:
            #        print(e0,d0)

        smalldosregion = np.array(smalldosregion)
        smalldosregionsum = np.sum(smalldosregion,axis=0)
        # =1 if both smalldos is 1
        smalldosregionsum = smalldosregionsum//len(dossumlist) # This is the same as AND operation
        #print(smalldosregionsum)

        #for ie,(ene,val) in enumerate(zip(dossumlist[0][:,0],smalldosregionsum)):
        #    print(ie,ene,val)
    
        # [01] list -> region list
        start = False
        istart = 0
        smallregion = []
        for i,t in enumerate(smalldosregionsum):
            if start and not t:
                # end
                smallregion.append([istart,i])
                start = False
            if t and start==False:
                start = True
                istart = i
        if start:
            smallregion.append([istart,len(smalldos)])
        #print(smallregion)

        dossum = dossumlist[0] # use only energy
        print("gap region")
        print(dossum[smallregion,0])
        self.smallregion = smallregion
    
    def calc_new_ewidth(self,ewidth):
        """
        ewidth : 現在のewidth
        
        ewidthはpositive, 下で見る下限はnegative
        """
        eth = 0.30
        ediff = 0.20        
        
        smallregion = self.smallregion
        
        # high energyからに変換
        smallregion = list(reversed(smallregion))

        newewidth = None  
        canusethisewidth = -1
        dossum = self.dossumlist[0]  # for energy
        
        # ewidthを動かさないで良いか？
        for iregion,ii in enumerate(smallregion):
            e1,e2 = dossum[ii[0],0],dossum[ii[1],0]
            #print(e1,e2,e2-e1)
            if e2-e1> eth and e1<-ewidth and -ewidth < e2 :
                if -ewidth < e2-ediff:
                    canusethisewidth = iregion
                    #print("can hold this ewidth in region",canusethisewidth)
                    break

        # ewidthを動かす場合
        #print("move",smallregion)
        newewidth_choice = []
        for ii in smallregion:
            e1,e2 = dossum[ii[0],0],dossum[ii[1],0]
            #print(e1,e2,e2-e1,eth, (e2-e1)>eth)
            if e2-e1> eth:
                newewidth = e2-ediff -0.01 # add additional 0.01Ry
                newewidth_choice.append( newewidth )
        #print("possible choice",newewidth_choice)
        newewidth_choice = -np.array(newewidth_choice)
        # newewidth_choice must be a list 
        newewidth_choice = list(newewidth_choice)
        
        if canusethisewidth>=0:
            #newewidth = -ewidth
            return "old",ewidth,newewidth_choice,canusethisewidth
        else:
            if len(newewidth_choice)==0:
                return "fail",None,None,None
            else:
                return "new",newewidth_choice[0],newewidth_choice,0
        """
        return flag, new_ewidth, candidates, which_candidate
        """


# In[105]:

class HEArun:
    def __init__(self,polytyp,heakey,ewidth,  *arg):
        # update dict
        #self.dic = { "bzqlty":5 , "maxitr":500 }
        self.dic = {}
        for x in arg:
            self.dic.update(x)


        self.polytyp = polytyp
        self.heakey = heakey
        self.ewidth = ewidth
        if "ewidth_dos" in self.dic.keys():
            self.ewidth_dos = self.dic["ewidth_dos"]
        else:
            #self.ewidth_dos = ewidth + 1.0
            self.ewidth_dos = 3.0
            # default energy range -> [-ewidth*0.75: ewith*0.25]

        
    def make_heainput(self,heakey,go="go",brvtyp="bcc",ewidth=1.7):
        """
        make akaikkr input 
        
        input heakey: .e.g, "32487481"
        """
        dic = self.dic
        dic.update( {"go":go,"a":0,"brvtyp":brvtyp, "ewidth":ewidth } )
        if go=="go":
            dic.update({"record":"2nd"} )
            inp = KKRinput_brv(dic)
        elif go=="dos":
            dic.update({"record":"2nd"} )
            inp = KKRinput_brv(dic)
        elif go == "j":
            dic.update({"record":"2nd"} )
            inp = KKRinput_brv(dic)

        zlist = heakey2zlist(heakey)
        typ = ["HEA"]
        heatyp = []
        for z in zlist:
            heatyp.append([z,25])
        heatyp = [heatyp]
        inp.set_type(typ,heatyp)

        atmname = ["HEA"]
        atmicx = [[0,0,0]]
        inp.set_atom(atmname,atmicx)
        inp.make_atom_card()

        return inp

    def make_and_changedir(self,dirname):
        """
        make directory if not exist and change directory
        """
        p = dirname
        if not os.path.exists(p):
            os.mkdir(p)
        os.chdir(p)
    
    def make_GOinput_and_run(self,prefix,force=False):
        """
        make GO input and run it
        
        """
        #print(os.getcwd())
        heakey = self.heakey
        polytyp = self.polytyp
        ewidth = self.ewidth
        
        specx = self.dic["specx"]
        
        # save input
        inputname = self.dic["inputcard_go"]
        inp = self.make_heainput(heakey,go="go",brvtyp=polytyp,ewidth=ewidth,)
        p = os.path.join(prefix,inputname)
        if not os.path.exists(p):
            with open(p,"w") as f:
                f.write(str(inp))

        outputname = self.dic["out_go_log"]
        p = os.path.join(prefix,outputname)
        #print(p,"exists?",os.path.exists(p))
        execute_it = False
        if not os.path.exists(p):
            execute_it = True
        if os.path.exists(p):
            output = OutputCommon(p)
            if not output.normal_exit:
                print("exist {}, but not normal exit. Run it again") 
                execute_it = True

        if execute_it or force:
            cmd = "cd {}; {} < {} > {}".format(prefix,specx,inputname,outputname)
            subprocess.call(cmd,shell=True)
        
    def make_DOSinput_and_run(self,prefix,force=False):
        """
        make DOS input and run it
        
        """
        #print(os.getcwd())
        
        heakey = self.heakey
        polytyp = self.polytyp
        ewidth = self.ewidth_dos
        specx = self.dic["specx"]

        inputname = self.dic["inputcard_dos"]
        inp = self.make_heainput(heakey,go="dos",ewidth=ewidth,brvtyp=polytyp)
        #print("make_DOSinput_and_run",inp)
        p = os.path.join(prefix,inputname)
        if not os.path.exists(p) or force:
            with open(p,"w") as f:
                f.write(str(inp))

        outputname = self.dic["out_dos_log"]
        p = os.path.join(prefix,outputname)

        execute_it = False
        if not os.path.exists(p):
            execute_it = True
        if os.path.exists(p):
            output = OutputCommon(p)
            if not output.normal_exit:
                print("exist {}, but not normal exit. Run it again") 
                execute_it = True

        if execute_it or force:
            cmd = "cd {}; {} < {} > {}".format(prefix,specx,inputname,outputname)
            subprocess.call(cmd,shell=True)
        
    def make_Jinput_and_run(self,prefix,force=False):
        """
        make J input and run it
        """
        #print(os.getcwd())
        
        heakey = self.heakey
        polytyp = self.polytyp
        ewidth = self.ewidth
        
        specx = self.dic["specx"]
                       
        inputname = self.dic[ "inputcard_j" ]
        inp = self.make_heainput(heakey,go="j",ewidth=ewidth,brvtyp=polytyp)
        p = os.path.join(prefix,inputname)
        if not os.path.exists(p):
            with open(p,"w") as f:
                f.write(str(inp))

        outputname = self.dic["out_j_log"]
        p = os.path.join(prefix,outputname)

        execute_it = False
        if not os.path.exists(p):
            execute_it = True
        if os.path.exists(p):
            output = OutputCommon(p)
            if not output.normal_exit:
                print("exist {}, but not normal exit. Run it again") 
                execute_it = True

        if execute_it or force:
            cmd = "cd {}; {} < {} > {}".format(prefix,specx,inputname,outputname)
            subprocess.call(cmd,shell=True)        

    def run_all(self, prefix, copyfiles_from = None, stopifnotconverged = False):
        """
        run GO, DOS and J
        """
        #self.make_and_changedir(prefix)
        p = prefix
        if not os.path.exists(p):
            os.mkdir(p)
            print(p,"is made.")

        if copyfiles_from is not None:
            print("copy {} from {}".format(self.dic["file"], copyfiles_from,) )
            shutil.copy2( os.path.join( copyfiles_from, self.dic["file"]) , os.path.join(p, self.dic["file"]))
            
        heakey = self.heakey
        self.make_GOinput_and_run(prefix)
        out = OutputGo(os.path.join(prefix,self.dic["out_go_log"]),heakey)
        print("CONVERGED?",out.dic["converged"])
        self.converged = out.dic["converged"]
        if stopifnotconverged :
            if self.converged:
                self.make_DOSinput_and_run(prefix)
                self.make_Jinput_and_run(prefix)
        else:
            self.make_DOSinput_and_run(prefix)
            self.make_Jinput_and_run(prefix)            

        return self.converged

# In[106]:


class HEADos:
    def __init__(self,gofile,dosfile,plotdos=False,plotcdos=False,*arg):
        figsize=(10,1.5)            

        heakey = ""
        polytyp = ""
        elements = []
        
        self.dic = {}
        
        if len(arg)>0:
            for x in arg:
                self.dic.update(x)
            if "heakey" in arg[0]:
                heakey = arg[0]["heakey"]
            if "polytyp" in arg[0]:
                polytyp = arg[0]["polytyp"]
            if len(heakey)>0:
                elements = heakey2elements(heakey)
            
        #outgo = OutFirstDOS(os.path.join(prefix,"out_go.log"))
        outgo = OutFirstDOS(gofile)
        self.ewidth = outgo.ewidth
        #print("ewidth_in_GO",self.ewidth)
        #outdos = OutFirstDOS(os.path.join(prefix,"out_dos.log"))
        outdos = OutFirstDOS(dosfile)
        #outvalence = OutputGo(os.path.join(prefix,"out_go.log"))
        outvalence = OutputGo(gofile)        
        df_nlconfig_list = make_nlconfig_list(outgo.zlist,outgo.zconfig,outvalence.finalcoreconf)
        
        # assume nspin==2
        dossum = outdos.totaldos[0] + outdos.totaldos[1]
        dossum[:,0] = dossum[:,0]*0.5  # energyは倍になってはいけいない。半分にする。
        self.dossum = dossum
        #以下はplot
        emin = dossum[:,0].min()
        emax = dossum[:,0].max()
        
        if plotdos:
            # linear DOS
            plt.figure(figsize=figsize)
            plt.plot(dossum[:,0],dossum[:,1],label="total DOS")
            plt.plot([emin,emax],[0,0], label="0")
            plt.xlim( (emin,emax))
            plt.title("{} {} {}".format(polytyp,heakey,"".join(elements)))
            plt.ylabel("DOS")
            plt.legend()
            plt.show()

            # log DOS
            plt.figure(figsize=figsize)
            plt.plot(dossum[:,0], np.log10(dossum[:,1]))
            plt.plot([emin,emax],[-2,-2], label="0")

            plt.xlim( (emin,emax))
            plt.title("{} {} {}".format(polytyp,heakey,"".join(elements)))
            plt.ylabel("log10(DOS)")
            plt.show()


        # spin sum
        cdos = outdos.componentdos[0]+outdos.componentdos[1]
        cdos[:,0] = cdos[:,0]*0.5
        
        if plotcdos:
            # log10 atomic PDOS
            for ic in range(outdos.ncmpx):
                plt.figure(figsize=figsize)
                plt.title("{} {} {}".format(polytyp,heakey,elements[ic]))
                for lm  in range(1,4):
                    plt.plot( cdos[ic,:,0],np.log10(cdos[ic,:,lm]) ,label=str(lm-1))
                plt.legend()
                plt.xlim( (emin,emax))
                plt.ylabel("log10(DOS)")

                plt.show()        


