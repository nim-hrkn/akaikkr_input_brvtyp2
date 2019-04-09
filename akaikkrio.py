import numpy as np
import os
import sys
import copy
import glob
import progressbar 
import fitequationofstate
import matplotlib.pyplot as plt
import pandas as pd

class CompleteCheck:
    def __init__(self):
        pass

    def execute(self,dirs):
        files = ["out_jij.log","out_disp.log", "out_equiv.log"]
        notcompletes = []
        for x in dirs:
            if self.isIncomplete(x):
                notcompletes.append(x)

        print(len(notcompletes))
        notcompletes2 = []
        for x in notcompletes:
            y = os.path.split(x)
            notcompletes2.append(y[-1])

        filename = "notcompleted"
        with open(filename,"w") as f:
            f.write("\n".join(notcompletes2))
        print("filename:{} is made".format(filename))

    def isIncomplete(self,dir_):
        files = ["out_jij.log","out_disp.log", "out_equiv.log"]
        notcompletes = []
        x = dir_
        for y in files:
            z = os.path.join(x,y)
            if not os.path.exists(z):
                notcompletes.append(x)
        notcompletes = set(notcompletes)
        if len(notcompletes)==0:
            return False
        else:
            return True



class OutJij: 
    def __init__(self):
        self.dic = {}
    def read(self,outputfile,key ):
        """
          analyze outjij.out
        """
        debug_ = False
        outputlist = open(outputfile).readlines()
        moment_comp=[]
        for output in outputlist:
                if 'total energy' in output:
                        str1=output.split() ; te=float(str1[2])
                if 'chr,spn' in output:
                        str1=output.split() ; moment=float(str1[6])
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

        self.dic = {"totalenergy":te, "moment":moment,"moment_comp": moment_comp, "Tc": Tc, "alen":alen,"vol":vol ,"key":key }

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

    def __init__(self):
        self.dic = {}

    def short(self):
        sss = []
        for d in self.dic["property"]:
            s = [ d["alen"] ,d["len_h_err"], d["h_err"][-1], d["h_moment"][-1],d["totalenergy"] ]
            s.extend(d["moment"]) 
            sss.append(" ".join(list(map(str,s))))
        return "\n".join(sss)

#    def test_dirs():
#        for id_,d in enumerate(dirs):
#            key = os.path.split(d)[-1]
#            p = d+'/out_jij.log'
#            print(read_jij(p,key,id_))
#            if id_>=3:
#               break

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

    def read(self,filename,key):
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
                elif 'chr,spn' in output:
                        str1=output.split() ; 
                        for i,x in enumerate(str1):
                            if x =="chr,spn":
                                moment = str1[i+2]
                                dic["totalmoment"] = float(moment)
                                break
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
            dic["len_h_err"] = len(errlist)
            dic["h_err"] = errlist[-1:]
            dic["h_moment"] = h_momentlist[-1:]
            dic["h_te"] = h_telist[-1:]

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

def test_jij_disp():
    # test code
    key = "24307274"
    filename = "/home/kino/kino/work/fukushima_HEA/result/result_4elements/{}/out_jij.log".format(key)

    disp = OutJij()
    disp.read(filename,key)
    print(disp)
    #
    #
    filename = "/home/kino/kino/work/fukushima_HEA/result/result_4elements/{}/out_disp.log".format(key)

    disp = OutDisp()
    disp.read(filename,key)
    print(disp.short())
#test_jij_disp()

class DispDB:
    def __init__(self,collection):
        self.collection = collection 
    def insert_all(self,dirmask):
        print("start insert_all",dirmask)
        use_progressbar = True
        collection = self.collection
        dirs = glob.glob(dirmask)
        nmax = len(dirs)
        if use_progressbar:
            p = progressbar.ProgressBar(0,nmax, widgets=[
                    progressbar.CurrentTime(), ':',
                        '(', progressbar.Counter(), ' of {}) '.format(nmax),
                         progressbar.Bar(), ' ', progressbar.ETA(), ])
        iadd = 0 
        for i,d in enumerate(dirs):
            if use_progressbar:
                p.update(i+1)
            key = os.path.split(d)[-1]
            one = collection.find_one({"key":key}) 
            if one is None:
                incomplete = CompleteCheck().isIncomplete(d) 
                if not incomplete:
                    disp = OutDisp()
                    disp.read(d+'/out_disp.log',key)
                    dic = disp.dic 
                    dic.update({"directory":d})
                    collection.insert_one(dic)
                    iadd +=1
        print("inserted {} items".format(iadd))

    def count(self):
        post = self.collection.count_documents({})
        print("n=",post)
        return post

    def show_all(self):
        for x in self.collection.find({}):
            print(x["key"])

    def make_eqos_fig(self,key,prefix="",action=["save"]):
        print("key",key)
        collection = self.collection
        prop = collection.find_one({"key":key})
        x = []
        y = []
        for dic  in prop["property"]:
            print( dic["alen"],len(dic["h_err"]), dic["h_err"][-1], dic["h_moment"][-1],dic["totalenergy"], dic["moment"],dic["valence_charge"] )
            x.append( dic["alen"] )
            y.append( dic["totalenergy"] )
        plt.figure()
        plt.plot(x,y,".")
        plt.title(key)
        savefile = os.path.join(os.path.join(prefix,str(key)),"aE.png")
        if "save" in action:
            print("savefile",savefile)
            plt.savefig(savefile)
        if "plot" in action:
            plt.show()


def test_dispDB():
    from pymongo import MongoClient
    client = MongoClient('mongodb://localhost:27017/')
    db = client['HEA4']
    collection = db['results']
    #collection.drop()
    #collection2 = db['eqos']

    dispDB = DispDB(collection)
    dispDB.insert_all("/home/kino/kino/work/fukushima_HEA/result/result_4elements/[0-9]*")
    dispDB.count()
    #dispDB.show_all()

    print("done")

def test_dispDB_short():
    from pymongo import MongoClient
    client = MongoClient('mongodb://localhost:27017/')
    db = client['HEA4']
    collection = db['results']
    #collection.drop()
    #collection2 = db['eqos']

    dispDB = DispDB(collection)
    dispDB.insert_all("/home/kino/kino/work/fukushima_HEA/result/result_4elements/[0-9]*")
    dispDB.count()
    #dispDB.show_all()

    print("done")


#test_dispDB()
#sys.exit(0)

            
class EqosDB:
    def __init__(self,collection,collection2):
        self.collection = collection
        self.collection2 = collection2

    def make(self):
        collection = self.collection
        collection2 = self.collection2
        n = collection.count({}) 
        nmax = n
        use_progressbar = True
        if use_progressbar:
            p = progressbar.ProgressBar(0,nmax, widgets=[
                    progressbar.CurrentTime(), ':',
                        '(', progressbar.Counter(), ' of {}) '.format(nmax),
                         progressbar.Bar(), ' ', progressbar.ETA(), ])
 
        iadd = 0
        for i,d in enumerate(collection.find()): 
            if use_progressbar:
                p.update(i+1)
            key = d["key"]
            conv = d["converged"]
            one = collection2.find_one({"key":key}) 
            if one is None:
                #print("add {} of {}, key:{}".format(i,n,key))
                data = []
                for dic in d["property"]:
                    data.append([ dic["alen"],len(dic["h_err"]), dic["h_err"][-1], dic["h_moment"][-1],dic["totalenergy"] ] ) 
                tmp_path = "temporary"
                if os.path.exists(tmp_path):
                    #print("tmp_path",tmp_path)
                    write_path = os.path.join(tmp_path,'xy')
                    with open(write_path, 'w') as f:
                        for d in data:
                            f.write(" ".join(map(str,d))+"\n")
                    try:
                        dic = fitequationofstate.fit_Murnaghan_eq(tmp_path,write_path)
                    except:
                        continue
                    dic.update({"key":key, "converged":conv})
                    dic["argmin"] = int(dic["argmin"])
                    collection2.insert_one(dic)
                    iadd += 1
            else:
               #print("skip {}".format(key))
               pass 
        print("inserted {} items".format(iadd))

    def show_all_rmse(self):
        collection2 = self.collection2
        keylist = []
        rrmselist = []
        iargminlist = []
        for dic in collection2.find():
            key = dic["key"]
            yrange = dic["yrange"]
            iargminlist.append( dic["argmin"] )
            rrmselist.append( dic["RMSE"]/yrange )
            keylist.append(key)
            
        if True:
            plt.figure()
            plt.xlabel("RMSE")
            plt.ylabel("occurence")
            plt.hist(rrmselist,range=(0,1),bins=100)
            plt.show()

            plt.figure()
            plt.xlabel("argmin")
            plt.ylabel("renormalized RMSE")
            plt.plot(iargminlist,rrmselist,".")
            plt.show()
 
def test_EqosDB():
    from pymongo import MongoClient
    client = MongoClient('mongodb://localhost:27017/')
    db = client['HEA4']
    collection = db['results']
    #collection.drop()
    collection2 = db['eqos']
    #collection2.drop()

    db2eqos = EqosDB(collection,collection2)
    db2eqos.make()
    #db2eqos.show_all_rmse()

#test_EqosDB()
#sys.exit(0)

def test_eqos_figure():
    from pymongo import MongoClient
    client = MongoClient('mongodb://localhost:27017/')
    db = client['HEA4']
    collection = db['results']
    #collection.drop()
    collection2 = db['eqos']
    #collection2.drop()

    dispDB = DispDB(collection)
    #db2eqos = EqosDB(collection,collection2)
 
    xlist = []
    alllist = collection2.find({})
    for x in alllist:
        key = x["key"]
        rrmse = x["RMSE"] / x["yrange"] 
        #if rrmse > 0.2:
        #     dispDB.make_eqos_fig(key,action=["plot"])
        xlist.append([int(key),float(rrmse),float(x["RMSE"]),float(x["yrange"])])
    xlist = np.array(xlist)
    print(xlist.shape)
    print(xlist)
    if True:
        hh,edge = np.histogram(xlist[:,1],bins=30)
        print(hh)
        print(edge)
    plt.figure()
    plt.hist(xlist[:,1],bins=30,range=(0,0.5))
    plt.show()

    alllist = collection2.find({})
    for x in alllist:
        key = x["key"]
        rrmse = x["RMSE"] / x["yrange"] 
        if rrmse > 0.1 and rrmse< 0.2  and x["converged"] :
             print("rrmse",rrmse,"converged",x["converged"])
             dispDB.make_eqos_fig(key,action=["plot"])
 

#test_eqos_figure()
#sys.exit(0)
def extract_shortresult(path1):

    outdisp = OutDisp()
    outjij = OutJij()
    key = os.path.split(path1)[-1]
    try:
        path = path1 + "/out_disp.log"
        outdisp.read(path,key)
        if outdisp.dic["converged"]:
            path = path1+"/out_equiv.log"
            outdisp.read(path,key)
            if outdisp.dic["converged"]:
                path = path1+"/out_jij.log"
                outjij.read(path,key)
                #print(outjij)
                return outjij.__str__(",")
    except:
        pass
    return key

def test_permally():

    path0 = "/home/kino/kino/work/fukushima_HEA/result/result.permalloy/*/result/*"
    dirs = glob.glob(path0)

    for path1 in dirs:
        s = extract_shortresult(path1)
        print(s)

def test_4element():

    path0 = "/home/kino/kino/work/fukushima_HEA/result/result_4elements/[0-9]???????"
    dirs = glob.glob(path0)

    nmax = len(dirs)
    p = progressbar.ProgressBar(0,nmax, widgets=[
            progressbar.CurrentTime(), ':',
            '(', progressbar.Counter(), ' of {}) '.format(nmax),
             progressbar.Bar(), ' ', progressbar.ETA(), ])

    f = open("result_4element.csv","w")
    for i,path1 in enumerate(dirs):
        p.update(i+1)
        s = extract_shortresult(path1)
        f.write(s+"\n")


if __name__ == "__main__":

    #test_dispDB()
    #sys.exit(0)

    #test_EqosDB()
    #sys.exit(0)
    #test_eqos_figure()
    #sys.exit(0)

    #path0 = "/home/kino/kino/work/fukushima_HEA/result/result_4elements/{}/out_disp.log".format(key)

    #test_permally()
    test_4element()

