
import copy
class KKRinput_brv:
    def __init__(self,dic):
        
        self.defaultdic = {"go":"go","file":"pot.dat",            "a":0, "c/a":1.0, "b/a": 1.0,              "alpha":90, "beta":90, "gamma":90,              "edelt":0.0001, "ewidth":1.7, "reltyp":"sra", "sdftyp":"pbe",              "magtyp":"mag", "record":"init",              "outtyp":"update", "bzqlty":10, "maxitr":1000,"pmix":0.005 }
        
        self.dic = copy.deepcopy(self.defaultdic)
        self.dic.update(dic)
        
    def get_defalut_param_typ(self):
        dic = {"rmt":1,"field":0.0,"mxl":2}
        return dic
    
    def set_type(self,typ,typconc):
#        for ityp,itypeconc in zip(typ,typconc):
#            print("TYP:",ityp,itypeconc)
        self.dic["ntyp"] = len(typ)
        self.dic["typ"] = typ
        self.dic["typconc"] = typconc
        
    def make_type_card(self):
        card = []
        for ityp,itypeconc in zip(self.dic["typ"],self.dic["typconc"]):
            card.append("#--- type ncmp")
            s = [ityp,len(itypeconc)]
            card.append(" ".join(list(map(str,s))))
            s = []
            dic = self.get_defalut_param_typ()
            head =["#-"]
            for key in [ "rmt"  ,  "field"  , "mxl"]:
                head.append(key)
                s.append(str(dic[key]))
            card.append(" ".join(head))
            card.append(" ".join(s))
            card.append("#- anclr conc")
            for eachtypeconc in itypeconc:
                card.append( " ".join(list(map(str,eachtypeconc))))
            
        return card

    def set_atom(self,atmname,atmicx,):
        self.dic["atmicx"] = copy.deepcopy(atmicx)
        self.dic["natm"] = len(atmicx)
        self.dic["atmtyp"] = atmname
#        for iatm,iatmname in zip(atmicx,atmname):
#            print("ATOM:",iatm,iatmname)
            
    def make_atom_card(self):
        s = ["#---","atmicx","atmtyp"]
        card = [" ".join(s)]
        for iatmicx,iatmname in zip(self.dic["atmicx"],self.dic["atmtyp"]):
            s = copy.deepcopy(iatmicx)
            s.append(iatmname)
            s = list(map(str,s))
            card.append(" ".join(s))
        #print(card)
        return card
    
    def __str__(self):
        card = []
        
        s = []
        head = ["#---"]
        for key in ["go","file"]:
            head.append(key)
            s.append(str(self.dic[key]))
        card.append(" ".join(head))
        card.append(" ".join(s))
        
        s = []
        head = ["#---"]
        for key in ["brvtyp",     "a",         "c/a",   "b/a",   "alpha",   "beta",   "gamma"]:
            head.append(key)
            s.append(str(self.dic[key]))
        card.append(" ".join(head))
        card.append(" ".join(s))
        
        s = []
        head = ["#---"]
        for key in ["edelt",    "ewidth",    "reltyp",   "sdftyp",   "magtyp",   "record"]:
            head.append(key)
            s.append(str(self.dic[key]))
        card.append(" ".join(head))
        card.append(" ".join(s))        
        
        s = []
        head = ["#---"]
        for key in ["outtyp",    "bzqlty",   "maxitr",   "pmix"]:
            head.append(key)
            s.append(str(self.dic[key]))
        card.append(" ".join(head))
        card.append(" ".join(s))
        
        s = []
        head = ["#---"]
        for key in ["ntyp"]:
            head.append(key)
            s.append(str(self.dic[key]))
        card.append(" ".join(head))
        card.append(" ".join(s))        
        
        card.extend(self.make_type_card())
        
        s = []
        head = ["#---"]
        for key in ["natm"]:
            head.append(key)
            s.append(str(self.dic[key]))
        card.append(" ".join(head))
        card.append(" ".join(s))      
        
        card.extend(self.make_atom_card())

        card.append("")
        
        return "\n".join(card)


# In[ ]:





# In[266]:
if __name__ == "__main__":

    inp = KKRinput_brv({"brvtyp":"fcc","c/a":10})


    # In[267]:


    typ = ["HEA"]
    heatyp = []
    for i in range(4):
        heatyp.append(["atom{}".format(i),25])
    heatyp = [heatyp]
    print("heatype",heatyp)

    inp.set_type(typ,heatyp)
#   print(inp.make_type_card())


    # In[270]:


    atmname = ["HEA","ATOM2"]
    atmicx = [[0,0,0],[0.5,0.5,0.5]]
    atmname = ["HEA"]
    atmicx = [[0,0,0]]
    inp.set_atom(atmname,atmicx)
    inp.make_atom_card()


    # In[271]:


    print(inp)


    # In[ ]:





    # In[ ]:





    # In[ ]:





    # In[ ]:




