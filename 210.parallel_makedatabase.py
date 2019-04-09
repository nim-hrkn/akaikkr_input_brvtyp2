
# coding: utf-8

# In[81]:

import numpy as np
import os
import subprocess
import itertools
import yaml
import sys

#os.chdir("/home/kino/kino/work/fukushima_HEA/script")

# define components of HEA
atom_2p=np.arange(13,15)
atom_3d=np.arange(22,33)
atom_4d=np.arange(40,51)
atom_5d=np.arange(72,84) # include Bi=83
atomlist=np.hstack([atom_3d,atom_2p,atom_4d,atom_5d])
print("atomlist",atomlist)
# check total number
for i,v in enumerate(itertools.combinations(atomlist, 4)):
    pass
print("allthecombinations = ",1+i)

# load control parameter

if len(sys.argv)!=2:
    print("argument error")
    sys.exit(100)

n1 = int(sys.argv[1])
controlparam = {"n1":n1, "n2":n1+1 }

print("control",controlparam)


resultdir = "result_p"
if not os.path.isdir(resultdir):
    os.mkdir(resultdir)


# In[82]:

dryrun = False

# number of points for optimization, range of search
ntime=10 ; arange=1.0

# AkaiKKR executable file
unlimit = ' '
#EXEPATH='/home/kino/work/fukushima_HEA/akaikkr_jij'
EXEPATH='/home/kino/work/fukushima_HEA/akaikkr_jij'
EXE='specx'

basedir = os.getcwd()
print("basedir",basedir)
print("resultdir",resultdir)


# In[83]:


# define function
def search_str(str1,str2):
        check=str1 in str2
        return check

# read inputcard.orig
#input_data_orig=open('inputcard.orig').read()
input_data_orig = """
#-------------------HEA--------------------------------------
     go.exe pot.dat
#------------------------------------------------------------
#   brvtyp     a         c/a   b/a   alpha   beta   gamma
     bcc    99999,  1.00 , 1.00  ,   90 ,   90  , 90  ,
#------------------------------------------------------------
#   edelt    ewidth    reltyp   sdftyp   magtyp   record
    0.0001    1.7         sra     pbe      mag    init
#------------------------------------------------------------
#   outtyp    bzqlty   maxitr   pmix
    update     10        1000     0.005
#------------------------------------------------------------
#    ntyp
      1
#------------------------------------------------------------
#   type     ncmp    rmt    field   mxl  anclr   conc
    HEA       4      1    0.000     2     atom1   25
                                          atom2   25
                                          atom3   25
                                          atom4   25
#------------------------------------------------------------
#   natm
     1
#------------------------------------------------------------
#   atmicx                            atmtyp
     0       ,   0       ,   0       , HEA
#------------------------------------------------------------
"""

#print("have read inputdata")


# In[84]:


# loop for components: atom1, atom2, atom3, atom4
for i,atoms in enumerate(itertools.combinations(atomlist, 4)):

    if controlparam["n1"]<=i and  i<controlparam["n2"]:
        
        atoms = sorted(atoms)
        n = i 
        dirname = basedir+"/"+resultdir+"/"+ "".join(map(str,atoms))

        if os.path.exists(dirname):
            print("\nskip ",i,atoms,dirname)
            continue
            
        print("\ngenerating ",i,atoms,dirname)

        try:
            os.mkdir(dirname)
        except FileExistsError:
            print("failed to make the directory ",dirname)
            continue
        
        os.chdir(dirname)


        # set exe and elements
        tmpfile = 'hoge1.tmp'
        atom1,atom2,atom3,atom4 = atoms
        input_data=input_data_orig.replace('go.exe','go')        .replace('atom1',str(atom1))        .replace('atom2',str(atom2))        .replace('atom3',str(atom3))        .replace('atom4',str(atom4))
        open(tmpfile,'w').write(input_data)

        if not dryrun: 

            inputfile = "inputcard_disp"
            outputfile = 'out_disp.log'

            # generate inputfile for AkaiKKR
            cmd= unlimit+ EXEPATH+'/make_inputfile.exe '+str(ntime)+' '+str(arange)        +' < {} > {}'.format(tmpfile,inputfile)
            print("execute ",cmd)
            subprocess.run(cmd,shell=True)
				

            # execute optimization by AkaiKKR
            cmd=unlimit+' '+EXEPATH+'/'+EXE+' < {} > {}'.format(inputfile,outputfile)
            print("execute" ,cmd)
            subprocess.run(cmd,shell=True)

            # fitting by lsm 
            data1_tmp=np.loadtxt('pot.dat.info')
            a_latt=data1_tmp[:,0] ; total_energy=data1_tmp[:,1]
            data2_tmp=np.polyfit(a_latt,total_energy,2)
            print("fitted coeffficient",data2_tmp)
            a_latt_opt=-data2_tmp[1]/(2.0*data2_tmp[0])
            print("fitted lattice length",a_latt_opt)

            # execute SCF for optimized lattice constant

            inputfile = "inputcard_equiv"
            outputfile = 'out_equiv.log'

            input_data=input_data.replace('99999',str(a_latt_opt))
            open(inputfile,'w').write(input_data)
            if os.path.exists(EXEPATH+'/'+EXE):
                cmd=unlimit+' '+EXEPATH+'/'+EXE+' < {} > {}'.format(inputfile,outputfile)
                print("execute",cmd)
                subprocess.run(cmd,shell=True)
            else:
                print("failed to find ",EXEPATH+'/'+EXE)
                raise

            # execute jij calculation

            inputfile = "inputcard_jij"
            outputfile = 'out_jij.log'

            input_data=input_data.replace('go','j').replace('init','2nd')
            open(inputfile,'w').write(input_data)
            cmd=unlimit+' '+EXEPATH+'/'+EXE+' < {} > {}'.format(inputfile,outputfile)
            print("execute",cmd)
            subprocess.run(cmd,shell=True)

            # gather results
            output=open(outputfile).readlines()
            moment_comp=[]
            for i in range(len(output)):
                    if search_str('total energy',output[i]):
                            str1=output[i].split() ; te=str1[2]
                    if search_str('chr,spn',output[i]):
                            str1=output[i].split() ; moment=str1[6]
                    elif search_str('spin moment',output[i]):
                            str1=output[i].split() ; moment_comp.append(str1[2])
                    elif search_str('mean field approximation',output[i]):
                            str1=output[i].split() ; Tc=str1[6]
                            Tc=Tc.split('K') ; Tc=Tc[0]

            # write component and data
            try:
                    open('component','a')                    .write('{:>8}{:>6d}{:>2d}{:>2d}{:>2d}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>18}{}'                    .format(n,atom1,atom2,atom3,atom4                    ,moment_comp[0],moment_comp[1],moment_comp[2],moment_comp[3]                    ,moment,Tc,te,'\n'))
            except IndexError:
                    open('component','a')                    .write('{:>8}{:>6d}{:>2d}{:>2d}{:>2d}   === error ===   {}'                    .format(n,atom1,atom2,atom3,atom4                    ,'\n'))
                    
            # delete unnecessary files
            unlinkfiles = [tmpfile]
            for x in unlinkfiles:
                os.unlink(x)
            print("delete",unlinkfiles)
                

        print("chdir",basedir)
        os.chdir(basedir)
        


# In[85]:

subprocess.run("ls")


# In[ ]:



