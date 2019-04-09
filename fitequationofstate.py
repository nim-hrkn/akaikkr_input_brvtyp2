import sys
import subprocess
import numpy as np
import os

def fit_Murnaghan_eq(tmpdir,datafile):
    p = datafile
    with open(p) as f:
       lines = f.readlines()

    x = []
    y = []
    for line in lines:
        s = line.split()
        if len(s)==0:
           break
        v = list(map(float,s))
        x.append(v[0])
        y.append(v[-1])

    x = np.array(x)
    y = np.array(y)
    n = x.shape[0]
    i = y.argmin()
    iargmin = i
    # estimtation of ymax = y[i+4] or y[i-4]
    i1 = i-4
    i2 = i+4
    if i1>=0:
        imax = i1
    elif i2<n:
        imax = i2
    yrange = y[imax]-y.min()
    narg = n
    E0 = str(y[i])
    V0 = str(x[i]**3)
    #print("ymin at ", i, "of ",n)
    if False:
        if i==n-1 or i==0:
            print("error on y.min range",datafile,i,n)

    script = """
    V0 = V0VALUE
    E0 = E0VALUE
    K0 = 0.0039
    K0d = 9.2
    Murnaghaneq(V) = E0 + K0*V0*( 1/( K0d*(K0d-1)) *(V/V0)**(1-K0d) + V/(K0d*V0) - 1/(K0d-1)  )

    n=3
    fit Murnaghaneq(x) "DATAFILE" u ($1)**n:5 via E0,V0,K0,K0d

    """
    script_plot = """
    plot "DATAFILE" u ($1)**n:5, Murnaghaneq(x) 
    pause -1
    """

    #script += script_plot

    script = script.replace("V0VALUE",V0)
    script = script.replace("E0VALUE",E0)
    script = script.replace("DATAFILE",datafile)

    with open(os.path.join(tmpdir,"fit.gnu"),"w") as f:
        f.write(script) 

    cmd = "gnuplot {}/fit.gnu".format(tmpdir)
    f1 = open(os.path.join(tmpdir,"stdout_data"),"wb")
    f2 = open(os.path.join(tmpdir,"stderr_data"),"wb")
    p = subprocess.Popen(cmd,shell=True, stdout=f1, stderr=f2)
    p.wait()
    #p = subprocess.call(cmd,shell=True)
    with open("fit.log") as f:
        lines = f.readlines()

    dic = {"argmin":iargmin,"narg":narg,"yrange":yrange}
    for line in lines:
        #print(line)
        if "fit converged" in line:    
           dic["fit_converged"] = True
        elif "rms of residuals" in line:
           s = line.split()
           dic["RMSE"] =  float(s[-1])

    istart = None
    for i,line in enumerate(lines):
        if line.startswith("Final set of parameters"):
           istart = i

    if istart is None:
        print("failed to setup istart")
        raise

    for line in lines[istart:]:
        s = line.split()
        if line.startswith("E0 "):
            x =  s[2]
            var = s[4]
            dic["E0"] = list(map(float,[x,var]))
        elif line.startswith("V0 "):
            x =  s[2]
            var = s[4]
            dic["V0"] = list(map(float,[x,var]))
        elif line.startswith("K0 "):
            x =  s[2]
            var = s[4]
            dic["K0"] = list(map(float,[x,var]))
        elif line.startswith("K0d"):
            x =  s[2]
            var = s[4]
            dic["K0d"] = list(map(float,[x,var]))
            break    

    #print(dic)
    return dic

if __name__ == "__main__":
    fit_Murnaghan_eq("xy")
    fit_Murnaghan_eq("xy2")
    fit_Murnaghan_eq("xy3") 



