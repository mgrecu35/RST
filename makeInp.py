from numpy import *
from rosen import *
# absair,abswv = gasabsr98(f,tk,rhowv,pa)

def makeInp(atmosFile):
    lout=[]
    lines=open("template.INP","r").readlines()
    for l in lines:
        if "atmosphere_file" in l:
            l='atmosphere_file %s'%atmosFile
        
        if "mol_tau_file" in l:
            lines2=open(atmosFile,"r").readlines()
            h=[]
            for l2 in lines2:
                if "#" not in l2:
                    l2s=l2.split()[0:7]
                    if len(l2s)>0:
                        #print l2s
                        z,p,t,wv=float(l2s[0]),float(l2s[1]),\
                                  float(l2s[2]),float(l2s[-1])/float(l2s[3])
                        h.append([z,p,t,wv])
            h=array(h)
    return lines2,h

l2,h=makeInp("../examples/MLS70.UVSPEC")
