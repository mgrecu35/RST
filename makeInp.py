from numpy import *
from grecuRT import *
import matplotlib.pyplot as plt
# absair,abswv = gasabsr98(f,tk,rhowv,pa)
from netCDF4 import Dataset
f=Dataset('iwc19990717.nc')
iwc=f['iwc'][:,:,:]*1.
x=f['x'][:]
z=f['z'][:]
import pickle
[resKext,resKsca,resG]=pickle.load(open('scaCoeff_340.pklz','rb'))

def makeInp(atmosFile,wc,freq):
    lout=[]
    lines=open("template.INP","r").readlines()
    f2=open("actual.INP","w")
    for l in lines:
        if "atmosphere_file" in l:
            l='atmosphere_file %s\n'%atmosFile
        
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
        if "wavelength" in l and "output" not in l:
            l='wavelength %8.2f \n'%(300./freq*1e6)
        f2.write(l)
    f=open('molecular_tau_zero','w')
    rho=h[:,1]/287./h[:,2]*1e2
    #print h[0:10,1]*100.
    for h1,i in zip(h[1:,0],range(1,h.shape[0])):
        tau=0.001

        rhowv=h[i,3]*rho[i]
        absair,abswv = gasabsr98(340.,h[i,2],rhowv,h[i,1]*100.)
        tau=(absair+abswv)
        if tau<0.000000001:
            tau=0.000000001
        #print '%7.2f %11.9f %10.6f'%(h[i,0],tau,rhowv)
        f.write('%7.2f %11.9f\n'%(h[i,0],tau))
    return lines2,h
import os
j=100
from numpy import *
tbL=[]
for i in range(10,210):
    print i
    freq=340.
    l2,h=makeInp("../examples/MLS70.UVSPEC_short",iwc[:,100,i],freq)
    iwcavg=zeros((20),float)
    cavg=zeros((20),float)
    for k in range(iwc.shape[0]):
        if iwc.mask[k,j,i]==False:
            k0=int(z[k]/1000.+0.5)
            iwcavg[k0]+=iwc[k,j,i]
            cavg[k0]+=1
            
    a=nonzero(cavg>0)
    iwcavg[a]=iwcavg[a]/cavg[a]
    ihyd=1
    if ihyd==1:
        f=open("moment_files","w")
        if len(a[0])>0:
            f.write("%6.0f null.layer\n"%(a[0][-1]+1))
        else:
            f.write("%6.0f null.layer\n"%(13))
            f.write("%6.0f null.layer\n"%(12))
        fopt=open("null.layer","w")
        output='%8.2f '%(300./freq*1e6)
        output+='%8.5f %6.4f %6.4f %6.4f %6.4f %6.4f'%(0,0.,1.,1,1**2,1**3)
        fopt.write(output)
        fopt.close()
        for k in range(20)[::-1]:
            if iwcavg[k]>0:
                kext=exp(polyval(resKext,log(iwcavg[k])))
                ssa=exp(polyval(resKsca,log(iwcavg[k])))
                g=(polyval(resG,log(iwcavg[k])))
                print k,iwcavg[k],kext*1e3,ssa,g
                f.write("%6.0f optical_%2.2i.layer\n"%(k*1.,k))
                fopt=open("optical_%2.2i.layer"%k,"w")
                output='%8.2f '%(300./freq*1e6)
                output+='%8.5f %6.4f %6.4f %6.4f %6.4f %6.4f'%(kext*1e3,ssa,1.,g,g**2,g**3)
                fopt.write(output)
                fopt.close()
        f.close()
    else:
        f=open("moment_files","w")
        f.write("%6.0f null.layer\n"%(a[0][-1]+1))
        fopt=open("null.layer","w")
        output='%8.2f '%(300./freq*1e6)
        output+='%8.5f %6.4f %6.4f %6.4f %6.4f %6.4f'%(0,0.,1.,1,1**2,1**3)
        fopt.write(output)
        fopt.close()
        ic=0
        for k in range(20)[::-1]:
            if iwcavg[k]>0 and ic<1:
                kext=exp(polyval(resKext,log(iwcavg[k])))
                ssa=exp(polyval(resKsca,log(iwcavg[k])))
                g=(polyval(resG,log(iwcavg[k])))
                print k,iwcavg[k],kext*1e3,ssa,g
                f.write("%6.0f optical_%2.2i.layer\n"%(k*1.,k))
                fopt=open("optical_%2.2i.layer"%k,"w")
                output='%8.2f '%(300./freq*1e6)
                output+='%8.5f %6.4f %6.4f %6.4f %6.4f %6.4f'%(kext*1e3,ssa,1.,g,g**2,g**3)
                fopt.write(output)
                fopt.close()
                ic+=1
        f.close()

    ilam=1
    umu=0.60181
    rho=h[:,1]/287./h[:,2]*1e2
    kextL=[]
    salbL=[]
    asymL=[]
    tempL=[]
    tempL.append(h[-1,2])
    for h1,i in zip(h[::-1,0],range(h.shape[0])[::-1]):
        tau=0.001
        rhowv=h[i,3]*rho[i]
        if h1<21.:
            absair,abswv = gasabsr98(340.,h[i,2],rhowv,h[i,1]*100.)
            print absair,abswv,h1
            k=24-i
            
            if k<20:
                #print h1, k
                if iwcavg[k]>0:
                    kext=exp(polyval(resKext,log(iwcavg[k])))*1000.
                    ssa=exp(polyval(resKsca,log(iwcavg[k])))
                    g=(polyval(resG,log(iwcavg[k])))
                    print 'kext=',kext,ssa,g, k, h1
                else:
                    kext=0.
                    ssa=0.
                    g=0.
            else:
                kext=0.
                ssa=0.
                g=0.
           
            kextold=kext
            kext+=(absair+abswv)
            kextL.append(kext)
            salbL.append(ssa)
            asymL.append(g)
            tempL.append(h[i:i+2,2].mean())
    nlyr=21
    dz=1.
    btemp=tempL[0]
    lyrhgt=arange(nlyr+1)*dz
    fisot=2.7
    albedo=0.027
    emis=1-albedo
    ebar=emis
    ilam=0
    tb = radtran(umu,nlyr,btemp,array(tempL),lyrhgt,array(kextL),\
                 array(salbL),array(asymL),fisot,emis,ebar)


    os.system("../bin/uvspec< actual.INP > out 2> out2")
    tb2=float(open("out","r").readlines()[-1].split()[1])
    print tb,tb2
    tbL.append([tb,tb2])
