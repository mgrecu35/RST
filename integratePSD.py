from numpy import *
N0=1e6
lambd=3000.
D=(0.5+arange(301)*0.05)*0.1#(in cm)
dD=0.05*1e-3 #(in m)
alpha=5.74*1e-3
mD=alpha*(D)**2.1 # in g
lambdL=[100.,150.,200.,250,300,350.,400,450.,500,550.,600,650,700,750.,800,\
        850.,900.,950.,1000.,1250.,1500.,1750.,2000,2250.,2500.,2750.,3000.]

import glob
from netCDF4 import Dataset

fh=Dataset('scatdb.hdf5')
scattdB=fh['scatdb']['floatMat'][:,:]
afreq=nonzero(abs(scattdB[:,0]-340)<0.1)

import matplotlib.pyplot as plt

#flaketype,frequencyghz,temperaturek,aeffum,max_dimension_mm,cabs(m^2),cbk(m^2),cext(m^2),csca (m^2),g,ar

massDB=pi*4/3.*(scattdB[afreq[0],2]*1e-6)**3*900.*1e3
#scattdB[afreq[0],6])
dL=[]
mDL=[]
kextL=[]
kscaL=[]
gscaL=[]
for i in range(D.shape[0]-1):
    a1=nonzero((massDB-mD[i])*(massDB-mD[i+1])<0)
    if len(a1[0])>0:
        dL.append(scattdB[afreq[0][a1],3].mean()*1e-1)
        mDL.append(massDB[a1].mean())
        kextL.append((scattdB[afreq[0][a1],6]+scattdB[afreq[0][a1],4]).mean())
        kscaL.append(scattdB[afreq[0][a1],7].mean())
        gscaL.append((scattdB[afreq[0][a1],7]*scattdB[afreq[0][a1],8]).mean())
dL=array(dL)
mDL=array(mDL)

plt.plot(mDL,kextL)
ax=plt.gca()
ax.set_yscale('log')
kextL=array(kextL)
iwcf=[]
kextf=[]
kscaf=[]
gscaf=[]
for lambd in lambdL:
    M=sum(N0*exp(-lambd*dL*1e-2)*mDL*dD)
    kext=sum(N0*exp(-lambd*dL*1e-2)*kextL*dD)
    ksca=sum(N0*exp(-lambd*dL*1e-2)*kscaL*dD)
    g=sum(N0*exp(-lambd*dL*1e-2)*gscaL*dD)
    iwcf.append(M)
    kextf.append(kext)
    kscaf.append(ksca)
    gscaf.append(g)
plt.figure()
kscaf=array(kscaf)
kextf=array(kextf)
gscaf=array(gscaf)
#plt.plot(iwcf,kextf)
plt.plot(iwcf,kscaf/kextf)
plt.plot(iwcf,gscaf/kscaf)

resKext=polyfit(log(iwcf),log(kextf),2)
resKsca=polyfit(log(iwcf),log(kscaf/kextf),2)
resG=polyfit(log(iwcf),gscaf/kscaf,2)

print resKext
print resKsca
print resG
import pickle
pickle.dump([resKext,resKsca,resG],open('scaCoeff_340.pklz','wb'))
ax=plt.gca()
#ax.set_yscale('log')
#ax.set_xscale('log')
