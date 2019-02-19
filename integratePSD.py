from numpy import *
N0=1e6
lambd=3000.
D=(1+arange(101)*0.1)*0.1 #(in cm)
dD=0.1*1e-3 #(in m)
alpha=5.74*1e-3
mD=alpha*(D)**2.1 # in g
lambdL=[300,350.,400,450.,500,550.,600,650,700,750.,800,\
        850.,900.,950.,1000.,1250.,1500.,1750.,2000,2250.,2500.,2750.,3000.]

import glob
from netCDF4 import Dataset

fh=Dataset('scatdb.hdf5')
scattdB=fh['scatdb']['floatMat'][:,:]
afreq=nonzero(abs(scattdB[:,0]-340)<0.1)

#flaketype,frequencyghz,temperaturek,aeffum,max_dimension_mm,cabs(m^2),cbk(m^2),cext(m^2),csca,g,ar
massDB=pi*4/3.*(scattdB[afreq[0],2]*1e-6)**3*900.*1e3
import matplotlib.pyplot as plt
fig=plt.figure()
ax=plt.gca()
plt.scatter(massDB,scattdB[afreq[0],6])
ax.set_yscale('log')
for i in range(D.shape[0]-1):
    a1=nonzero((massDB-mD[i])*(massDB-mD[i+1])<0)
    if len(a1[0])>0:
        print scattdB[afreq[0][a1],2].mean(),massDB[a1].mean(), i
for lambd in lambdL:
    M=sum(N0*exp(-lambd*D*1e-2)*mD*dD)
    print M
