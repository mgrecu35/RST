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
import glob
fs=glob.glob('/media/grecu/ExtraDrive1/SSD/TotallyRandom/Ice/Aggregates/Pristine/'+\
             'EvansSnowAggregate_Id1/TotallyRandom/*')

from netCDF4 import Dataset
m_eqL=[]
cscaL=[]
#g=1/2 \integral P(cos(theta) cos(theta) d(cos(theta))
# \integral P(cos(theta) d(cos(theta))=2
theta=arange(361)*0.5
singleCrys=['5bulletRosette_Id2','6BulletRosette_Id6','BlockColumn_Id12'\
            'ColumnType1_Id7','Flat3BulletRosette_Id5',\
            'Flat4BulletRosette_Id11',\
            'GemCloudIce_Id31','IceSphere_Id24',\
            'IconCloudIce_Id27','LongColumn_Id14',\
            'Perpendicular3BulletRosette_Id4','Perpendicular4BulletRosette_Id10',
            'PlateType1_Id9','SectorSnowflake_Id3',\
            'ShortColumn_Id13','ThickPlate_Id15',\
            'ThinPlate_Id16']
aggregates=['8ColumnAggregate_Id8','EvansSnowAggregate_Id1',
            'GemHail_Id34','GemSnow_Id32',
            'IconSnow_Id28', 'LargeBlockAggregate_Id22',
            'LargeColumnAggregate_Id18',  'LargePlateAggregate_Id20',
            'SmallBlockAggregate_Id21',  'SmallColumnAggregate_Id17',
            'SmallPlateAggregate_Id19', 'TyynelaDendriteAggregate_Id26']

for s in aggregates:
    fs=glob.glob('/media/grecu/ExtraDrive1/SSD/TotallyRandom/Ice/SingleCrystals/Pristine/'+\
                 '%s/TotallyRandom/*'%s)
    fs=glob.glob('/media/grecu/ExtraDrive1/SSD/TotallyRandom/Ice/Aggregates/Pristine/'+\
                 '%s/TotallyRandom/*'%s)
    print s, len(cscaL)
    for f in fs:
        try:
            fh=Dataset(f,'r')['Freq0336.100GHz_T230.00K']['SingleScatteringData']
            d_eq=Dataset(f,'r')['Freq0336.100GHz_T230.00K']['ShapeData']['diameter_vol_eq'][:]
            theta=fh['za_scat'][:]
            sint=sin(theta/180.*pi)*pi*(2*pi)
            m_eq=pi*4/3.*(d_eq/2)**3*900.*1e3
            m_eqL.append(m_eq)
            phaseF=fh['phaMat_data'][0,0,0,0,:]
            cscaSSD=mean(phaseF*sint)/fh['extMat_data'][0,0,0]
            cscaL.append(cscaSSD)
        except:
            pass

#flaketype,frequencyghz,temperaturek,aeffum,max_dimension_mm,cabs(m^2),cbk(m^2),cext(m^2),csca(m^2),g,ar

massDB=pi*4/3.*(scattdB[afreq[0],2]*1e-6)**3*900.*1e3
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 12})
fig=plt.figure()
ax=plt.gca()
plt.scatter(massDB,scattdB[afreq[0],7]/scattdB[afreq[0],6])
#ax.set_yscale('log')
plt.ylim(1e-10,2e-4)
ax.set_xscale('log')
plt.xlim(0.000001,0.008)
#fig=plt.figure()
#ax=plt.gca()
plt.scatter(m_eqL,cscaL)
#ax.set_yscale('log')
ax.set_xscale('log')
plt.ylim(2e-1,1.)
plt.xlim(0.0000001,0.008)
plt.xlabel('Particle mass (g)')
plt.ylabel('Scattering Albedo')
plt.title('340 GHz Scattering Albedo')
plt.legend(['Liu 2008','SSD Aggregates'])
plt.savefig('scattering340_Aggregates.png')
for i in range(D.shape[0]-1):
    a1=nonzero((massDB-mD[i])*(massDB-mD[i+1])<0)
    if len(a1[0])>0:
        print scattdB[afreq[0][a1],2].mean(),massDB[a1].mean(),  i
for lambd in lambdL:
    M=sum(N0*exp(-lambd*D*1e-2)*mD*dD)
    print M
