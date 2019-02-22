from numpy import *
dbins=[50.0,   70.0,    90.0,   112.5,   137.5,   175.0,   225.0,   275.0,   325.0,   375.0,   437.5,   512.5,   587.5,   662.5,   750.0,   850.0,   950.0,  1100.0,  1300.0,  1500.0,  1700.0,  2000.0,  2400.0,  2800.0,  3200.0,  3600.0,  4000.0,  4400.0,  4800.0,  5500.0,  6500.0,  7500.0,  8500.0,  9500.0, 11000.0, 13000.0, 15000.0, 17000.0, 19000.0, 22500.0,  27500.0]
dbins=array(dbins)
dint=0.5*(dbins[0:-1]+dbins[1:])
dintL=[]
dintL.append(2*dbins[0]-dint[0])
for d in dint:
    dintL.append(d)
dintL.append(2*dbins[-1]-dint[-1])
a=zeros((42,42),float)
b=zeros((42),float)
a[0,0]=1
b[0]=40
for i in range(41):
    a[i+1,i]=0.5
    a[i+1,i+1]=0.5
    b[i+1]=dbins[i]
dint=dot(linalg.pinv(a),b)
ds=dint[1:]-dint[0:-1]

mass=0.0061*(1e-4*dbins)**2.05
massInt=0.0061*(1e-4*dint)**2.05
import pickle
luo_Tables=pickle.load(open('liu_snow_Tables_Extended.pklz','rb'))
luo_Tables2=pickle.load(open('liu_snow_Tables_Extended.pklz','rb'))


dbinsM=(mass*6/pi/0.92)**(1./3.)*10.
dbinsI=(mass*6/pi)**(1./3.)*10.

import pytmatrix.refractive
wl=[pytmatrix.refractive.wl_Ku,pytmatrix.refractive.wl_Ka,pytmatrix.refractive.wl_W]

sback95=0.5*(10**luo_Tables[0][0]*wl[2]**4/pi**5+10**luo_Tables2[0][0]*wl[2]**4/pi**5)
sback37=0.5*(10**luo_Tables[1][0]*wl[1]**4/pi**5+10**luo_Tables2[1][0]*wl[1]**4/pi**5)
sback13=0.5*(10**luo_Tables[2][0]*wl[0]**4/pi**5+10**luo_Tables2[2][0]*wl[0]**4/pi**5)
ext95=0.5*(10**luo_Tables[0][1]+10**luo_Tables2[0][1])
ext37=0.5*(10**luo_Tables[1][1]+10**luo_Tables2[1][1])
ext13=0.5*(10**luo_Tables[2][1]+10**luo_Tables2[2][1])

sca95=0.5*(10**luo_Tables[0][2]+10**luo_Tables2[0][2])
sca37=0.5*(10**luo_Tables[1][2]+10**luo_Tables2[1][2])
sca13=0.5*(10**luo_Tables[2][2]+10**luo_Tables2[2][2])

gsca95=0.5*(luo_Tables[0][3]+luo_Tables2[0][3])
gsca37=0.5*(luo_Tables[1][3]+luo_Tables2[1][3])
gsca13=0.5*(luo_Tables[2][3]+luo_Tables2[2][3])

#luo_Tables=pickle.load(open('luo_Graupel_Tables.pklz','rb'))
#luo_Tables2=pickle.load(open('luo_Graupel_Tables.pklz','rb'))

mi1= pytmatrix.refractive.mi(wl[0],0.92)
mi2= pytmatrix.refractive.mi(wl[1],0.92)
mi3= pytmatrix.refractive.mi(wl[2],0.92)
m=array([mi1,mi2,mi3])
K2r=abs((m**2-1)/(m**2+2))**2
K2s=((m**2-1)/(m**2+2))**2
#print K2r
K2=0.93

dbins1=0.5*dbinsI*1e3+0.5*dbins
xW1=4*3.142/wl[0]*(0.30*dbins1*1e-3)
xf1=(1+0.159*xW1**2)/(1+(0.159+1./3)*xW1**2+0.164*xW1**4);
xW2=4*3.142/wl[1]*(0.30*dbins1*1e-3)
xf2=(1+0.159*xW2**2)/(1+(0.159+1./3)*xW2**2+0.164*xW2**4);
xW3=4*3.142/wl[2]*(0.30*dbins1*1e-3)

xf3=(1+0.159*xW3**2)/(1+(0.159+1./3)*xW3**2+0.164*xW3**4);
#print sca37/ext37
#print gsca37
#stop
def readAH(fname):
    lines=open(fname,'r').readlines()
    datL=[]
    datL2=[]
    for l in lines[:]:
        #print l
        try:
            dat1=[float(v) for v in l.split()]
        except:
            dat1=[]
            pass
        if len(dat1)==44:
            if dat1[2]>0.0001 and dat1[2]<9.999e3:
                Nbins=array(dat1[3:])
                dm=sum(Nbins*ds*mass*dbinsM)/sum(Nbins*ds*mass)
                Z1=log10(sum(Nbins*ds*1e-6*(dbinsI)**6*xf1)*K2r[0]/K2)*10.
                Z2=log10(sum(Nbins*ds*1e-6*(dbinsI)**6*xf2)*K2r[1]/K2)*10.
                Z3=log10(sum(Nbins*ds*1e-6*(dbinsI)**6*xf3)*K2r[2]/K2)*10.
                Z3=log10(sum(Nbins*ds*1e-6*sback95/K2))*10.
                Z1=log10(sum(Nbins*ds*1e-6*sback13/K2))*10.
                Z2=log10(sum(Nbins*ds*1e-6*sback37/K2))*10.
                kext35=4.34*sum(Nbins*ds*1e-6*ext37)*1e-3
                kext95=4.34*sum(Nbins*ds*1e-6*ext95)*1e-3
                kext13=4.34*sum(Nbins*ds*1e-6*ext13)*1e-3
                ksca13=4.34*sum(Nbins*ds*1e-6*sca13)*1e-3
                ksca35=4.34*sum(Nbins*ds*1e-6*sca37)*1e-3
                ksca95=4.34*sum(Nbins*ds*1e-6*sca95)*1e-3
                _gsca13=4.34*sum(Nbins*ds*1e-6*sca13*gsca13)*1e-3
                _gsca35=4.34*sum(Nbins*ds*1e-6*sca37*gsca37)*1e-3
                _gsca95=4.34*sum(Nbins*ds*1e-6*sca95*gsca95)*1e-3
                #if Z1>250:
                #    print Z3, kext95, kext35
                Nw=4**4/pi*dat1[2]/(0.1*dm)**4/1e6
                #if abs(dm-1)<1e-3:
                    #print sum(Nbins*ds*mass*dbinsM)
                    #print sum(Nbins*ds*mass), dat1[2]
                dat2=[kext13,ksca13/kext13,ksca35/kext35,ksca95/kext95,\
                      _gsca13/ksca13,_gsca35/ksca35,_gsca95/ksca95]
                #print dat2
            else:
                Z1=-99
                Z2=-99
                Z3=-99
                kext35=0
                kext95=0
                kext13=0
                dm=0
                Nw=-99
                dat2=[kext13,0.,0.,0.,0.,0.,0.]
                #print [Nw,kext35,kext95,Z1,Z2,Z3,dm]
            dat1.extend([Nw,kext35,kext95,Z1,Z2,Z3,dm])
            
            datL.append(dat1)
            datL2.append(dat2)
            #print dat1
    #stop
    return array(datL),array(datL2)

def readAH2(fname):
    lines=open(fname,'r').readlines()
    datL=[]
    nbinsL=[]
    for l in lines[:]:
        #print l
        try:
            dat1=[float(v) for v in l.split()]
        except:
            dat1=[]
            pass
        if len(dat1)==44:
            if dat1[2]>0.0001 and dat1[2]<9.999e3:
                Nbins=array(dat1[3:])
                dm=sum(Nbins*ds*mass*dbinsM)/sum(Nbins*ds*mass)
                Z1=log10(sum(Nbins*ds*1e-6*(dbinsI)**6*xf1)*K2r[0]/K2)*10.
                Z2=log10(sum(Nbins*ds*1e-6*(dbinsI)**6*xf2)*K2r[1]/K2)*10.
                Z3=log10(sum(Nbins*ds*1e-6*(dbinsI)**6*xf3)*K2r[2]/K2)*10.
                Z3=log10(sum(Nbins*ds*1e-6*sback95/K2))*10.
                Z1=log10(sum(Nbins*ds*1e-6*sback13/K2))*10.
                Z2=log10(sum(Nbins*ds*1e-6*sback37/K2))*10.
                kext35=4.34*sum(Nbins*ds*1e-6*ext37)*1e-3
                kext95=4.34*sum(Nbins*ds*1e-6*ext95)*1e-3
                if Z1>250:
                    print Z3, kext95, kext35
                Nw=4**4/pi*dat1[2]/(0.1*dm)**4/1e6
                #if abs(dm-1)<1e-3:
                    #print sum(Nbins*ds*mass*dbinsM)
                    #print sum(Nbins*ds*mass), dat1[2]
            else:
                Z1=-99
                Z2=-99
                Z3=-99
                kext35=0
                kext95=0
                dm=0
                Nw=-99
                Nbins=array(dat1[3:])*0-99
                #print [Nw,kext35,kext95,Z1,Z2,Z3,dm]
            dat1.extend([Nw,kext35,kext95,Z1,Z2,Z3,dm])

            datL.append(dat1)
            nbinsL.append(Nbins)
            #print dat1
    return array(datL), array(nbinsL), sback95, sback13, sback37, K2, ds

#datAH=readAH()
