import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

dat10 = np.loadtxt('cal-C12K')
#dat10 = np.loadtxt('cal-k14')
#dat10 = np.loadtxt('cal-k16')
#dat10 = np.loadtxt('cal-k18')
temp10 = dat10[:,0]
m010   = dat10[:,1]
h010   = dat10[:,2]
h110   = dat10[:,3]
l010   = dat10[:,4]

data = np.loadtxt('C12K')
#data = np.loadtxt('k14')
#data = np.loadtxt('k16')
#data = np.loadtxt('k18')
temp = data[:,0]
m0   = data[:,1]
h0   = data[:,2]
h1   = data[:,3]
l0   = data[:,4]

plt.figure()
plt.ion()

plt.plot(m010,temp10,'ks',markersize=4,label='Cal.')
plt.plot(h010,temp10,'ks',markersize=4)
plt.plot(h110,temp10,'ks',markersize=4)
plt.plot(l010,temp10,'ks',markersize=4)

plt.plot(m0,temp,'ko',mfc='none',mec='k',markersize=9,label='Exp.')
plt.plot(h0,temp,'ko',mfc='none',mec='k',markersize=9)
plt.plot(h1,temp,'ko',mfc='none',mec='k',markersize=9)
plt.plot(l0,temp,'ko',mfc='none',mec='k',markersize=9)

#plt.plot(27.9725,45.1195,'k^',markersize=16,mfc='none',mec='k')
#plt.plot(55.8763,45.1195,'k^',markersize=16,mfc='none',mec='k')
#plt.plot(33.1959,104.096,'k^',markersize=16,mfc='none',mec='k')
#plt.plot(57.2509,104.096,'k^',markersize=16,mfc='none',mec='k')
#plt.plot(56.8504,85.8824,'k^',markersize=16,mfc='none',mec='k')
#plt.plot(31.3386,85.8824,'k^',markersize=16,mfc='none',mec='k')


m0par = np.polyfit(temp,m0,2)
h0par = np.polyfit(temp,h0,2)
h1par = np.polyfit(temp,h1,2)
l0par = np.polyfit(temp,l0,2)

m0fit = np.poly1d(m0par)
h0fit = np.poly1d(h0par)
h1fit = np.poly1d(h1par)
l0fit = np.poly1d(l0par)

tempx = np.linspace(30,180,501)

plt.plot(m0fit(tempx),tempx,'k-')
plt.plot(h0fit(tempx),tempx,'k-')
plt.plot(h1fit(tempx),tempx,'k-')
plt.plot(l0fit(tempx),tempx,'k-')

plt.ylabel('Temperature ($^\circ$C)',fontsize=12)
plt.xlabel('Weight Per Cent',fontsize=12)
#plt.title('Potassium Myristate C14K Water',fontsize=13)
plt.axis([100,0,30,190])

plt.text(25,100,'$L_1$',fontsize=11)
plt.text(50,100,'$E$',fontsize=11)
plt.text(70,100,'$D$',fontsize=11)

plt.legend()
plt.show()

input('press any key to continue ...')
plt.close()
