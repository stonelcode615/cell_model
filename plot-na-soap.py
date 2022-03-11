import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

#dat10 = np.loadtxt('cal-na12')
#dat10 = np.loadtxt('cal-na14')
#dat10 = np.loadtxt('cal-na16')
dat10 = np.loadtxt('cal-C18Na')
temp10 = dat10[:,0]
m010   = dat10[:,1]
h010   = dat10[:,2]
h110   = dat10[:,3]
l010   = dat10[:,4]

#data = np.loadtxt('na12')
#data = np.loadtxt('na14')
#data = np.loadtxt('na16')
data = np.loadtxt('C18Na')
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

m0par = np.polyfit(temp,m0,2)
h0par = np.polyfit(temp,h0,2)
h1par = np.polyfit(temp,h1,2)
l0par = np.polyfit(temp,l0,2)

m0fit = np.poly1d(m0par)
h0fit = np.poly1d(h0par)
h1fit = np.poly1d(h1par)
l0fit = np.poly1d(l0par)

tempx = np.linspace(86,137,501)

plt.plot(m0fit(tempx),tempx,'k-')
plt.plot(h0fit(tempx),tempx,'k-')
plt.plot(h1fit(tempx),tempx,'k-')
plt.plot(l0fit(tempx),tempx,'k-')

plt.ylabel('Temperature ($^\circ$C)',fontsize=12)
plt.xlabel('Weight Per Cent',fontsize=12)
plt.title('Sodium Myristate C14Na Water',fontsize=13)
plt.axis([100,0,60,190])

plt.text(25,100,'$L_1$',fontsize=11)
plt.text(42,100,'$E$',fontsize=11)
plt.text(70,100,'$D$',fontsize=11)

plt.legend()
plt.show()

input('press any key to continue ...')
plt.close()
