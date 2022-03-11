import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

data = np.loadtxt('k12')
temp = data[:,0]
m0   = data[:,1]
h0   = data[:,2]
h1   = data[:,3]
l0   = data[:,4]

plt.figure()
plt.ion()

m0par = np.polyfit(temp,m0,2)
h0par = np.polyfit(temp,h0,2)
h1par = np.polyfit(temp,h1,2)
l0par = np.polyfit(temp,l0,2)

m0fit = np.poly1d(m0par)
h0fit = np.poly1d(h0par)
h1fit = np.poly1d(h1par)
l0fit = np.poly1d(l0par)

tempx = np.linspace(80,175,501)

plt.plot(m0fit(tempx),tempx,'k-')
plt.plot(h0fit(tempx),tempx,'k-')
plt.plot(h1fit(tempx),tempx,'k-')
plt.plot(l0fit(tempx),tempx,'k-')

plt.plot(m0,temp,'ko',mfc='none',mec='k',markersize=9,label='Exp.')
plt.plot(h0,temp,'ko',mfc='none',mec='k',markersize=9)
plt.plot(h1,temp,'ko',mfc='none',mec='k',markersize=9)
plt.plot(l0,temp,'ko',mfc='none',mec='k',markersize=9)

plt.ylabel('Temperature ($^\circ$C)',fontsize=12)
plt.xlabel('Weight Per Cent',fontsize=12)
plt.title('Potassium Myristate C14K Water',fontsize=13)
plt.axis([100,0,10,200])

plt.text(25,100,'$L_1$',fontsize=11)
plt.text(50,100,'$E$',fontsize=11)
plt.text(70,100,'$D$',fontsize=11)

plt.legend()
plt.show()

input('press any key to continue ...')
plt.close()
