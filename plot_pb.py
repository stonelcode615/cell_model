import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

p1 = np.loadtxt('./log1')
#p2 = np.loadtxt('./log2')
#p3 = np.loadtxt('./log3')


plt.ion()
plt.figure()

plt.plot(p1[:,3],p1[:,2],'k-' ,mfc='none',mec='k',label='lamellar')
#plt.plot(p2[:,3],p2[:,2],'b--',mfc='none',mec='b',label='hexagonal')
#plt.plot(p3[:,3],p3[:,2],'g:' ,mfc='none',mec='g',label='micelle')

#plt.axis((-10.0,-0.0,-9.5,-7.5))
plt.legend()
plt.show()
input('enter any key to close...')
plt.close()
