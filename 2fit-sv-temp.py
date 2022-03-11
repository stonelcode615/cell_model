import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

def func(x,a,b,c):
    return a/x+b*x+c
    
k10_sv  = np.array([[45,65,86],[0.898,0.920,0.943]])
k14_sv  = np.array([[45,65,86,104],[0.957,0.976,0.996,1.014]])
k18_sv  = np.array([[86,104],[1.031,1.047]])

tem_x   = np.linspace(30+273,180+273,1001)

plt.ion()

print(k10_sv[0,:]+273)
print(k10_sv[1,:])

plt.figure()
plt.plot(k10_sv[0,:]+273,k10_sv[1,:],'ks',mfc='none',mec='k',label='k10',markersize=9)
plt.plot(k14_sv[0,:]+273,k14_sv[1,:],'ko',mfc='none',mec='k',label='k14',markersize=9)
plt.plot(k18_sv[0,:]+273,k18_sv[1,:],'k^',mfc='none',mec='k',label='k18',markersize=9)


popt_14,pcov_14 = curve_fit(func,k14_sv[0,:]+273,k14_sv[1,:])
print('k14 Temp: %8.4f %8.4f %8.4f'%(popt_14[0],popt_14[1],popt_14[2]))

k14_sv_86=0.996;  # reference
k18_sv_86=1.031;  diff_k18=k18_sv_86-k14_sv_86
k16_sv_86=1.015;  diff_k16=k16_sv_86-k14_sv_86
k12_sv_86=0.973;  diff_k12=k12_sv_86-k14_sv_86

na18_sv_86=1.055; diff_na18=na18_sv_86-k14_sv_86
na16_sv_86=1.041; diff_na16=na16_sv_86-k14_sv_86
na14_sv_86=1.022; diff_na14=na14_sv_86-k14_sv_86
na12_sv_86=1.002; diff_na12=na12_sv_86-k14_sv_86

print('k18')
print('a=%-10.6e; b=%-10.6f; c=%-10.6f;'%(popt_14[0], popt_14[1], popt_14[2]+diff_k18))
print('k16')
print('a=%-10.6e; b=%-10.6f; c=%-10.6f;'%(popt_14[0], popt_14[1], popt_14[2]+diff_k16))
print('k12')
print('a=%-10.6e; b=%-10.6f; c=%-10.6f;'%(popt_14[0], popt_14[1], popt_14[2]+diff_k12))
print('na18')
print('a=%-10.6e; b=%-10.6f; c=%-10.6f;'%(popt_14[0], popt_14[1], popt_14[2]+diff_na18))
print('na16')
print('a=%-10.6e; b=%-10.6f; c=%-10.6f;'%(popt_14[0], popt_14[1], popt_14[2]+diff_na16))
print('na14')
print('a=%-10.6e; b=%-10.6f; c=%-10.6f;'%(popt_14[0], popt_14[1], popt_14[2]+diff_na14))
print('na12')
print('a=%-10.6e; b=%-10.6f; c=%-10.6f;'%(popt_14[0], popt_14[1], popt_14[2]+diff_na12))

plt.plot(tem_x,func(tem_x,popt_14[0],popt_14[1],popt_14[2]-0.056),'k--')
plt.plot(tem_x,func(tem_x,*popt_14),'k--')
plt.plot(tem_x,func(tem_x,popt_14[0],popt_14[1],popt_14[2]+0.035),'k--')

plt.legend()
plt.xlabel('Temperature $^\circ$C',fontsize=12)
plt.ylabel('Specific Volume (ml/g)',fontsize=12)
plt.title('Specific volume vs Temperature',fontsize=13)


plt.show()
input('press any key to continue...')
plt.close()
