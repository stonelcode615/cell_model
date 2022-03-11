import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

def func_ln(x,a,b):
    return a*np.log(x)+b

def func(x,a,b,c):
    return a/x+b*x+c

    

# 86
Na_data = np.array([[8,12,16,18],[0.936,1.002,1.041,1.055]])
K_data  = np.array([[8,10,12,14,16,18,20,22],[0.904,0.943,0.973,0.996,1.015,1.031,1.044,1.055]])

carbon_x   = np.linspace(8,22,501)

plt.ion()

plt.figure()
plt.plot(Na_data[0,:],Na_data[1,:],'ks',mfc='none',mec='k',label='Na',markersize=9)
plt.plot(K_data[0,:],  K_data[1,:],'ko',mfc='none',mec='k',label='K', markersize=9)

#popt_na,pcov_na = curve_fit(func_ln,Na_data[0,:],Na_data[1,:])
#popt_k, pcov_k  = curve_fit(func_ln, K_data[0,:], K_data[1,:])
popt_na,pcov_na = curve_fit(func,Na_data[0,:],Na_data[1,:])
popt_k, pcov_k  = curve_fit(func, K_data[0,:], K_data[1,:])

#plt.plot(carbon_x,func_ln(carbon_x,popt_k[0],popt_k[1]+0.026),'k--')
#plt.plot(carbon_x,func_ln(carbon_x,*popt_k), 'k--')

plt.plot(carbon_x,func(carbon_x,popt_k[0],popt_k[1],popt_k[2]+0.026),'k--')
plt.plot(carbon_x,func(carbon_x,*popt_k), 'k--')

print('Na: %10.6f %10.6f %10.6f'%(popt_na[0],popt_na[1],popt_na[2]))
print('K:  %10.6f %10.6f %10.6f'%(popt_k[0], popt_k[1],popt_k[2]))

print('C10Na sv: %5.3f'%func(10,popt_k[0],popt_k[1],popt_k[2]+0.026))
plt.plot(10,func(10,popt_k[0],popt_k[1],popt_k[2]+0.026),'ks')
print('C14Na sv: %5.3f'%func(14,popt_k[0],popt_k[1],popt_k[2]+0.026))
plt.plot(14,func(14,popt_k[0],popt_k[1],popt_k[2]+0.026),'ks')

#Na_para = np.polyfit(Na_data[0,:],Na_data[1,:],2)
#Na_fit  = np.poly1d(Na_para)
#K_para = np.polyfit(K_data[0,:],K_data[1,:],2)
#K_fit  = np.poly1d(K_para)
#plt.plot(carbon_x,Na_fit(carbon_x),'k-')
#plt.plot(carbon_x,K_fit(carbon_x),'k-')
#print(Na_fit(14))

plt.legend()
plt.xlabel('Carbon Number',fontsize=12)
plt.ylabel('Specific Volume (ml/g)',fontsize=12)
plt.title('Specific volume vs Carbon Number',fontsize=13)


plt.show()
input('press any key to continue...')
plt.close()
