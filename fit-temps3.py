import os
import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

fit_mode = 1 # 1, parabola

def linear(x,a,b):
    return a*x+b

def parabola(x,a,b,c):
    return a*x**2+b*x+c

def func(x,a,b,c):
    return a/x+b*x+c

def func2(x,a,b,c,d):
    return a/x+b*x+c+d*np.log(x)

def main():
    #f_data  = np.loadtxt('./na12/temps3.dat')
    #f_data  = np.loadtxt('./na14/temps3.dat')
    f_data  = np.loadtxt('./na16/temps3.dat')
    #f_data  = np.loadtxt('./na18/temps3.dat')
    temp    = f_data[:,0]+273
    mu20    = f_data[:,1]
    mu30    = f_data[:,2]
    muw20   = f_data[:,3]*1000
    muw30   = f_data[:,4]*1000

    print(temp)

    plt.ion()

    temps   = np.linspace(temp[0]-5,temp[2]+5,501)

    if fit_mode == 1:
        popt22,pcov22 = curve_fit(func,temp[:3],mu20[:3])
        popt32,pcov32 = curve_fit(func,temp[:3],mu30[:3])
        popt42,pcov42 = curve_fit(func,temp[:3],muw20[:3])
        popt52,pcov52 = curve_fit(func,temp[:3],muw30[:3])
        print('mu20  a=%-10.6e; b=%-10.6f; c=%-10.6f;'%(popt22[0], popt22[1], popt22[2]))
        print('mu30  a=%-10.6e; b=%-10.6f; c=%-10.6f;'%(popt32[0], popt32[1], popt32[2]))
        print('muw20 a=%-10.8e; b=%-10.8f; c=%-10.8f;'%(popt42[0], popt42[1], popt42[2]))
        print('muw30 a=%-10.8e; b=%-10.8f; c=%-10.8f;'%(popt52[0], popt52[1], popt52[2]))
        #plt.plot(temps,func(temps,*popt22),'k-',label='quadratic')
        #plt.plot(temps,func(temps,*popt32),'k-')

    plt.plot(temp,mu20,'ko',mfc='none',mec='k',label='$\delta\mu^{0,2}_{mon}$',markersize=9)
    plt.plot(temp,mu30,'ks',mfc='none',mec='k',label='$\delta\mu^{0,2}_{mon}$',markersize=9)
    plt.xlabel('Temperature ($^\circ$C)',fontsize=12)
    plt.ylabel('$\Delta\mu^0_{correct}$ (kT)',fontsize=12)
    plt.title('$\Delta\mu^0_{correct}$ vs Temperature',fontsize=15)



    plt.figure()
    plt.plot(temp,muw20,'k+',label='$\delta\mu^{0,2}_{wat}$',markersize=9)
    plt.plot(temp,muw30,'kx',label='$\delta\mu^{0,3}_{wat}$',markersize=9)
    plt.plot(temps,func(temps,*popt42),'k-')
    plt.plot(temps,func(temps,*popt52),'k-')
    
    #a=4.72376733e+01; b=0.00024658; c=-0.21710243
    tps=[359., 393., 417., 363., 368., 373., 378., 383., 388., 398., 403., 408., 413.]
    for tp in tps:
        plt.plot(tp,func(tp,*popt42),'ko',mfc='none',mec='k')
#        plt.plot(tp,func(tp,a,b,c),'ko',mfc='none',mec='k')
    
    plt.legend()
    plt.show()

    input('press any key to continue ...')
    plt.close()

if __name__ == '__main__':
    main()

