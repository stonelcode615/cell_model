import os
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import sys


bopt_name = {1:'out01',   2:'out02',    3:'out03'}
bfix_name = {1:'out11',   2:'out12',    3:'out13'}
logF_name = {1:'log1',    2:'log2',     3:'log3'}


for num in (1,):
    bopt = bopt_name[num]
    bfix = bfix_name[num]
    logF = logF_name[num]
    
    FOPT = 0
    FFIX = 0
    if os.path.exists(bopt) and os.path.getsize(bopt) > 0:
        FOPT   = 1
    if os.path.exists(bfix) and os.path.getsize(bfix) > 0 :
        FFIX   = 1
    if FOPT*FFIX == 0:
        #print(bopt, ' or ', bfix, ' does not exit or empty\n')
        continue
    else:
        f_bopt = np.loadtxt(bopt,usecols=(0,1,4,5,2,3,6),ndmin=2)
        f_bfix = np.loadtxt(bfix,usecols=(0,1,4,5,2,3,6),ndmin=2)
        bmax   = f_bfix[0,1]

        for index in range(len(f_bopt[:,1])):
            if f_bopt[index,1] < bmax:
                break
#        print(num,index,f_bopt[index,1])             # index is the very first where b is smaller than bmax in opt
                                              
        if index == 0:
            logD = f_bopt[:]
            #print('all b in opt is smaller than bmax ')

        if index == len(f_bopt[:,1])-1:
            logD = f_bfix[:]                         # all b in opt is larger than bmax 
            #print('all b in opt is larger than bmax ')

        if index > 0 and index < len(f_bopt[:,1])-1:
            logD2 = f_bopt[index:]    
            Lopt  = logD2[0,0]                       # the first L in opt file
            #print('the 1st row in opt, b<bmax      L = %5.2f,b = %8.4f'%(Lopt,logD2[0,1]))
            for index2 in range(len(f_bfix[:,0])):
                if f_bfix[index2,0] == Lopt:          
                    break
            if index2 < len(f_bfix[:,0])-1:          
                logD1 = f_bfix[:index2]              
                logD  = np.row_stack((logD1,logD2))  # join opt fix togther
                #print('the row in fix to join with opt L = %5.2f,b = %8.4f'%(f_bfix[index2-1,0],f_bfix[index2-1,1])) 
            if index2 == len(f_bfix[:,0])-1:         # all L in fix is smaller than 1st L in opt
                logD  = logD2[:]
                #print('all L in fix is smaller than 1st L in opt')

        with open(logF,'w') as fsave:
            np.savetxt(fsave,logD,fmt='%-12.8f')




