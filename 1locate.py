import os
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import sys

DAccuracy=1e-3

logF = {1:'log1', 2:'log2', 3:'log3'}

def intersect(n1,n2):
    #n1 fixed, n2 move
    # or in other words, move points on n2 curve.

    log1 = np.loadtxt(logF[n1])
    log2 = np.loadtxt(logF[n2])

    concentration1 = log1[:,6]
    concentration2 = log2[:,6]

    ChemPoten_A1   = np.asarray(log1[:,2])
    ChemPoten_W1   = np.asarray(log1[:,3])
    ChemPoten_A2   = np.asarray(log2[:,2])
    ChemPoten_W2   = np.asarray(log2[:,3])

    len1 = len(ChemPoten_A1)

    def nearest_idx(x2):
        # x2,y2 is from Hexagon for Hexagon-Lamelle case
        # search point along curve x_in y_in, whose x >= x2 and 
        # have smallest distance
        # return the index of the point in log1
        for idx in range(len1):
            if ChemPoten_W1[idx] >= x2:
                break
        return(idx)

    def distance(idx2,idx1):
        # calculate distance between (x2,y2)   on Hexagon curve 
        # and the line through (x3,y3),(x1,y1) on Lamella curve
        # return the distance
        idx3 = idx1-1
        x2   = ChemPoten_A2[idx2]
        y2   = ChemPoten_W2[idx2]
        x1   = ChemPoten_A1[idx1]
        y1   = ChemPoten_W1[idx1]
        x3   = ChemPoten_A1[idx3]
        y3   = ChemPoten_W1[idx3]
        #print('x3=%8.4f x1=%8.4f'%(x3,x1))
        
        k = (y3-y1)/(x3-x1)
        m = (x3*y1-x1*y3)/(x3-x1)
        d = abs(-k*x2+y2-m)/(1+k**2)**0.5
        return(round(d,5))


    idx2_initial = 1     # starting search from the index 1, the first point on Hexagon curve
    dist_initial = 10    # initial guess of the distance between the first Hexagon point and lamella curve

    Accuracy = 1e-6

    A_dist = []
    L_idxs = []
    Get_idx2 = []
    Get_idx1 = []

    for idx2 in range(idx2_initial,len(ChemPoten_A2)):
        x2 = ChemPoten_W2[idx2]
        y2 = ChemPoten_A2[idx2]
        idx1 = nearest_idx(x2)
        dist = distance(idx2,idx1)
        #print('%2d %8.4f %8.4f %8.4f'%(idx2,dist,x2,y2))
        A_dist = np.append(A_dist,dist)
        L_idxs.append([idx2,idx1])       # the two lines record the distance and index 

        if dist <= dist_initial:
            dist_initial = dist          # swap and save the smallest distance

        #if dist - dist_initial >= 1e-3:  # this the tricky accuracy, to make sure that it stops when the distance begin to increase
        #    break


    sorted_index = np.argsort(A_dist)
#    print(sorted_index)

    if len(sorted_index) >= 1:
        Get_idx = sorted_index[:1]
    else:
        Get_idx = sorted_index

#    print(Get_idx)
    index = Get_idx[0]
    Get_idx2.append(L_idxs[index][0])
    Get_idx1.append(L_idxs[index][1])
    Get_idx1.append(L_idxs[index][1]-1)

    Gotten_delta = []
    Gotten_index = []
        
    for idx2 in Get_idx2:
        for idx1 in Get_idx1:
            Delta_ChemPoten_A = abs(ChemPoten_A1[idx1]-ChemPoten_A2[idx2])
            #print(idx2,idx1,ChemPoten_A2[idx2],ChemPoten_A1[idx1])
            Gotten_delta.append(Delta_ChemPoten_A)
            Gotten_index.append((idx2,idx1))
    idx       = Gotten_delta.index(min(Gotten_delta))

    idx2,idx1 = Gotten_index[idx]
    #print('H %5.2f %8.4f %8.4f %8.4f %8.4f'%(log2[idx2,0],log2[idx2,1],log2[idx2,2],log2[idx2,3],log2[idx2,6]))
    #print('L %5.2f %8.4f %8.4f %8.4f %8.4f'%(log1[idx1,0],log1[idx1,1],log1[idx1,2],log1[idx1,3],log1[idx1,6]))

    w2 = log2[idx2,6]
    w1 = log1[idx1,6]
    #print('%5.1f %5.1f'%(w2,w1))
    return(w2,w1)

hw2, lw = intersect(1,2) 
#mw, hw1 = intersect(2,3)

#print('%5.1f %5.1f %5.1f %5.1f'%(mw,hw1,hw2,lw))
print('%5.1f %5.1f'%(hw2,lw))
#print('%5.1f %5.1f'%(mw,hw1))

