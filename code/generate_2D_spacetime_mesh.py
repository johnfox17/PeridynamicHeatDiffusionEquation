######################################################################################################################################
## This will be a python code that solves example problem 6.1 from Professor Madenci's book "Peridynamic DO for Numerical Analysis".##
## pg 118 - 120. I will be using a fortran code given to me by the Professor as reference.                                          ##
##                                                                                                                                  ##
## OBJECTIVES:                                                                                                                      ##
## 1. Learn to use and solve problems by means of Peridynamic Differential Operators                                                ##
## 2. Solve diffusion equation by using PDDO for future use in applications such as image noise removal.                            ##
## 3. Finish my PhD                                                                                                                 ##
######################################################################################################################################

#TODO
# 1. Generate space-time mesh.
# 2. Generate node family.
# 3. Generate g polynomials 
# 4.
# 5. 

import os
import numpy as np
import math as mt
from scipy import spatial, integrate
import matplotlib.pyplot as plt

TOL = 1.000001

num_node_x = 41
num_node_y = 41
dxy=0.025

if num_node_x % 2 == 0 & num_node_y % 2 == 0:
    x = np.arange(dxy/2, 1, dxy)
    y = np.arange(dxy/2, 1, dxy)
else:
    x = np.arange(0, 1+dxy, dxy)
    y = np.arange(0, 1+dxy, dxy)


indexing = 'xy'
asymmetric_time = False

xx, yy = np.meshgrid(x, y, indexing=indexing)
xx = xx.reshape(-1, 1)
yy = yy.reshape(-1, 1)

CRD = np.concatenate((xx, yy), axis=-1)

horiz = 3
num_feature = int(horiz*2 + 1)**2

pts, ids, exi = [], [], []
for I, crdI in enumerate(CRD):
    if(I<num_node_x*num_node_y):
        IX, IY = np.unravel_index(I, (num_node_x, num_node_y))
        xref = np.linspace(-horiz, horiz, horiz*2+1, dtype=int)
        if np.any(xref+IX<0):
            xref -= np.min(xref+IX)
        elif np.any(xref+IX>=num_node_x):
            xref -= np.max(xref+IX)-num_node_x+1
        yref = np.linspace(-horiz, horiz, horiz*2+1, dtype=int)
        if np.any(yref+IY<0):
            yref -= np.min(yref+IY)
        elif np.any(yref+IY>=num_node_y):
            yref -= np.max(yref+IY)-num_node_y+1
       # tref = np.where(tref>0,-999,tref)
        xref, yref = np.meshgrid(xref, yref, indexing=indexing)
        exi_ref = np.concatenate((xref.reshape((-1,1)), yref.reshape((-1,1))), axis=-1)
        distI = np.sqrt(np.sum(exi_ref**2, axis=-1))
        orderI = np.argsort(distI, axis=0)
        exi_ref = exi_ref[orderI, :]
        exiI = exi_ref
        ptsI = np.array([[IX, IY]]) + exi_ref
        ptsI = np.ravel_multi_index(ptsI.T, (num_node_x, num_node_y))
        
        #The bottom line assures assymetric family for time coordinate
        if asymmetric_time:
            mask = np.where(ptsI>=(ptsI[0]+num_node_x)//num_node_x*num_node_x)
            ptsI = np.delete(ptsI, mask)
            #print(ptsI)
            #a = input('').split(" ")[0]
        
        pts.append(ptsI)
        exi.append(exiI)
        ids.append(I)

area = np.ones_like(xx) * (x.max()-x.min())/(num_node_x-1) * (y.max()-y.min())/(num_node_y-1)





path = '../data';
#if os.path.isfile(path) == False:
#    os.mkdir(path)
if not os.path.exists(path):
    os.mkdir(path)

np.savetxt('../data/uniform_horizon2d_std.txt', np.array(pts)[ids, :], fmt='%6d', delimiter=', ')
np.savetxt('../data/uniform_xi2d_x_std.txt', np.array(exi)[ids, :, 0], fmt='%.6f', delimiter=',')
np.savetxt('../data/uniform_xi2d_y_std.txt', np.array(exi)[ids, :, 1], fmt='%.6f', delimiter=',')
np.savetxt('../data/uniform_mesh2d_std.txt', CRD, fmt='%.6f', delimiter=',')
np.savetxt('../data/uniform_area2d_std.txt', area.reshape(-1, 1), fmt='%.6f', delimiter=',')





