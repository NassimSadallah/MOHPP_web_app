'''
Created on Jul 8, 2021

@author: nassim
'''
from math import floor, sqrt
import numpy as np
import matplotlib.pyplot as plt

UNDETECTED_OBSTACLE =-1 
NO_OBSTACLE = 0
KNOWN_OBSTACLE = 1
FORBIDDEN = 9999999999.0 
NEW_FORBIDDEN = 9999999999.0 
INFINI=99999999999.0

def coordinatesToIndex(current_coordinates=[-1, -1], d_ = [-1, -1]):
    
    idx = current_coordinates[0]
    idx += current_coordinates[1] * d_[0]
    return idx

def indexToCoordinates(idx, d_ = [-1, -1]):
    current_coordinates = [-1, -1]
    current_coordinates[1] = int(floor(idx/d_[0]))
    current_coordinates[0] = (idx-current_coordinates[1]*d_[0])      
    return current_coordinates

def getmaxcost(l):
    m = 0
    for i in l:
        if m < i.cost and not i.cost == INFINI:
            m = i.cost
        
    return m

def getmaxVel(l):
    m = 0
    for i in l:
        if m < i.v:# and not i.v == INFINI:
            m = i.v
        
    return m

def getNorth_East_Down(current, nextStep, down, default_altitude):
    
    dNorth = round((-(nextStep[1])+(current[1])),5)
    dEast =round((nextStep[0])-(current[0]),5)
    dDown = round((default_altitude-down) ,5)   
    
    return dNorth, dEast, dDown

def sqrt_dist(a,b,c = 0):
    return sqrt(a**2+b**2+c**2)

def wavePlot(r, r2,nodes):
    zVelocity = [[0 for _ in range(r)] for _ in range(r2)]
    xvec = np.linspace(0.,r,r)
    yvec = np.linspace(0.,r2,r2)                               
    x,y = np.meshgrid(xvec, yvec)
    
    ind = 0
    cou2 = 0
    
    o = getmaxcost(nodes)    
    
    while not cou2==r2-1:
        for c in range(0,r):
            if nodes[ind].cost >o:
                nodes[ind].cost = o+10
            zVelocity[cou2][c] = nodes[ind].cost
            ind+=1 
        cou2+=1                
    
    
    
    plt.contourf(x, y, zVelocity ,200)                             
    plt.colorbar() 

    
    plt.imshow(zVelocity, cmap='Greys',  interpolation='nearest')
    plt.show()    