'''
Created on Jul 8, 2021

@author: nassim
'''
from math import floor, sqrt
import numpy as np
import matplotlib.pyplot as plt
from heapq import heappop
from UavAndSensors.Sensors import Sensors


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

def isEmptyList(l):
    
    for i in l:
        if i !=[]:
            return False    
    return True

def locateBestIdx(l):
    
    b =(INFINI,INFINI)
    ids = 0
    for co,i in enumerate(l):
        if i !=[] and i[0]<b:
            b = i[0]
            ids = co
    #return the popped best element in the hosting list
    return heappop(l[ids])

def DetectUnexpectedObs(curIdx, nodes, extendedObs, safety_margin, sensRange, d_):
    
    isDetected, brake = False, False
    
    sensArea = Sensors()
    sensors = sensArea.sensArea()
    
    for s in sensors:#8 sensors embedded static oriented to the map
        
        if s['ID']=='s1':#the north of the map
            
            if s['VALUE'] <= sensRange:
                obs = int(floor(curIdx - d_[0]*s['VALUE']))
    
                if ((obs>=0) and(int(floor(obs/d_[1]))==int(floor(curIdx/d_[1])))):
                    isDetected = True
                    
                    if s['VALUE'] <=safety_margin:
                        brake = True 
                        
                    curobs = nodes[obs]
                    curobs.OBSTACLE = KNOWN_OBSTACLE
                    curobs.TAG = NEW_FORBIDDEN
                    curobs.cost = INFINI
                    extendedObs.append(curobs)
        
        elif s['ID']=='s2':#north east
            
            if s['VALUE'] <= sensRange:
                obs = int(floor(curIdx - d_[0]*s['VALUE'] + 1))
    
                if ((obs>=0) and(int(floor(obs/d_[1]))==int(floor(curIdx/d_[1])))):
                    isDetected = True
                    
                    if s['VALUE'] <=safety_margin:
                        brake = True 
                        
                    curobs = nodes[obs]
                    curobs.OBSTACLE = KNOWN_OBSTACLE
                    curobs.TAG = NEW_FORBIDDEN
                    curobs.cost = INFINI
                    extendedObs.append(curobs)            
        
        elif s['ID']=='s3': #east

            if s['VALUE'] <= sensRange:
                obs = int(floor(curIdx + s['VALUE']))
    
                if ((obs>=0) and(int(floor(obs/d_[1]))==int(floor(curIdx/d_[1])))):
                    isDetected = True
                    
                    if s['VALUE'] <=safety_margin:
                        brake = True 
                        
                    curobs = nodes[obs]
                    curobs.OBSTACLE = KNOWN_OBSTACLE
                    curobs.TAG = NEW_FORBIDDEN
                    curobs.cost = INFINI
                    extendedObs.append(curobs)            
            
        elif s['ID']=='s4': #south east

            if s['VALUE'] <= sensRange:
                obs = int(floor(curIdx + d_[0]*s['VALUE']+s['VALUE']+s['VALUE']))
    
                if ((obs>=0) and(int(floor(obs/d_[1]))==int(floor(curIdx/d_[1])))):
                    isDetected = True
                    
                    if s['VALUE'] <=safety_margin:
                        brake = True 
                        
                    curobs = nodes[obs]
                    curobs.OBSTACLE = KNOWN_OBSTACLE
                    curobs.TAG = NEW_FORBIDDEN
                    curobs.cost = INFINI
                    extendedObs.append(curobs)            
        
        elif s['ID']=='s5':#south
            
            if s['VALUE'] <= sensRange:
                obs = int(floor(curIdx + d_[0]*s['VALUE']))
    
                if ((obs>=0) and(int(floor(obs/d_[1]))==int(floor(curIdx/d_[1])))):
                    isDetected = True
                    
                    if s['VALUE'] <=safety_margin:
                        brake = True 
                        
                    curobs = nodes[obs]
                    curobs.OBSTACLE = KNOWN_OBSTACLE
                    curobs.TAG = NEW_FORBIDDEN
                    curobs.cost = INFINI
                    extendedObs.append(curobs)            
        
        elif s['ID']=='s6':#south west
            
            if s['VALUE'] <= sensRange:
                obs = int(floor(curIdx + d_[0]*s['VALUE'] - s['VALUE']))
    
                if ((obs>=0) and(int(floor(obs/d_[1]))==int(floor(curIdx/d_[1])))):
                    isDetected = True
                    
                    if s['VALUE'] <=safety_margin:
                        brake = True 
                        
                    curobs = nodes[obs]
                    curobs.OBSTACLE = KNOWN_OBSTACLE
                    curobs.TAG = NEW_FORBIDDEN
                    curobs.cost = INFINI
                    extendedObs.append(curobs)            
        
        elif s['ID']=='s7':#west
            
            if s['VALUE'] <= sensRange:
                obs = int(floor(curIdx - s['VALUE']))
    
                if ((obs>=0) and(int(floor(obs/d_[1]))==int(floor(curIdx/d_[1])))):
                    isDetected = True
                    
                    if s['VALUE'] <=safety_margin:
                        brake = True 
                        
                    curobs = nodes[obs]
                    curobs.OBSTACLE = KNOWN_OBSTACLE
                    curobs.TAG = NEW_FORBIDDEN
                    curobs.cost = INFINI
                    extendedObs.append(curobs)            
            
        elif s['ID']=='s8':#north west
            
            if s['VALUE'] <= sensRange:
                obs = int(floor(curIdx - d_[0]*s['VALUE']- s['VALUE']))
    
                if ((obs>=0) and(int(floor(obs/d_[1]))==int(floor(curIdx/d_[1])))):
                    isDetected = True
                    
                    if s['VALUE'] <=safety_margin:
                        brake = True 
                        
                    curobs = nodes[obs]
                    curobs.OBSTACLE = KNOWN_OBSTACLE
                    curobs.TAG = NEW_FORBIDDEN
                    curobs.cost = INFINI
                    extendedObs.append(curobs)            
            
    return extendedObs, isDetected, brake 
    



