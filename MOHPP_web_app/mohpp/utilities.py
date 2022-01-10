'''
Created on Jul 8, 2021

@author: nassim
'''
from math import floor, cos, sin, pi, sqrt, log10
import numpy as np
import matplotlib.pyplot as plt
from heapq import heappop
from Tkinter import *
import time
from autobahn.twisted.util import sleep
risks = 0
UNDETECTED_OBSTACLE =-1 
NO_OBSTACLE = 0
KNOWN_OBSTACLE = 1
FORBIDDEN = -1 
NEW_FORBIDDEN = -2 
INFINI=9999999.0
width, height = 200, 150
echelle = 2
FAR, KNOWN, FROZEN = -1,0,1
globalPath = []
d_ = [width, width*height]

def coordinatesToIndex(current_coordinates, d_ = [-1, -1]):
    
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
        if m < i.cost and i.cost < INFINI:
            m = i.cost  
    return m

def getRisk(l):
    global risks    
    for i in l:
        if i.cost != INFINI and risks< i.risk:      
            risks = i.risk
    return risks 

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

def wavePlot(r= -1, r2= -1,nodes = [],path = [], start = [], goal = []):
    zVelocity = [[0 for _ in range(r)] for _ in range(r2)]
    xvec = np.linspace(0.,r,r)
    yvec = np.linspace(0.,r2,r2)                               
    x,y = np.meshgrid(xvec, yvec)
    
    ind = 0
    cou2 = r2-1
    
    o = getmaxcost(nodes)    
    
    while not cou2==0:
        for c in range(0,r):
            if nodes[ind].cost >o:
                nodes[ind].cost = o+10
            zVelocity[cou2][c] = nodes[ind].cost
            ind+=1 
        cou2-=1                
    
    plt.imshow(zVelocity, cmap='Greys',  interpolation='nearest')
    
    for i in path:
       
        i[1] = r2-i[1]    
        plt.plot(i[0],i[1],'m+')
    
    plt.plot(start.abscice,r2-start.colonne,'^k:')
    plt.plot(goal.abscice,r2-goal.colonne,'go')
    plt.contourf(x, y, zVelocity ,200)                             
    plt.colorbar() 
    
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

def UavHeading(heading):
    return heading*pi/180

def degToCartesian(deg):
    return deg*pi/180

def DetectUnexpectedObs(sensType, sensors, UavOrient, curIdx, nodes, extendedObs, safety_margin, sensRange, d_):
    
    isDetected, brake = False, False
    
    cur = [nodes[curIdx].abscice,nodes[curIdx].colonne] # save the x, y coordinates of the current position of the UAV    
    
    if sensType == 'hcsr04':
        sensorsValues=sensors#call the method which read the sensors' values sent from the arduino nano
        sN, sE, sS, sW, sNE, sSE, sSW, sNW = [0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]

        if sensorsValues !=[]: 
            
            if sensorsValues[0] >=.8 and sensorsValues[0] < sensRange:
    
                #transforms the radian values to cartesian one
                sN= [int(round(int(sensorsValues[0])*cos(2*pi-UavOrient+pi/2),0)),-int(round(int(sensorsValues[0])*sin(2*pi-UavOrient+pi/2),0))]   
                obs = [cur[0]+sN[0],cur[1]+sN[1]]
                    
                if ((obs>=[0,0]) and(obs<=[width, height])):
                      
                    curobs = nodes[coordinatesToIndex(obs, d_)]
                
                    if curobs.OBSTACLE != KNOWN_OBSTACLE and curobs.TAG != FORBIDDEN and sN !=[0,0]:            
                        
                        isDetected = True
                        curobs.OBSTACLE = KNOWN_OBSTACLE
                        curobs.TAG = NEW_FORBIDDEN
                        curobs.cost = INFINI
                        extendedObs.append(curobs)
                                    
                        if sensorsValues[0]<=safety_margin:
                            
                            brake = True
            
            if sensorsValues[1] >=.8 and sensorsValues[1] < safety_margin:
                
                sE= [int(round(int(sensorsValues[1])*cos(2*pi-UavOrient),0)),int(round(int(sensorsValues[1])*sin(2*pi-UavOrient),0))]
                
                obs = [cur[0]+sE[0],cur[1]+sE[1]]
                    
                if ((obs>=[0,0]) and(obs<=[width, height])):
                    curobs = nodes[coordinatesToIndex(obs, d_)]
                    
                    if curobs.OBSTACLE != KNOWN_OBSTACLE and curobs.TAG != FORBIDDEN and sE !=[0,0]:            
                        
                        isDetected = True
                        curobs.OBSTACLE = KNOWN_OBSTACLE
                        curobs.TAG = NEW_FORBIDDEN
                        curobs.cost = INFINI
                        extendedObs.append(curobs)
                                    
                        if sensorsValues[1]<=safety_margin:
                            
                            brake = True
            
            if sensorsValues[2] >=.8 and sensorsValues[2] < sensRange:
                
                sS= [int(round(int(sensorsValues[2])*cos(2*pi-UavOrient+pi*3/2),0)),-int(round(int(sensorsValues[2])*sin(2*pi-UavOrient+pi*3/2),0))]          
                obs = [cur[0]+sS[0],cur[1]+sS[1]]
                    
                if ((obs>=[0,0]) and(obs<=[width, height])):
                
                    curobs = nodes[coordinatesToIndex(obs, d_)]
                    
                    if curobs.OBSTACLE != KNOWN_OBSTACLE and curobs.TAG != FORBIDDEN and sS !=[0,0]:            
                        
                        isDetected = True
                        curobs.OBSTACLE = KNOWN_OBSTACLE
                        curobs.TAG = NEW_FORBIDDEN
                        curobs.cost = INFINI
                        extendedObs.append(curobs)
                                    
                        if sensorsValues[2]<=safety_margin:
                      
                            brake = True
            
            if sensorsValues[3] >=.8 and sensorsValues[3] < sensRange:
                
                sW= [int(round(int(sensorsValues[3])*cos(2*pi-UavOrient+pi),0)),int(round(int(sensorsValues[3])*sin(2*pi-UavOrient+pi),0))]            
                obs = [cur[0]+sW[0],cur[1]+sW[1]]
                    
                if ((obs>=[0,0]) and(obs<=[width, height])):
                    
                    curobs = nodes[coordinatesToIndex(obs, d_)]
                    
                    if curobs.OBSTACLE != KNOWN_OBSTACLE and curobs.TAG != FORBIDDEN and sW !=[0,0]:            
                        
                        isDetected = True
                        curobs.OBSTACLE = KNOWN_OBSTACLE
                        curobs.TAG = NEW_FORBIDDEN
                        curobs.cost = INFINI
                        extendedObs.append(curobs)
                                    
                        if sensorsValues[3]<=safety_margin:
                          
                            brake = True
                            
            if sensorsValues[4] >=.8 and sensorsValues[4] < sensRange:
                
                sNE=[int(round(int(sensorsValues[4])*cos(2*pi-UavOrient+pi/4),0)),-int(round(int(sensorsValues[4])*sin(2*pi-UavOrient+pi/4),0))] 
                
                obs = [cur[0]+sNE[0],cur[1]+sNE[1]]
                    
                if ((obs>=[0,0]) and(obs<=[width, height])):            
                    
                    curobs = nodes[coordinatesToIndex(obs, d_)]
                    
                    if curobs.OBSTACLE != KNOWN_OBSTACLE and curobs.TAG != FORBIDDEN and sNE !=[0,0]:            
                        
                        isDetected = True
                        curobs.OBSTACLE = KNOWN_OBSTACLE
                        curobs.TAG = NEW_FORBIDDEN
                        curobs.cost = INFINI
                        extendedObs.append(curobs)
                                    
                        if sensorsValues[4]<=safety_margin:
                            
                            brake = True
                    
            if sensorsValues[5] >=.8 and sensorsValues[5] < sensRange:
                
                sSE=[int(round(int(sensorsValues[5])*cos(2*pi-UavOrient+7*pi/4),0)),-int(round(int(sensorsValues[5])*sin(2*pi-UavOrient+7*pi/4),0))]
                obs = [cur[0]+sSE[0],cur[1]+sSE[1]]
                    
                if ((obs>=[0,0]) and(obs<=[width, height])):            
                
                    curobs = nodes[coordinatesToIndex(obs, d_)]
                    
                    if curobs.OBSTACLE != KNOWN_OBSTACLE and curobs.TAG != FORBIDDEN and sSE !=[0,0]:            
                        
                        isDetected = True
                        curobs.OBSTACLE = KNOWN_OBSTACLE
                        curobs.TAG = NEW_FORBIDDEN
                        curobs.cost = INFINI
                        extendedObs.append(curobs)
                                    
                        if sensorsValues[5]<=safety_margin:
                          
                            brake = True
                
            if sensorsValues[6] >=.8 and sensorsValues[6] < sensRange:
                
                sSW=[int(round(int(sensorsValues[6])*cos(2*pi-UavOrient+5*pi/4),0)),-int(round(int(sensorsValues[6])*sin(2*pi-UavOrient+5*pi/4),0))] 
                obs = [cur[0]+sSW[0],cur[1]+sSW[1]]
                    
                if ((obs>=[0,0]) and(obs<=[width, height])):            
                                
                    curobs = nodes[coordinatesToIndex(obs, d_)]
                    
                    if curobs.OBSTACLE != KNOWN_OBSTACLE and curobs.TAG != FORBIDDEN and sSW !=[0,0]:            
                        
                        isDetected = True
                        curobs.OBSTACLE = KNOWN_OBSTACLE
                        curobs.TAG = NEW_FORBIDDEN
                        curobs.cost = INFINI
                        extendedObs.append(curobs)
                                    
                        if sensorsValues[6]<=safety_margin:
                          
                            brake = True
                        
            if sensorsValues[7] >=.8 and sensorsValues[7] < sensRange:
    
                sNW=[int(round(int(sensorsValues[7])*cos(2*pi-UavOrient+3*pi/4),0)),-int(round(int(sensorsValues[7])*sin(2*pi-UavOrient+3*pi/4),0))]      
                obs = [cur[0]+sNW[0],cur[1]+sNW[1]]
                    
                if ((obs>=[0,0]) and(obs<=[width, height])):             
                    
                    curobs = nodes[coordinatesToIndex(obs, d_)]
                    
                    if curobs.OBSTACLE != KNOWN_OBSTACLE and curobs.TAG != FORBIDDEN and sNW !=[0,0]:            
                        
                        isDetected = True
                        curobs.OBSTACLE = KNOWN_OBSTACLE
                        curobs.TAG = NEW_FORBIDDEN
                        curobs.cost = INFINI
                        extendedObs.append(curobs)
                                    
                        if sensorsValues[7]<=safety_margin:
                           
                            brake = True
                        
            #print  sensorsValues, sN,sE,sS,sW,sNE, sSE, sSW, sNW   
            sensorsValues = []  
        
    elif sensType=='lidar':#if lidar is the integrated sensor

        risks = getRisk(nodes)
        sensorsValues = sensors.items()        
        cartesval = set()
            
        for s in sensorsValues:
            
            if s[1]<4.0 and s[1] >=0.50:
                
                cartesianS = [int(round(int(s[1])*cos(2*pi-UavOrient+degToCartesian(s[0])),0)),-int(round(int(s[1])*sin(2*pi-UavOrient+degToCartesian(s[0])),0))]
                obs = (cur[0]+cartesianS[0],cur[1]+cartesianS[1])
                cartesval.add(obs)
                
        for obs in cartesval:        
                
                curobs = nodes[coordinatesToIndex(obs, d_)]

                       
                
                if (curobs.OBSTACLE != KNOWN_OBSTACLE and curobs.TAG != FORBIDDEN and cartesianS !=[0,0]):            
                    
                    isDetected = True
                    curobs.OBSTACLE = KNOWN_OBSTACLE
                    curobs.TAG = NEW_FORBIDDEN
                    curobs.cost = INFINI
                    extendedObs.append(curobs)
                                
                    if s[1]<=safety_margin:
                       
                        brake = True
                
                    #print obs, curobs.indice

    #extendedObs = list(set(extendedObs))        
    return extendedObs, isDetected, brake#cartesval#
          
def dessiner(noeud,col, grille):
           
            grille.create_rectangle(noeud.abscice*echelle,noeud.colonne*echelle,noeud.abscice*echelle+echelle,
                                        noeud.colonne*echelle+echelle,outline = col,fill=col)
def dessinerPath(x,y,co, grille, fenetre):
 
            
            grille.create_line(x[0]*echelle,x[1]*echelle,y[0]*echelle,y[1]*echelle,width= 2,fill=co)
            fenetre.update()    
   
def drawTK(Nodes, path = globalPath):
    fenetre =Tk()
    fenetre.title("MAP MOPP")
    grille = Canvas(fenetre,width=width*2,height=height*2,bg='grey')
    grille.pack()
    grille.delete(ALL)
    for i in Nodes:
        #ffcccb ffd5d3 ffdddc ffe6e5 ffeeed fff6f6 ffffff
        pourcent = round(i.v,2)*100
        if pourcent ==0:
            dessiner(i, 'black', grille)
        elif pourcent ==1:
            dessiner(i, 'white', grille)
        else: 
            coleur = "grey"+str(int(pourcent))
            dessiner(i, coleur, grille)
        if i.OBSTACLE ==UNDETECTED_OBSTACLE:
            dessiner(i, 'grey', grille)
    if path !=[]:
        x = path[0]
        path.remove(path[0])
        for i in path:
            dessinerPath(x, i, 'magenta', grille, fenetre)
            dessiner(Nodes[coordinatesToIndex([100,74], d_)], 'green', grille)
            dessiner(Nodes[coordinatesToIndex([20,50], d_)], 'red', grille)
            fenetre.update()
            x= i  
    fenetre.mainloop()     
    
