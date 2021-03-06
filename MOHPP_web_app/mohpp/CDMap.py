'''
Created on Jul 8, 2021

@author: nassim
'''

from mohpp import FastMarching
from utilities import INFINI, getmaxcost, FAR
import numpy as np
import time

def get_Vel_Cost(dim = -1, nodes =[], srcObs=[], Start=[], Goal=-1, alpha=0, Velocity=0, d_=[], seq = 0, block = 0):
    #working on 3 dimension
    #if dim==3 and seq ==10:
    print '--> Velocity map3Dparallel ...', block
    
    
    
    
    #vmap3D = FastMarching.ParallelMSFM( srcObs, -1 , nodes, d_, cross =True, seq = 1, block = block)

    '''
    Scaling the velocity to the MAX velocity of the UAV
    '''
    maxvalue = 0
    for i in range(block):
        maxdistance = alpha
        if maxvalue<getmaxcost(nodes[i]):
            maxvalue = getmaxcost(nodes[i])
        maxvelocity = maxdistance/1.0
    print(maxvelocity, maxvalue)
        #time.sleep(10)
    for j in range(block):    
        for i in nodes[j]:
    
            vel = i.cost/maxvalue
            #i.risk = (1/0.001+vel)
            
            if vel<maxvelocity:
                i.v = (vel/maxvelocity)*Velocity
            else :
                i.v = 1*Velocity
            
            i.cost = INFINI
            i.type = FAR  
    

    '''
    assigns the final travel time cost to each node based on its scaled velocity 
    '''
    
    #print '--> Cost computation...'
    #CDmap3D = FastMarching.ParallelMSFM(Goal, Start, vmap3D.nodesList, d_, cross = True,heuristicAct=True, seq =1, block = block)
    
    return nodes#CDmap3D.nodesList


     

        
    x_nodes = np.asarray(nodes)
    print '--> Velocity map...'
    vmap = FastMarching.fm2D(srcObs, -1 , nodes, d_, cross =False, seq = 1, block = block)
    #vmap.processFMM()
    '''
    Scaling the velocity to the MAX velocity of the UAV
    '''
    maxdistance = alpha
    maxvalue = getmaxcost(nodes)
    maxvelocity = maxdistance/1.0
  
    for i in nodes:

        vel = i.cost/maxvalue
        #i.risk = (1/0.001+vel)
        
        if vel<maxvelocity:
            i.v = (vel/maxvelocity)*Velocity
        else :
            i.v = 1*Velocity
        
        i.cost = INFINI
        i.type = FAR  
        

    '''
    assigns the final travel time cost to each node based on its scaled velocity 
    '''
    
    print '--> CDMap...'
    CDmap = FastMarching.heuristicFM2D(Goal, Start, nodes, d_, cross = True, seq =0, block = block)
    
    return CDmap.nodesList
