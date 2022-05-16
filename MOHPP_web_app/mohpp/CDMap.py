'''
Created on Jul 8, 2021

@author: nassim
'''

from mohpp import FastMarching
from utilities import INFINI, getmaxcost, FAR
import numpy as np

def get_Vel_Cost(dim, nodes, srcObs, Start, Goal, alpha, Velocity, d_, seq = 0, block = 0):
    #working on 3 dimension
    if dim==3:
        print '--> Velocity map3D...', block
        vmap3D = FastMarching.VBLHGS3D(block, srcObs, -1 , nodes, d_, cross =True, seq = 1, block = block)
        vmap3D.processMSFM3D()
        '''
        Scaling the velocity to the MAX velocity of the UAV
        '''
        maxdistance = alpha
        maxvalue = getmaxcost(vmap3D.nodesList)
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
        CDmap3D = FastMarching.MSfm3D(Goal, Start, nodes, d_, cross = True,heuristicAct=True, seq =0, block = block)
        
        return CDmap3D.nodesList
    



        
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
