'''
Created on Jul 8, 2021

@author: nassim
'''
from UavAndSensors.Sensors import Sensors
from utilities import indexToCoordinates, coordinatesToIndex, NEW_FORBIDDEN, DetectUnexpectedObs
from MainMOHPP import  d_
import AStar

def ObsInPath(l, nodes):

    for i in l:
        n = indexToCoordinates(i,d_)
        if nodes[n].TAG == NEW_FORBIDDEN:
            return True
    return False
        
    
def processONPS(start, goal, extendedObs, is_detected, nodes, d_):
    
    current = start
    replannedPath = []

    while extendedObs !=[]:
        
        if replannedPath !=[]:
            
            coords = replannedPath[0]
            curIdx = coordinatesToIndex(coords, d_)
            replannedPath.remove(coords)
            
            extendedObs, is_detected, brake = DetectUnexpectedObs(curIdx, nodes, extendedObs, 2, 4, d_)
            crossObsPath = ObsInPath(replannedPath, nodes)
            
            if is_detected:
                
                if crossObsPath:
                    replannedPath, length = AStar(current, goal, extendedObs, nodes)
                
                else:
                    continue
        else:
            replannedPath, length = AStar(current, goal, extendedObs, nodes,  mode = 1)
            if length ==0 and not brake:
                extendedObs = []
    return replannedPath[-1]
            
                
                
                
            
            
            
            
    