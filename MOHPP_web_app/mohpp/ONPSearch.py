
'''
Created on Jul 8, 2021

@author: nassim
'''

from mohpp.utilities import indexToCoordinates, coordinatesToIndex, NEW_FORBIDDEN, DetectUnexpectedObs,sqrt_dist,getNorth_East_Down, globalPath
import AStar
from UavAndSensors import VehiclesMethods as VeMeth
import time
from math import floor



def ObsInPath(l, nodes, d_):

    for i in l:
        n = coordinatesToIndex(i,d_)
        if nodes[n].TAG == NEW_FORBIDDEN:
            return True
    return False
        
    
def processONPS(start, goal, heading,extendedObs, is_detected, nodes, d_,sensors, sensType, UAV, default_alt):
   
    current,curCoord = start, indexToCoordinates(start, d_)
    
    replannedPath = []
    onpsPath = []
    while extendedObs !=[]:
        
        if replannedPath !=[] and replannedPath !=None:
            print "replanned ",len(replannedPath), replannedPath
            nextStep = replannedPath[0]
            replannedPath.remove(nextStep)
            current = coordinatesToIndex(nextStep, d_)
            nodeVel = nodes[current].v
            n, e, d = getNorth_East_Down(curCoord, nextStep, UAV.location.local_frame.down, default_alt)
            print n, e, d,sqrt_dist(n, e, d)
            #set the appropriate speed at which the UAV should travel through the point
            UAV.airspeed = nodeVel
            #send command with the north, east, down( -z) distance to move on 
            VeMeth.UAV().send_NED_velocity(n, e,d, UAV)
            globalPath.append(nextStep)
            #pause the script for the corresponding travel time
            time.sleep(sqrt_dist(n, e, d))
            #nex = input('nextONPS')
            onpsPath.append(nextStep)
            
            curCoord = nextStep
            
            extendedObs, is_detected, brake = DetectUnexpectedObs(sensType,sensors, heading, current, nodes, extendedObs, 2, 4, d_)
            crossObsPath = ObsInPath(replannedPath, nodes, d_)
            
            if is_detected:
                
                if crossObsPath:
                    print 'now to replan'
                    astar= AStar.Astar(nodes[current], goal, extendedObs, nodes)
                    replannedPath = astar.computePath()
                
                else:
                    continue
        else:
            print 'into Astar empty list'
            astar = AStar.Astar(nodes[current], goal, extendedObs, nodes,  mode = 1)
            replannedPath = astar.computePath()
            
        if replannedPath ==[] and not brake:
            extendedObs = []
                
    return onpsPath[-1]
            
                
"""
'''
Created on Jul 8, 2021

@author: nassim
'''
from UavAndSensors.Sensors import Sensors
from utilities import indexToCoordinates, coordinatesToIndex, FORBIDDEN, INFINI, NEW_FORBIDDEN, DetectUnexpectedObs,width, height, d_,sqrt_dist
from heapq import heappush, heappop
from math import floor
import AStar

#from AStar import Astart

def ObsInPath(l, nodes, d_):

    for i in l:
        print i
        n = coordinatesToIndex(i,d_)
        if nodes[n].TAG == NEW_FORBIDDEN:
            return True
    return False
        
    
def processONPS(start, goal, heading, extendedObs, is_detected, nodes, d_,sensors):
    brake = False
    current = nodes[start]
    replannedPath = []

    while extendedObs !=[]:
        
        if replannedPath !=[]:
            print 'replanned '
            conextStep replannedPath[0]
            print coords
            curIdx = coordinatesToIndex(coords, d_)
            replannedPath.remove(coords)
            
            extendedObs, is_detected, brake = DetectUnexpectedObs(sensors,heading, curIdx, nodes, extendedObs, .5, 4, d_)
            crossObsPath = ObsInPath(replannedPath, nodes, d_)
            
            if is_detected:
                
                if crossObsPath:
            #ast = AStar.Astar(current, goal, extendedObs, nodes)
                    
                    replannedPath = AStar.Astar(current, nodes[goal], extendedObs, nodes)
                
                else:
                    continue
        else:
            
            #ast = AStar.Astar(current, goal, extendedObs, nodes,  mode = 1)

            replannedPath = AStar.Astar(current, nodes[goal], extendedObs, nodes, mode = 1)
            print 'length ', (replannedPath), 'brake ', brake
            
            if replannedPath ==[] and not brake:
                extendedObs = []
                
    print 'replanned path',replannedPath[-1],'mode 1' 
    return replannedPath[-1]



def  Astart(current, goal, extendedObs, nodes, mode = None):
        opened, closed, neighbors, path = [],[],[],[]
        g = current.G = 0
        h = current.H = 0
        current.F = g + h 
        current.parent = current
        heappush(opened, (current.F, current))
    
        while opened !=[]:
        
            current = heappop(opened)[1]
            closed.append(current)
            
            if mode==1:
                intersect = IndetectionArea(3, current, nodes)
        
            else:
            
            intersect = checkLineIntersection(current, goal, extendedObs)
    
            if current==goal or not intersect:#.abscice \and self.current.colonne==self.nodes[self.goal].colonne
                
                Path = came_from(current)
                coords = [-1,-1]
                i = len(Path)
                Paths = []
                while i>0:
                    coords = indexToCoordinates(Path[i-1].indice,d_)
                    Paths.append([coords[0],coords[1]])
                    i-=1
                return Paths
            
            else:
                        neighborList = Neighbours_of(current, nodes)
                
                        if neighborList==[]:
                                return BaseException
                
                        for n in neighborList:
                            
                            n.G = n.parent.G + sqrt_dist((n.parent.colonne-n.colonne),(n.parent.abscice-n.abscice))
                        n.H =(n.cost)*199.5#+dist_between(n, goal)*99)#*0.3+(dist_between(n, goal)*0.7)
                        n.F = n.G + n.H
                        m =InTheList(n,closed,0)
                        
                        if m !=False :  
                            if n.F < m.F :
                                    closed.remove(m)
                        else:  
            
                            m = InTheList(n,opened,1) 
                            if m !=False: 
                                if n.F < m.F :
    
                                        opened.remove(m)
                                else:                           
                                    n.parent = current
                                    heappush(opened, (n.F,n))

                                    
def InTheList(c,cl,chif):
        if chif==0:
            for k in cl:
                if c.colonne == k.colonne and c.abscice == k.abscice:
                    return c
        else : 
            for k in cl:
                if c.colonne == k[1].colonne and c.abscice == k[1].abscice:
                    return c
        return False 

def Neighbours_of(c, nodes):
    v = [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]
    nei =[]
    for i in v:
        coordSuivante = [c.abscice + i[0],c.colonne + i[1]]
    #print coordSuivante, c.abscice, c.colonne
        
        if coordSuivante[0]>=0 and coordSuivante[0]< width and\
        coordSuivante[1]>=0 and coordSuivante[1]<height:           
            
            n = nodes[coordinatesToIndex(coordSuivante, d_)]        
            
            if n.TAG==FORBIDDEN or n.TAG==NEW_FORBIDDEN or n.cost==INFINI :#si le noeud n'est pas un mur

                continue
            else:    
                nei.append(n)
    return nei
              
                                

def came_from(n):
    Path = []
    
    while n.parent != n:
        Path.append(n)
        n = n.parent
    #print 'came from',Path
    return Path    
            
def checkLineIntersection(co, goa, zo):       
       
        result = False

        if goa.abscice==co.abscice:
           
            for obstacle in zo:
                
                if obstacle.abscice ==goa.abscice and\
                abs(goa.colonne - co.colonne)== (abs(obstacle.colonne-goa.colonne)):
                    return True
            
        elif goa.colonne==co.colonne:
            
            for obstacle in zo:
               
                if obstacle.colonne == co.colonne and\
                abs(goa.abscice-co.abscice) > (abs(obstacle.abscice-goa.abscice)):
                    return True     
        else:    
            
            m = (goa.colonne - co.colonne)/(goa.abscice-co.abscice)# coefficient directeur
            p = co.colonne - (co.abscice*m)
            dCuGo = ((co.abscice-goa.abscice)**2)+((co.colonne-goa.colonne)**2)
            
            for obstacle in zo:
                
                dObsGo =(((goa.abscice-obstacle.abscice)**2)+((goa.colonne-obstacle.colonne)**2))
                y = m*obstacle.abscice + p            
            
                if obstacle.colonne<=floor(y) and obstacle.colonne >= floor(y) and dCuGo > dObsGo:
                        
                        return True                        
        return False

def IndetectionArea(Rang,curElem, nodes):
        
    for i in range(-Rang,Rang): 
        if i+curElem.colonne >=0 and i+curElem.colonne<height : 
            for j in range(-Rang,Rang):
                if j+curElem.abscice>=0 and j+curElem.abscice<width:
                    if curElem.colonne+i== curElem.colonne and curElem.abscice ==j+curElem.abscice:
                        continue
                    else:
                  
                        if (((((curElem.colonne+i) -curElem.colonne)**2)+(((curElem.abscice+j) - curElem.abscice)**2))<=Rang**2):
                                if nodes[coordinatesToIndex([curElem.abscice+i,curElem.colonne+j],d_)].TAG == NEW_FORBIDDEN:
                                    return True
    return False

                
                
            
            
"""            
            
    