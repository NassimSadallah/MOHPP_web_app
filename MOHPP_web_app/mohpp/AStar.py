'''
Created on Jul 8, 2021

@author: nassim
'''
from heapq import heappush, heappop
from math import floor
from utilities import NEW_FORBIDDEN
from MainMOHPP import width, height, d_
from MAP import GridMap
from mohpp.utilities import indexToCoordinates

class Astar(object):


    def __init__(self, start, goal, extendedObs, nodes, mode = None ):
 
        '''
        Constructor
        '''
        self.current = start
        self.goal = goal
        self.extendedObs = extendedObs
        self.mode = mode
        self.nodes = nodes
        self.opened, self.closed, self.neighbors, self.path = [],[],[],[]
        
        self.computePath()
        
def came_from(self):

    n = self.closed[-1]
    while n.parent != n:

        self.Path.append(n)
        n = n.parent
    return self.Path
     
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

def IndetectionArea(Rang,curElem):
        
    for i in range(-Rang,Rang): 
        if i+curElem.colonne >=0 and i+curElem.colonne<height : 
            for j in range(-Rang,Rang):
                if j+curElem.abscice>=0 and j+curElem.abscice<width:
                    if curElem.colonne+i== curElem.colonne and curElem.abscice ==j+curElem.abscice:
                        continue
                    else:
                        if (((((curElem.colonne+i) -curElem.colonne)**2)+(((curElem.abscice+j) - curElem.abscice)**2))<=Rang**2):
                                if GridMap[curElem.colonne+i][curElem.abscice+j].TAG == NEW_FORBIDDEN:
                                    return True
    return False

def  computePath(self):
    
    g = self.nodes[self.current].G = 0
    h = self.nodes[self.current].H = 0
    self.nodes[self.current].F = g + h 
    self.nodes[self.current].parent = self.start
    heappush(self.opened, (self.nodes[self.current].F, self.nodes[self.current]))
    
    while self.opened !=[]:
        
        self.current = heappop(self.opened)[1]
        self.closed.append(self.current)
        
        if self.mode==1:
            intersect = IndetectionArea(4, self.current)
        
        else:
            intersect = checkLineIntersection(self.current, self.goal, self.extendedObs)
    
        if self.current==self.nodes[self.goal] or not intersect:#.abscice \and self.current.colonne==self.nodes[self.goal].colonne
            came_from()
            coords = [-1,-1]
            i = len(self.Path)
            Paths = []
            while i>0:
                coords = indexToCoordinates(self.Path[i-1].indice,d_)
                Paths.append([coords[0],coords[1]])
                i-=1
    
            return Paths,len(Paths)                
                