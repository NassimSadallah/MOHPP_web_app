'''
Created on Jul 8, 2021

@author: nassim
'''
from heapq import heappush, heappop
from math import floor
from utilities import NEW_FORBIDDEN, FORBIDDEN, INFINI, width, height, d_

from mohpp.utilities import indexToCoordinates, coordinatesToIndex, sqrt_dist

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
        self.opened, self.closed, self.neighbors, self.Path = [],[],[],[]
        
        
        
    def came_from(self):
    
        n = self.closed[-1]
        while n.parent != n:
    
            self.Path.append(n)
            n = n.parent
        return self.Path

    def  computePath(self):
        
        g =self.current.G = 0
        h = self.current.H = 0
        self.current.F = g + h 
        self.current.parent = self.current
        heappush(self.opened, (self.current.F, self.current))
        
        while self.opened !=[]:
            
            self.current = heappop(self.opened)[1]
            self.closed.append(self.current)
           
            if self.mode==1:
                intersect = IndetectionArea(4, self.current, self.nodes)
            
            else:
                
                intersect = checkLineIntersection(self.current, self.nodes[self.goal], self.extendedObs)
        
            if self.current==self.nodes[self.goal] or not intersect:#.abscice \and self.current.colonne==self.nodes[self.goal].colonne
                self.came_from()
                coords = [-1,-1]
                i = len(self.Path)
                Paths = []
                while i>0:
                    coords = indexToCoordinates(self.Path[i-1].indice,d_)
                    Paths.append(coords)
                    i-=1
                print 'path by A*',len(Paths)
                return Paths 
            
            else:
                
                neighborList = Neighbours_of(self.current, self.nodes)
                
                if neighborList==[]:
                    return BaseException
                
                for n in neighborList:
                            
                    n.G = n.parent.G + sqrt_dist((n.parent.colonne-n.colonne),(n.parent.abscice-n.abscice))

                    n.H =(n.cost)*199.5#+dist_between(n, goal)*99)#*0.3+(dist_between(n, goal)*0.7)
                    n.F = n.G + n.H
                    m =InTheList(n,self.closed,0)
                    if m !=False :  
                        if n.F < m.F :
                            self.closed.remove(m)
                    else:  
            
                        m = InTheList(n,self.opened,1) 
                        if m !=False: 
                            if n.F < m.F :
    
                                self.opened.remove(m)
                        else:                           
                                n.parent = self.current
                                heappush(self.opened, (n.F,n))


                                            
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
    
        if coordSuivante[0]>=0 and coordSuivante[0]< height and\
        coordSuivante[1]>=0 and coordSuivante[1]<width:           
            
            n = nodes[coordinatesToIndex(coordSuivante, d_)]        
            
            if n.TAG==FORBIDDEN or n.TAG==NEW_FORBIDDEN or n.cost==INFINI :#si le noeud n'est pas un mur

                continue
            else:    
                nei.append(n)
    return nei
                               
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
                                if nodes[coordinatesToIndex([curElem.abscice+j,curElem.colonne+i],d_)].TAG == NEW_FORBIDDEN:
                                    return True
    return False             