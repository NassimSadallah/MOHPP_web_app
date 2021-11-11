'''

Created on Jul 8, 2021

@author: nassim
'''
from __builtin__ import False
import time
from math import floor, sqrt
from heapq import _siftdown, heappop, heappush
from utilities import isEmptyList, locateBestIdx, FORBIDDEN, INFINI, NEW_FORBIDDEN
from mohpp.utilities import sqrt_dist


def getneighbors(self,idx):  
          
    self.nbr_neighbors = 0
    self.nbr_neighbors = getNeighborsInDim(self, idx,self.neighborsList, -1, False)
    return self.neighborsList

def getNeighborsInDim(self,idx,ne,dim , crossing):
    
    c1, c2, c3, c4, cc1, cc2, cc3, cc4 = -1, -1, -1, -1, -1, -1, -1, -1
    
    if ne==self.neighborsList:
        # neighbors over the X and Y axes
        c1 = idx-1
        c2 = idx+1
        c3 = int(floor(idx-self.d_[0]))
        c4 = int(floor(idx+self.d_[0]))

        if ((c1>=0) and(int(floor(c1/self.d_[1]))==int(floor(idx/self.d_[1])))):
            self.nodesList[c1].cros = False
            self.neighborsList[self.nbr_neighbors] = int(c1)
            self.nbr_neighbors+=1
        
        if(c2>=0)and (int(floor(c2/self.d_[1]))==int(floor(idx/self.d_[1]))):
            self.nodesList[c2].cros = False
            self.neighborsList[self.nbr_neighbors]= int(c2)
            self.nbr_neighbors+=1
            
        if ((c3>=0) and(int(floor(c3/self.d_[1]))==int(floor(idx/self.d_[1])))):
            self.nodesList[c3].cros = False
            self.neighborsList[self.nbr_neighbors] = int(c3)
            self.nbr_neighbors+=1
        
        if(c4>=0)and (int(floor(c4/self.d_[1]))==int(floor(idx/self.d_[1]))):
            self.nodesList[c4].cros = False
            self.neighborsList[self.nbr_neighbors]= int(c4)
            self.nbr_neighbors+=1
            
        '''
        the crossing neighbors if activated
        '''
            
        if self.cross:
            # neighbors indices over the crossing stencils
            cc1 = int(floor(idx+self.d_[0]+1))
            cc2 = int(floor(idx-self.d_[0]-1))
            cc3 = int(floor(idx+self.d_[0]-1))
            cc4 = int(floor(idx-self.d_[0]+1))
            
            if ((cc1>=0) and (cc1<self.d_[1])and (int(floor(cc1/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                self.nodesList[int(cc1)].cros = True
                self.neighborsList[self.nbr_neighbors] = int(cc1)
                self.nbr_neighbors+=1
                
            if((cc2>=0) and (cc2<self.d_[1]) and (int(floor(cc2/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                self.nodesList[int(cc2)].cros = True
                self.neighborsList[self.nbr_neighbors] = int(cc2)
                self.nbr_neighbors+=1
                
                
            if((cc3>=0) and (cc3<self.d_[1]) and (int(floor(cc3/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                self.nodesList[int(cc3)].cros = True
                self.neighborsList[self.nbr_neighbors] = int(cc3)
                self.nbr_neighbors+=1
                
            if((cc4>=0) and (cc4<self.d_[1]) and (int(floor(cc4/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                self.nodesList[int(cc4)].cros = True
                self.neighborsList[self.nbr_neighbors] = int(cc4)
                self.nbr_neighbors+=1

        return self.nbr_neighbors

    elif ne==self.minneighborsList:
        
        '''
        looks for the best neighbors over the same axe
        '''          
        if crossing == False:
            c1 = idx-1
            c2 = idx+1
            c3 = int(floor(idx-self.d_[0]))
            c4 = int(floor(idx+self.d_[0]))
            
            if dim == 0:
                    
                if ((c1>=0) and(int(floor(c1/self.d_[1]))==int(floor(idx/self.d_[1])))):
                    self.minneighborsList[self.nbr_neighbors] = int(c1)
                    self.nbr_neighbors+=1
                
                if(c2>=0)and (int(floor(c2/self.d_[1]))==int(floor(idx/self.d_[1]))):
                    self.minneighborsList[self.nbr_neighbors]= int(c2)
                    self.nbr_neighbors+=1
            
            elif dim==1:    
            
                if ((c3>=0) and(int(floor(c3/self.d_[1]))==int(floor(idx/self.d_[1])))):
                    self.minneighborsList[self.nbr_neighbors] = int(c3)
                    self.nbr_neighbors+=1
                
                if(c4>=0)and (int(floor(c4/self.d_[1]))==int(floor(idx/self.d_[1]))):
                    self.minneighborsList[self.nbr_neighbors]= int(c4)
                    self.nbr_neighbors+=1
            return self.nbr_neighbors

        else:
            
            cc1 = int(floor(idx+self.d_[0]+1))
            cc2 = int(floor(idx-self.d_[0]-1))
            cc3 = int(floor(idx+self.d_[0]-1))
            cc4 = int(floor(idx-self.d_[0]+1))
            
            if dim == 0:
                
                if ((cc1>=0) and (cc1<self.d_[1])and (int(floor(cc1/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                    self.minneighborsList[self.nbr_neighbors] = int(cc1)
                    self.nbr_neighbors+=1
                    
                if((cc2>=0) and (cc2<self.d_[1]) and (int(floor(cc2/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                    self.minneighborsList[self.nbr_neighbors] = int(cc2)
                    self.nbr_neighbors+=1
            
            elif dim == 1:
                      
                if((cc3>=0) and (cc3<self.d_[1]) and (int(floor(cc3/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                    self.minneighborsList[self.nbr_neighbors] = int(cc3)
                    self.nbr_neighbors+=1
                    
                if((cc4>=0) and (cc4<self.d_[1]) and (int(floor(cc4/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                    self.minneighborsList[self.nbr_neighbors] = int(cc4)
                    self.nbr_neighbors+=1
            
            return self.nbr_neighbors

def getMinValueInDim(self,idx, dim, crossing):

    self.nbr_neighbors = 0      
    
    if crossing == False: 
        getNeighborsInDim(self, idx, self.minneighborsList, dim, False)
       
        if (self.nbr_neighbors==0):
            return self.nodesList[self.minneighborsList[0]].cost
        
        elif self.nodesList[self.minneighborsList[0]].cost < self.nodesList[self.minneighborsList[1]].cost:
                return self.nodesList[self.minneighborsList[0]].cost
        else:
                return self.nodesList[self.minneighborsList[1]].cost  
    
    else:
        getNeighborsInDim(self, idx, self.minneighborsList, dim, True)
       
        if (self.nbr_neighbors==0):
            return self.nodesList[self.minneighborsList[0]].cost
        
        elif self.nodesList[self.minneighborsList[0]].cost < self.nodesList[self.minneighborsList[1]].cost:
            return self.nodesList[self.minneighborsList[0]].cost
        
        else:
            return self.nodesList[self.minneighborsList[1]].cost             
   
def solveEikonalNDims(self,idx,dim):
    print dim
    curnode = self.nodesList[idx]
    
    if dim==1:
   
        if curnode.cros:
         
            return self.Tvalue[0]+sqrt(2)/curnode.v
            
        else:
            return self.Tvalue[0]+1/curnode.v
    
    sumT = 0
    sumTT = 0    
    
    for i in range(0,dim):
        sumT +=self.Tvalue[i]
        sumTT +=self.Tvalue[i]*self.Tvalue[i]
    
    a=dim
    b = -2*sumT
    c=  0 
    if curnode.cros:
        c = sumTT-sqrt(2)*sqrt(2)/(curnode.v * curnode.v)
    else:
        c = sumTT-1*1/(curnode.v * curnode.v)
    q = b*b - 4*a*c
    if q<0:
        return INFINI
    else:
        return ((-b+sqrt(q))/(2*a))

def SolveEikonal(self,idx, crossing):
    nodecur = self.nodesList[idx]
    if crossing ==False:
        a = 2
                 
        self.Tvalue = []
        
        for i in range(a):
            minT = 0
            minT = getMinValueInDim(self, idx,i, False)
                       
            if minT != INFINI and minT< nodecur.cost:
                heappush(self.Tvalue,minT) 
            
            else:
                a-=1
        
        if a==0:
            return INFINI
        
        self.Tvalue = sorted(self.Tvalue)
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = solveEikonalNDims(self, idx,i)
            
            if i==a or updatedT-self.Tvalue[i]<=0:
                break
        
        return updatedT
    
    if crossing == True:
        
        a = 2
                 
        self.Tvalue = []
        
        for i in range(a):
            minT = 0
            minT = getMinValueInDim(self, idx,i, True)
                       
            if minT != INFINI and minT< nodecur.cost:
                heappush(self.Tvalue,minT) 
            
            else:
                a-=1
        
        if a==0:
            return INFINI
        
        self.Tvalue = sorted(self.Tvalue)
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = solveEikonalNDims(self, idx,i)
            
            if i==a or updatedT-self.Tvalue[i]<=0:
                break
        
        return updatedT

class fmm(object):
    '''
    Fast Marching method with (8 and/or 16 neighbors )
    '''
    
    def __init__(self, start_points=[], goal=None, nodesList=[], d_ =[], cross=False, heuristicAct = False, seq = 0, block = 0):
        
        self.d_ = d_
        self.narrowBand = []
        self.stopWave = False
        self.nbr_neighbors = 0
        self.Tvalue = []
        self.start_points = start_points
        self.goal = goal
        self.nodesList = nodesList
        self.cross = cross
        self.block = block
        self.heuristicAct = heuristicAct
        if self.cross:
            if len(self.d_)==2:
                self.neighborsList = [-1 for _ in range(8)]
            else:
                self.neighborsList = [-1 for _ in range(26)]
        else:
            self.neighborsList = [-1 for _ in range(2**(len(d_)-1)+2)]
        self.minneighborsList = [-1 for _ in range(2)]

        if seq == 0:
            self.processFMM()
        elif seq == 1:
            self.processSequentialFMM()

    def processFMM(self):

        itera = 0
        
        '''
        initialization of the starting nodes
        '''      
        
        for i in self.start_points:
            self.nodesList[i].cost = 0
            heappush(self.narrowBand, (self.nodesList[i].cost, i))
        
        '''
        main loop over the node list
        '''
        initime = time.time()
        while self.narrowBand !=[] and not self.stopWave:

            itera+=1
            best = heappop(self.narrowBand)
            idxMin = best[1]
            neighbors = self.getneighbors(idxMin)
            self.nodesList[idxMin].type = 'A'
            
            for j in neighbors:
                
                if j==-1:
                    continue
                 
                current = self.nodesList[j]
             
                if current.type == 'A' or current.TAG ==FORBIDDEN or current.TAG == NEW_FORBIDDEN:
                    continue 
               
                cost = self.SolveEikonal(j, current.cros)
                
                if self.narrowBand.__contains__((current.cost,j)):
                    
                    if cost<current.cost:
                        indexe = self.narrowBand.index((current.cost,j))
                        current.cost = cost 
                        _siftdown(self.narrowBand,0,indexe)
                else:            
                    current.cost = cost
                    current.type = 'N'
                    heappush(self.narrowBand,(cost,j))
        
            if idxMin == self.goal:
                self.stopWave = True                    
        print 'seq 0:',time.time()-initime,'seconds '    
        return self.nodesList    
    
    def processSequentialFMM(self):

        blk_nodes = [[] for _ in range(self.block)]
        '''
        initialization of the starting nodes over the blocks
        '''      
        
        for i in self.start_points:
            self.nodesList[i].cost = 0
            heappush(blk_nodes[self.nodesList[i].block - 1], (self.nodesList[i].cost, i))
        EmptyList = isEmptyList(blk_nodes)
        '''
        main loop over the node list
        '''
        initime = time.time()
        while not EmptyList and not self.stopWave:

            best = locateBestIdx(blk_nodes)
            idxMin = best[1]
            neighbors = getneighbors(self, idxMin)
            self.nodesList[idxMin].type = 'A'
            
            for j in neighbors:
                
                if j==-1:
                    continue
                 
                current = self.nodesList[j]
             
                if current.type == 'A' or current.TAG ==FORBIDDEN or current.TAG == NEW_FORBIDDEN:
                    continue 
               
                cost = SolveEikonal(self, j, current.cros)
                
                if blk_nodes[current.block - 1].__contains__((current.cost,j)):
                    
                    if cost<current.cost:
                        indexe = blk_nodes[current.block - 1].index((current.cost,j))
                        current.cost = cost 
                        _siftdown(blk_nodes[current.block - 1],0,indexe)            
                else:            
                    current.cost = cost
                    current.type = 'N'
                    heappush(blk_nodes[current.block - 1],(cost,j))
            
            EmptyList = isEmptyList(blk_nodes)
            
            if idxMin == self.goal:
                self.stopWave = True                    
        print time.time() - initime, 'seconds'    
        return self.nodesList                

class heuristicFMM(object):
    '''
    Fast Marching heuristic activated method with (8 and/or 16 neighbors (multi-stencils FMM))
    '''
    
    def __init__(self, start_points=[], goal=None, nodesList=[], d_ =[], cross=False, heuristicAct=False, seq = 0, block = 0):
        
        self.d_ = d_
        self.narrowBand = []
        self.stopWave = False
        self.nbr_neighbors = 0
        self.Tvalue = []
        self.start_points = start_points
        self.goal = goal
        self.nodesList = nodesList
        self.cross = cross
        self.heuristicActivate = heuristicAct
        self.block = block
        if self.cross:
            self.neighborsList = [-1 for _ in range(8)]
            
        else:
            self.neighborsList = [-1 for _ in range(4)]
        self.minneighborsList = [-1 for _ in range(2)]
        
        if seq == 0:
            self.processHFMM()
        elif seq == 1:
            self.processSequentialHFMM()

    def processHFMM(self):

        itera = 0
        
        '''
        initialization of the starting nodes
        '''      
        
        for i in self.start_points:
            self.nodesList[i].cost = 0
            self.nodesList[i].hCost = sqrt_dist((self.nodesList[i].abscice-self.nodesList[self.goal].abscice),
                                            (self.nodesList[i].colonne-self.nodesList[self.goal].colonne))
            self.nodesList[i].full =self.nodesList[i].cost + (self.nodesList[i].hCost/self.nodesList[i].v -
                                                          self.nodesList[i].hCost/1.0) #1.0 is the max VELOCITY            
            heappush(self.narrowBand, (self.nodesList[i].full, i))
        
        '''
        main loop over the node list
        '''
        initime = time.time()
        while self.narrowBand !=[] and not self.stopWave:

            itera+=1
            best = heappop(self.narrowBand)
            idxMin = best[1]
            neighbors = getneighbors(self, idxMin)
            self.nodesList[idxMin].type = 'A'
            
            for j in neighbors:
                
                if j==-1:
                    continue
                 
                current = self.nodesList[j]
             
                if current.type == 'A' or current.TAG ==FORBIDDEN or current.TAG == NEW_FORBIDDEN:
                    continue 
               
                cost = SolveEikonal(self, j, current.cros)
                
                hcost = sqrt_dist((current.abscice-self.nodesList[self.goal].abscice),
                                            (current.colonne-self.nodesList[self.goal].colonne))
                fcost = cost + (hcost/current.v - hcost/1.0) #1.0 is the max VELOCITY 
                
                if self.narrowBand.__contains__((current.full,j)):
                    
                    if fcost<current.full:
                        indexe = self.narrowBand.index((current.full,j))
                        current.cost = cost
                        current.hCost = hcost
                        current.full = fcost 
                        _siftdown(self.narrowBand,0,indexe)
                else:            
                    current.cost = cost
                    current.hCost = hcost
                    current.full = fcost                    
                    current.type = 'N'
                    heappush(self.narrowBand,(fcost,j))
        
            if idxMin == self.goal:
                self.stopWave = True                    
        print 'seq 0:',time.time()-initime,'seconds '    
        return self.nodesList    
    
    def processSequentialHFMM(self):

        blk_nodes = [[] for _ in range(self.block)]
        '''
        initialization of the starting nodes over the blocks
        '''      
        
        for i in self.start_points:
            self.nodesList[i].cost = 0
            self.nodesList[i].hCost = sqrt_dist((self.nodesList[i].abscice-self.nodesList[self.goal].abscice),
                                            (self.nodesList[i].colonne-self.nodesList[self.goal].colonne))
            self.nodesList[i].full =self.nodesList[i].cost + (self.nodesList[i].hCost/self.nodesList[i].v -
                                                          self.nodesList[i].hCost/1.0) #1.0 is the max VELOCITY
               
            heappush(blk_nodes[self.nodesList[i].block - 1], (self.nodesList[i].full, i))
        EmptyList = isEmptyList(blk_nodes)
        '''
        main loop over the node list
        '''
        initime = time.time()
        while not EmptyList and not self.stopWave:

            best = locateBestIdx(blk_nodes)
            idxMin = best[1]
            neighbors = getneighbors(self, idxMin)
            self.nodesList[idxMin].type = 'A'
            
            for j in neighbors:
                
                if j==-1:
                    continue
                 
                current = self.nodesList[j]
             
                if current.type == 'A' or current.TAG ==FORBIDDEN or current.TAG == NEW_FORBIDDEN:
                    continue 
               
                cost = SolveEikonal(self, j, current.cros)
                hcost = sqrt_dist((current.abscice-self.nodesList[self.goal].abscice),
                                            (current.colonne-self.nodesList[self.goal].colonne))
                fcost = cost + (hcost/current.v - hcost/1.0) #1.0 is the max VELOCITY                 
                
                if blk_nodes[current.block - 1].__contains__((current.full,j)):
                    
                    if cost<current.cost:
                        indexe = blk_nodes[current.block - 1].index((current.full,j))
                        current.cost = cost 
                        current.hCost = hcost
                        current.full = fcost                         
                        _siftdown(blk_nodes[current.block - 1],0,indexe)            
                else:            
                    current.cost = cost
                    current.hCost = hcost
                    current.full = fcost                     
                    current.type = 'N'
                    heappush(blk_nodes[current.block - 1],(fcost,j))
            
            EmptyList = isEmptyList(blk_nodes)
            
            if idxMin == self.goal:
                self.stopWave = True                    
        print time.time() - initime, 'seconds'    
        return self.nodesList                

'''
3 dimension version
'''
class fmm3D(object):
    
    
    def __init__(self, start_points=[], goal=None, nodesList=[], d_ =[], cross=False, heuristicAct = False, seq = 0, block = 0):
        
        self.d_ = d_
        self.narrowBand = []
        self.stopWave = False
        self.nbr_neighbors = 0
        self.Tvalue = []
        self.start_points = start_points
        self.goal = goal
        self.nodesList = nodesList
        self.cross = cross
        self.block = block
        self.heuristicAct = heuristicAct
        if self.cross:

            self.neighborsList = [-1 for _ in range(26)]
        else:
            self.neighborsList = [-1 for _ in range(2**(len(d_)-1)+2)]
        self.minneighborsList = [-1 for _ in range(2)]
    
    def getneighbors(self, idx,ne,dims):
        
        global neighbors,neighbors_cross,neighbors_cross1,neighbors_cross2,neighbors_cross3,\
        neighbors_cross4
        
        if ne == neighbors:
            for i in range(dims): #x,y,z dimension
                nei =  getNeighborsInDim(idx,neighbors,i)
                n_neighs = nei     
            return n_neighs   
     
        elif ne ==neighbors_cross :            
            # cross dimension
            n_neighs_cross=self.getNeighborsInDimcross(idx,neighbors_cross)
            return n_neighs_cross
        
        elif ne ==neighbors_cross1 :            
            # cross dimension
            n_neighs_cross1=self.getNeighborsInDimcross(idx,neighbors_cross1)
            return n_neighs_cross1
            
        elif ne ==neighbors_cross2 :            
            # cross dimension
            n_neighs_cross2=self.getNeighborsInDimcross(idx,neighbors_cross2)
            return n_neighs_cross2
        
        elif ne ==neighbors_cross3 :            
            # cross dimension
            n_neighs_cross3=self.getNeighborsInDimcross(idx,neighbors_cross3)
            return n_neighs_cross3
        
        elif ne ==neighbors_cross4 :            
            # cross dimension
            n_neighs_cross4=self.getNeighborsInDimcross(idx,neighbors_cross4)
            return n_neighs_cross4
       
    def getNeighborsInDim(self, idx,ne,i):
          
        global c1,c2,n_neighs,neighbors
        
        if ne==neighbors:
         
            if i==0:
                c1 = idx-1
                c2 = idx+1
              
            elif i==1:
                c1 = int(floor(idx-self.d_[0]))
                c2 = int(floor(idx+self.d_[0]))
            
            elif i==2:
                c1 = int(floor(idx-self.d_[1]))
                c2 = int(floor(idx+self.d_[1]))
            
                
            if ((c1>=0) and (int(floor(c1/self.d_[i]))==int(floor(idx/self.d_[i])))):
                ne[n_neighs] = int(c1)
                n_neighs+=1
            
            if((c2>=0) and (int(floor(c2/self.d_[i]))==int(floor(idx/self.d_[i])))):
                ne[n_neighs]= int(c2)
                n_neighs+=1
            
            
            return n_neighs
        
        elif ne==minneighbors:
            n = 0
            if i==0:
                c1 = idx-1
                c2 = idx+1
            else:
                c1 = int(floor(idx-self.d_[i-1]))
                c2 = int(floor(idx+self.d_[i-1]))
    
            if ((c1>=0) and (int(floor(c1/self.d_[i]))==int(floor(idx/self.d_[i])))):
                ne[n] = int(c1)
                n+=1
                
            
            if((c2>=0) and(int(floor(c2/self.d_[i]))==int(floor(idx/self.d_[i])))):
                ne[n]= int(c2)
                n+=1
            
        return n

    def getNeighborsInDimcross(self, idx,nec):
        
        global  cc1,cc2,cc3,cc4,tcc1,tcc2,tcc3,tcc4,tcc5,tcc6,tcc7,tcc8,\
        bcc1,bcc2,bcc3,bcc4,bcc5,bcc6,bcc7,bcc8,\
        n_neighs_cross,n_neighs_cross1,n_neighs_cross2,n_neighs_cross3,n_neighs_cross4,\
        neighbors_cross,neighbors_cross1,neighbors_cross2,neighbors_cross3,neighbors_cross4
        n_neighs_cross,n_neighs_cross1,n_neighs_cross2,n_neighs_cross3,n_neighs_cross4 = 0,0,0,0,0
        
        
        if nec==neighbors_cross:
        
            cc1 = int(floor(idx+self.d_[0]+1))
            cc2 = int(floor(idx-self.d_[0]-1))
            cc3 = int(floor(idx+self.d_[0]-1))
            cc4 = int(floor(idx-self.d_[0]+1))
            
            if ((cc1>=0) and (cc1<self.d[1]) and (int(floor(cc1/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                nec[n_neighs_cross] = int(cc1)       
                n_neighs_cross+=1
            
            if((cc2>=0) and (cc2<self.d[1]) and (int(floor(cc2/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                nec[n_neighs_cross]= int(cc2)    
                n_neighs_cross+=1
            
            if((cc3>=0) and (cc3<self.d[1]) and (int(floor(cc3/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                nec[n_neighs_cross]= int(cc3)  
                n_neighs_cross+=1
            
            if((cc4>=0) and (cc4<self.d[1]) and (int(floor(cc4/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                nec[n_neighs_cross]= int(cc4)
                n_neighs_cross+=1

            return n_neighs_cross
        
        elif nec==neighbors_cross1:
            
            tcc1 = int(floor(idx+self.d_[1]-1))
            tcc2 = int(floor(idx-self.d_[1]+1))
            tcc3 = int(floor(idx-self.d_[1]-1))
            tcc4 = int(floor(idx+self.d_[1]+1))
            #Top diagonals
            if ((tcc1>=0) and (tcc1<self.d[2])and (int(floor(tcc1/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross1] = int(tcc1)       
                n_neighs_cross1+=1
            
            if((tcc2>=0) and (tcc2<self.d[2]) and (int(floor(tcc2/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross1]= int(tcc2)    
                n_neighs_cross1+=1
            
            if((tcc3>=0) and (tcc3<self.d[2]) and (int(floor(tcc3/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross1]= int(tcc3)  
                n_neighs_cross1+=1
            
            if((tcc4>=0) and (cc4<self.d[2]) and (int(floor(tcc4/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross1]= int(tcc4)
                n_neighs_cross1+=1

                  
            return n_neighs_cross1

        elif nec==neighbors_cross2:        
            
            tcc5 = int(floor(idx+self.d_[1]-self.d_[0]))
            tcc6 = int(floor(idx-self.d_[1]+self.d_[0]))
            tcc7 = int(floor(idx-self.d_[1]-self.d_[0]))
            tcc8 = int(floor(idx+self.d_[1]+self.d_[0])) 
                                   
            if ((tcc5>=0) and (tcc5<self.d[2])and (int(floor(tcc5/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross2] = int(tcc5)       
                n_neighs_cross2+=1
            
            if((tcc6>=0) and (tcc6<self.d[2]) and (int(floor(tcc6/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross2]= int(tcc6)    
                n_neighs_cross2+=1
            
            if((tcc7>=0) and (tcc7<self.d[2]) and (int(floor(tcc7/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross2]= int(tcc7)  
                n_neighs_cross2+=1
            
            if((tcc8>=0) and (tcc8<self.d[2]) and (int(floor(tcc8/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross2]= int(tcc8)
                n_neighs_cross2+=1
        
            return n_neighs_cross2

        elif nec==neighbors_cross3:
        
            bcc1 = int(floor(idx-self.d_[1]+self.d_[0]-1))
            bcc2 = int(floor(idx+self.d_[1]-self.d_[0]+1))
            bcc3 = int(floor(idx+self.d_[1]+self.d_[0]-1))
            bcc4 = int(floor(idx-self.d_[1]-self.d_[0]+1))   
            #Bottom diagonals   
            if ((bcc1>=0) and (bcc1<self.d[2])and (int(floor(bcc1/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross3] = int(bcc1)       
                n_neighs_cross3+=1
            
            if((bcc2>=0) and (bcc2<self.d[2]) and (int(floor(bcc2/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross3]= int(bcc2)    
                n_neighs_cross3+=1
            
            if((bcc3>=0) and (bcc3<self.d[2]) and (int(floor(bcc3/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross3]= int(bcc3)  
                n_neighs_cross3+=1
            
            if((bcc4>=0) and (bcc4<self.d[2]) and (int(floor(bcc4/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross3]= int(bcc4)
                n_neighs_cross3+=1
            
            return n_neighs_cross3
    
        elif nec==neighbors_cross4:
            
            bcc5 = int(floor(idx+self.d_[1]-self.d_[0]-1))
            bcc6 = int(floor(idx-self.d_[1]+self.d_[0]+1))
            bcc7 = int(floor(idx-self.d_[1]-self.d_[0]-1))
            bcc8 = int(floor(idx+self.d_[1]+self.d_[0]+1))
            
            if ((bcc5>=0) and (bcc5<self.d[2])and (int(floor(bcc5/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross4] = int(bcc5)       
                n_neighs_cross4+=1
            
            if((bcc6>=0) and (bcc6<self.d[2]) and (int(floor(bcc6/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross4]= int(bcc6)    
                n_neighs_cross4+=1
            
            if((bcc7>=0) and (bcc7<self.d[2]) and (int(floor(bcc7/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross4]= int(bcc7)  
                n_neighs_cross4+=1
            
            if((bcc8>=0) and (bcc8<self.d[2]) and (int(floor(bcc8/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross4]= int(bcc8)
                n_neighs_cross4+=1
            
            return n_neighs_cross4
    
    def getMinNeighborsInDimcross(self, idx, ne, i):    
        
        global c1,c2,c3,c4
        c1,c2,c3,c4=-1,-1,-1,-1
        n = 0
        if i==0:
            c1 = int(floor(idx+self.d_[0]+1))
            c2 = int(floor(idx-self.d_[0]-1))
        
            if ((c1>=0) and (c1<self.d[1]) and (int(floor(c1/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[n] = int(c1)
                n+=1
            
            if((c2>=0) and (c2<self.d[1]) and(int(floor(c2/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[n]= int(c2)
                n+=1
                
        elif i==1:
            c3 = int(floor(idx+self.d_[0]-1))
            c4 = int(floor(idx-self.d_[0]+1))
         
            if ((c3>=0) and (c3<self.d[1]) and (int(floor(c3/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[n] = int(c3)
                n+=1
            
            if((c4>=0) and (c4<self.d[1])  and(int(floor(c4/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[n]= int(c4)
                n+=1
        # cross Neighbors +
        elif i==2:
            c1 = int(floor(idx-self.d_[1]-1))
            c2 = int(floor(idx+self.d_[1]+1))
         
            if ((c1>=0) and (c1<self.d[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n] = int(c1)
                n+=1
            
            if((c2>=0) and (c2<self.d[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n]= int(c2)
                n+=1
       
        elif i==3:
            c1 = int(floor(idx-self.d_[1]+1))
            c2 = int(floor(idx+self.d_[1]-1))
         
            if ((c1>=0) and (c1<self.d[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n] = int(c1)
                n+=1
            
            if((c2>=0) and (c2<self.d[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n]= int(c2)
                n+=1
                
        # cross neighbors +        
        elif i==4:
            c1 = int(floor(idx-self.d_[1]-self.d_[0]))
            c2 = int(floor(idx+self.d_[1]+self.d_[0]))
         
            if ((c1>=0) and (c1<self.d[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n] = int(c1)
                n+=1
            
            if((c2>=0) and (c2<self.d[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n]= int(c2)
                n+=1
                
        elif i==5:
            c1 = int(floor(idx-self.d_[1]+self.d_[0]))
            c2 = int(floor(idx+self.d_[1]-self.d_[0]))
         
            if ((c1>=0) and (c1<self.d[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n] = int(c1)
                n+=1
            
            if((c2>=0) and (c2<self.d[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n]= int(c2)
                n+=1
        #cross neighbors x        
        elif i==6:
            c1 = int(floor(idx-self.d_[1]-self.d_[0]-1))
            c2 = int(floor(idx+self.d_[1]+self.d_[0]+1))
         
            if ((c1>=0) and (c1<self.d[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n] = int(c1)
                n+=1
            
            if((c2>=0) and (c2<self.d[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n]= int(c2)
                n+=1
        
        elif i==7:
            c1 = int(floor(idx+self.d_[1]-self.d_[0]-1))
            c2 = int(floor(idx-self.d_[1]+self.d_[0]+1))
         
            if ((c1>=0) and (c1<self.d[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n] = int(c1)
                n+=1
            
            if((c2>=0) and (c2<self.d[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n]= int(c2)
                n+=1
                
        elif i==8:
            c1 = int(floor(idx-self.d_[1]+self.d_[0]-1))
            c2 = int(floor(idx+self.d_[1]-self.d_[0]+1))
         
            if ((c1>=0) and (c1<self.d[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n] = int(c1)
                n+=1
            
            if((c2>=0) and (c2<self.d[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n]= int(c2)
                n+=1
                
        elif i==9:
            c1 = int(floor(idx+self.d_[1]+self.d_[0]-1))
            c2 = int(floor(idx-self.d_[1]-self.d_[0]+1))
         
            if ((c1>=0) and (c1<self.d[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n] = int(c1)
                n+=1
            
            if((c2>=0) and (c2<self.d[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n]= int(c2)
                n+=1
                
        return n
                    
    def getMinValueInDim(self, idx,dim):
        global n,minneighbors
        n=0
        n =  getNeighborsInDim(idx, minneighbors, dim)
       
        if (n==0) or self.nodesList[minneighbors[0]].cost < self.nodesList[minneighbors[1]].cost:
            return self.nodesList[minneighbors[0]].cost
        else:
            return self.nodesList[minneighbors[1]].cost
    
    def getMinValueInDimcross(self, idx,dim):
        global minneighbors
        n = 0
        
        n=self.getMinNeighborsInDimcross(idx, minneighbors, dim)
       
        if (n==0) or self.nodesList[minneighbors[0]].cost < self.nodesList[minneighbors[1]].cost:
            return self.nodesList[minneighbors[0]].cost
        else:
            return self.nodesList[minneighbors[1]].cost

    def solveEikonalNDims(self, idx,dim):
               
            if dim==1:
                return Tvalues[0]+(1/self.nodesList[idx].v)
            
            sumT = 0
            sumTT = 0    
        
            for i in range(0,dim):
                sumT +=Tvalues[i]
                sumTT +=Tvalues[i]*Tvalues[i]
        
            a=dim
            b = -2*sumT
            c = sumTT-(1*1/pow(self.nodesList[idx].v,2))# *getNode(idx).v)
            q = b*b - 4*a*c
            
            if q<0:
                return INFINI
            else:
                return ((-b+sqrt(q))/(2*a))

    def solveEikonalNDimscross(self, idx,dim):
               
            if dim==1:
                return Tvalues[0]+(sqrt(2)/self.nodesList[idx].v)
            
            sumT = 0
            sumTT = 0    
        
            for i in range(0,dim):
                sumT +=Tvalues[i]
                sumTT +=Tvalues[i]*Tvalues[i]
        
            a=dim
            b = -2*sumT
            c = sumTT-(sqrt(2)*sqrt(2)/pow(self.nodesList[idx].v,2))# *getNode(idx).v)
            q = b*b - 4*a*c
            
            if q<0:
                return INFINI
            else:
                return ((-b+sqrt(q))/(2*a))
              
    def SolveEikonal(self, idx):
        
        global Tvalues       
        a = 3         
        Tvalues = []
        for dim in range(a):
            minT = 0
            minT = getMinValueInDim(idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[idx].cost:
                heappush(Tvalues,minT) 
            
            else:
                a-=1
        
        if a==0:
            return INFINI
        
        Tvalues = sorted(Tvalues)
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDims(idx,i)
            if i==a or updatedT-Tvalues[i]<=0:
                break
        
        return updatedT
                
    def SolveEikonalcross(self, idx):
        
        global Tvalues       
        a = 2#dimension x,y         
        Tvalues = []
       
        for dim in range(a):
            minT = 0
            minT = self.getMinValueInDimcross(idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[idx].cost:
                heappush(Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(idx,i)
            if i==a or updatedT-Tvalues[i]<0:
                break
        
        return updatedT
         
    def SolveEikonalcross1(self, idx):
        
        global Tvalues       
        a = 2#dimension x,y         
        Tvalues = []
       
        for dim in range(2,a+2):
            minT = 0
            minT = self.getMinValueInDimcross(idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[idx].cost:
                heappush(Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(idx,i)
            if i==a or updatedT-Tvalues[i]<0:
                break
        
        return updatedT
            

    def SolveEikonalcross2(self, idx):
        
        global Tvalues       
        a = 2#dimension x,y         
        Tvalues = []
       
        for dim in range(4,4+a):
            minT = 0
            minT = self.getMinValueInDimcross(idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[idx].cost:
                heappush(Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(idx,i)
            if i==a or updatedT-Tvalues[i]<0:
                break
        
        return updatedT
            
    
    
        
    ########################################################################        

    def SolveEikonalcross3(self, idx):
        
        global Tvalues       
        a = 2#dimension x,y         
        Tvalues = []
       
        for dim in range(6,6+a):
            minT = 0
            minT = self.getMinValueInDimcross(idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[idx].cost:
                heappush(Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(idx,i)
            if i==a or updatedT-Tvalues[i]<0:
                break
        
        return updatedT     

    def SolveEikonalcross4(self, idx):
        
        global Tvalues       
        a = 2#dimension x,y         
        Tvalues = []
       
        for dim in range(8,8+a):
            minT = 0
            minT = self.getMinValueInDimcross(idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[idx].cost:
                heappush(Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(idx,i)
            if i==a or updatedT-Tvalues[i]<0:
                break
        
        return updatedT 

    def processMSFMM(self):
        
        for i in self.start_points:
            initCur = self.nodesList[i]
            
            initCur.cost = 0
            heappush(self.NarrowBand,(initCur.cost,i))
       
        while self.NarrowBand != [] and not self.stopWaveProp:#the main loop of FMM
    
            neighbors = [-1,-1,-1,-1,-1,-1]
            neighbors_cross = [-1,-1,-1,-1]
            neighbors_cross1 = [-1,-1,-1,-1]  
            neighbors_cross2 = [-1,-1,-1,-1]  
            neighbors_cross3 = [-1,-1,-1,-1]  
            neighbors_cross4 = [-1,-1,-1,-1]  
            n_neighs,n_neighs_cross,n_neighs_cross1,n_neighs_cross2,n_neighs_cross3,n_neighs_cross4 = 0,0,0,0,0,0
           
            best = heappop(self.NarrowBand)   
            idxMin = best[1]
                    
            n_neighs        = getneighbors(idxMin,neighbors,3)
            n_neighs_cross  = getneighbors(idxMin,neighbors_cross,3)
            n_neighs_cross1 = getneighbors(idxMin,neighbors_cross1,3)
            n_neighs_cross2 = getneighbors(idxMin,neighbors_cross2,3)
            n_neighs_cross3 = getneighbors(idxMin,neighbors_cross3,3)
            n_neighs_cross4 = getneighbors(idxMin,neighbors_cross4,3)
            
            self.nodesList[idxMin].type ='A'
            
           
            #1
            for s in range(n_neighs):#MANHATAN NEIGHBORS
                
                j = neighbors[s]
                curMN = self.nodesList[j]
                if  (j==-1 or curMN.type =='A'  or curMN.TAG==FORBIDDEN or curMN.TAG==NEW_FORBIDDEN):
                    
                    continue
                
                else:
                   
                    cost =SolveEikonal(j)
                    
                    if self.NarrowBand.__contains__((curMN.cost,j)):
                       
                        if cost<curMN.cost:
                            
                            indexe = self.NarrowBand.index((curMN.cost,j))
                            curMN.cost = cost
                          
                            _siftdown(self.NarrowBand,0,indexe)
                    else:
                        
                        curMN.cost = cost
                        curMN.type = 'N'
                        heappush(self.NarrowBand,(cost,j))
                
                #print 'cur ',curMN.type,' ',curMN.cost,' pos',curMN.x,curMN.y,curMN.z,'  Obst',curMN.OBSTACLE 
            #2
            for s in range(n_neighs_cross):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross[s]
                ccur = self.nodesList[j]
                     
                if  (j==-1 or ccur.type =='A'  or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN):
                
                    continue
                
                else:
                   
                    cost =self.SolveEikonalcross(j)
                    
                    if self.NarrowBand.__contains__((ccur.cost,j)):
                       
                        if cost<ccur.cost:
                            
                            indexe = self.NarrowBand.index((ccur.cost,j))
                            ccur.cost = cost
                            _siftdown(self.NarrowBand,0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.type = 'N'
                        heappush(self.NarrowBand,(cost,j))
                
                #print 'ccur ',ccur.type,' ',ccur.cost,' pos',ccur.x,ccur.y,ccur.z,'  Obst',ccur.OBSTACLE 
            #3    
            for s in range(n_neighs_cross1):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross1[s]
                ccur = self.nodesList[j]
                     
                if  (j==-1 or ccur.type =='A'  or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN):
                
                    continue
                
                else:
                   
                    cost =self.SolveEikonalcross1(j)
                    
                    if self.NarrowBand.__contains__((ccur.cost,j)):
                       
                        if cost<ccur.cost:
                            
                            indexe = self.NarrowBand.index((ccur.cost,j))
                            ccur.cost = cost
                            _siftdown(self.NarrowBand,0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.type = 'N'
                        heappush(self.NarrowBand,(cost,j))
                
                #print 'ccur ',ccur.type,' ',ccur.cost,' pos',ccur.x,ccur.y,ccur.z,'  Obst',ccur.OBSTACLE 
            #4
            for s in range(n_neighs_cross2):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross2[s]
                ccur = self.nodesList[j]
                     
                if  (j==-1 or ccur.type =='A'  or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN):
                
                    continue
                
                else:
                   
                    cost =self.SolveEikonalcross2(j)
                    
                    if self.NarrowBand.__contains__((ccur.cost,j)):
                       
                        if cost<ccur.cost:
                            
                            indexe = self.NarrowBand.index((ccur.cost,j))
                            ccur.cost = cost
                            _siftdown(self.NarrowBand,0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.type = 'N'
                        heappush(self.NarrowBand,(cost,j))
                
                #print 'ccur ',ccur.type,' ',ccur.cost,' pos',ccur.x,ccur.y,ccur.z,'  Obst',ccur.OBSTACLE 
            #5            
            for s in range(n_neighs_cross3):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross3[s]
                ccur = self.nodesList[j]
                     
                if  (j==-1 or ccur.type =='A'  or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN):
                
                    continue
                
                else:
                   
                    cost =self.SolveEikonalcross3(j)
                    
                    if self.NarrowBand.__contains__((ccur.cost,j)):
                       
                        if cost<ccur.cost:
                            
                            indexe = self.NarrowBand.index((ccur.cost,j))
                            ccur.cost = cost
                            _siftdown(self.NarrowBand,0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.type = 'N'
                        heappush(self.NarrowBand,(cost,j))
                
                #print 'ccur ',ccur.type,' ',ccur.cost,' pos',ccur.x,ccur.y,ccur.z,'  Obst',ccur.OBSTACLE 
                
            #6           
            for s in range(n_neighs_cross4):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross4[s]
                ccur = self.nodesList[j]
                     
                if  (j==-1 or ccur.type =='A'  or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN):
                
                    continue
                
                else:
                   
                    cost =self.SolveEikonalcross4(j)
                    
                    if self.NarrowBand.__contains__((ccur.cost,j)):
                       
                        if cost<ccur.cost:
                            
                            indexe = self.NarrowBand.index((ccur.cost,j))
                            ccur.cost = cost
                            _siftdown(self.NarrowBand,0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.type = 'N'
                        heappush(self.NarrowBand,(cost,j))
                
                #print 'ccur ',ccur.type,' ',ccur.cost,' pos',ccur.x,ccur.y,ccur.z,'  Obst',ccur.OBSTACLE 
                
            if idxMin == self.goal:
                self.stopWaveProp = True
    
