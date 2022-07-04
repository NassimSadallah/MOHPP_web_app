'''

Created on Jul 8, 2021

@author: nassim
'''
from __builtin__ import False
import time
from math import floor, sqrt
from heapq import _siftdown, heappop, heappush
from utilities import isEmptyList, locateBestIdx, FORBIDDEN, INFINI, NEW_FORBIDDEN, KNOWN, FROZEN_FIX,FROZEN,sqrt_dist,d_, d_parallel, d_parallels
import multiprocessing as mp

c1,c2,n_neighs,neighbors=0,0,0,0
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

class fm2D(object):
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
            heappush(self.narrowBand, (0, i))
       
        '''
        main loop over the node list
        '''
        initime = time.time()
        while self.narrowBand !=[] and not self.stopWave:
            itera+=1
            best = heappop(self.narrowBand)
            idxMin = best[1]
            neighbors = getneighbors(self,idxMin)
            self.nodesList[idxMin].type = FROZEN
            
            for j in neighbors:
                
                if j==-1:
                    continue
                 
                current = self.nodesList[j]
             
                if current.type == FROZEN or current.TAG ==FORBIDDEN or current.TAG == NEW_FORBIDDEN:
                    
                    continue 
               
                cost = SolveEikonal(self,j, current.cros)
                
                if self.narrowBand.__contains__((current.cost,j)):
                    
                    if cost<current.cost:
                        indexe = self.narrowBand.index((current.cost,j))
                        current.cost = cost 
                        _siftdown(self.narrowBand,0,indexe)
                else:            
                    current.cost = cost
                    current.type = KNOWN
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
            self.nodesList[idxMin].type = FROZEN
            
            for j in neighbors:
                
                if j==-1:
                    continue
                 
                current = self.nodesList[j]
             
                if current.type == FROZEN or current.TAG ==FORBIDDEN or current.TAG == NEW_FORBIDDEN:
                    continue 
               
                cost = SolveEikonal(self, j, current.cros)
                
                if blk_nodes[current.block - 1].__contains__((current.cost,j)):
                    
                    if cost<current.cost:
                        indexe = blk_nodes[current.block - 1].index((current.cost,j))
                        current.cost = cost 
                        _siftdown(blk_nodes[current.block - 1],0,indexe)            
                else:            
                    current.cost = cost
                    current.type = KNOWN
                    heappush(blk_nodes[current.block - 1],(cost,j))
            
            EmptyList = isEmptyList(blk_nodes)
            
            if idxMin == self.goal:
                self.stopWave = True                    
        print time.time() - initime, 'seconds'    
        return self.nodesList                

class heuristicFM2D(object):
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
            self.nodesList[i].cost = 0+sqrt_dist((self.nodesList[i].abscice-self.nodesList[self.goal].abscice),
                                                 (self.nodesList[i].colonne-self.nodesList[self.goal].colonne))/self.nodesList[i].v 
            
            #self.nodesList[i].full =self.nodesList[i].cost + (self.nodesList[i].hCost/self.nodesList[i].v -
            #                                              self.nodesList[i].hCost/1.0) #1.0 is the max VELOCITY            
            heappush(self.narrowBand, (self.nodesList[i].cost, i))
        
        '''
        main loop over the node list
        '''
        initime = time.time()
        while self.narrowBand !=[] and not self.stopWave:

            itera+=1
            best = heappop(self.narrowBand)
            idxMin = best[1]
            neighbors = getneighbors(self, idxMin)
            self.nodesList[idxMin].type = FROZEN
            
            for j in neighbors:
                
                if j==-1:
                    continue
                 
                current = self.nodesList[j]
             
                if current.type == FROZEN or current.TAG ==FORBIDDEN or current.TAG == NEW_FORBIDDEN:
                    continue 
               
                cost = SolveEikonal(self, j, current.cros)+sqrt_dist((current.abscice-self.nodesList[self.goal].abscice),
                                            (current.colonne-self.nodesList[self.goal].colonne))/current.v
                
                #hcost = sqrt_dist((current.abscice-self.nodesList[self.goal].abscice),
                #                            (current.colonne-self.nodesList[self.goal].colonne))/current.v
                #fcost = cost + (hcost/current.v - hcost/1.0) #1.0 is the max VELOCITY 
                
                if self.narrowBand.__contains__((current.cost,j)):
                    
                    if cost<current.cost:
                        indexe = self.narrowBand.index((current.full,j))
                        current.cost = cost
                        #current.hCost = hcost
                        #current.full = fcost 
                        _siftdown(self.narrowBand,0,indexe)
                else:            
                    current.cost = cost
                    #current.hCost = hcost
                    #current.full = fcost                    
                    current.type = KNOWN
                    heappush(self.narrowBand,(cost,j))
        
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
            self.nodesList[i].cost = 0+sqrt_dist((self.nodesList[i].abscice-self.nodesList[self.goal].abscice),
                                            (self.nodesList[i].colonne-self.nodesList[self.goal].colonne))/self.nodesList[i].v 
            #self.nodesList[i].hCost = sqrt_dist((self.nodesList[i].abscice-self.nodesList[self.goal].abscice),
            #                                (self.nodesList[i].colonne-self.nodesList[self.goal].colonne))
            #self.nodesList[i].full =self.nodesList[i].cost + (self.nodesList[i].hCost/self.nodesList[i].v -
            #                                              self.nodesList[i].hCost/1.0) #1.0 is the max VELOCITY
               
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
            self.nodesList[idxMin].type = FROZEN
            
            for j in neighbors:
                
                if j==-1:
                    continue
                 
                current = self.nodesList[j]
             
                if current.type == FROZEN or current.TAG ==FORBIDDEN or current.TAG == NEW_FORBIDDEN:
                    continue 
               
                cost = SolveEikonal(self, j, current.cros)+sqrt_dist((current.abscice-self.nodesList[self.goal].abscice),
                                            (current.colonne-self.nodesList[self.goal].colonne))/current.v
                #hcost = sqrt_dist((current.abscice-self.nodesList[self.goal].abscice),
                #                            (current.colonne-self.nodesList[self.goal].colonne))
                #fcost = cost + (hcost/current.v - hcost/1.0) #1.0 is the max VELOCITY                 
                
                if blk_nodes[current.block - 1].__contains__((current.cost,j)):
                    
                    if cost<current.cost:
                        indexe = blk_nodes[current.block - 1].index((current.cost,j))
                        current.cost = cost 
                        #current.hCost = hcost
                        #current.full = fcost                         
                        _siftdown(blk_nodes[current.block - 1],0,indexe)            
                else:            
                    current.cost = cost
                    #current.hCost = hcost
                    #current.full = fcost                     
                    current.type = KNOWN
                    heappush(blk_nodes[current.block - 1],(cost,j))
            
            EmptyList = isEmptyList(blk_nodes)
            
            if idxMin == self.goal:
                self.stopWave = True                    
        print time.time() - initime, 'seconds'    
        return self.nodesList                

'''
3 dimension version Multi Stencils Fm
'''
class MSfm3D(object):

    def __init__(self,start_points=[], goal=None, nodesList=[], d_ =[], cross=False, heuristicAct = False, seq = 0, block = 0):
        
        self.d_ = d_
        self.NarrowBand = []
        self.stopWaveProp = False
        self.nbr_neighbors = 0
        self.Tvalue = []
        self.start_points = start_points
        self.goal = goal
        self.nodesList = nodesList
        self.cross = cross
        self.block = block 
        self.heuristicAct = heuristicAct
        self.Hcost = sqrt_dist
            
        if self.cross:
            self.neighborsList = [-1 for _ in range(26)]
        else:
            self.neighborsList = [-1 for _ in range(2**(len(d_)-1)+2)]
        self.minneighborsList = [-1 for _ in range(2)]
    
    def getneighbors(self, idx,ne,nume,dims):

        self.nbr_neighbors = 0
        if nume == 0:#standard neighbors
            for i in range(dims): #x,y,z dimension
                self.nbr_neighbors =  self.getNeighborsInDim(idx,ne,0,i)    
            return self.nbr_neighbors   
     
        elif nume ==1 :            
            # cross neighbors 1
            n_neighs_cross=self.getNeighborsInDimcross(idx,1,ne)
            return n_neighs_cross
        
        elif nume ==2 :            
            # cross neighbors 2
            n_neighs_cross1=self.getNeighborsInDimcross(idx,2,ne)
            return n_neighs_cross1
            
        elif nume ==3:            
            # cross neighbors 3
            n_neighs_cross2=self.getNeighborsInDimcross(idx,3,ne)
            return n_neighs_cross2
        
        elif nume ==4:            
            # cross neighbors 4
            n_neighs_cross3=self.getNeighborsInDimcross(idx,4,ne)
            return n_neighs_cross3
        
        elif nume ==5 :            
            # cross neighbors 5
            n_neighs_cross4=self.getNeighborsInDimcross(idx,5,ne)
            return n_neighs_cross4
       
    def getNeighborsInDim(self, idx,ne,nume,i):
          
        c1,c2= 0,0
     
        if i==0:
            c1 = idx-1
            c2 = idx+1
          
        else:
            c1 = int(floor(idx-self.d_[i-1]))
            c2 = int(floor(idx+self.d_[i-1]))
            
        
        if nume==0:#neighbors
            
            if ((c1>=0) and (int(floor(c1/self.d_[i]))==int(floor(idx/self.d_[i])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (int(floor(c2/self.d_[i]))==int(floor(idx/self.d_[i])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
            print self.nbr_neighbors
            return self.nbr_neighbors
        
        else:#min_neighbors
            
            
            if ((c1>=0) and (int(floor(c1/self.d_[i]))==int(floor(idx/self.d_[i])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1

            if((c2>=0) and(int(floor(c2/self.d_[i]))==int(floor(idx/self.d_[i])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
            print self.nbr_neighbors
            return self.nbr_neighbors

    def getNeighborsInDimcross(self, idx,nbr,nec):
        
        n_neighs_cross,n_neighs_cross1,n_neighs_cross2,n_neighs_cross3,n_neighs_cross4 = 0,0,0,0,0
        
        
        if nbr==1:
        
            cc1 = int(floor(idx+self.d_[0]+1))
            cc2 = int(floor(idx-self.d_[0]-1))
            cc3 = int(floor(idx+self.d_[0]-1))
            cc4 = int(floor(idx-self.d_[0]+1))
            
            if ((cc1>=0) and (cc1<self.d_[1]) and (int(floor(cc1/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                nec[n_neighs_cross] = int(cc1)       
                n_neighs_cross+=1
            
            if((cc2>=0) and (cc2<self.d_[1]) and (int(floor(cc2/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                nec[n_neighs_cross]= int(cc2)    
                n_neighs_cross+=1
            
            if((cc3>=0) and (cc3<self.d_[1]) and (int(floor(cc3/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                nec[n_neighs_cross]= int(cc3)  
                n_neighs_cross+=1
            
            if((cc4>=0) and (cc4<self.d_[1]) and (int(floor(cc4/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                nec[n_neighs_cross]= int(cc4)
                n_neighs_cross+=1

            return n_neighs_cross
        
        elif nbr==2:
            
            tcc1 = int(floor(idx+self.d_[1]-1))
            tcc2 = int(floor(idx-self.d_[1]+1))
            tcc3 = int(floor(idx-self.d_[1]-1))
            tcc4 = int(floor(idx+self.d_[1]+1))
            #Top diagonals
            if ((tcc1>=0) and (tcc1<self.d_[2])and (int(floor(tcc1/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross1] = int(tcc1)       
                n_neighs_cross1+=1
            
            if((tcc2>=0) and (tcc2<self.d_[2]) and (int(floor(tcc2/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross1]= int(tcc2)    
                n_neighs_cross1+=1
            
            if((tcc3>=0) and (tcc3<self.d_[2]) and (int(floor(tcc3/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross1]= int(tcc3)  
                n_neighs_cross1+=1
            
            if((tcc4>=0) and (tcc4<self.d_[2]) and (int(floor(tcc4/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross1]= int(tcc4)
                n_neighs_cross1+=1

                  
            return n_neighs_cross1

        elif nbr==3:        
            
            tcc5 = int(floor(idx+self.d_[1]-self.d_[0]))
            tcc6 = int(floor(idx-self.d_[1]+self.d_[0]))
            tcc7 = int(floor(idx-self.d_[1]-self.d_[0]))
            tcc8 = int(floor(idx+self.d_[1]+self.d_[0])) 
                                   
            if ((tcc5>=0) and (tcc5<self.d_[2])and (int(floor(tcc5/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross2] = int(tcc5)       
                n_neighs_cross2+=1
            
            if((tcc6>=0) and (tcc6<self.d_[2]) and (int(floor(tcc6/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross2]= int(tcc6)    
                n_neighs_cross2+=1
            
            if((tcc7>=0) and (tcc7<self.d_[2]) and (int(floor(tcc7/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross2]= int(tcc7)  
                n_neighs_cross2+=1
            
            if((tcc8>=0) and (tcc8<self.d_[2]) and (int(floor(tcc8/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross2]= int(tcc8)
                n_neighs_cross2+=1
        
            return n_neighs_cross2

        elif nbr==4:
        
            bcc1 = int(floor(idx-self.d_[1]+self.d_[0]-1))
            bcc2 = int(floor(idx+self.d_[1]-self.d_[0]+1))
            bcc3 = int(floor(idx+self.d_[1]+self.d_[0]-1))
            bcc4 = int(floor(idx-self.d_[1]-self.d_[0]+1))   
            #Bottom diagonals   
            if ((bcc1>=0) and (bcc1<self.d_[2])and (int(floor(bcc1/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross3] = int(bcc1)       
                n_neighs_cross3+=1
            
            if((bcc2>=0) and (bcc2<self.d_[2]) and (int(floor(bcc2/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross3]= int(bcc2)    
                n_neighs_cross3+=1
            
            if((bcc3>=0) and (bcc3<self.d_[2]) and (int(floor(bcc3/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross3]= int(bcc3)  
                n_neighs_cross3+=1
            
            if((bcc4>=0) and (bcc4<self.d_[2]) and (int(floor(bcc4/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross3]= int(bcc4)
                n_neighs_cross3+=1
            
            return n_neighs_cross3
    
        elif nbr==5:
            
            bcc5 = int(floor(idx+self.d_[1]-self.d_[0]-1))
            bcc6 = int(floor(idx-self.d_[1]+self.d_[0]+1))
            bcc7 = int(floor(idx-self.d_[1]-self.d_[0]-1))
            bcc8 = int(floor(idx+self.d_[1]+self.d_[0]+1))
            
            if ((bcc5>=0) and (bcc5<self.d_[2])and (int(floor(bcc5/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross4] = int(bcc5)       
                n_neighs_cross4+=1
            
            if((bcc6>=0) and (bcc6<self.d_[2]) and (int(floor(bcc6/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross4]= int(bcc6)    
                n_neighs_cross4+=1
            
            if((bcc7>=0) and (bcc7<self.d_[2]) and (int(floor(bcc7/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[n_neighs_cross4]= int(bcc7)  
                n_neighs_cross4+=1
            
            if((bcc8>=0) and (bcc8<self.d_[2]) and (int(floor(bcc8/self.d_[2]))==int(floor((idx)/self.d_[2])))):
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
        
            if ((c1>=0) and (c1<self.d_[1]) and (int(floor(c1/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[n] = int(c1)
                n+=1
            
            if((c2>=0) and (c2<self.d_[1]) and(int(floor(c2/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[n]= int(c2)
                n+=1
                
        elif i==1:
            c3 = int(floor(idx+self.d_[0]-1))
            c4 = int(floor(idx-self.d_[0]+1))
         
            if ((c3>=0) and (c3<self.d_[1]) and (int(floor(c3/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[n] = int(c3)
                n+=1
            
            if((c4>=0) and (c4<self.d_[1])  and(int(floor(c4/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[n]= int(c4)
                n+=1
        # cross Neighbors +
        elif i==2:
            c1 = int(floor(idx-self.d_[1]-1))
            c2 = int(floor(idx+self.d_[1]+1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n] = int(c1)
                n+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n]= int(c2)
                n+=1
       
        elif i==3:
            c1 = int(floor(idx-self.d_[1]+1))
            c2 = int(floor(idx+self.d_[1]-1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n] = int(c1)
                n+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n]= int(c2)
                n+=1
                
        # cross neighbors +        
        elif i==4:
            c1 = int(floor(idx-self.d_[1]-self.d_[0]))
            c2 = int(floor(idx+self.d_[1]+self.d_[0]))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n] = int(c1)
                n+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n]= int(c2)
                n+=1
                
        elif i==5:
            c1 = int(floor(idx-self.d_[1]+self.d_[0]))
            c2 = int(floor(idx+self.d_[1]-self.d_[0]))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n] = int(c1)
                n+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n]= int(c2)
                n+=1
        #cross neighbors x        
        elif i==6:
            c1 = int(floor(idx-self.d_[1]-self.d_[0]-1))
            c2 = int(floor(idx+self.d_[1]+self.d_[0]+1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n] = int(c1)
                n+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n]= int(c2)
                n+=1
        
        elif i==7:
            c1 = int(floor(idx+self.d_[1]-self.d_[0]-1))
            c2 = int(floor(idx-self.d_[1]+self.d_[0]+1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n] = int(c1)
                n+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n]= int(c2)
                n+=1
                
        elif i==8:
            c1 = int(floor(idx-self.d_[1]+self.d_[0]-1))
            c2 = int(floor(idx+self.d_[1]-self.d_[0]+1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n] = int(c1)
                n+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n]= int(c2)
                n+=1
                
        elif i==9:
            c1 = int(floor(idx+self.d_[1]+self.d_[0]-1))
            c2 = int(floor(idx-self.d_[1]-self.d_[0]+1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n] = int(c1)
                n+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[n]= int(c2)
                n+=1
                
        return n
                    
    def getMinValueInDim(self, idx,dim):
        #global n,minneighbors
        n=0
        n =  self.getNeighborsInDim(idx, self.minneighborsList,1, dim)
       
        if (n==0) or self.nodesList[self.minneighborsList[0]].cost < self.nodesList[self.minneighborsList[1]].cost:
            return self.nodesList[self.minneighborsList[0]].cost
        else:
            return self.nodesList[self.minneighborsList[1]].cost
    
    def getMinValueInDimcross(self, idx,dim):
        global minneighbors
        n = 0
        
        n=self.getMinNeighborsInDimcross(idx, self.minneighborsList, dim)
       
        if (n==0) or self.nodesList[self.minneighborsList[0]].cost < self.nodesList[self.minneighborsList[1]].cost:
            return self.nodesList[self.minneighborsList[0]].cost
        else:
            return self.nodesList[self.minneighborsList[1]].cost

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
        
        #global Tvalues       
        a = 3         
        self.Tvalues = []
        for dim in range(a):
            minT = 0
            minT = self.getMinValueInDim(idx,dim)
                       
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

    def processMSFM3D(self):
        
        if self.goal !=None:
            g = self.nodesList[self.goal]
        
        print int(self.heuristicAct)    
        
        for i in self.start_points:
            
            initCur = self.nodesList[i]
            initCur.cost = 0+self.Hcost((initCur.x-g.x),(initCur.y-g.y),(initCur.z-g.z))/initCur.v
            heappush(self.NarrowBand,(initCur.cost,i))
        initime = time.time()
        ite = 0
        while self.NarrowBand != [] and not self.stopWaveProp:#the main loop of FMM
            print('iteration ',ite)
            neighbors = [-1,-1,-1,-1,-1,-1]
            neighbors_cross = [-1,-1,-1,-1]
            neighbors_cross1 = [-1,-1,-1,-1]  
            neighbors_cross2 = [-1,-1,-1,-1]  
            neighbors_cross3 = [-1,-1,-1,-1]  
            neighbors_cross4 = [-1,-1,-1,-1]  
            n_neighs,n_neighs_cross,n_neighs_cross1,n_neighs_cross2,n_neighs_cross3,n_neighs_cross4 = 0,0,0,0,0,0
           
            best = heappop(self.NarrowBand)   
            idxMin = best[1]
                    
            n_neighs        = self.getneighbors(idxMin,neighbors,0,3)
            n_neighs_cross  = self.getneighbors(idxMin,neighbors_cross,1,3)
            n_neighs_cross1 = self.getneighbors(idxMin,neighbors_cross1,2,3)
            n_neighs_cross2 = self.getneighbors(idxMin,neighbors_cross2,3,3)
            n_neighs_cross3 = self.getneighbors(idxMin,neighbors_cross3,4,3)
            n_neighs_cross4 = self.getneighbors(idxMin,neighbors_cross4,5,3)
            
            self.nodesList[idxMin].type =FROZEN
            
           
            #1
            for s in range(n_neighs):#MANHATAN NEIGHBORS
                
                j = neighbors[s]
                curMN = self.nodesList[j]
                if  (j==-1 or curMN.type ==FROZEN  or curMN.TAG==FORBIDDEN or curMN.TAG==NEW_FORBIDDEN):
                    
                    continue
                
                else:
                   
                    cost =self.SolveEikonal(j)+self.Hcost((curMN.x-g.x),(curMN.y-g.y),(curMN.z-g.z))/curMN.v
                    
                    if self.NarrowBand.__contains__((curMN.cost,j)):
                       
                        if cost<curMN.cost:
                            
                            indexe = self.NarrowBand.index((curMN.cost,j))
                            curMN.cost = cost
                          
                            _siftdown(self.NarrowBand,0,indexe)
                    else:
                        
                        curMN.cost = cost
                        curMN.type = KNOWN
                        heappush(self.NarrowBand,(cost,j))
                
                #print 'cur ',curMN.type,' ',curMN.cost,' pos',curMN.x,curMN.y,curMN.z,'  Obst',curMN.OBSTACLE 
            #2
            for s in range(n_neighs_cross):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross[s]
                ccur = self.nodesList[j]
                     
                if  (j==-1 or ccur.type ==FROZEN  or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN):
                
                    continue
                
                else:
                   
                    cost =self.SolveEikonalcross(j)+self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))#/ccur.v
                    
                    if self.NarrowBand.__contains__((ccur.cost,j)):
                       
                        if cost<ccur.cost:
                            
                            indexe = self.NarrowBand.index((ccur.cost,j))
                            ccur.cost = cost
                            _siftdown(self.NarrowBand,0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.type = KNOWN
                        heappush(self.NarrowBand,(cost,j))
                
                #print 'ccur ',ccur.type,' ',ccur.cost,' pos',ccur.x,ccur.y,ccur.z,'  Obst',ccur.OBSTACLE 
            #3    
            for s in range(n_neighs_cross1):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross1[s]
                ccur = self.nodesList[j]
                     
                if  (j==-1 or ccur.type ==FROZEN  or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN):
                
                    continue
                
                else:
                   
                    cost =self.SolveEikonalcross1(j)+self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))/ccur.v
                    
                    if self.NarrowBand.__contains__((ccur.cost,j)):
                       
                        if cost<ccur.cost:
                            
                            indexe = self.NarrowBand.index((ccur.cost,j))
                            ccur.cost = cost
                            _siftdown(self.NarrowBand,0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.type = KNOWN
                        heappush(self.NarrowBand,(cost,j))
                
                #print 'ccur ',ccur.type,' ',ccur.cost,' pos',ccur.x,ccur.y,ccur.z,'  Obst',ccur.OBSTACLE 
            #4
            for s in range(n_neighs_cross2):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross2[s]
                ccur = self.nodesList[j]
                     
                if  (j==-1 or ccur.type ==FROZEN  or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN):
                
                    continue
                
                else:
                   
                    cost =self.SolveEikonalcross2(j)+self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))/ccur.v
                    
                    if self.NarrowBand.__contains__((ccur.cost,j)):
                       
                        if cost<ccur.cost:
                            
                            indexe = self.NarrowBand.index((ccur.cost,j))
                            ccur.cost = cost
                            _siftdown(self.NarrowBand,0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.type = KNOWN
                        heappush(self.NarrowBand,(cost,j))
                
                #print 'ccur ',ccur.type,' ',ccur.cost,' pos',ccur.x,ccur.y,ccur.z,'  Obst',ccur.OBSTACLE 
            #5            
            for s in range(n_neighs_cross3):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross3[s]
                ccur = self.nodesList[j]
                     
                if  (j==-1 or ccur.type ==FROZEN  or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN):
                
                    continue
                
                else:
                   
                    cost =self.SolveEikonalcross3(j)+self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))/ccur.v
                    
                    if self.NarrowBand.__contains__((ccur.cost,j)):
                       
                        if cost<ccur.cost:
                            
                            indexe = self.NarrowBand.index((ccur.cost,j))
                            ccur.cost = cost
                            _siftdown(self.NarrowBand,0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.type = KNOWN
                        heappush(self.NarrowBand,(cost,j))
                
                #print 'ccur ',ccur.type,' ',ccur.cost,' pos',ccur.x,ccur.y,ccur.z,'  Obst',ccur.OBSTACLE 
                
            #6           
            for s in range(n_neighs_cross4):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross4[s]
                ccur = self.nodesList[j]
                     
                if  (j==-1 or ccur.type ==FROZEN  or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN):
                
                    continue
                
                else:
                   
                    cost =self.SolveEikonalcross4(j)+self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))/ccur.v
                    
                    if self.NarrowBand.__contains__((ccur.cost,j)):
                       
                        if cost<ccur.cost:
                            
                            indexe = self.NarrowBand.index((ccur.cost,j))
                            ccur.cost = cost
                            _siftdown(self.NarrowBand,0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.type = KNOWN
                        heappush(self.NarrowBand,(cost,j))
                
                #print 'ccur ',ccur.type,' ',ccur.cost,' pos',ccur.x,ccur.y,ccur.z,'  Obst',ccur.OBSTACLE 
            ite+=1    
            if idxMin == self.goal:
                self.stopWaveProp = True
        print time.time() - initime, 'seconds'
        return self.nodesList

class VBLHGS3D(object):

    def __init__(self, start_points=[], goal=None, nodesList=[], d_ =[], cross=False, heuristicAct = False, seq = 0, block = 0):

        self.block = block 
        print(self.block)
        self.d_ = d_
        self.NarrowList = [[] for _ in range(block)]
        self.stopWaveProp = False
        self.Tvalue = []
        self.start_points = start_points
        self.goal = goal
        self.nodesList = nodesList
        self.nbr_neighbors = 0
        self.cross = cross
        
        self.heuristicAct = heuristicAct
        self.Hcost = sqrt_dist
            
        if self.cross:
            self.neighborsList = [-1 for _ in range(26)]
        else:
            self.neighborsList = [-1 for _ in range(2**(len(d_)-1)+2)]
        self.minneighborsList = [-1 for _ in range(2)]
        self.processMSFM3D()

    def getneighbors(self, idx,ne,nume,dims):
        
        self.nbr_neighbors = 0
        if nume == 0:#standard neighbors
            #for i in range(dims): #x,y,z dimension
            self.nbr_neighbors =  self.getNeighborsInDim(idx,ne,0,-1)
                    
            return self.nbr_neighbors   
     
        elif nume ==1 :            
            # cross neighbors 1
            self.nbr_neighbors=self.getNeighborsInDimcross(idx,1,ne)
            return self.nbr_neighbors
        
        elif nume ==2 :            
            # cross neighbors 2
            self.nbr_neighbors=self.getNeighborsInDimcross(idx,2,ne)
            return self.nbr_neighbors
            
        elif nume ==3:            
            # cross neighbors 3
            self.nbr_neighbors=self.getNeighborsInDimcross(idx,3,ne)
            return self.nbr_neighbors
        
        elif nume ==4:            
            # cross neighbors 4
            self.nbr_neighbors=self.getNeighborsInDimcross(idx,4,ne)
            return self.nbr_neighbors
        
        elif nume ==5 :            
            # cross neighbors 5
            self.nbr_neighbors=self.getNeighborsInDimcross(idx,5,ne)
            return self.nbr_neighbors
       
    def getNeighborsInDim(self, idx,ne,nume, i):
          
        c1,c2,c3,c4,c5,c6 = -1,-1,-1,-1,-1,-1
     
        if nume==0:#neighbors

            c1 = idx-1
            c2 = idx+1
            c3 = int(floor(idx-self.d_[0]))
            c4 = int(floor(idx+self.d_[0]))
            c5= int(floor(idx-self.d_[1]))
            c6= int(floor(idx+self.d_[1]))   

            if ((c1>=0) and (int(floor(c1/self.d_[0]))==int(floor(idx/self.d_[0])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (int(floor(c2/self.d_[0]))==int(floor(idx/self.d_[0])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
            
            if ((c3>=0) and (int(floor(c3/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[self.nbr_neighbors] = int(c3)
                self.nbr_neighbors+=1
            
            if((c4>=0) and (int(floor(c4/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[self.nbr_neighbors]= int(c4)
                self.nbr_neighbors+=1
                
            if ((c5>=0) and (int(floor(c5/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c5)
                self.nbr_neighbors+=1
            
            if((c6>=0) and (int(floor(c6/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c6)
                self.nbr_neighbors+=1
                            
            return self.nbr_neighbors
        
        else:#min_neighbors
    
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
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1

            if((c2>=0) and(int(floor(c2/self.d_[i]))==int(floor(idx/self.d_[i])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
            
            return self.nbr_neighbors

    def getNeighborsInDimcross(self, idx,nbr,nec):
        cc1,cc2,cc3,cc4 = -1,-1,-1,-1
        if nbr==1:
        
            cc1 = int(floor(idx+self.d_[0]+1))
            cc2 = int(floor(idx-self.d_[0]-1))
            cc3 = int(floor(idx+self.d_[0]-1))
            cc4 = int(floor(idx-self.d_[0]+1))
            
            if ((cc1>=0) and (cc1<self.d_[1]) and (int(floor(cc1/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                nec[self.nbr_neighbors] = int(cc1)       
                self.nbr_neighbors+=1
            
            if((cc2>=0) and (cc2<self.d_[1]) and (int(floor(cc2/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                nec[self.nbr_neighbors]= int(cc2)    
                self.nbr_neighbors+=1
            
            if((cc3>=0) and (cc3<self.d_[1]) and (int(floor(cc3/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                nec[self.nbr_neighbors]= int(cc3)  
                self.nbr_neighbors+=1
            
            if((cc4>=0) and (cc4<self.d_[1]) and (int(floor(cc4/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                nec[self.nbr_neighbors]= int(cc4)
                self.nbr_neighbors+=1

            return self.nbr_neighbors
        
        elif nbr==2:
            tcc1,tcc2,tcc3,tcc4 = -1,-1,-1,-1
            
            tcc1 = int(floor(idx+self.d_[1]-1))
            tcc2 = int(floor(idx-self.d_[1]+1))
            tcc3 = int(floor(idx-self.d_[1]-1))
            tcc4 = int(floor(idx+self.d_[1]+1))
            #Top diagonals
            if ((tcc1>=0) and (tcc1<self.d_[2])and (int(floor(tcc1/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors] = int(tcc1)       
                self.nbr_neighbors+=1
            
            if((tcc2>=0) and (tcc2<self.d_[2]) and (int(floor(tcc2/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(tcc2)    
                self.nbr_neighbors+=1
            
            if((tcc3>=0) and (tcc3<self.d_[2]) and (int(floor(tcc3/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(tcc3)  
                self.nbr_neighbors+=1
            
            if((tcc4>=0) and (tcc4<self.d_[2]) and (int(floor(tcc4/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(tcc4)
                self.nbr_neighbors+=1

                  
            return self.nbr_neighbors

        elif nbr==3:        
            tcc5,tcc6,tcc7,tcc8 = -1,-1,-1,-1
            tcc5 = int(floor(idx+self.d_[1]-self.d_[0]))
            tcc6 = int(floor(idx-self.d_[1]+self.d_[0]))
            tcc7 = int(floor(idx-self.d_[1]-self.d_[0]))
            tcc8 = int(floor(idx+self.d_[1]+self.d_[0])) 
                                   
            if ((tcc5>=0) and (tcc5<self.d_[2])and (int(floor(tcc5/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors] = int(tcc5)       
                self.nbr_neighbors+=1
            
            if((tcc6>=0) and (tcc6<self.d_[2]) and (int(floor(tcc6/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(tcc6)    
                self.nbr_neighbors+=1
            
            if((tcc7>=0) and (tcc7<self.d_[2]) and (int(floor(tcc7/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(tcc7)  
                self.nbr_neighbors+=1
            
            if((tcc8>=0) and (tcc8<self.d_[2]) and (int(floor(tcc8/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(tcc8)
                self.nbr_neighbors+=1
        
            return self.nbr_neighbors

        elif nbr==4:
            bcc1,bcc2,bcc3,bcc4 = -1,-1,-1,-1
            bcc1 = int(floor(idx-self.d_[1]+self.d_[0]-1))
            bcc2 = int(floor(idx+self.d_[1]-self.d_[0]+1))
            bcc3 = int(floor(idx+self.d_[1]+self.d_[0]-1))
            bcc4 = int(floor(idx-self.d_[1]-self.d_[0]+1))   
            #Bottom diagonals   
            if ((bcc1>=0) and (bcc1<self.d_[2])and (int(floor(bcc1/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors] = int(bcc1)       
                self.nbr_neighbors+=1
            
            if((bcc2>=0) and (bcc2<self.d_[2]) and (int(floor(bcc2/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(bcc2)    
                self.nbr_neighbors+=1
            
            if((bcc3>=0) and (bcc3<self.d_[2]) and (int(floor(bcc3/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(bcc3)  
                self.nbr_neighbors+=1
            
            if((bcc4>=0) and (bcc4<self.d_[2]) and (int(floor(bcc4/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(bcc4)
                self.nbr_neighbors+=1
            
            return self.nbr_neighbors
    
        elif nbr==5:
            bcc5,bcc6,bcc7,bcc8 = -1,-1,-1,-1
            bcc5 = int(floor(idx+self.d_[1]-self.d_[0]-1))
            bcc6 = int(floor(idx-self.d_[1]+self.d_[0]+1))
            bcc7 = int(floor(idx-self.d_[1]-self.d_[0]-1))
            bcc8 = int(floor(idx+self.d_[1]+self.d_[0]+1))
            
            if ((bcc5>=0) and (bcc5<self.d_[2])and (int(floor(bcc5/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors] = int(bcc5)       
                self.nbr_neighbors+=1
            
            if((bcc6>=0) and (bcc6<self.d_[2]) and (int(floor(bcc6/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(bcc6)    
                self.nbr_neighbors+=1
            
            if((bcc7>=0) and (bcc7<self.d_[2]) and (int(floor(bcc7/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(bcc7)  
                self.nbr_neighbors+=1
            
            if((bcc8>=0) and (bcc8<self.d_[2]) and (int(floor(bcc8/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(bcc8)
                self.nbr_neighbors+=1
            
            return self.nbr_neighbors
    
    def getMinNeighborsInDimcross(self, idx, ne, i):    
        
        #global c1,c2,c3,c4
        c1,c2,c3,c4=-1,-1,-1,-1
        self.nbr_neighbors = 0
        if i==0:
            c1 = int(floor(idx+self.d_[0]+1))
            c2 = int(floor(idx-self.d_[0]-1))
        
            if ((c1>=0) and (c1<self.d_[1]) and (int(floor(c1/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[1]) and(int(floor(c2/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
                
        elif i==1:
            c3 = int(floor(idx+self.d_[0]-1))
            c4 = int(floor(idx-self.d_[0]+1))
         
            if ((c3>=0) and (c3<self.d_[1]) and (int(floor(c3/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[self.nbr_neighbors] = int(c3)
                self.nbr_neighbors+=1
            
            if((c4>=0) and (c4<self.d_[1])  and(int(floor(c4/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[self.nbr_neighbors]= int(c4)
                self.nbr_neighbors+=1
        # cross Neighbors +
        elif i==2:
            c1 = int(floor(idx-self.d_[1]-1))
            c2 = int(floor(idx+self.d_[1]+1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
       
        elif i==3:
            c1 = int(floor(idx-self.d_[1]+1))
            c2 = int(floor(idx+self.d_[1]-1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
                
        # cross neighbors +        
        elif i==4:
            c1 = int(floor(idx-self.d_[1]-self.d_[0]))
            c2 = int(floor(idx+self.d_[1]+self.d_[0]))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
                
        elif i==5:
            c1 = int(floor(idx-self.d_[1]+self.d_[0]))
            c2 = int(floor(idx+self.d_[1]-self.d_[0]))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
        #cross neighbors x        
        elif i==6:
            c1 = int(floor(idx-self.d_[1]-self.d_[0]-1))
            c2 = int(floor(idx+self.d_[1]+self.d_[0]+1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
        
        elif i==7:
            c1 = int(floor(idx+self.d_[1]-self.d_[0]-1))
            c2 = int(floor(idx-self.d_[1]+self.d_[0]+1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
                
        elif i==8:
            c1 = int(floor(idx-self.d_[1]+self.d_[0]-1))
            c2 = int(floor(idx+self.d_[1]-self.d_[0]+1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
                
        elif i==9:
            c1 = int(floor(idx+self.d_[1]+self.d_[0]-1))
            c2 = int(floor(idx-self.d_[1]-self.d_[0]+1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
                
        return self.nbr_neighbors
                    
    def getMinValueInDim(self, idx,dim):
        
        self.nbr_neighbors = 0        
        self.nbr_neighbors =  self.getNeighborsInDim(idx, self.minneighborsList,1, dim)
       
        if (self.nbr_neighbors==0) or self.nodesList[self.minneighborsList[0]].cost < self.nodesList[self.minneighborsList[1]].cost:
            return self.nodesList[self.minneighborsList[0]].cost
        else:
            return self.nodesList[self.minneighborsList[1]].cost
    
    def getMinValueInDimcross(self, idx,dim):
        #global minneighbors
        self.nbr_neighbors = 0
        
        self.nbr_neighbors=self.getMinNeighborsInDimcross(idx, self.minneighborsList, dim)
       
        if (self.nbr_neighbors==0) or self.nodesList[self.minneighborsList[0]].cost < self.nodesList[self.minneighborsList[1]].cost:
            return self.nodesList[self.minneighborsList[0]].cost
        else:
            return self.nodesList[self.minneighborsList[1]].cost

    def solveEikonalNDims(self, idx,dim):
            #global Tvalues   
            if dim==1:
                return self.Tvalues[0]+(1/self.nodesList[idx].v)
            
            sumT = 0
            sumTT = 0    
        
            for i in range(0,dim):
                sumT +=self.Tvalues[i]
                sumTT +=self.Tvalues[i]*self.Tvalues[i]
        
            a=dim
            b = -2*sumT
            c = sumTT-(1*1/pow(self.nodesList[idx].v,2))# *getNode(idx).v)
            q = b*b - 4*a*c
            
            if q<0:
                return INFINI
            else:
                return ((-b+sqrt(q))/(2*a))

    def solveEikonalNDimscross(self, curnode,dim):
               
            if dim==1:
              
                return self.Tvalues[0]+(sqrt(2)/curnode.v)
            
            sumT = 0
            sumTT = 0    
        
            for i in range(0,dim):
                sumT +=self.Tvalues[i]
                sumTT +=self.Tvalues[i]*self.Tvalues[i]
        
            a=dim
            b = -2*sumT
            c = sumTT-(sqrt(2)*sqrt(2)/pow(curnode.v,2))# *getNode(idx).v)
            q = b*b - 4*a*c
            
            if q<0:
                return INFINI
            else:
                return ((-b+sqrt(q))/(2*a))
              
    def SolveEikonal(self, idx):
     
        a = 3         
        self.Tvalues = []
        for dim in range(a):
            minT = 0
            minT = self.getMinValueInDim(idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
        
        if a==0:
            return INFINI
        
        self.Tvalues = sorted(self.Tvalues)
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDims(idx,i)
            if i==a or updatedT-self.Tvalues[i]<=0:
                break
        
        return updatedT
                
    def SolveEikonalcross(self, idx):
             
        a = 2#dimension x,y         
        self.Tvalues = []
        curn = self.nodesList[idx]
        for dim in range(a):
            minT = 0
            minT = self.getMinValueInDimcross(idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(curn,i)
            if i==a or updatedT-self.Tvalues[i]<0:
                break
        
        return updatedT
         
    def SolveEikonalcross1(self, idx):
        
        #global Tvalues       
        a = 2#dimension x,y         
        self.Tvalues = []
        curn = self.nodesList[idx]
       
        for dim in range(2,a+2):
            minT = 0
            minT = self.getMinValueInDimcross(idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(curn,i)
            if i==a or updatedT-self.Tvalues[i]<0:
                break
        
        return updatedT
            
    def SolveEikonalcross2(self, idx):
        
        #global Tvalues       
        a = 2#dimension x,y    
        curn = self.nodesList[idx]     
        self.Tvalues = []
       
        for dim in range(4,4+a):
            minT = 0
            minT = self.getMinValueInDimcross(idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(curn,i)
            if i==a or updatedT-self.Tvalues[i]<0:
                break
        
        return updatedT
            
    
    
        
    ########################################################################        

    def SolveEikonalcross3(self, idx):
        
        #global Tvalues       
        a = 2#dimension x,y 
        curn = self.nodesList[idx]        
        self.Tvalues = []
       
        for dim in range(6,6+a):
            minT = 0
            minT = self.getMinValueInDimcross(idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(curn,i)
            if i==a or updatedT-self.Tvalues[i]<0:
                break
        
        return updatedT     

    def SolveEikonalcross4(self, idx):
        
        #global Tvalues       
        curn = self.nodesList[idx]
        a = 2#dimension x,y         
        self.Tvalues = []
       
        for dim in range(8,8+a):
            minT = 0
            minT = self.getMinValueInDimcross(idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(curn,i)
            if i==a or updatedT-self.Tvalues[i]<0:
                break
        
        return updatedT 

    def processMSFM3D(self):
        print("VBLHGS3D begins, PLEASE Wait! ...")
        if self.goal !=None:
            g = self.nodesList[self.goal]

        h=0
        for i in self.start_points:

            initCur = self.nodesList[i]
            initCur.cost = 0
            if self.heuristicAct:        
                h = self.Hcost((initCur.x-g.x),(initCur.y-g.y),(initCur.z-g.z))#/initCur.v
            
            initCur.type = FROZEN
            initCur.full = 0+h
            heappush(self.NarrowList[initCur.block-1],(initCur.full,i))
            
        initime = time.time()
        ite = 0
        isEmptylist =isEmptyList(self.NarrowList)

        while not isEmptylist and not self.stopWaveProp:#the main loop of FMM
                       
            neighbors = [-1,-1,-1,-1,-1,-1]
            neighbors_cross = [-1,-1,-1,-1]
            neighbors_cross1 = [-1,-1,-1,-1]  
            neighbors_cross2 = [-1,-1,-1,-1]  
            neighbors_cross3 = [-1,-1,-1,-1]  
            neighbors_cross4 = [-1,-1,-1,-1]  
            #n_neighs,n_neighs_cross,n_neighs_cross1,n_neighs_cross2,n_neighs_cross3,n_neighs_cross4 = 0,0,0,0,0,0
           
            best = locateBestIdx(self.NarrowList)
            #print('best ',best)
            if best ==-1:
                break    
            idxMin = best[1]
            self.nodesList[idxMin].type=FROZEN
                    
            n_neighs        = self.getneighbors(idxMin,neighbors,0,3)
            n_neighs_cross  = self.getneighbors(idxMin,neighbors_cross,1,3)
            n_neighs_cross1 = self.getneighbors(idxMin,neighbors_cross1,2,3)
            n_neighs_cross2 = self.getneighbors(idxMin,neighbors_cross2,3,3)
            n_neighs_cross3 = self.getneighbors(idxMin,neighbors_cross3,4,3)
            n_neighs_cross4 = self.getneighbors(idxMin,neighbors_cross4,5,3)

            #1
            for s in range(n_neighs):#MANHATAN NEIGHBORS
                
                
                j = neighbors[s]
                
                if  (j==-1):
                    continue
                curMN = self.nodesList[j]
                
                if (curMN.type ==FROZEN  or curMN.TAG==FORBIDDEN or curMN.TAG==NEW_FORBIDDEN or curMN.v==0.0):
                    continue
                
                else:
                   
                    cost =self.SolveEikonal(j)
                    if self.heuristicAct:
                        h=self.Hcost((curMN.x-g.x),(curMN.y-g.y),(curMN.z-g.z))#/curMN.v

                    
                    if self.NarrowList[curMN.block-1].__contains__((curMN.full,j)):

                        if cost+h<curMN.full:

                            indexe = self.NarrowList[curMN.block-1].index((curMN.full,j))
                            curMN.cost = cost
                            curMN.full = cost+h
                            _siftdown(self.NarrowList[curMN.block-1],0,indexe)
                    else:
                       
                        curMN.cost = cost
                        curMN.full = cost+h
                        curMN.type = KNOWN
                        heappush(self.NarrowList[curMN.block-1],(curMN.full,j))
            #2
            for s in range(n_neighs_cross):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross[s]    
                if  (j==-1):
                    continue
                
                ccur = self.nodesList[j]
                     
                if  (ccur.type ==FROZEN  or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN or ccur.v==0.0):
                
                    continue
                
                else:
                    if self.heuristicAct:
                        h=self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))#/ccur.v
                    cost =self.SolveEikonalcross(j)
                    
                    
                    if self.NarrowList[ccur.block-1].__contains__((ccur.full,j)):
                       
                        if cost+h<ccur.full:

                            indexe = self.NarrowList[ccur.block-1].index((ccur.full,j))
                            ccur.cost = cost
                            ccur.full = cost+h
                            _siftdown(self.NarrowList[ccur.block-1],0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.full = cost+h
                        ccur.type = KNOWN
                        heappush(self.NarrowList[ccur.block-1],(ccur.full,j))
            #3    
            for s in range(n_neighs_cross1):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross1[s]
                if  (j==-1):
                    continue
                
                ccur = self.nodesList[j]
                     
                if  (ccur.type ==FROZEN  or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN or ccur.v==0.0):
                
                    continue
                
                else:
                    
                    if self.heuristicAct:
                        h=self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))#/ccur.v
                    cost =self.SolveEikonalcross1(j)
                    print('iteration', ite,h, ccur.TAG, ccur.v,ccur.cost)
                    if self.NarrowList[ccur.block-1].__contains__((ccur.full,j)):
                       
                        if cost+h<ccur.full:

                            indexe = self.NarrowList[ccur.block-1].index((ccur.full,j))
                            ccur.cost = cost
                            ccur.full = cost+h
                            _siftdown(self.NarrowList[ccur.block-1],0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.full = cost+h
                        ccur.type = KNOWN
                        heappush(self.NarrowList[ccur.block-1],(ccur.full,j))
            #4
            for s in range(n_neighs_cross2):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross2[s]
                if  (j==-1):
                    continue
                
                ccur = self.nodesList[j]
                     
                if  (ccur.type ==FROZEN  or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN or ccur.v==0.0):
                
                    continue
                
                else:
                    if self.heuristicAct:
                        h=self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))#/ccur.v
                    cost =self.SolveEikonalcross2(j)
                    
                    if self.NarrowList[ccur.block-1].__contains__((ccur.full,j)):
                       
                        if cost+h<ccur.full:
                            indexe = self.NarrowList[ccur.block-1].index((ccur.full,j))
                            ccur.cost = cost
                            ccur.full = cost+h
                            _siftdown(self.NarrowList[ccur.block-1],0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.full = cost+h
                        ccur.type = KNOWN
                        heappush(self.NarrowList[ccur.block-1],(ccur.full,j))
            #5            
            for s in range(n_neighs_cross3):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross3[s]
                if  (j==-1):
                    continue
                                
                ccur = self.nodesList[j]
                     
                if  (ccur.type ==FROZEN  or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN or ccur.v==0.0):
                
                    continue
                
                else:
                    if self.heuristicAct:
                        h=self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))#/ccur.v
                    cost =self.SolveEikonalcross3(j)
                    
                    if self.NarrowList[ccur.block-1].__contains__((ccur.full,j)):
                       
                        if cost+h<ccur.full:

                            indexe = self.NarrowList[ccur.block-1].index((ccur.full,j))
                            ccur.cost = cost
                            ccur.full = cost+h
                            _siftdown(self.NarrowList[ccur.block-1],0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.full = cost+h
                        ccur.type = KNOWN
                        heappush(self.NarrowList[ccur.block-1],(ccur.full,j))
            #6           
            for s in range(n_neighs_cross4):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross4[s]
                if  (j==-1):
                    continue
                                
                ccur = self.nodesList[j]
                     
                if  (ccur.type ==FROZEN  or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN or ccur.v==0.0):
                
                    continue
                
                else:
                    if self.heuristicAct:
                        h=self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))#/ccur.v
                    cost =self.SolveEikonalcross4(j)
                    
                    
                    if self.NarrowList[ccur.block-1].__contains__((ccur.full,j)):
                       
                        if cost+h<ccur.full:

                            indexe = self.NarrowList[ccur.block-1].index((ccur.full,j))
                            ccur.cost = cost
                            ccur.full = cost+h
                            _siftdown(self.NarrowList[ccur.block-1],0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.full = cost+h
                        ccur.type = KNOWN
                        heappush(self.NarrowList[ccur.block-1],(ccur.full,j))
            ite+=1
                
            if idxMin == self.goal:
                self.stopWaveProp = True
            
            isEmptylist =isEmptyList(self.NarrowList)
            
        print time.time() - initime, 'seconds'
        return self.nodesList

class ParallelMSFM(object):
    
    def __init__(self,state = mp.Event(),queue =None,start_points=[], goal=None, nodeset = None,nodesList=[], d_ =[], cross=False, heuristicAct = False, seq = 0, block = 0, procID=-1):
        
        self.state = state
  
        self.block = block 
        self.d_ = d_

        self.procID = procID
        self.stopWaveProp = False
        self.Tvalue = []
        self.start_points = start_points
        self.goal = goal
        self.nodesList = nodesList
        self.queue = queue
        self.Narrowbnd = []
        self.nbr_neighbors = 0
        self.cross = cross
        self.nodeset = nodeset
        self.heuristicAct = heuristicAct
        self.Hcost = sqrt_dist
        
       
        self.minneighborsList = [-1 for _ in range(2)]

        
        print("waiting")
        self.state.wait()
        print("running") 
        
        self.processMSFM3D()

    def getneighbors(self, idx,ne,nume,dims):
        
        self.nbr_neighbors = 0
        if nume == 0:#standard neighbors
            #for i in range(dims): #x,y,z dimension
            self.nbr_neighbors =  self.getNeighborsInDim(idx,ne,0,-1)
                    
            return self.nbr_neighbors   
     
        elif nume ==1 :            
            # cross neighbors 1
            self.nbr_neighbors=self.getNeighborsInDimcross(idx,1,ne)
            return self.nbr_neighbors
        
        elif nume ==2 :            
            # cross neighbors 2
            self.nbr_neighbors=self.getNeighborsInDimcross(idx,2,ne)
            return self.nbr_neighbors
            
        elif nume ==3:            
            # cross neighbors 3
            self.nbr_neighbors=self.getNeighborsInDimcross(idx,3,ne)
            return self.nbr_neighbors
        
        elif nume ==4:            
            # cross neighbors 4
            self.nbr_neighbors=self.getNeighborsInDimcross(idx,4,ne)
            return self.nbr_neighbors
        
        elif nume ==5 :            
            # cross neighbors 5
            self.nbr_neighbors=self.getNeighborsInDimcross(idx,5,ne)
            return self.nbr_neighbors
       
    def getNeighborsInDim(self, idx,ne,nume, i):
          
        c1,c2,c3,c4,c5,c6 = -1,-1,-1,-1,-1,-1
     
        if nume==0:#neighbors

            c1 = idx-1
            c2 = idx+1
            c3 = int(floor(idx-self.d_[0]))
            c4 = int(floor(idx+self.d_[0]))
            c5= int(floor(idx-self.d_[1]))
            c6= int(floor(idx+self.d_[1]))   

            if ((c1>=0) and (int(floor(c1/self.d_[0]))==int(floor(idx/self.d_[0])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (int(floor(c2/self.d_[0]))==int(floor(idx/self.d_[0])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
            
            if ((c3>=0) and (int(floor(c3/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[self.nbr_neighbors] = int(c3)
                self.nbr_neighbors+=1
            
            if((c4>=0) and (int(floor(c4/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[self.nbr_neighbors]= int(c4)
                self.nbr_neighbors+=1
                
            if ((c5>=0) and (int(floor(c5/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c5)
                self.nbr_neighbors+=1
            
            if((c6>=0) and (int(floor(c6/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c6)
                self.nbr_neighbors+=1
                            
            return self.nbr_neighbors
        
        else:#min_neighbors
    
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
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1

            if((c2>=0) and(int(floor(c2/self.d_[i]))==int(floor(idx/self.d_[i])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
            
            return self.nbr_neighbors

    def getNeighborsInDimcross(self, idx,nbr,nec):
        cc1,cc2,cc3,cc4 = -1,-1,-1,-1
        if nbr==1:
        
            cc1 = int(floor(idx+self.d_[0]+1))
            cc2 = int(floor(idx-self.d_[0]-1))
            cc3 = int(floor(idx+self.d_[0]-1))
            cc4 = int(floor(idx-self.d_[0]+1))
            
            if ((cc1>=0) and (cc1<self.d_[1]) and (int(floor(cc1/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                nec[self.nbr_neighbors] = int(cc1)       
                self.nbr_neighbors+=1
            
            if((cc2>=0) and (cc2<self.d_[1]) and (int(floor(cc2/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                nec[self.nbr_neighbors]= int(cc2)    
                self.nbr_neighbors+=1
            
            if((cc3>=0) and (cc3<self.d_[1]) and (int(floor(cc3/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                nec[self.nbr_neighbors]= int(cc3)  
                self.nbr_neighbors+=1
            
            if((cc4>=0) and (cc4<self.d_[1]) and (int(floor(cc4/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                nec[self.nbr_neighbors]= int(cc4)
                self.nbr_neighbors+=1

            return self.nbr_neighbors
        
        elif nbr==2:
            tcc1,tcc2,tcc3,tcc4 = -1,-1,-1,-1
            
            tcc1 = int(floor(idx+self.d_[1]-1))
            tcc2 = int(floor(idx-self.d_[1]+1))
            tcc3 = int(floor(idx-self.d_[1]-1))
            tcc4 = int(floor(idx+self.d_[1]+1))
            #Top diagonals
            if ((tcc1>=0) and (tcc1<self.d_[2])and (int(floor(tcc1/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors] = int(tcc1)       
                self.nbr_neighbors+=1
            
            if((tcc2>=0) and (tcc2<self.d_[2]) and (int(floor(tcc2/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(tcc2)    
                self.nbr_neighbors+=1
            
            if((tcc3>=0) and (tcc3<self.d_[2]) and (int(floor(tcc3/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(tcc3)  
                self.nbr_neighbors+=1
            
            if((tcc4>=0) and (tcc4<self.d_[2]) and (int(floor(tcc4/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(tcc4)
                self.nbr_neighbors+=1

                  
            return self.nbr_neighbors

        elif nbr==3:        
            tcc5,tcc6,tcc7,tcc8 = -1,-1,-1,-1
            tcc5 = int(floor(idx+self.d_[1]-self.d_[0]))
            tcc6 = int(floor(idx-self.d_[1]+self.d_[0]))
            tcc7 = int(floor(idx-self.d_[1]-self.d_[0]))
            tcc8 = int(floor(idx+self.d_[1]+self.d_[0])) 
                                   
            if ((tcc5>=0) and (tcc5<self.d_[2])and (int(floor(tcc5/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors] = int(tcc5)       
                self.nbr_neighbors+=1
            
            if((tcc6>=0) and (tcc6<self.d_[2]) and (int(floor(tcc6/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(tcc6)    
                self.nbr_neighbors+=1
            
            if((tcc7>=0) and (tcc7<self.d_[2]) and (int(floor(tcc7/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(tcc7)  
                self.nbr_neighbors+=1
            
            if((tcc8>=0) and (tcc8<self.d_[2]) and (int(floor(tcc8/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(tcc8)
                self.nbr_neighbors+=1
        
            return self.nbr_neighbors

        elif nbr==4:
            bcc1,bcc2,bcc3,bcc4 = -1,-1,-1,-1
            bcc1 = int(floor(idx-self.d_[1]+self.d_[0]-1))
            bcc2 = int(floor(idx+self.d_[1]-self.d_[0]+1))
            bcc3 = int(floor(idx+self.d_[1]+self.d_[0]-1))
            bcc4 = int(floor(idx-self.d_[1]-self.d_[0]+1))   
            #Bottom diagonals   
            if ((bcc1>=0) and (bcc1<self.d_[2])and (int(floor(bcc1/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors] = int(bcc1)       
                self.nbr_neighbors+=1
            
            if((bcc2>=0) and (bcc2<self.d_[2]) and (int(floor(bcc2/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(bcc2)    
                self.nbr_neighbors+=1
            
            if((bcc3>=0) and (bcc3<self.d_[2]) and (int(floor(bcc3/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(bcc3)  
                self.nbr_neighbors+=1
            
            if((bcc4>=0) and (bcc4<self.d_[2]) and (int(floor(bcc4/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(bcc4)
                self.nbr_neighbors+=1
            
            return self.nbr_neighbors
    
        elif nbr==5:
            bcc5,bcc6,bcc7,bcc8 = -1,-1,-1,-1
            bcc5 = int(floor(idx+self.d_[1]-self.d_[0]-1))
            bcc6 = int(floor(idx-self.d_[1]+self.d_[0]+1))
            bcc7 = int(floor(idx-self.d_[1]-self.d_[0]-1))
            bcc8 = int(floor(idx+self.d_[1]+self.d_[0]+1))
            
            if ((bcc5>=0) and (bcc5<self.d_[2])and (int(floor(bcc5/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors] = int(bcc5)       
                self.nbr_neighbors+=1
            
            if((bcc6>=0) and (bcc6<self.d_[2]) and (int(floor(bcc6/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(bcc6)    
                self.nbr_neighbors+=1
            
            if((bcc7>=0) and (bcc7<self.d_[2]) and (int(floor(bcc7/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(bcc7)  
                self.nbr_neighbors+=1
            
            if((bcc8>=0) and (bcc8<self.d_[2]) and (int(floor(bcc8/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(bcc8)
                self.nbr_neighbors+=1
            
            return self.nbr_neighbors
    
    def getMinNeighborsInDimcross(self, idx, ne, i):    
        
        #global c1,c2,c3,c4
        c1,c2,c3,c4=-1,-1,-1,-1
        self.nbr_neighbors = 0
        if i==0:
            c1 = int(floor(idx+self.d_[0]+1))
            c2 = int(floor(idx-self.d_[0]-1))
        
            if ((c1>=0) and (c1<self.d_[1]) and (int(floor(c1/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[1]) and(int(floor(c2/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
                
        elif i==1:
            c3 = int(floor(idx+self.d_[0]-1))
            c4 = int(floor(idx-self.d_[0]+1))
         
            if ((c3>=0) and (c3<self.d_[1]) and (int(floor(c3/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[self.nbr_neighbors] = int(c3)
                self.nbr_neighbors+=1
            
            if((c4>=0) and (c4<self.d_[1])  and(int(floor(c4/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[self.nbr_neighbors]= int(c4)
                self.nbr_neighbors+=1
        # cross Neighbors +
        elif i==2:
            c1 = int(floor(idx-self.d_[1]-1))
            c2 = int(floor(idx+self.d_[1]+1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
       
        elif i==3:
            c1 = int(floor(idx-self.d_[1]+1))
            c2 = int(floor(idx+self.d_[1]-1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
                
        # cross neighbors +        
        elif i==4:
            c1 = int(floor(idx-self.d_[1]-self.d_[0]))
            c2 = int(floor(idx+self.d_[1]+self.d_[0]))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
                
        elif i==5:
            c1 = int(floor(idx-self.d_[1]+self.d_[0]))
            c2 = int(floor(idx+self.d_[1]-self.d_[0]))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
        #cross neighbors x        
        elif i==6:
            c1 = int(floor(idx-self.d_[1]-self.d_[0]-1))
            c2 = int(floor(idx+self.d_[1]+self.d_[0]+1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
        
        elif i==7:
            c1 = int(floor(idx+self.d_[1]-self.d_[0]-1))
            c2 = int(floor(idx-self.d_[1]+self.d_[0]+1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
                
        elif i==8:
            c1 = int(floor(idx-self.d_[1]+self.d_[0]-1))
            c2 = int(floor(idx+self.d_[1]-self.d_[0]+1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
                
        elif i==9:
            c1 = int(floor(idx+self.d_[1]+self.d_[0]-1))
            c2 = int(floor(idx-self.d_[1]-self.d_[0]+1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
                
        return self.nbr_neighbors
                    
    def getMinValueInDim(self, idx,dim):
        
        self.nbr_neighbors = 0        
        self.nbr_neighbors =  self.getNeighborsInDim(idx, self.minneighborsList,1, dim)
       
        if (self.nbr_neighbors==0) or self.nodesList[self.minneighborsList[0]].cost < self.nodesList[self.minneighborsList[1]].cost:
            return self.nodesList[self.minneighborsList[0]].cost
        else:
            return self.nodesList[self.minneighborsList[1]].cost
    
    def getMinValueInDimcross(self, idx,dim):
        #global minneighbors
        self.nbr_neighbors = 0
        
        self.nbr_neighbors=self.getMinNeighborsInDimcross(idx, self.minneighborsList, dim)
        #print idx,'min neighbors',self.minneighborsList
        if (self.nbr_neighbors==0) or self.nodesList[self.minneighborsList[0]].cost < self.nodesList[self.minneighborsList[1]].cost:
            return self.nodesList[self.minneighborsList[0]].cost
        else:
            return self.nodesList[self.minneighborsList[1]].cost

    def solveEikonalNDims(self, idx,dim):
            #global Tvalues   
            if dim==1:
                return self.Tvalues[0]+(1/self.nodesList[idx].v)
            
            sumT = 0
            sumTT = 0    
        
            for i in range(0,dim):
                sumT +=self.Tvalues[i]
                sumTT +=self.Tvalues[i]*self.Tvalues[i]
        
            a=dim
            b = -2*sumT
            c = sumTT-(1*1/pow(self.nodesList[idx].v,2))# *getNode(idx).v)
            q = b*b - 4*a*c
            
            if q<0:
                return INFINI
            else:
                return ((-b+sqrt(q))/(2*a))

    def solveEikonalNDimscross(self, curnode,dim):
               
            if dim==1:
              
                return self.Tvalues[0]+(sqrt(2)/curnode.v)
            
            sumT = 0
            sumTT = 0    
        
            for i in range(0,dim):
                sumT +=self.Tvalues[i]
                sumTT +=self.Tvalues[i]*self.Tvalues[i]
        
            a=dim
            b = -2*sumT
            c = sumTT-(sqrt(2)*sqrt(2)/pow(curnode.v,2))# *getNode(idx).v)
            q = b*b - 4*a*c
            
            if q<0:
                return INFINI
            else:
                return ((-b+sqrt(q))/(2*a))
              
    def SolveEikonal(self, idx):
     
        a = 3         
        self.Tvalues = []
        for dim in range(a):
            minT = 0
            minT = self.getMinValueInDim(idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
        
        if a==0:
            return INFINI
        
        self.Tvalues = sorted(self.Tvalues)
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDims(idx,i)
            if i==a or updatedT-self.Tvalues[i]<=0:
                break
        
        return updatedT
                
    def SolveEikonalcross(self, idx):
             
        a = 2#dimension x,y         
        self.Tvalues = []
        curn = self.nodesList[idx]
        for dim in range(a):
            minT = 0
            minT = self.getMinValueInDimcross(idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(curn,i)
            if i==a or updatedT-self.Tvalues[i]<0:
                break
        
        return updatedT
         
    def SolveEikonalcross1(self, idx):
        
        #global Tvalues       
        a = 2#dimension x,y         
        self.Tvalues = []
        curn = self.nodesList[idx]
       
        for dim in range(2,a+2):
            minT = 0
            minT = self.getMinValueInDimcross(idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(curn,i)
            if i==a or updatedT-self.Tvalues[i]<0:
                break
        
        return updatedT
            
    def SolveEikonalcross2(self, idx):
        
        #global Tvalues       
        a = 2#dimension x,y    
        curn = self.nodesList[idx]     
        self.Tvalues = []
       
        for dim in range(4,4+a):
            minT = 0
            minT = self.getMinValueInDimcross(idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(curn,i)
            if i==a or updatedT-self.Tvalues[i]<0:
                break
        
        return updatedT
            
    
    
        
    ########################################################################        

    def SolveEikonalcross3(self, idx):
        
        #global Tvalues       
        a = 2#dimension x,y 
        curn = self.nodesList[idx]        
        self.Tvalues = []
       
        for dim in range(6,6+a):
            minT = 0
            minT = self.getMinValueInDimcross(idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(curn,i)
            if i==a or updatedT-self.Tvalues[i]<0:
                break
        
        return updatedT     

    def SolveEikonalcross4(self, idx):
        
        #global Tvalues       
        curn = self.nodesList[idx]
        a = 2#dimension x,y         
        self.Tvalues = []
       
        for dim in range(8,8+a):
            minT = 0
            minT = self.getMinValueInDimcross(idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(curn,i)
            if i==a or updatedT-self.Tvalues[i]<0:
                break
        
        return updatedT 

    def processMSFM3D(self):
        leng=len(self.nodesList)
        print("Parallel begins, PLEASE Wait! ...")
        if self.goal !=None:
            g = self.nodesList[self.goal]#-leng*self.procID]
        h=0
        
        for i in self.start_points:

            initCur = self.nodesList[i]
            initCur.cost = 0
            if self.heuristicAct:        
                h = self.Hcost((initCur.x-g.x),(initCur.y-g.y),(initCur.z-g.z))#/initCur.v
            
            initCur.type = KNOWN
            initCur.full = 0+h
            #self.nodeset[initCur.indice] = initCur
            heappush(self.Narrowbnd,(initCur.full,i))
            #print((initCur.full,i))
                    
        print '===> PROCESS N ',self.procID,' start len',len(self.start_points)
        #time.sleep(50000)
        h,ite=0,0
        t = time.time()
        
        #print len(self.Narrowbnd )
        
        id2, idxMin = -2,-1
        
        while self.Narrowbnd  !=[] and not self.stopWaveProp:#the main loop of FMM
            id2 = idxMin
            neighbors = [-1,-1,-1,-1,-1,-1]
            neighbors_cross = [-1,-1,-1,-1]
            neighbors_cross1 = [-1,-1,-1,-1]  
            neighbors_cross2 = [-1,-1,-1,-1]  
            neighbors_cross3 = [-1,-1,-1,-1]  
            neighbors_cross4 = [-1,-1,-1,-1]  
            #n_neighs,n_neighs_cross,n_neighs_cross1,n_neighs_cross2,n_neighs_cross3,n_neighs_cross4 = 0,0,0,0,0,0
            ite+=1
            best = heappop(self.Narrowbnd)
            
            #print('best ',best)
            if best ==-1:
                break    
            idxMin = best[1]
 
            self.nodesList[idxMin].type=FROZEN
            #self.nodeset[self.nodesList[idxMin].indice] =FROZEN
            #print "Pr.",self.procID,"itr ",ite, 'Lst',len(self.nodesList)," i ",idxMin, 'nrd',len(self.Narrowbnd)#,'blk', self.nodesList[idxMin-leng*self.procID].block         
            n_neighs        = self.getneighbors(idxMin,neighbors,0,3)
            n_neighs_cross  = self.getneighbors(idxMin,neighbors_cross,1,3)
            n_neighs_cross1 = self.getneighbors(idxMin,neighbors_cross1,2,3)
            n_neighs_cross2 = self.getneighbors(idxMin,neighbors_cross2,3,3)
            n_neighs_cross3 = self.getneighbors(idxMin,neighbors_cross3,4,3)
            n_neighs_cross4 = self.getneighbors(idxMin,neighbors_cross4,5,3)

            #1
            for s in range(n_neighs):#MANHATAN NEIGHBORS
                
                
                j = neighbors[s]
                
                if  (j==-1):
                    continue
                curMN = self.nodesList[j]#-leng*self.procID]
                
                if (curMN.type == FROZEN  or curMN.TAG==FORBIDDEN or curMN.TAG==NEW_FORBIDDEN or curMN.v==0.0):
                    continue
                
                else:
                   
                    cost =self.SolveEikonal(j)
                    if self.heuristicAct:
                        h=self.Hcost((curMN.x-g.x),(curMN.y-g.y),(curMN.z-g.z))#/curMN.v

                    #print curMN.cost
                    if self.Narrowbnd .__contains__((curMN.full,j)):

                        if cost+h<curMN.full:

                            indexe = self.Narrowbnd .index((curMN.full,j))
                            curMN.cost = cost
                            curMN.full = cost+h
                            _siftdown(self.Narrowbnd ,0,indexe)
                    else:
                       
                        curMN.cost = cost
                        curMN.full = cost+h
                        curMN.type = KNOWN
                        heappush(self.Narrowbnd ,(curMN.full,j))
                    #self.nodeset[curMN.indice] = curMN
            #2
            for s in range(n_neighs_cross):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross[s]    
                if  (j==-1):
                    continue
                
                ccur = self.nodesList[j]#-leng*self.procID]
                     
                if  (ccur.type ==FROZEN  or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN or ccur.v==0.0 ):
                
                    continue
                
                else:
                    if self.heuristicAct:
                        h=self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))#/ccur.v
                    cost =self.SolveEikonalcross(j)
                    
                    
                    if self.Narrowbnd .__contains__((ccur.full,j)):
                       
                        if cost+h<ccur.full:

                            indexe = self.Narrowbnd .index((ccur.full,j))
                            ccur.cost = cost
                            ccur.full = cost+h
                            _siftdown(self.Narrowbnd ,0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.full = cost+h
                        ccur.type = KNOWN
                        heappush(self.Narrowbnd ,(ccur.full,j))
                    #self.nodeset[ccur.indice] = ccur
            #3    
            for s in range(n_neighs_cross1):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross1[s]
                if  (j==-1):
                    continue
                
                ccur = self.nodesList[j]#-leng*self.procID]
                     
                if  (ccur.type ==FROZEN  or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN or ccur.v==0.0 ):
                
                    continue
                
                else:
                    
                    if self.heuristicAct:
                        h=self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))#/ccur.v
                    cost =self.SolveEikonalcross1(j)
                    #print('iteration',h, ccur.TAG, ccur.v,ccur.cost)
                    if self.Narrowbnd .__contains__((ccur.full,j)):
                       
                        if cost+h<ccur.full:

                            indexe = self.Narrowbnd .index((ccur.full,j))
                            ccur.cost = cost
                            ccur.full = cost+h
                            _siftdown(self.Narrowbnd ,0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.full = cost+h
                        ccur.type = KNOWN
                        heappush(self.Narrowbnd ,(ccur.full,j))
                    #self.nodeset[ccur.indice] = ccur
            #4
            for s in range(n_neighs_cross2):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross2[s]
                if  (j==-1):
                    continue
                
                ccur = self.nodesList[j]#-leng*self.procID]
                     
                if  (ccur.type ==FROZEN  or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN or ccur.v==0.0 ):
                
                    continue
                
                else:
                    if self.heuristicAct:
                        h=self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))#/ccur.v
                    cost =self.SolveEikonalcross2(j)
                    
                    if self.Narrowbnd .__contains__((ccur.full,j)):
                       
                        if cost+h<ccur.full:
                            indexe = self.Narrowbnd .index((ccur.full,j))
                            ccur.cost = cost
                            ccur.full = cost+h
                            _siftdown(self.Narrowbnd ,0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.full = cost+h
                        ccur.type = KNOWN
                        heappush(self.Narrowbnd ,(ccur.full,j))
                    #self.nodeset[ccur.indice] = ccur
            #5            
            for s in range(n_neighs_cross3):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross3[s]
                if  (j==-1):
                    continue
                                
                ccur = self.nodesList[j]#-leng*self.procID]
                     
                if  (ccur.type ==FROZEN  or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN or ccur.v==0.0 ):
                
                    continue
                
                else:
                    if self.heuristicAct:
                        h=self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))#/ccur.v
                    cost =self.SolveEikonalcross3(j)
                    
                    if self.Narrowbnd .__contains__((ccur.full,j)):
                       
                        if cost+h<ccur.full:

                            indexe = self.Narrowbnd .index((ccur.full,j))
                            ccur.cost = cost
                            ccur.full = cost+h
                            _siftdown(self.Narrowbnd ,0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.full = cost+h
                        ccur.type = KNOWN
                        heappush(self.Narrowbnd ,(ccur.full,j))
                    #self.nodeset[ccur.indice] = ccur
            #6           
            for s in range(n_neighs_cross4):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross4[s]
                if  (j==-1):
                    continue
                                
                ccur = self.nodesList[j]#-leng*self.procID]
                     
                if  (ccur.type ==FROZEN  or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN or ccur.v==0.0  ):
                
                    continue
                
                else:
                    if self.heuristicAct:
                        h=self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))#/ccur.v
                    cost =self.SolveEikonalcross4(j)
                    
                    
                    if self.Narrowbnd .__contains__((ccur.full,j)):
                       
                        if cost+h<ccur.full:

                            indexe = self.Narrowbnd .index((ccur.full,j))
                            ccur.cost = cost
                            ccur.full = cost+h
                            _siftdown(self.Narrowbnd ,0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.full = cost+h
                        ccur.type = KNOWN
                        heappush(self.Narrowbnd ,(ccur.full,j))
                    #self.nodeset[ccur.indice] = ccur

            if idxMin == self.goal:
                self.stopWaveProp = True
        del self.Narrowbnd
        del self.start_points
            
        #self.queue.put((time.time()-t))
        return self.nodesList

class MSfm3D_SeqPar(object):

    def __init__(self,start_points=[],sblk = -1, goal=None,gblk = -1, nodesList=[], ghost = [], d_ =[],cross=False, heuristicAct = False, active = 0, block = 0):

        self.d_ = d_
        self.stopWaveProp = False
        self.nbr_neighbors = 0
        self.Tvalue = []
        self.start_points = start_points
        if goal !=None and goal !=-1: 
            self.goal = goal.idx
        else:
            self.goal = goal
        self.gblk = gblk
        self.sblk = sblk
        self.nodesList = nodesList
        self.lenlist = len(self.nodesList)
        self.ghost = ghost
        self.cross = cross
        self.block = block 
        #print self.block  
        self.boo = False     
        self.NarrowBand = [[] for _ in range(self.block)]
        self.iter = [0 for _ in range(self.block)]
        
        self.heuristicAct = heuristicAct
        self.Hcost = sqrt_dist
        self.actives = active
        self.exchanges =[False for _ in range(self.block)]    
        self.exchangesDATA = [[] for _ in range(self.block)]
        
        self.minneighborsList = [-1 for _ in range(2)]
        
        self.processMSFM3D()
    
    def getneighbors(self, blk, idx,ne,nume,dims):
        
        self.nbr_neighbors = 0
        if nume == 0:#standard neighbors
            #for i in range(dims): #x,y,z dimension
            self.nbr_neighbors =  self.getNeighborsInDim(blk, idx,ne,0,-1)
                    
            return self.nbr_neighbors   
     
        elif nume ==1 :            
            # cross neighbors 1
            self.nbr_neighbors=self.getNeighborsInDimcross(blk,idx,1,ne)
            return self.nbr_neighbors
        
        elif nume ==2 :            
            # cross neighbors 2
            self.nbr_neighbors=self.getNeighborsInDimcross(blk,idx,2,ne)
            return self.nbr_neighbors
            
        elif nume ==3:            
            # cross neighbors 3
            self.nbr_neighbors=self.getNeighborsInDimcross(blk,idx,3,ne)
            return self.nbr_neighbors
        
        elif nume ==4:            
            # cross neighbors 4
            self.nbr_neighbors=self.getNeighborsInDimcross(blk,idx,4,ne)
            return self.nbr_neighbors
        
        elif nume ==5 :            
            # cross neighbors 5
            self.nbr_neighbors=self.getNeighborsInDimcross(blk,idx,5,ne)
            return self.nbr_neighbors
       
    def getNeighborsInDim(self, blk,idx,ne,nume, i):
          
        c1,c2,c3,c4,c5,c6 = -1,-1,-1,-1,-1,-1
     
        if nume==0:#neighbors

            c1 = idx-1
            c2 = idx+1
            c3 = int(floor(idx-self.d_[blk][0]))
            c4 = int(floor(idx+self.d_[blk][0]))
            c5= int(floor(idx-self.d_[blk][1]))
            c6= int(floor(idx+self.d_[blk][1]))   
            print self.d_[blk], c1,c2,c3,c4,c5,c6
            if ((c1>=0) and (c1<self.d_[blk][2]) and (int(floor(c1/self.d_[blk][0]))==int(floor(idx/self.d_[blk][0])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[blk][2]) and (int(floor(c2/self.d_[blk][0]))==int(floor(idx/self.d_[blk][0])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
            
            if ((c3>=0) and (c3<self.d_[blk][2]) and (int(floor(c3/self.d_[blk][1]))==int(floor(idx/self.d_[blk][1])))):
                ne[self.nbr_neighbors] = int(c3)
                self.nbr_neighbors+=1
            
            if((c4>=0) and (c4<self.d_[blk][2])and (int(floor(c4/self.d_[blk][1]))==int(floor(idx/self.d_[blk][1])))):
                ne[self.nbr_neighbors]= int(c4)
                self.nbr_neighbors+=1
                
            if ((c5>=0) and (c5<self.d_[blk][2]) and (int(floor(c5/self.d_[blk][2]))==int(floor(idx/self.d_[blk][2])))):
                ne[self.nbr_neighbors] = int(c5)
                self.nbr_neighbors+=1
            
            if((c6>=0) and (c6<self.d_[blk][2]) and (int(floor(c6/self.d_[blk][2]))==int(floor(idx/self.d_[blk][2])))):
                ne[self.nbr_neighbors]= int(c6)
                self.nbr_neighbors+=1
                            
            return self.nbr_neighbors
        
        else:#min_neighbors
    
            if i==0:
                c1 = idx-1
                c2 = idx+1
            elif i==1:
                c1 = int(floor(idx-self.d_[blk][0]))
                c2 = int(floor(idx+self.d_[blk][0]))
            elif i==2:    
                c1 = int(floor(idx-self.d_[blk][1]))
                c2 = int(floor(idx+self.d_[blk][1]))
                            
            if ((c1>=0) and (c1<self.d_[blk][2]) and (int(floor(c1/self.d_[blk][i]))==int(floor(idx/self.d_[blk][i])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1

            if((c2>=0) and (c2<self.d_[blk][2]) and(int(floor(c2/self.d_[blk][i]))==int(floor(idx/self.d_[blk][i])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
            
            return self.nbr_neighbors

    def getNeighborsInDimcross(self, blk,idx,nbr,nec):
        cc1,cc2,cc3,cc4 = -1,-1,-1,-1
        if nbr==1:
        
            cc1 = int(floor(idx+self.d_[blk][0]+1))
            cc2 = int(floor(idx-self.d_[blk][0]-1))
            cc3 = int(floor(idx+self.d_[blk][0]-1))
            cc4 = int(floor(idx-self.d_[blk][0]+1))
            
            if ((cc1>=0) and (cc1<self.d_[blk][1]) and (int(floor(cc1/self.d_[blk][1]))==int(floor((idx)/self.d_[blk][1])))):
                nec[self.nbr_neighbors] = int(cc1)       
                self.nbr_neighbors+=1
            
            if((cc2>=0) and (cc2<self.d_[blk][1]) and (int(floor(cc2/self.d_[blk][1]))==int(floor((idx)/self.d_[blk][1])))):
                nec[self.nbr_neighbors]= int(cc2)    
                self.nbr_neighbors+=1
            
            if((cc3>=0) and (cc3<self.d_[blk][1]) and (int(floor(cc3/self.d_[blk][1]))==int(floor((idx)/self.d_[blk][1])))):
                nec[self.nbr_neighbors]= int(cc3)  
                self.nbr_neighbors+=1
            
            if((cc4>=0) and (cc4<self.d_[blk][1]) and (int(floor(cc4/self.d_[blk][1]))==int(floor((idx)/self.d_[blk][1])))):
                nec[self.nbr_neighbors]= int(cc4)
                self.nbr_neighbors+=1

            return self.nbr_neighbors
        
        elif nbr==2:
            tcc1,tcc2,tcc3,tcc4 = -1,-1,-1,-1
            
            tcc1 = int(floor(idx+self.d_[blk][1]-1))
            tcc2 = int(floor(idx-self.d_[blk][1]+1))
            tcc3 = int(floor(idx-self.d_[blk][1]-1))
            tcc4 = int(floor(idx+self.d_[blk][1]+1))
            #Top diagonals
            if ((tcc1>=0) and (tcc1<self.d_[blk][2])and (int(floor(tcc1/self.d_[blk][2]))==int(floor((idx)/self.d_[blk][2])))):
                nec[self.nbr_neighbors] = int(tcc1)       
                self.nbr_neighbors+=1
            
            if((tcc2>=0) and (tcc2<self.d_[blk][2]) and (int(floor(tcc2/self.d_[blk][2]))==int(floor((idx)/self.d_[blk][2])))):
                nec[self.nbr_neighbors]= int(tcc2)    
                self.nbr_neighbors+=1
            
            if((tcc3>=0) and (tcc3<self.d_[blk][2]) and (int(floor(tcc3/self.d_[blk][2]))==int(floor((idx)/self.d_[blk][2])))):
                nec[self.nbr_neighbors]= int(tcc3)  
                self.nbr_neighbors+=1
            
            if((tcc4>=0) and (tcc4<self.d_[blk][2]) and (int(floor(tcc4/self.d_[blk][2]))==int(floor((idx)/self.d_[blk][2])))):
                nec[self.nbr_neighbors]= int(tcc4)
                self.nbr_neighbors+=1

                  
            return self.nbr_neighbors

        elif nbr==3:        
            tcc5,tcc6,tcc7,tcc8 = -1,-1,-1,-1
            tcc5 = int(floor(idx+self.d_[blk][1]-self.d_[blk][0]))
            tcc6 = int(floor(idx-self.d_[blk][1]+self.d_[blk][0]))
            tcc7 = int(floor(idx-self.d_[blk][1]-self.d_[blk][0]))
            tcc8 = int(floor(idx+self.d_[blk][1]+self.d_[blk][0])) 
                                   
            if ((tcc5>=0) and (tcc5<self.d_[blk][2])and (int(floor(tcc5/self.d_[blk][2]))==int(floor((idx)/self.d_[blk][2])))):
                nec[self.nbr_neighbors] = int(tcc5)       
                self.nbr_neighbors+=1
            
            if((tcc6>=0) and (tcc6<self.d_[blk][2]) and (int(floor(tcc6/self.d_[blk][2]))==int(floor((idx)/self.d_[blk][2])))):
                nec[self.nbr_neighbors]= int(tcc6)    
                self.nbr_neighbors+=1
            
            if((tcc7>=0) and (tcc7<self.d_[blk][2]) and (int(floor(tcc7/self.d_[blk][2]))==int(floor((idx)/self.d_[blk][2])))):
                nec[self.nbr_neighbors]= int(tcc7)  
                self.nbr_neighbors+=1
            
            if((tcc8>=0) and (tcc8<self.d_[blk][2]) and (int(floor(tcc8/self.d_[blk][2]))==int(floor((idx)/self.d_[blk][2])))):
                nec[self.nbr_neighbors]= int(tcc8)
                self.nbr_neighbors+=1
        
            return self.nbr_neighbors

        elif nbr==4:
            bcc1,bcc2,bcc3,bcc4 = -1,-1,-1,-1
            bcc1 = int(floor(idx-self.d_[blk][1]+self.d_[blk][0]-1))
            bcc2 = int(floor(idx+self.d_[blk][1]-self.d_[blk][0]+1))
            bcc3 = int(floor(idx+self.d_[blk][1]+self.d_[blk][0]-1))
            bcc4 = int(floor(idx-self.d_[blk][1]-self.d_[blk][0]+1))   
            #Bottom diagonals   
            if ((bcc1>=0) and (bcc1<self.d_[blk][2])and (int(floor(bcc1/self.d_[blk][2]))==int(floor((idx)/self.d_[blk][2])))):
                nec[self.nbr_neighbors] = int(bcc1)       
                self.nbr_neighbors+=1
            
            if((bcc2>=0) and (bcc2<self.d_[blk][2]) and (int(floor(bcc2/self.d_[blk][2]))==int(floor((idx)/self.d_[blk][2])))):
                nec[self.nbr_neighbors]= int(bcc2)    
                self.nbr_neighbors+=1
            
            if((bcc3>=0) and (bcc3<self.d_[blk][2]) and (int(floor(bcc3/self.d_[blk][2]))==int(floor((idx)/self.d_[blk][2])))):
                nec[self.nbr_neighbors]= int(bcc3)  
                self.nbr_neighbors+=1
            
            if((bcc4>=0) and (bcc4<self.d_[blk][2]) and (int(floor(bcc4/self.d_[blk][2]))==int(floor((idx)/self.d_[blk][2])))):
                nec[self.nbr_neighbors]= int(bcc4)
                self.nbr_neighbors+=1
            
            return self.nbr_neighbors
    
        elif nbr==5:
            bcc5,bcc6,bcc7,bcc8 = -1,-1,-1,-1
            bcc5 = int(floor(idx+self.d_[blk][1]-self.d_[blk][0]-1))
            bcc6 = int(floor(idx-self.d_[blk][1]+self.d_[blk][0]+1))
            bcc7 = int(floor(idx-self.d_[blk][1]-self.d_[blk][0]-1))
            bcc8 = int(floor(idx+self.d_[blk][1]+self.d_[blk][0]+1))
            
            if ((bcc5>=0) and (bcc5<self.d_[blk][2])and (int(floor(bcc5/self.d_[blk][2]))==int(floor((idx)/self.d_[blk][2])))):
                nec[self.nbr_neighbors] = int(bcc5)       
                self.nbr_neighbors+=1
            
            if((bcc6>=0) and (bcc6<self.d_[blk][2]) and (int(floor(bcc6/self.d_[blk][2]))==int(floor((idx)/self.d_[blk][2])))):
                nec[self.nbr_neighbors]= int(bcc6)    
                self.nbr_neighbors+=1
            
            if((bcc7>=0) and (bcc7<self.d_[blk][2]) and (int(floor(bcc7/self.d_[blk][2]))==int(floor((idx)/self.d_[blk][2])))):
                nec[self.nbr_neighbors]= int(bcc7)  
                self.nbr_neighbors+=1
            
            if((bcc8>=0) and (bcc8<self.d_[blk][2]) and (int(floor(bcc8/self.d_[blk][2]))==int(floor((idx)/self.d_[blk][2])))):
                nec[self.nbr_neighbors]= int(bcc8)
                self.nbr_neighbors+=1
            
            return self.nbr_neighbors
    
    def getMinNeighborsInDimcross(self, blk,idx, ne, i):    
        
        #global c1,c2,c3,c4
        c1,c2,c3,c4=-1,-1,-1,-1
        self.nbr_neighbors = 0
        if i==0:
            c1 = int(floor(idx+self.d_[blk][0]+1))
            c2 = int(floor(idx-self.d_[blk][0]-1))
        
            if ((c1>=0) and (c1<self.d_[blk][2]) and (int(floor(c1/self.d_[blk][1]))==int(floor(idx/self.d_[blk][1])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[blk][2]) and (int(floor(c2/self.d_[blk][1]))==int(floor(idx/self.d_[blk][1])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
            print self.d_[blk], c1,c2
                
        elif i==1:
            c3 = int(floor(idx+self.d_[blk][0]-1))
            c4 = int(floor(idx-self.d_[blk][0]+1))
         
            if ((c3>=0) and (c3<self.d_[blk][2]) and (int(floor(c3/self.d_[blk][1]))==int(floor(idx/self.d_[blk][1])))):
                ne[self.nbr_neighbors] = int(c3)
                self.nbr_neighbors+=1
            
            if((c4>=0) and (c4<self.d_[blk][2])  and(int(floor(c4/self.d_[blk][1]))==int(floor(idx/self.d_[blk][1])))):
                ne[self.nbr_neighbors]= int(c4)
                self.nbr_neighbors+=1
            print self.d_[blk],c3,c4
        # cross Neighbors +
        elif i==2:
            c1 = int(floor(idx-self.d_[blk][1]-1))
            c2 = int(floor(idx+self.d_[blk][1]+1))
         
            if ((c1>=0) and (c1<self.d_[blk][2]) and (int(floor(c1/self.d_[blk][2]))==int(floor(idx/self.d_[blk][2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[blk][2])  and(int(floor(c2/self.d_[blk][2]))==int(floor(idx/self.d_[blk][2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
            print self.d_[blk], c1,c2
        elif i==3:
            c1 = int(floor(idx-self.d_[blk][1]+1))
            c2 = int(floor(idx+self.d_[blk][1]-1))
         
            if ((c1>=0) and (c1<self.d_[blk][2]) and (int(floor(c1/self.d_[blk][2]))==int(floor(idx/self.d_[blk][2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[blk][2])  and(int(floor(c2/self.d_[blk][2]))==int(floor(idx/self.d_[blk][2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
            print self.d_[blk], c1,c2    
        # cross neighbors +        
        elif i==4:
            c1 = int(floor(idx-self.d_[blk][1]-self.d_[blk][0]))
            c2 = int(floor(idx+self.d_[blk][1]+self.d_[blk][0]))
         
            if ((c1>=0) and (c1<self.d_[blk][2]) and (int(floor(c1/self.d_[blk][2]))==int(floor(idx/self.d_[blk][2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[blk][2])  and(int(floor(c2/self.d_[blk][2]))==int(floor(idx/self.d_[blk][2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
            print self.d_[blk], c1,c2    
        elif i==5:
            c1 = int(floor(idx-self.d_[blk][1]+self.d_[blk][0]))
            c2 = int(floor(idx+self.d_[blk][1]-self.d_[blk][0]))
         
            if ((c1>=0) and (c1<self.d_[blk][2]) and (int(floor(c1/self.d_[blk][2]))==int(floor(idx/self.d_[blk][2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[blk][2])  and(int(floor(c2/self.d_[blk][2]))==int(floor(idx/self.d_[blk][2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
            print self.d_[blk], c1,c2
        #cross neighbors x        
        elif i==6:
            c1 = int(floor(idx-self.d_[blk][1]-self.d_[blk][0]-1))
            c2 = int(floor(idx+self.d_[blk][1]+self.d_[blk][0]+1))
         
            if ((c1>=0) and (c1<self.d_[blk][2]) and (int(floor(c1/self.d_[blk][2]))==int(floor(idx/self.d_[blk][2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[blk][2])  and(int(floor(c2/self.d_[blk][2]))==int(floor(idx/self.d_[blk][2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
            print self.d_[blk], c1,c2
        elif i==7:
            c1 = int(floor(idx+self.d_[blk][1]-self.d_[blk][0]-1))
            c2 = int(floor(idx-self.d_[blk][1]+self.d_[blk][0]+1))
         
            if ((c1>=0) and (c1<self.d_[blk][2]) and (int(floor(c1/self.d_[blk][2]))==int(floor(idx/self.d_[blk][2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            
            if((c2>=0) and (c2<self.d_[blk][2])  and(int(floor(c2/self.d_[blk][2]))==int(floor(idx/self.d_[blk][2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
            print self.d_[blk], c1,c2    
        elif i==8:
            c1 = int(floor(idx-self.d_[blk][1]+self.d_[blk][0]-1))
            c2 = int(floor(idx+self.d_[blk][1]-self.d_[blk][0]+1))
         
            if ((c1>=0) and (c1<self.d_[blk][2]) and (int(floor(c1/self.d_[blk][2]))==int(floor(idx/self.d_[blk][2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[blk][2])  and(int(floor(c2/self.d_[blk][2]))==int(floor(idx/self.d_[blk][2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
            print self.d_[blk], c1,c2    
        elif i==9:
            c1 = int(floor(idx+self.d_[blk][1]+self.d_[blk][0]-1))
            c2 = int(floor(idx-self.d_[blk][1]-self.d_[blk][0]+1))
         
            if ((c1>=0) and (c1<self.d_[blk][2]) and (int(floor(c1/self.d_[blk][2]))==int(floor(idx/self.d_[blk][2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[blk][2])  and(int(floor(c2/self.d_[blk][2]))==int(floor(idx/self.d_[blk][2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
            print self.d_[blk], c1,c2    
        return self.nbr_neighbors
                    
    def getMinValueInDim(self,blk, idx,dim):
        
        self.nbr_neighbors = 0        
        self.nbr_neighbors =  self.getNeighborsInDim(blk, idx, self.minneighborsList,1, dim)

        if (self.nbr_neighbors==0) or self.nodesList[blk][self.minneighborsList[0]].cost < self.nodesList[blk][self.minneighborsList[1]].cost:
            #print '                 ',dim, len(self.nodesList[blk]), self.d_[blk], self.minneighborsList[0]
            return self.nodesList[blk][self.minneighborsList[0]].cost
        else:
            return self.nodesList[blk][self.minneighborsList[1]].cost
    
    def getMinValueInDimcross(self, blk, idx,dim):
        #global minneighbors
        self.nbr_neighbors = 0        
        self.nbr_neighbors=self.getMinNeighborsInDimcross(blk,idx, self.minneighborsList, dim)
        print 'list ',len(self.nodesList[blk]),self.minneighborsList
        if (self.nbr_neighbors==0) or self.nodesList[blk][self.minneighborsList[0]].cost < self.nodesList[blk][self.minneighborsList[1]].cost:
            return self.nodesList[blk][self.minneighborsList[0]].cost
        else:
            return self.nodesList[blk][self.minneighborsList[1]].cost

    def solveEikonalNDims(self, blk, idx,dim):
            #global Tvalues   
            if dim==1:
                return self.Tvalues[0]+(1/self.nodesList[blk][idx].v)
            
            sumT = 0
            sumTT = 0    
        
            for i in range(0,dim):
                sumT +=self.Tvalues[i]
                sumTT +=self.Tvalues[i]*self.Tvalues[i]
        
            a=dim
            b = -2*sumT
            c = sumTT-(1*1/pow(self.nodesList[blk][idx].v,2))# *getNode(idx).v)
            q = b*b - 4*a*c
            
            if q<0:
                return INFINI
            else:
                return ((-b+sqrt(q))/(2*a))

    def solveEikonalNDimscross(self, curnode,dim):
               
            if dim==1:
              
                return self.Tvalues[0]+(sqrt(2)/curnode.v)
            
            sumT = 0
            sumTT = 0    
        
            for i in range(0,dim):
                sumT +=self.Tvalues[i]
                sumTT +=self.Tvalues[i]*self.Tvalues[i]
        
            a=dim
            b = -2*sumT
            c = sumTT-(sqrt(2)*sqrt(2)/pow(curnode.v,2))# *getNode(idx).v)
            q = b*b - 4*a*c
            
            if q<0:
                return INFINI
            else:
                return ((-b+sqrt(q))/(2*a))
              
    def SolveEikonal(self,blk, idx):
     
        a = 3         
        self.Tvalues = []
        for dim in range(a):
            minT = 0
            minT = self.getMinValueInDim(blk,idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[blk][idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
        
        if a==0:
            return INFINI
        
        self.Tvalues = sorted(self.Tvalues)
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDims(blk,idx,i)
            if i==a or updatedT-self.Tvalues[i]<=0:
                break
        
        return updatedT
                
    def SolveEikonalcross(self,blk, idx):
             
        a = 2#dimension x,y         
        self.Tvalues = []
        curn = self.nodesList[blk][idx]
        for dim in range(a):
            minT = 0
            minT = self.getMinValueInDimcross(blk,idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[blk][idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(curn,i)
            if i==a or updatedT-self.Tvalues[i]<0:
                break
        
        return updatedT
         
    def SolveEikonalcross1(self, blk,idx):
        
        #global Tvalues       
        a = 2#dimension x,y         
        self.Tvalues = []
        curn = self.nodesList[blk][idx]
       
        for dim in range(2,a+2):
            minT = 0
            minT = self.getMinValueInDimcross(blk,idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[blk][idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(curn,i)
            if i==a or updatedT-self.Tvalues[i]<0:
                break
        
        return updatedT
            
    def SolveEikonalcross2(self, blk,idx):
        
        #global Tvalues       
        a = 2#dimension x,y    
        curn = self.nodesList[blk][idx]     
        self.Tvalues = []
       
        for dim in range(4,4+a):
            minT = 0
            minT = self.getMinValueInDimcross(blk,idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[blk][idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(curn,i)
            if i==a or updatedT-self.Tvalues[i]<0:
                break
        
        return updatedT
            
    
    
        
    ########################################################################        

    def SolveEikonalcross3(self, blk,idx):
        
        #global Tvalues       
        a = 2#dimension x,y 
        curn = self.nodesList[blk][idx]        
        self.Tvalues = []
       
        for dim in range(6,6+a):
            minT = 0
            minT = self.getMinValueInDimcross(blk,idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[blk][idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(curn,i)
            if i==a or updatedT-self.Tvalues[i]<0:
                break
        
        return updatedT     

    def SolveEikonalcross4(self, blk,idx):
        
        #global Tvalues       
        curn = self.nodesList[blk][idx]
        
        a = 2#dimension x,y         
        self.Tvalues = []
       
        for dim in range(8,8+a):
            minT = 0
            minT = self.getMinValueInDimcross(blk,idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[blk][idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(curn,i)
            if i==a or updatedT-self.Tvalues[i]<0:
                break
        
        return updatedT 
    
    def init_heap(self, idblk, s_p, g, nrbnd):
        
        h=0
        for i in s_p:
            initCur = self.nodesList[idblk][i]
            #print 'blk init',idblk,initCur.block
            if self.heuristicAct:        
                h = self.Hcost((initCur.x-g.x),(initCur.y-g.y),(initCur.z-g.z))#/initCur.v
                           
            initCur.cost = 0+h
            initCur.type = FROZEN_FIX
            heappush(self.NarrowBand[idblk],(initCur.cost,initCur.idx))
            self.actives[idblk] = True
            self.exchanges[idblk] = False  

    def marchNarrowBnd(self, blk, g):
        
        ite = 0

        print "Block number ", blk
        if g != 0:
            print len(self.nodesList[blk]),len(self.NarrowBand[blk]) 
            
        while self.NarrowBand[blk] != [] and not self.stopWaveProp:#the main loop of FMM

            neighbors = [-1,-1,-1,-1,-1,-1]
            neighbors_cross = [-1,-1,-1,-1]
            neighbors_cross1 = [-1,-1,-1,-1]  
            neighbors_cross2 = [-1,-1,-1,-1]  
            neighbors_cross3 = [-1,-1,-1,-1]  
            neighbors_cross4 = [-1,-1,-1,-1]  
            n_neighs,n_neighs_cross,n_neighs_cross1,n_neighs_cross2,n_neighs_cross3,n_neighs_cross4 = 0,0,0,0,0,0
           
            best = heappop(self.NarrowBand[blk])   
            idxMin = best[1]
                    
            n_neighs        = self.getneighbors(blk,idxMin,neighbors,0,3)
            n_neighs_cross  = self.getneighbors(blk,idxMin,neighbors_cross,1,3)
            n_neighs_cross1 = self.getneighbors(blk,idxMin,neighbors_cross1,2,3)
            n_neighs_cross2 = self.getneighbors(blk,idxMin,neighbors_cross2,3,3)
            n_neighs_cross3 = self.getneighbors(blk,idxMin,neighbors_cross3,4,3)
            n_neighs_cross4 = self.getneighbors(blk,idxMin,neighbors_cross4,5,3)
            #print('block, iteration ',blk, ite,'blk nds', nodesList[idxMin].block,'2 ?', nodesList[idxMin].type)
            if self.nodesList[blk][idxMin].type != FROZEN_FIX:
                self.nodesList[blk][idxMin].type =FROZEN
                
            if self.nodesList[blk][idxMin].block != blk:
                #if self.iter[self.nodesList[blk][idxMin].block] ==0 :
                self.exchangesDATA[self.nodesList[blk][idxMin].block].append(self.nodesList[blk][idxMin])
                self.exchanges[self.nodesList[blk][idxMin].block] = True
                if g != 0: 
                    print 'ghost ',self.nodesList[blk][idxMin].block,self.exchanges[self.nodesList[blk][idxMin].block]                     
                    
            if g!=0:
                print 'itera', self.d_[blk],ite,self.nodesList[blk][idxMin].TAG, self.nodesList[blk][idxMin].OBSTACLE,\
                self.nodesList[blk][idxMin].x, self.nodesList[blk][idxMin].y, self.nodesList[blk][idxMin].z,\
                len(self.nodesList[blk]),len(self.NarrowBand[blk]), idxMin,\
                 n_neighs, n_neighs_cross, n_neighs_cross1, n_neighs_cross2,n_neighs_cross3, n_neighs_cross4            
            h =0
            #1
            for s in range(n_neighs):#MANHATAN NEIGHBORS
                
                j = neighbors[s]
                curMN = self.nodesList[blk][j]
                
                if  (j==-1):
                    continue
                
                elif curMN.type == FROZEN_FIX or curMN.type ==FROZEN  or curMN.TAG==FORBIDDEN or curMN.TAG==NEW_FORBIDDEN:
                    
                    continue
                
                else:
                    if self.heuristicAct:
                        h=self.Hcost((curMN.x-g.x),(curMN.y-g.y),(curMN.z-g.z))#/curMN.v
                    
                    cost =self.SolveEikonal(blk, j)+h#self.Hcost((curMN.x-g.x),(curMN.y-g.y),(curMN.z-g.z))/curMN.v
                    
                    if self.NarrowBand[blk].__contains__((curMN.cost,j)):
                       
                        if cost<curMN.cost:
                            
                            indexe = self.NarrowBand[blk].index((curMN.cost,j))
                            curMN.cost = cost
                          
                            _siftdown(self.NarrowBand[blk],0,indexe)
                    else:
                        
                        curMN.cost = cost
                        curMN.type = KNOWN
                        heappush(self.NarrowBand[blk],(cost,j))
                    """
                    if curMN.block!=blk:
                        #print  blk,curMN.block
                        self.exchangesDATA[curMN.block].append(curMN)
                        self.exchanges[curMN.block] = True
                    """
                #print 'cur ',curMN.type,' ',curMN.cost,' pos',curMN.x,curMN.y,curMN.z,'  Obst',curMN.OBSTACLE 
            #2
            for s in range(n_neighs_cross):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross[s]
                ccur = self.nodesList[blk][j]
                     
                if  (j==-1 or ccur.type ==FROZEN or ccur.type ==FROZEN_FIX or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN):
                
                    continue
                
                else:
                    if self.heuristicAct:
                        h=self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))#/curMN.v
                                       
                    cost =self.SolveEikonalcross(blk,j)+h
                    
                    if self.NarrowBand[blk].__contains__((ccur.cost,j)):
                       
                        if cost<ccur.cost:
                            
                            indexe = self.NarrowBand[blk].index((ccur.cost,j))
                            ccur.cost = cost
                            _siftdown(self.NarrowBand[blk],0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.type = KNOWN
                        heappush(self.NarrowBand[blk],(cost,j))
                    """
                    if ccur.block!=blk:
                        #print  blk,ccur.block
                        self.exchangesDATA[ccur.block].append(ccur)
                        self.exchanges[ccur.block] = True                
                    """
            #3    
            for s in range(n_neighs_cross1):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross1[s]
                ccur = self.nodesList[blk][j]
                     
                if  (j==-1 or ccur.type ==FROZEN  or ccur.type ==FROZEN_FIX or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN):
                
                    continue
                
                else:
                    if self.heuristicAct:
                        h=self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))#/curMN.v                   
                    cost =self.SolveEikonalcross1(blk,j)+h
                    
                    if self.NarrowBand[blk].__contains__((ccur.cost,j)):
                       
                        if cost<ccur.cost:
                            
                            indexe = self.NarrowBand[blk].index((ccur.cost,j))
                            ccur.cost = cost
                            _siftdown(self.NarrowBand[blk],0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.type = KNOWN
                        heappush(self.NarrowBand[blk],(cost,j))
                    """    
                    if ccur.block!=blk:
                        #print  blk,ccur.block
                        self.exchangesDATA[ccur.block].append(ccur)
                        self.exchanges[ccur.block] = True                
                    """    
                #print 'ccur ',ccur.type,' ',ccur.cost,' pos',ccur.x,ccur.y,ccur.z,'  Obst',ccur.OBSTACLE 
            #4
            for s in range(n_neighs_cross2):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross2[s]
                ccur = self.nodesList[blk][j]
                     
                if  (j==-1 or ccur.type ==FROZEN or ccur.type ==FROZEN_FIX or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN):
                
                    continue
                
                else:
                    if self.heuristicAct:
                        h=self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))#/curMN.v                   
                    cost =self.SolveEikonalcross2(blk,j)+h
                    
                    if self.NarrowBand[blk].__contains__((ccur.cost,j)):
                       
                        if cost<ccur.cost:
                            
                            indexe = self.NarrowBand[blk].index((ccur.cost,j))
                            ccur.cost = cost
                            _siftdown(self.NarrowBand[blk],0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.type = KNOWN
                        heappush(self.NarrowBand[blk],(cost,j))
                    """
                    if ccur.block!=blk:
                        #print  blk,ccur.block
                        self.exchangesDATA[ccur.block].append(ccur)
                        self.exchanges[ccur.block] = True
                    """                

                #print 'ccur ',ccur.type,' ',ccur.cost,' pos',ccur.x,ccur.y,ccur.z,'  Obst',ccur.OBSTACLE 
            #5            
            for s in range(n_neighs_cross3):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross3[s]
                ccur = self.nodesList[blk][j]
                     
                if  (j==-1 or ccur.type ==FROZEN or ccur.type ==FROZEN_FIX or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN):
                
                    continue
                
                else:
                    if self.heuristicAct:
                        h=self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))#/curMN.v                   
                    cost =self.SolveEikonalcross3(blk,j)+h
                    
                    if self.NarrowBand[blk].__contains__((ccur.cost,j)):
                       
                        if cost<ccur.cost:
                            
                            indexe = self.NarrowBand[blk].index((ccur.cost,j))
                            ccur.cost = cost
                            _siftdown(self.NarrowBand[blk],0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.type = KNOWN
                        heappush(self.NarrowBand[blk],(cost,j))
                    """
                    if ccur.block!=blk:
                        #print  blk,ccur.block
                        
                        self.exchangesDATA[ccur.block].append(ccur)
                        self.exchanges[ccur.block] = True
                    """                

                #print 'ccur ',ccur.type,' ',ccur.cost,' pos',ccur.x,ccur.y,ccur.z,'  Obst',ccur.OBSTACLE 
                
            #6           
            for s in range(n_neighs_cross4):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross4[s]
                ccur = self.nodesList[blk][j]
                     
                if  (j==-1 or ccur.type ==FROZEN or ccur.type ==FROZEN_FIX or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN):
                
                    continue
                
                else:
                    if self.heuristicAct:
                        h=self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))#/curMN.v                   
                    cost =self.SolveEikonalcross4(blk,j)+h
                    
                    if self.NarrowBand[blk].__contains__((ccur.cost,j)):
                       
                        if cost<ccur.cost:
                            
                            indexe = self.NarrowBand[blk].index((ccur.cost,j))
                            ccur.cost = cost
                            _siftdown(self.NarrowBand[blk],0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.type = KNOWN
                        heappush(self.NarrowBand[blk],(cost,j))
                    """
                    if ccur.block!=blk:
                        #print  blk,ccur.block
                        
                        self.exchangesDATA[ccur.block].append(ccur)
                        self.exchanges[ccur.block] = True
                    """
                        
                #print 'ccur ',ccur.type,' ',ccur.cost,' pos',ccur.x,ccur.y,ccur.z,'  Obst',ccur.OBSTACLE 

            ite+=1    
            if idxMin == self.goal:
                print 'goalll'
                self.stopWaveProp = True
                    
             
        self.iter[blk]=1        
        self.actives[blk]= False

    def GetGhostDataNeighs(self, idBlk):
        
        for  ghst in self.exchangesDATA[idBlk]:
            
            cur=self.nodesList[idBlk][ghst.sndidx]

            if (ghst.cost<cur.cost) and  cur.type != FROZEN_FIX:
                #print('nrbd ',idBlk,'idx ', ghst.sndidx,'cidx', cur.idx,'-',ghst.type,cur.type, '--c g',cur.cost, ghst.cost)
                if self.NarrowBand[idBlk].__contains__((cur.cost,cur.idx)):

                    indexe = self.NarrowBand[idBlk].index((cur.cost,cur.idx))
                    cur.cost = ghst.cost
                    cur.type = ghst.type
 
                    _siftdown(self.NarrowBand[idBlk],0,indexe)
               
                else:
                    cur.cost = ghst.cost
                    cur.type = ghst.type
                    heappush(self.NarrowBand[idBlk],(cur.cost,cur.idx))
            else:    
                cur.type = FROZEN_FIX
               
            
        self.exchangesDATA[idBlk] = []   
        self.exchanges[idBlk] = False  

        if self.NarrowBand[idBlk] !=[]:

            if self.iter[idBlk]==0: 
                self.actives[idBlk] =  True
                print 'exchange ',idBlk, 'activated ' 
                
                self.boo = True
        else:
            self.actives[idBlk] =  False

    def processMSFM3D(self):
        g=0
        if self.goal !=None and self.goal !=-1:
            g = self.nodesList[self.gblk][self.goal] 

        initime = time.time()

        #initialize the local heaps
        print('initializing local Heaps ...', self.block)
        
        if self.sblk !=-1:
            
            self.init_heap(self.sblk, self.start_points,g, self.NarrowBand[self.sblk])
            print 'block',self.sblk, self.start_points, g, self.gblk, self.NarrowBand[self.sblk]
        else:
                
            for i in range(0,self.block):#
                #print i
                if self.start_points[i] != []:
                    self.init_heap(i, self.start_points[i],g, self.NarrowBand[i])


        self.boo = True
        while self.boo:#marching steps
            
            
            print('Marching step ...')
            for i in range(self.block):#marching blocks
                if self.actives[i]:                  
                    
                    self.marchNarrowBnd(i, g)
                else:
                    print 'blocked ', i
                    
            self.boo = False
            print('Exchange step ...')
            for i in range(self.block):
                
                if self.exchanges[i]:
                        
                    self.GetGhostDataNeighs(i)
                    #del self.exchangesDATA
                          

        print time.time() - initime, 'seconds'
        
        return self.nodesList 

"""

class MSfm3D_SeqPar(object):

    def __init__(self,start_points=[], goal=None, nodesList=[], ghost = [], d_ =[], cross=False, heuristicAct = False, active = 0, block = 0):
        
        self.d_[blk] = d_
        
        self.stopWaveProp = False
        self.nbr_neighbors = 0
        self.Tvalue = []
        self.start_points = start_points
        self.goal = goal
        self.nodesList = nodesList
        self.ghost = ghost
        self.cross = cross
        self.block = block 
        self.NarrowBand = [[] for _ in range(self.block)]
        self.heuristicAct = heuristicAct
        self.Hcost = sqrt_dist
        self.actives = active
        self.exchanges =[False for _ in range(self.block)]    
        if self.cross:
            self.neighborsList = [-1 for _ in range(26)]
        else:
            self.neighborsList = [-1 for _ in range(2**(len(d_)-1)+2)]
        self.minneighborsList = [-1 for _ in range(2)]
        self.processMSFM3D()
    
    def getneighbors(self, idx,ne,nume,dims):
        
        self.nbr_neighbors = 0
        if nume == 0:#standard neighbors
            #for i in range(dims): #x,y,z dimension
            self.nbr_neighbors =  self.getNeighborsInDim(idx,ne,0,-1)
                    
            return self.nbr_neighbors   
     
        elif nume ==1 :            
            # cross neighbors 1
            self.nbr_neighbors=self.getNeighborsInDimcross(idx,1,ne)
            return self.nbr_neighbors
        
        elif nume ==2 :            
            # cross neighbors 2
            self.nbr_neighbors=self.getNeighborsInDimcross(idx,2,ne)
            return self.nbr_neighbors
            
        elif nume ==3:            
            # cross neighbors 3
            self.nbr_neighbors=self.getNeighborsInDimcross(idx,3,ne)
            return self.nbr_neighbors
        
        elif nume ==4:            
            # cross neighbors 4
            self.nbr_neighbors=self.getNeighborsInDimcross(idx,4,ne)
            return self.nbr_neighbors
        
        elif nume ==5 :            
            # cross neighbors 5
            self.nbr_neighbors=self.getNeighborsInDimcross(idx,5,ne)
            return self.nbr_neighbors
       
    def getNeighborsInDim(self, idx,ne,nume, i):
          
        c1,c2,c3,c4,c5,c6 = -1,-1,-1,-1,-1,-1
     
        if nume==0:#neighbors

            c1 = idx-1
            c2 = idx+1
            c3 = int(floor(idx-self.d_[0]))
            c4 = int(floor(idx+self.d_[0]))
            c5= int(floor(idx-self.d_[1]))
            c6= int(floor(idx+self.d_[1]))   

            if ((c1>=0) and (int(floor(c1/self.d_[0]))==int(floor(idx/self.d_[0])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (int(floor(c2/self.d_[0]))==int(floor(idx/self.d_[0])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
            
            if ((c3>=0) and (int(floor(c3/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[self.nbr_neighbors] = int(c3)
                self.nbr_neighbors+=1
            
            if((c4>=0) and (int(floor(c4/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[self.nbr_neighbors]= int(c4)
                self.nbr_neighbors+=1
                
            if ((c5>=0) and (int(floor(c5/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c5)
                self.nbr_neighbors+=1
            
            if((c6>=0) and (int(floor(c6/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c6)
                self.nbr_neighbors+=1
                            
            return self.nbr_neighbors
        
        else:#min_neighbors
    
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
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1

            if((c2>=0) and(int(floor(c2/self.d_[i]))==int(floor(idx/self.d_[i])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
            
            return self.nbr_neighbors

    def getNeighborsInDimcross(self, idx,nbr,nec):
        cc1,cc2,cc3,cc4 = -1,-1,-1,-1
        if nbr==1:
        
            cc1 = int(floor(idx+self.d_[0]+1))
            cc2 = int(floor(idx-self.d_[0]-1))
            cc3 = int(floor(idx+self.d_[0]-1))
            cc4 = int(floor(idx-self.d_[0]+1))
            
            if ((cc1>=0) and (cc1<self.d_[1]) and (int(floor(cc1/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                nec[self.nbr_neighbors] = int(cc1)       
                self.nbr_neighbors+=1
            
            if((cc2>=0) and (cc2<self.d_[1]) and (int(floor(cc2/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                nec[self.nbr_neighbors]= int(cc2)    
                self.nbr_neighbors+=1
            
            if((cc3>=0) and (cc3<self.d_[1]) and (int(floor(cc3/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                nec[self.nbr_neighbors]= int(cc3)  
                self.nbr_neighbors+=1
            
            if((cc4>=0) and (cc4<self.d_[1]) and (int(floor(cc4/self.d_[1]))==int(floor((idx)/self.d_[1])))):
                nec[self.nbr_neighbors]= int(cc4)
                self.nbr_neighbors+=1

            return self.nbr_neighbors
        
        elif nbr==2:
            tcc1,tcc2,tcc3,tcc4 = -1,-1,-1,-1
            
            tcc1 = int(floor(idx+self.d_[1]-1))
            tcc2 = int(floor(idx-self.d_[1]+1))
            tcc3 = int(floor(idx-self.d_[1]-1))
            tcc4 = int(floor(idx+self.d_[1]+1))
            #Top diagonals
            if ((tcc1>=0) and (tcc1<self.d_[2])and (int(floor(tcc1/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors] = int(tcc1)       
                self.nbr_neighbors+=1
            
            if((tcc2>=0) and (tcc2<self.d_[2]) and (int(floor(tcc2/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(tcc2)    
                self.nbr_neighbors+=1
            
            if((tcc3>=0) and (tcc3<self.d_[2]) and (int(floor(tcc3/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(tcc3)  
                self.nbr_neighbors+=1
            
            if((tcc4>=0) and (tcc4<self.d_[2]) and (int(floor(tcc4/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(tcc4)
                self.nbr_neighbors+=1

                  
            return self.nbr_neighbors

        elif nbr==3:        
            tcc5,tcc6,tcc7,tcc8 = -1,-1,-1,-1
            tcc5 = int(floor(idx+self.d_[1]-self.d_[0]))
            tcc6 = int(floor(idx-self.d_[1]+self.d_[0]))
            tcc7 = int(floor(idx-self.d_[1]-self.d_[0]))
            tcc8 = int(floor(idx+self.d_[1]+self.d_[0])) 
                                   
            if ((tcc5>=0) and (tcc5<self.d_[2])and (int(floor(tcc5/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors] = int(tcc5)       
                self.nbr_neighbors+=1
            
            if((tcc6>=0) and (tcc6<self.d_[2]) and (int(floor(tcc6/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(tcc6)    
                self.nbr_neighbors+=1
            
            if((tcc7>=0) and (tcc7<self.d_[2]) and (int(floor(tcc7/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(tcc7)  
                self.nbr_neighbors+=1
            
            if((tcc8>=0) and (tcc8<self.d_[2]) and (int(floor(tcc8/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(tcc8)
                self.nbr_neighbors+=1
        
            return self.nbr_neighbors

        elif nbr==4:
            bcc1,bcc2,bcc3,bcc4 = -1,-1,-1,-1
            bcc1 = int(floor(idx-self.d_[1]+self.d_[0]-1))
            bcc2 = int(floor(idx+self.d_[1]-self.d_[0]+1))
            bcc3 = int(floor(idx+self.d_[1]+self.d_[0]-1))
            bcc4 = int(floor(idx-self.d_[1]-self.d_[0]+1))   
            #Bottom diagonals   
            if ((bcc1>=0) and (bcc1<self.d_[2])and (int(floor(bcc1/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors] = int(bcc1)       
                self.nbr_neighbors+=1
            
            if((bcc2>=0) and (bcc2<self.d_[2]) and (int(floor(bcc2/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(bcc2)    
                self.nbr_neighbors+=1
            
            if((bcc3>=0) and (bcc3<self.d_[2]) and (int(floor(bcc3/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(bcc3)  
                self.nbr_neighbors+=1
            
            if((bcc4>=0) and (bcc4<self.d_[2]) and (int(floor(bcc4/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(bcc4)
                self.nbr_neighbors+=1
            
            return self.nbr_neighbors
    
        elif nbr==5:
            bcc5,bcc6,bcc7,bcc8 = -1,-1,-1,-1
            bcc5 = int(floor(idx+self.d_[1]-self.d_[0]-1))
            bcc6 = int(floor(idx-self.d_[1]+self.d_[0]+1))
            bcc7 = int(floor(idx-self.d_[1]-self.d_[0]-1))
            bcc8 = int(floor(idx+self.d_[1]+self.d_[0]+1))
            
            if ((bcc5>=0) and (bcc5<self.d_[2])and (int(floor(bcc5/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors] = int(bcc5)       
                self.nbr_neighbors+=1
            
            if((bcc6>=0) and (bcc6<self.d_[2]) and (int(floor(bcc6/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(bcc6)    
                self.nbr_neighbors+=1
            
            if((bcc7>=0) and (bcc7<self.d_[2]) and (int(floor(bcc7/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(bcc7)  
                self.nbr_neighbors+=1
            
            if((bcc8>=0) and (bcc8<self.d_[2]) and (int(floor(bcc8/self.d_[2]))==int(floor((idx)/self.d_[2])))):
                nec[self.nbr_neighbors]= int(bcc8)
                self.nbr_neighbors+=1
            
            return self.nbr_neighbors
    
    def getMinNeighborsInDimcross(self, idx, ne, i):    
        
        #global c1,c2,c3,c4
        c1,c2,c3,c4=-1,-1,-1,-1
        self.nbr_neighbors = 0
        if i==0:
            c1 = int(floor(idx+self.d_[0]+1))
            c2 = int(floor(idx-self.d_[0]-1))
        
            if ((c1>=0) and (c1<self.d_[1]) and (int(floor(c1/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[1]) and(int(floor(c2/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
                
        elif i==1:
            c3 = int(floor(idx+self.d_[0]-1))
            c4 = int(floor(idx-self.d_[0]+1))
         
            if ((c3>=0) and (c3<self.d_[1]) and (int(floor(c3/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[self.nbr_neighbors] = int(c3)
                self.nbr_neighbors+=1
            
            if((c4>=0) and (c4<self.d_[1])  and(int(floor(c4/self.d_[1]))==int(floor(idx/self.d_[1])))):
                ne[self.nbr_neighbors]= int(c4)
                self.nbr_neighbors+=1
        # cross Neighbors +
        elif i==2:
            c1 = int(floor(idx-self.d_[1]-1))
            c2 = int(floor(idx+self.d_[1]+1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
       
        elif i==3:
            c1 = int(floor(idx-self.d_[1]+1))
            c2 = int(floor(idx+self.d_[1]-1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
                
        # cross neighbors +        
        elif i==4:
            c1 = int(floor(idx-self.d_[1]-self.d_[0]))
            c2 = int(floor(idx+self.d_[1]+self.d_[0]))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
                
        elif i==5:
            c1 = int(floor(idx-self.d_[1]+self.d_[0]))
            c2 = int(floor(idx+self.d_[1]-self.d_[0]))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
        #cross neighbors x        
        elif i==6:
            c1 = int(floor(idx-self.d_[1]-self.d_[0]-1))
            c2 = int(floor(idx+self.d_[1]+self.d_[0]+1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
        
        elif i==7:
            c1 = int(floor(idx+self.d_[1]-self.d_[0]-1))
            c2 = int(floor(idx-self.d_[1]+self.d_[0]+1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
                
        elif i==8:
            c1 = int(floor(idx-self.d_[1]+self.d_[0]-1))
            c2 = int(floor(idx+self.d_[1]-self.d_[0]+1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
                
        elif i==9:
            c1 = int(floor(idx+self.d_[1]+self.d_[0]-1))
            c2 = int(floor(idx-self.d_[1]-self.d_[0]+1))
         
            if ((c1>=0) and (c1<self.d_[2]) and (int(floor(c1/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors] = int(c1)
                self.nbr_neighbors+=1
            
            if((c2>=0) and (c2<self.d_[2])  and(int(floor(c2/self.d_[2]))==int(floor(idx/self.d_[2])))):
                ne[self.nbr_neighbors]= int(c2)
                self.nbr_neighbors+=1
                
        return self.nbr_neighbors
                    
    def getMinValueInDim(self,blk, idx,dim):
        
        self.nbr_neighbors = 0        
        self.nbr_neighbors =  self.getNeighborsInDim(idx, self.minneighborsList,1, dim)
       
        if (self.nbr_neighbors==0) or self.nodesList[blk][self.minneighborsList[0]].cost < self.nodesList[blk][self.minneighborsList[1]].cost:
            return self.nodesList[blk][self.minneighborsList[0]].cost
        else:
            return self.nodesList[blk][self.minneighborsList[1]].cost
    
    def getMinValueInDimcross(self, blk, idx,dim):
        #global minneighbors
        self.nbr_neighbors = 0
        
        self.nbr_neighbors=self.getMinNeighborsInDimcross(idx, self.minneighborsList, dim)
        #print self.minneighborsList,' - ', self.nodesList[blk][self.minneighborsList[0]]
        if (self.nbr_neighbors==0) or self.nodesList[blk][self.minneighborsList[0]].cost < self.nodesList[blk][self.minneighborsList[1]].cost:
            return self.nodesList[blk][self.minneighborsList[0]].cost
        else:
            return self.nodesList[blk][self.minneighborsList[1]].cost

    def solveEikonalNDims(self, blk, idx,dim):
            #global Tvalues   
            if dim==1:
                return self.Tvalues[0]+(1/self.nodesList[blk][idx].v)
            
            sumT = 0
            sumTT = 0    
        
            for i in range(0,dim):
                sumT +=self.Tvalues[i]
                sumTT +=self.Tvalues[i]*self.Tvalues[i]
        
            a=dim
            b = -2*sumT
            c = sumTT-(1*1/pow(self.nodesList[blk][idx].v,2))# *getNode(idx).v)
            q = b*b - 4*a*c
            
            if q<0:
                return INFINI
            else:
                return ((-b+sqrt(q))/(2*a))

    def solveEikonalNDimscross(self, curnode,dim):
               
            if dim==1:
              
                return self.Tvalues[0]+(sqrt(2)/curnode.v)
            
            sumT = 0
            sumTT = 0    
        
            for i in range(0,dim):
                sumT +=self.Tvalues[i]
                sumTT +=self.Tvalues[i]*self.Tvalues[i]
        
            a=dim
            b = -2*sumT
            c = sumTT-(sqrt(2)*sqrt(2)/pow(curnode.v,2))# *getNode(idx).v)
            q = b*b - 4*a*c
            
            if q<0:
                return INFINI
            else:
                return ((-b+sqrt(q))/(2*a))
              
    def SolveEikonal(self,blk, idx):
     
        a = 3         
        self.Tvalues = []
        for dim in range(a):
            minT = 0
            minT = self.getMinValueInDim(blk,idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[blk][idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
        
        if a==0:
            return INFINI
        
        self.Tvalues = sorted(self.Tvalues)
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDims(blk,idx,i)
            if i==a or updatedT-self.Tvalues[i]<=0:
                break
        
        return updatedT
                
    def SolveEikonalcross(self,blk, idx):
             
        a = 2#dimension x,y         
        self.Tvalues = []
        curn = self.nodesList[blk][idx]
        for dim in range(a):
            minT = 0
            minT = self.getMinValueInDimcross(blk,idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[blk][idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(curn,i)
            if i==a or updatedT-self.Tvalues[i]<0:
                break
        
        return updatedT
         
    def SolveEikonalcross1(self, blk,idx):
        
        #global Tvalues       
        a = 2#dimension x,y         
        self.Tvalues = []
        curn = self.nodesList[blk][idx]
       
        for dim in range(2,a+2):
            minT = 0
            minT = self.getMinValueInDimcross(blk,idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[blk][idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(curn,i)
            if i==a or updatedT-self.Tvalues[i]<0:
                break
        
        return updatedT
            
    def SolveEikonalcross2(self, blk,idx):
        
        #global Tvalues       
        a = 2#dimension x,y    
        curn = self.nodesList[blk][idx]     
        self.Tvalues = []
       
        for dim in range(4,4+a):
            minT = 0
            minT = self.getMinValueInDimcross(blk,idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[blk][idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(curn,i)
            if i==a or updatedT-self.Tvalues[i]<0:
                break
        
        return updatedT
            
    
    
        
    ########################################################################        

    def SolveEikonalcross3(self, blk,idx):
        
        #global Tvalues       
        a = 2#dimension x,y 
        curn = self.nodesList[blk][idx]        
        self.Tvalues = []
       
        for dim in range(6,6+a):
            minT = 0
            minT = self.getMinValueInDimcross(blk,idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[blk][idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(curn,i)
            if i==a or updatedT-self.Tvalues[i]<0:
                break
        
        return updatedT     

    def SolveEikonalcross4(self, blk,idx):
        
        #global Tvalues       
        curn = self.nodesList[blk][idx]
        
        a = 2#dimension x,y         
        self.Tvalues = []
       
        for dim in range(8,8+a):
            minT = 0
            minT = self.getMinValueInDimcross(blk,idx,dim)
                       
            if minT != INFINI and minT< self.nodesList[blk][idx].cost:
                heappush(self.Tvalues,minT) 
            
            else:
                a-=1
             
        if a==0:
            return INFINI
        
        
        updatedT = -1
        
        for i in range(1,1+a):
        
            updatedT = self.solveEikonalNDimscross(curn,i)
            if i==a or updatedT-self.Tvalues[i]<0:
                break
        
        return updatedT 
    
    def init_heap(self, idblk, s_p, g, nrbnd):
        h=0
        for i in s_p:
            initCur = self.nodesList[idblk][i]
            if self.heuristicAct:        
                h = self.Hcost((initCur.x-g.x),(initCur.y-g.y),(initCur.z-g.z))#/initCur.v
                           
            initCur.cost = 0+h
            initCur.type = FROZEN_FIX
            heappush(nrbnd,(initCur.cost,initCur.indice))
            self.actives[idblk] = True  

    def init_ghosts(self, idblk):
        
        
        
        
        return -1


    def marchNarrowBnd(self, blk, NarrowBand, nodesList, ghostnodes, g):
        
        ite = 0
        
        
        if NarrowBand==[]:
                
                self.actives[blk] = False
                print('            ****False block, Deactive ',blk,self.actives[blk])
                return False
        
        while NarrowBand != [] and not self.stopWaveProp:#the main loop of FMM
            
            print('block, iteration ',blk, ite)
            neighbors = [-1,-1,-1,-1,-1,-1]
            neighbors_cross = [-1,-1,-1,-1]
            neighbors_cross1 = [-1,-1,-1,-1]  
            neighbors_cross2 = [-1,-1,-1,-1]  
            neighbors_cross3 = [-1,-1,-1,-1]  
            neighbors_cross4 = [-1,-1,-1,-1]  
            n_neighs,n_neighs_cross,n_neighs_cross1,n_neighs_cross2,n_neighs_cross3,n_neighs_cross4 = 0,0,0,0,0,0
           
            best = heappop(NarrowBand)   
            idxMin = best[1]
                    
            n_neighs        = self.getneighbors(idxMin,neighbors,0,3)
            n_neighs_cross  = self.getneighbors(idxMin,neighbors_cross,1,3)
            n_neighs_cross1 = self.getneighbors(idxMin,neighbors_cross1,2,3)
            n_neighs_cross2 = self.getneighbors(idxMin,neighbors_cross2,3,3)
            n_neighs_cross3 = self.getneighbors(idxMin,neighbors_cross3,4,3)
            n_neighs_cross4 = self.getneighbors(idxMin,neighbors_cross4,5,3)
            
            if nodesList[idxMin].type != FROZEN_FIX:
                nodesList[idxMin].type =FROZEN
            
            h =0
            #1
            for s in range(n_neighs):#MANHATAN NEIGHBORS
                
                j = neighbors[s]
                curMN = nodesList[j]
                if  (j==-1 or curMN.type ==FROZEN_FIX or curMN.type ==FROZEN  or curMN.TAG==FORBIDDEN or curMN.TAG==NEW_FORBIDDEN):
                    
                    continue
                
                else:
                    if self.heuristicAct:
                        h=self.Hcost((curMN.x-g.x),(curMN.y-g.y),(curMN.z-g.z))#/curMN.v
                    
                    cost =self.SolveEikonal(blk, j)+h#self.Hcost((curMN.x-g.x),(curMN.y-g.y),(curMN.z-g.z))/curMN.v
                    
                    if NarrowBand.__contains__((curMN.cost,j)):
                       
                        if cost<curMN.cost:
                            
                            indexe = NarrowBand.index((curMN.cost,j))
                            curMN.cost = cost
                          
                            _siftdown(NarrowBand,0,indexe)
                    else:
                        
                        curMN.cost = cost
                        curMN.type = KNOWN
                        heappush(NarrowBand,(cost,j))
                
                #print 'cur ',curMN.type,' ',curMN.cost,' pos',curMN.x,curMN.y,curMN.z,'  Obst',curMN.OBSTACLE 
            #2
            for s in range(n_neighs_cross):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross[s]
                ccur = nodesList[j]
                     
                if  (j==-1 or ccur.type ==FROZEN or ccur.type ==FROZEN_FIX or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN):
                
                    continue
                
                else:
                    if self.heuristicAct:
                        h=self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))#/curMN.v
                                       
                    cost =self.SolveEikonalcross(blk,j)+h
                    
                    if NarrowBand.__contains__((ccur.cost,j)):
                       
                        if cost<ccur.cost:
                            
                            indexe = NarrowBand.index((ccur.cost,j))
                            ccur.cost = cost
                            _siftdown(NarrowBand,0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.type = KNOWN
                        heappush(NarrowBand,(cost,j))
                
                #print 'ccur ',ccur.type,' ',ccur.cost,' pos',ccur.x,ccur.y,ccur.z,'  Obst',ccur.OBSTACLE 
            #3    
            for s in range(n_neighs_cross1):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross1[s]
                ccur = nodesList[j]
                     
                if  (j==-1 or ccur.type ==FROZEN  or ccur.type ==FROZEN_FIX or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN):
                
                    continue
                
                else:
                    if self.heuristicAct:
                        h=self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))#/curMN.v                   
                    cost =self.SolveEikonalcross1(blk,j)+h
                    
                    if NarrowBand.__contains__((ccur.cost,j)):
                       
                        if cost<ccur.cost:
                            
                            indexe = NarrowBand.index((ccur.cost,j))
                            ccur.cost = cost
                            _siftdown(NarrowBand,0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.type = KNOWN
                        heappush(NarrowBand,(cost,j))
                
                #print 'ccur ',ccur.type,' ',ccur.cost,' pos',ccur.x,ccur.y,ccur.z,'  Obst',ccur.OBSTACLE 
            #4
            for s in range(n_neighs_cross2):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross2[s]
                ccur = nodesList[j]
                     
                if  (j==-1 or ccur.type ==FROZEN or ccur.type ==FROZEN_FIX or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN):
                
                    continue
                
                else:
                    if self.heuristicAct:
                        h=self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))#/curMN.v                   
                    cost =self.SolveEikonalcross2(blk,j)+h
                    
                    if NarrowBand.__contains__((ccur.cost,j)):
                       
                        if cost<ccur.cost:
                            
                            indexe = NarrowBand.index((ccur.cost,j))
                            ccur.cost = cost
                            _siftdown(NarrowBand,0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.type = KNOWN
                        heappush(NarrowBand,(cost,j))
                
                #print 'ccur ',ccur.type,' ',ccur.cost,' pos',ccur.x,ccur.y,ccur.z,'  Obst',ccur.OBSTACLE 
            #5            
            for s in range(n_neighs_cross3):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross3[s]
                ccur = nodesList[j]
                     
                if  (j==-1 or ccur.type ==FROZEN or ccur.type ==FROZEN_FIX or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN):
                
                    continue
                
                else:
                    if self.heuristicAct:
                        h=self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))#/curMN.v                   
                    cost =self.SolveEikonalcross3(blk,j)+h
                    
                    if NarrowBand.__contains__((ccur.cost,j)):
                       
                        if cost<ccur.cost:
                            
                            indexe = NarrowBand.index((ccur.cost,j))
                            ccur.cost = cost
                            _siftdown(NarrowBand,0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.type = KNOWN
                        heappush(NarrowBand,(cost,j))
                
                #print 'ccur ',ccur.type,' ',ccur.cost,' pos',ccur.x,ccur.y,ccur.z,'  Obst',ccur.OBSTACLE 
                
            #6           
            for s in range(n_neighs_cross4):#xy EUCLIDIAN NEIGHBORS 
                
                j = neighbors_cross4[s]
                ccur = nodesList[j]
                     
                if  (j==-1 or ccur.type ==FROZEN or ccur.type ==FROZEN_FIX or ccur.TAG==FORBIDDEN or ccur.TAG==NEW_FORBIDDEN):
                
                    continue
                
                else:
                    if self.heuristicAct:
                        h=self.Hcost((ccur.x-g.x),(ccur.y-g.y),(ccur.z-g.z))#/curMN.v                   
                    cost =self.SolveEikonalcross4(blk,j)+h
                    
                    if NarrowBand.__contains__((ccur.cost,j)):
                       
                        if cost<ccur.cost:
                            
                            indexe = NarrowBand.index((ccur.cost,j))
                            ccur.cost = cost
                            _siftdown(NarrowBand,0,indexe)
                    else:
                       
                        ccur.cost = cost
                        ccur.type = KNOWN
                        heappush(NarrowBand,(cost,j))
                
                #print 'ccur ',ccur.type,' ',ccur.cost,' pos',ccur.x,ccur.y,ccur.z,'  Obst',ccur.OBSTACLE 
            ite+=1    
            if idxMin == self.goal:
                self.stopWaveProp = True
        self.actives[blk]= False
        print self.actives[blk]
        return self.actives[blk]

    def GetGhostDataNeighs(self, idblk):
        print('ok') 
        
        
        
        
        
        

    def processMSFM3D(self):
        
        if self.goal !=None:
            g = self.nodesList[self.goal] 
        #print('start ',self.start_points[2])
        initime = time.time()
        #initialize the local heaps
        for i in range(0,self.block):#
            self.actives[i] = False
            self.init_heap(i, self.start_points[i],g, self.NarrowBand[i])
        
        con = True    
        
        
        while con:#marching steps
            con =False
            
            for i in range(self.block):#marching blocks
                if self.actives[i]:
                    self.actives[i] = self.marchNarrowBnd(i, self.NarrowBand[i], self.nodesList[i], self.ghost[i], g)
                    con = True
            
            for i in range(self.block):
                if self.exchanges[i]:
                    self.exchanges[i] = self.GetGhostDataNeighs(i)
                           

        print time.time() - initime, 'seconds'
        return self.nodesList 

"""