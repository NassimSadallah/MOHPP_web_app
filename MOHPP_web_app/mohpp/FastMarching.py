'''
Created on Jul 8, 2021

@author: nassim
'''
from __builtin__ import False
import time
from math import floor, sqrt
from heapq import _siftdown, heappop, heappush
from utilities import FORBIDDEN, INFINI, NEW_FORBIDDEN
from utilities import isEmptyList, locateBestIdx
class fmm(object):
    '''
    Fast Marching method with (8 and/or 16 neighbors )
    '''
    
    def __init__(self, start_points = [], goal = None, nodesList = [], d_ = [], cross = False, heuristicActivate = False, seq = 0, block = 0):
        
        self.d_ = d_
        self.narrowBand = []
        self.stopWave = False
        self.nbr_neighbors = 0
        self.Tvalue = []
        self.start_points = start_points
        self.goal = goal
        self.nodesList = nodesList
        self.cross = cross
        self.heuristicActivate = heuristicActivate
        #self.seq = seq
        self.block = block
        if self.cross:
            self.neighborsList = [-1 for _ in range(8)]
            
        else:
            self.neighborsList = [-1 for _ in range(4)]
        self.minneighborsList = [-1 for _ in range(2)]
        
        if seq == 0:
            self.processFMM()
        elif seq == 1:
            self.processSequentialFMM()
        
    
    def getneighbors(self,idx):        
        self.nbr_neighbors = 0

        self.nbr_neighbors = self.getNeighborsInDim(idx,self.neighborsList, -1, False)
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
                
                elif dim>0:    
                
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
            self.getNeighborsInDim(idx, self.minneighborsList, dim, False)
           
            if (self.nbr_neighbors==0):
                return self.nodesList[self.minneighborsList[0]].cost
            
            elif self.nodesList[self.minneighborsList[0]].cost < self.nodesList[self.minneighborsList[1]].cost:
                    return self.nodesList[self.minneighborsList[0]].cost
            else:
                    return self.nodesList[self.minneighborsList[1]].cost  
        
        else:
            self.getNeighborsInDim(idx, self.minneighborsList, dim, True)
           
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
            c = sumTT-sqrt(2)*sqrt(2)/(self.nodesList[idx].v * self.nodesList[idx].v)
        else:
            c = sumTT-1*1/(self.nodesList[idx].v * self.nodesList[idx].v)
        q = b*b - 4*a*c
        if q<0:
            return INFINI
        else:
            return ((-b+sqrt(q))/(2*a))

    def SolveEikonal(self,idx, crossing):
        
        if crossing ==False:
            a = 2
                     
            self.Tvalue = []
            
            for i in range(a):
                minT = 0
                minT = self.getMinValueInDim(idx,i, False)
                           
                if minT != INFINI and minT< self.nodesList[idx].cost:
                    heappush(self.Tvalue,minT) 
                
                else:
                    a-=1
            
            if a==0:
                return INFINI
            
            self.Tvalue = sorted(self.Tvalue)
            updatedT = -1
            
            for i in range(1,1+a):
            
                updatedT = self.solveEikonalNDims(idx,i)
                
                if i==a or updatedT-self.Tvalue[i]<=0:
                    break
            
            return updatedT
        
        if crossing == True:
            
            a = 2
                     
            self.Tvalue = []
            
            for i in range(a):
                minT = 0
                minT = self.getMinValueInDim(idx,i, True)
                           
                if minT != INFINI and minT< self.nodesList[idx].cost:
                    heappush(self.Tvalue,minT) 
                
                else:
                    a-=1
            
            if a==0:
                return INFINI
            
            self.Tvalue = sorted(self.Tvalue)
            updatedT = -1
            
            for i in range(1,1+a):
            
                updatedT = self.solveEikonalNDims(idx,i)
                
                if i==a or updatedT-self.Tvalue[i]<=0:
                    break
            
            return updatedT
    
        
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
            neighbors = self.getneighbors(idxMin)
            self.nodesList[idxMin].type = 'A'
            
            for j in neighbors:
                
                if j==-1:
                    continue
                 
                current = self.nodesList[j]
             
                if current.type == 'A' or current.TAG ==FORBIDDEN or current.TAG == NEW_FORBIDDEN:
                    continue 
               
                cost = self.SolveEikonal(j, current.cros)
                
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
