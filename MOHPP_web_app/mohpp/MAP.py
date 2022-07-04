'''
Created on Jul 8, 2021

@author: nassim
'''
import time
from PIL import Image
from Node import node, node3
from utilities import height, d_, width, depth, FORBIDDEN, CONTOUR,NO_OBSTACLE, KNOWN_OBSTACLE, UNDETECTED_OBSTACLE,d_parallel
from math import floor
from copy import deepcopy

tabGrid, GridMap, nodeList,Nodelist = [],[],[],[]

class MAP_2D(object):
        
    '''
    Read a binary map with static obstacles
    The image should be at the same package as the current script
    r, r2 are, respectively, width and the height of the matrix
    '''
    
    #r = width, r2 = height
    def readMap(self,r, r2,binaryMap):
    
        image = Image.open(binaryMap)
        
        for i in range (r):
            tabGrid[0][i] = -1
            tabGrid[r2-1][i] = -1
        for i in range (r2):
            tabGrid[i][0] = -1
            tabGrid[i][r-1] = -1
                
        for x in range(r2):# height
            for y in range(r): #width
                coordinate = y,x
            
                if image.getpixel(coordinate)==(0, 0, 0):
                    tabGrid[x][y] = -1
               
    '''   
    define which pixels are referring to the obstacles
    ''' 
    def defineObstacles(self, r, r2): #defines the obstacle nodes
        
        for l in range(r2):
            for c in range(r):
                if tabGrid[l][c]==-1:
                    GridMap[l][c].OBSTACLE = KNOWN_OBSTACLE
                    
                elif tabGrid[l][c]==-2:
                    GridMap[l][c].OBSTACLE = UNDETECTED_OBSTACLE
                    
                else:
                    GridMap[l][c].OBSTACLE = NO_OBSTACLE       
    
    '''
    Split up the GridMap into n blocks,
    and defines the corresponding block number for each node 
    '''
    def blockSegmentation(self, r, r2, nbrBlock):
        
        blkPerCoord = nbrBlock**(1.0/2)
        subResolu, subResolu2 = [], []
        blkLenght_r = int(r/blkPerCoord)
        blkLenght_r2 = int(r2/blkPerCoord)
        
        for i in range(int(float(str(blkPerCoord)))):
            subResolu.append((blkLenght_r*i, blkLenght_r*(i+1)))
            subResolu2.append((blkLenght_r2*i, blkLenght_r2*(i+1))) 
        
        blk_nbr = 0
        
        for i in subResolu:
            for j in subResolu2:
                for x in range(i[0], i[1]):
                    for y in range(j[0], j[1]):
                        GridMap[y][x].block = blk_nbr
                blk_nbr+=1    
        print 'blk per axe:', len(subResolu), len(subResolu2), 'total blk:', blk_nbr
        return blk_nbr
                  
    '''
        SEQ is a variable which defines if the algorithm is parallelized or not:
            seq=0: one block execution (default value)
            seq=1: use sequenced list of nodes (shorter lenghth to reduce the computation time)
            seq=2: use multi threading (not tested yet)
            seq=3: use multiprocesses on a CPU (not finished)
            seq=4: use multiprocesses on a GPU (not finished)
            nbr_blocks>=0: is the number of blocks that you want to use to parallelize the computation (default: do not split up the map)
                                                                    
    '''
    def processMap(self, r, r2, binaryMap, seq = 0, nbr_blocks=0):
        global tabGrid, GridMap, nodeList
        tabGrid = [[0 for _ in range(r)] for _ in range(r2)]
        GridMap = [[0 for _ in range(r)] for _ in range(r2)]
        num, nodeList, srcObstacles = 0, [], []
        print('--> MAP processing...')
    
        '''
        assign node for each cell 
        '''
        for i in range(r2):
            for j in range(r):
                GridMap[i][j] = node(i, j)
        
        self.readMap(r,r2,binaryMap)
        
        self.defineObstacles(r,r2)
        
        blk_nbr = 0 
        if seq == 1:
            blk_nbr = self.blockSegmentation(r, r2, nbr_blocks)
        '''
        define the obstacle nodes and index them & return the indexed list of all map's nodes
        '''
        for i in range(r2):
            for j in range(r):
                GridMap[i][j].indice = num
                
                if GridMap[i][j].OBSTACLE ==KNOWN_OBSTACLE:
                    GridMap[i][j].TAG = FORBIDDEN
                    srcObstacles.append(GridMap[i][j].indice)
                nodeList.append(GridMap[i][j])
                num+=1
        
        print('                    MAP ready !')
        
        return nodeList, srcObstacles, blk_nbr
         
class MAP_3D(object):
    
    def MiniMAp3D(self):
        
        #obs1
        for i in range(0,10):
            for j in range(10,15):
                for k in range(5,15):
                    tabGrid[i][j][k] =-1
        #obs1
        for i in range(0,30):
            for j in range(20,30):
                for k in range(20,30):
                    tabGrid[i][j][k] =-1
                #obs1
        for i in range(0,38):
            for j in range(5,15):
                for k in range(35,40):
                    tabGrid[i][j][k] =-1       
        #contours
        for i in range(0,depth):
            for j in range(0,height):
                    tabGrid[i][j][0] =-2
        for i in range(0,depth):
            for j in range(width):
                    tabGrid[i][0][j] =-2
        for i in range(0,height):
            for j in range(0,width):
                    tabGrid[0][i][j] =-2
        #contours
        for i in range(0,depth):
            for j in range(0,height):
                    tabGrid[i][j][width-1] =-2
        for i in range(0,depth):
            for j in range(width):
                    tabGrid[i][height-1][j] =-2
        for i in range(0,height):
            for j in range(0,width):
                    tabGrid[depth-1][i][j] =-2
    
    def lesObstacles3D(self):
        
        #obs1
        for i in range(0,40):
            for j in range(20,35):
                for k in range(25,35):
                    tabGrid[i][j][k] =-1
        
        #contours
        for i in range (depth):
            for j in range(height):
                tabGrid[i][j][0] = -1
                #tabGrid[i][j][width-1] = -1
        
        for i in range (height):
            for j in range (width):
                tabGrid[0][i][j] = -1
                #tabGrid[depth-1][i][j] = -1 
                
        for i in range (depth):
            for j in range (width):
                tabGrid[i][0][j] = -1
                #tabGrid[i][height-1][j] = -1 
                                         
        
        
        for i in range(0,35):
            for j in range(45,55):
                for k in range(25,35):
                    tabGrid[i][j][k] =-1#ligne(axe Y), colonnes(axe X), profondeur(axe Z)
        
        for i in range(0,15):
            for j in range(35,55):
                for k in range(45,65):
                    tabGrid[i][j][k] =-1
                    
        for i in range(0,21):
            for j in range(65,75):
                for k in range(25,55):
                    tabGrid[i][j][k] =-1#ligne(axe Y), colonnes(axe X), profondeur(axe Z)

        for i in range(0,45):
            for j in range(85,95):
                for k in range(25,35):
                    tabGrid[i][j][k] =-1#ligne(axe Y), colonnes(axe X), profondeur(axe Z)
        
        #obs2

        for i in range(0,25):
            for j in range(25,40):
                for k in range(75,95):
                    tabGrid[i][j][k] =-1#ligne(axe Y), colonnes(axe X), profondeur(axe Z)
        
        
        for i in range(0,30):
            for j in range(50,6):
                for k in range(80,90):
                    tabGrid[i][j][k] =-1#ligne(axe Y), colonnes(axe X), profondeur(axe Z)
        
        for i in range(0,25):
            for j in range(65,75):
                for k in range(70,80):
                    tabGrid[i][j][k] =-1#ligne(axe Y), colonnes(axe X), profondeur(axe Z)
        
        for i in range(0,20):
            for j in range(85,95):
                for k in range(85,95):
                    tabGrid[i][j][k] =-1#ligne(axe Y), colonnes(axe X), profondeur(axe Z)
        #obs3
        
        for i in range(0,35):
            for j in range(30,40):
                for k in range(20,30):
                    tabGrid[i][j][k] =-1#ligne(axe Y), colonnes(axe X), profondeur(axe Z)
        
        
        for i in range(0,15):
            for j in range(75,80):
                for k in range(35,45):
                    tabGrid[i][j][k] =-1#ligne(axe Y), colonnes(axe X), profondeur(axe Z)
        
        for i in range(0,35):
            for j in range(85,95):
                for k in range(50,65):
                    tabGrid[i][j][k] =-1#ligne(axe Y), colonnes(axe X), profondeur(axe Z)

        for i in range(0,29):
            for j in range(85,95):
                for k in range(70,85):
                    tabGrid[i][j][k] =-1#ligne(axe height), colonnes(axe X), profondeur(axe Z)

    def defineObstacles3D(self):
     
        for l in range(depth):
            for c in range(height):
                for d in range(width):
                    
                    if tabGrid[l][c][d]==-1:
                        
                        GridMap[l][c][d].OBSTACLE = KNOWN_OBSTACLE
                        
                    elif tabGrid[l][c][d]==-2:
                        GridMap[l][c][d].OBSTACLE = UNDETECTED_OBSTACLE
                    
                    elif tabGrid[l][c][d]==-5:
                        GridMap[l][c][d].OBSTACLE = CONTOUR                        
                                
                    else:
                        GridMap[l][c][d].OBSTACLE = NO_OBSTACLE  
    
    def segmention(self, x, y, z,nbrBlk):
        
        blkx, blky, blkz = nbrBlk/(x), nbrBlk/(y/2), nbrBlk/(z/2)       
        #blkx, blky, blkz = nbrBlk/(x), nbrBlk/(y), nbrBlk/(z)#1000
        #blkx, blky, blkz = 1,2,2
        #length of blocks
        subx, suby, subz = x/blkx,y/blky,z/blkz
        X,Y,Z = [],[],[]
        Nodelist = [[] for _ in range(nbrBlk)]
        print( blkx, blky, blkz ,subx, suby, subz )

        #time.sleep(5000)
        for i in range(blkx):
            X.append((subx*i,subx*(1+i)))
        #print X
        for i in range(blky):
            Y.append((suby*i,suby*(1+i)))
        #print Y
        for i in range(blkz):    
            Z.append((subz*i,subz*(1+i)))
        nbr = 0            
        for i in X:
            for j in Y:
                for k in Z:
                    #co = 0
                    for a in range(i[0],i[1]):
                        for b in range(j[0],j[1]):
                            for c in range(k[0],k[1]):

                                
                                if i[0] > 0 :
                                    GridMap[i[0]][b][c].OBSTACLE = CONTOUR
                                    GridMap[i[0]][b][c].ghost = True
                                if i[1]<x:
                                    GridMap[i[1]-1][b][c].OBSTACLE = CONTOUR
                                    GridMap[i[1]-1][b][c].ghost = True
                                if j[0] > 0 :
                                    GridMap[a][j[0]][c].OBSTACLE = CONTOUR
                                    GridMap[a][j[0]][c].ghost = True
                                if j[1]<y:
                                    GridMap[a][j[1]-1][c].OBSTACLE = CONTOUR
                                    GridMap[a][j[1]-1][c].ghost = True
                                if k[0] > 0 :
                                    GridMap[a][b][k[0]].OBSTACLE = CONTOUR
                                    GridMap[a][b][k[0]].ghost = True
                                if k[1]<y:    
                                    GridMap[a][b][k[1]-1].OBSTACLE = CONTOUR
                                    GridMap[a][b][k[1]-1].ghost = True

                                #print('abc ',a,b,c,' nbr ', nbr)  
                                GridMap[a][b][c].block = nbr
                                
                                #GridMap[a][b][c].idx,GridMap[a][b][c].sndidx  = co,co
                                #GridMap[a][b][c].sndidx  = co
                                
                                #co+=1
                    #print(nbr, co)
                    nbr+=1
        
        directionblock = [[0,0,0] for _ in range(nbrBlk)]
        
        ghostnodes = [[] for _ in range(nbrBlk)]
   
        nbr = 0 
        Blft, Brgt, Bfrwd, Bbckrd, Bup, Bdwn = False, False,False, False,False, False
        for i in X:
            for j in Y:
                for k in Z: 
                    co = 0
                    lft, rgt, frwd, bckrd, up, dwn = nbr-1,nbr+1,nbr-d_parallel[0], nbr+d_parallel[0], nbr-d_parallel[1],nbr+d_parallel[1]
                    Blft, Brgt, Bfrwd, Bbckrd, Bup, Bdwn = False, False,False, False,False, False
                    for a in range(i[0],i[1]):
                        for b in range(j[0],j[1]):
                            for c in range(k[0],k[1]):
                                
                                #left right
                                
                                if lft>=0 and int(floor(lft/d_parallel[0]))==int(floor(nbr/d_parallel[0])):
            
                                    if k[0] > 0 and c == k[0]:
 
                                        if GridMap[a][b][c-1].ghost:
                                            #print GridMap[a][b][c-1].idx 
                                            ghost = deepcopy(GridMap[a][b][c-1])#
                                            ghost.idx = co
                                            
                                            #print 'test Gridmap and ghosts and ghost ',nbr,id(GridMap[a][b][c-1]), id(ghost)
     
                                            ghostnodes[nbr].append(ghost) 
                                            Nodelist[nbr].append(ghost)
                                            Blft = True
                                            co+=1
                        
                                if rgt<d_parallel[2] and (int(floor(rgt/d_parallel[0]))==int(floor(nbr/d_parallel[0]))):
                        
                                    if k[1] < d_parallel[1] and c +1 ==k[1]:
                         
                                        if GridMap[a][b][k[1]].ghost:

                                            
                                            ghost = deepcopy(GridMap[a][b][c+1])
                                            ghost.idx = co
                                            ghostnodes[nbr].append(ghost) 
                                            co+=1
                                            Nodelist[nbr].append(ghost)
                                            Brgt = True
                                                                                        
                                # front back side
                                if frwd>=0 and (int(floor(frwd/d_parallel[1]))==int(floor(nbr/d_parallel[1]))):
                      
                                    if j[0] > 0 and b==j[0]:
                             
                                        if GridMap[a][j[0]-1][c].ghost:
                                            
                                            ghost = deepcopy(GridMap[a][b-1][c])
                                            ghost.idx = co
                                            ghostnodes[nbr].append(ghost)
                                            co+=1
                                            #print ' frwd  <<',ghost.indice,GridMap[a][b][c].indice,'>>'
                                            Nodelist[nbr].append(ghost)
                                                    
                                            Bfrwd =True                                    
                                            
                                if bckrd<d_parallel[2] and (int(floor(bckrd/d_parallel[1]))==int(floor(nbr/d_parallel[1]))):            
                    
                                    if j[1] < d_parallel[1] and b+1==j[1]:

                                        if GridMap[a][j[1]][c].ghost:

                                            
                                            ghost = deepcopy(GridMap[a][b+1][c])
                                            ghost.idx = co
                                            ghostnodes[nbr].append(ghost)
                                            co+=1
                                            #print ' bck  <<',ghost.indice,GridMap[a][b+1][c].indice,'>>', ghost.x, ghost.y, ghost.z,' --',GridMap[a][b+1][c].x,GridMap[a][b+1][c].y,GridMap[a][b+1][c].z
                                            Nodelist[nbr].append(ghost)
                                            
                                            Bbckrd = True                                            
                                            
                                #up down neighbor                    
                                if up>=0 and (int(floor(up/d_parallel[2]))==int(floor(nbr/d_parallel[2]))): 
                                    if i[0] > 0 and a ==i[0]:
                            
                                        if GridMap[i[0]-1][b][c].ghost:
                                            ghost = deepcopy(GridMap[a-1][b][c])
                                            ghost.idx = co                                            
                                            ghostnodes[nbr].append(ghost)
                                            co+=1
                                            #print ' up  <<',ghost.indice,GridMap[a-1][b][c].indice,'>>', ghost.x, ghost.y, ghost.z,' --',GridMap[a-1][b][c].x,GridMap[a-1][b][c].y,GridMap[a-1][b][c].z
                                            Nodelist[nbr].append(ghost)
                                            
                                            Bup = True

                                if dwn <=d_parallel[2] and (int(floor(dwn/d_parallel[2]))==int(floor(nbr/d_parallel[2]))):             
                                    
                                    if i[1] < d_parallel[1] and a+1==i[1]:

                                        if GridMap[i[1]][b][c].ghost:
                                            ghost = deepcopy(GridMap[a+1][b][c])
                                            ghost.idx = co
                                            ghostnodes[nbr].append(ghost)                    
                                            #print'i1 ',co
                                            #print ' dwn  <<',ghost.indice,GridMap[a][b][c].indice,'>>', ghost.x, ghost.y, ghost.z,' --',GridMap[a][b][c].x,GridMap[a][b][c].y,GridMap[a][b][c].z 
                                            Nodelist[nbr].append(ghost)
                                            co+=1
                                            Bdwn = True
                                
                                #print'nbr c ',nbr, co
                                GridMap[a][b][c].idx = co
                                co+=1
                                
                                Nodelist[nbr].append(GridMap[a][b][c]) 
                                #GridMap[a][b][c].sndidx = len(Nodelist[nbr]-1)          
                    #ghostnodes.append(directionblock)
                    if Blft:
                        directionblock[nbr][0]+=1
                    if Brgt:
                        directionblock[nbr][0]+=1
                    if Bup:
                        directionblock[nbr][2]+=1
                    if Bdwn:
                        directionblock[nbr][2]+=1
                    if Bbckrd:
                        directionblock[nbr][1]+=1
                    if Bfrwd:
                        directionblock[nbr][1]+=1
                        
                    nbr+=1
        for i in range(nbrBlk):
            for nds in ghostnodes[i]:
                nds.sndidx = GridMap[nds.x][nds.y][nds.z].idx
            
           
        del  X,Y,Z

        return  ghostnodes, directionblock, Nodelist,subz, suby, subx
       
    def processMap_3D(self, seq = 0, nbrBlk=0):       
        global tabGrid, GridMap, nodeList        
        tabGrid = [[[0 for _ in range(width)] for _ in range(height)] for _ in range(depth)]  
        GridMap = [[[0 for _ in range(width)] for _ in range(height)] for _ in range(depth)]
        
        num,num2, nodeList, srcObstacles = 0,0, [], []
        print('--> MAP processing...')
    
        '''
        assign node for each cell 
        '''
        for d in range(depth):
            for h in range(height):
                for w in range(width):
                    GridMap[d][h][w] = node3(d, h, w)
        
        self.lesObstacles3D()
        self.defineObstacles3D()
        
        if seq==1:
            #import multiprocessing as mp
            nbrBlk =500#mp.cpu_count()
                    
            
       
        for d in range(depth):
            for h in range(height):
                for w in range(width):
   
                    GridMap[d][h][w].indice = num
                    
                    if GridMap[d][h][w].OBSTACLE == KNOWN_OBSTACLE:
                        GridMap[d][h][w].TAG = FORBIDDEN
                        srcObstacles.append(GridMap[d][h][w].indice)
                        
                    nodeList.append(GridMap[d][h][w])
                        
                    #Nodelist[GridMap[d][h][w].block].append(GridMap[d][h][w])
        
                    num+=1
        ghostzones, directionBlk,Nodelist,dp,he,wi = self.segmention(depth, height, width, nbrBlk)
        #GhostNodes = self.ghostSets(nbrBlk,Nodelist,x,y,z)            
        del GridMap, tabGrid
        print('               MAP ready')#, len(srcObstacles))

        return ghostzones,directionBlk,Nodelist, nodeList, srcObstacles,nbrBlk,dp,he,wi


    """
    def segmention(self, x, y, z,nbrBlk):
        
        blkx, blky, blkz = nbrBlk/(x), nbrBlk/(y/2), nbrBlk/(z/2)       
        #blkx, blky, blkz = nbrBlk/(x), nbrBlk/(y), nbrBlk/(z)#1000
        #blkx, blky, blkz = 1,2,2
        #length of blocks
        subx, suby, subz = x/blkx,y/blky,z/blkz
        X,Y,Z = [],[],[]
        nbr = 0
        print( blkx, blky, blkz ,subx, suby, subz )

        #time.sleep(5000)
        for i in range(blkx):
            X.append((subx*i,subx*(1+i)))
        #print X
        for i in range(blky):
            Y.append((suby*i,suby*(1+i)))
        #print Y
        for i in range(blkz):    
            Z.append((subz*i,subz*(1+i)))
            
        for i in X:
            for j in Y:
                for k in Z:
                    co = 0
                    for a in range(i[0],i[1]):
                        for b in range(j[0],j[1]):
                            for c in range(k[0],k[1]):
                                if i[0] > 0 :
                                    GridMap[i[0]][b][c].OBSTACLE = CONTOUR
                                    GridMap[i[0]][b][c].ghost = True
                                if i[1]<x:
                                    GridMap[i[1]-1][b][c].OBSTACLE = CONTOUR
                                    GridMap[i[1]-1][b][c].ghost = True
                                if j[0] > 0 :
                                    GridMap[a][j[0]][c].OBSTACLE = CONTOUR
                                    GridMap[a][j[0]][c].ghost = True
                                if j[1]<y:
                                    GridMap[a][j[1]-1][c].OBSTACLE = CONTOUR
                                    GridMap[a][j[1]-1][c].ghost = True
                                if k[0] > 0 :
                                    GridMap[a][b][k[0]].OBSTACLE = CONTOUR
                                    GridMap[a][b][k[0]].ghost = True
                                if k[1]<y:    
                                    GridMap[a][b][k[1]-1].OBSTACLE = CONTOUR
                                    GridMap[a][b][k[1]-1].ghost = True

                                #print('abc ',a,b,c,' nbr ', nbr)  
                                GridMap[a][b][c].block = nbr
                                GridMap[a][b][c].idx,GridMap[a][b][c].sndidx  = co,co
                                
                                co+=1
                    #print(nbr, co)
                    nbr+=1
        
        directionblock = [[0,0,0] for _ in range(nbrBlk)]
        
        ghostnodes = [[] for _ in range(nbrBlk)]
   
        nbr = 0 

        for i in X:
            for j in Y:
                for k in Z: 
                    GridMap[a][b][c].idx = co
                    lft, rgt, frwd, bckrd, up, dwn = nbr-1,nbr+1,nbr-d_parallel[0], nbr+d_parallel[0], nbr-d_parallel[1],nbr+d_parallel[1]
                    
                    #left right
                    if lft>=0 and int(floor(lft/d_parallel[0]))==int(floor(nbr/d_parallel[0])):
            
                        if k[0] > 0 :

                            for a in range(i[0],i[1]):
                                for c in range(j[0],j[1]):
                                    
                                    if GridMap[a][c][k[0]-1].ghost:
                                        #print('<<nbr',nbr,' blk', GridMap[a][c][k[0]-1].block,' a b c', GridMap[a][c][k[0]-1].x, GridMap[a][c][k[0]-1].y, GridMap[a][c][k[0]-1].z,'(',GridMap[a][c][k[0]-1].indice,GridMap[a][c][k[0]].indice,')',GridMap[a][c][k[0]].block, GridMap[a][c][k[0]].x, GridMap[a][c][k[0]].y, GridMap[a][c][k[0]].z)
                                        ghostnodes[nbr].append(GridMap[a][c][k[0]-1]) 
                            directionblock[nbr][0]+=1
                            
                    if rgt<d_parallel[2] and (int(floor(rgt/d_parallel[0]))==int(floor(nbr/d_parallel[0]))):
                        
                        if k[1] < d_parallel[1]:
                            
                            for a in range(i[0],i[1]):
                                for c in range(j[0],j[1]):
                                                                
                                    if GridMap[a][c][k[1]].ghost:
                                        #print('>> nbr',nbr,' blk, a b c-', GridMap[a][c][k[1]].block, GridMap[a][c][k[1]].x, GridMap[a][c][k[1]].y, GridMap[a][c][k[1]].z,'(',GridMap[a][c][k[1]].indice,GridMap[a][c][k[1]-1].indice,')',GridMap[a][c][k[1]-1].block, GridMap[a][c][k[1]-1].x, GridMap[a][c][k[1]-1].y, GridMap[a][c][k[1]-1].z)
                                        ghostnodes[nbr].append(GridMap[a][c][k[1]])  
                            directionblock[nbr][0]+=1
                    # front back side
                    if frwd>=0 and (int(floor(frwd/d_parallel[1]))==int(floor(nbr/d_parallel[1]))):
                      
                        if j[0] > 0 :
                            for a in range(i[0],i[1]):
                                for c in range(k[0],k[1]):                            
                                    if GridMap[a][j[0]-1][c].ghost:
                                        #print('frnt nbr blk, a b c',nbr,'-', GridMap[a][j[0]-1][c].block, GridMap[a][j[0]-1][c].x, GridMap[a][j[0]-1][c].y, GridMap[a][j[0]-1][c].z,'(',GridMap[a][j[0]-1][c].indice,GridMap[a][j[0]][c].indice,')',GridMap[a][j[0]][c].block, GridMap[a][j[0]][c].x, GridMap[a][j[0]][c].y, GridMap[a][j[0]][c].z)
                                        ghostnodes[nbr].append(GridMap[a][j[0]-1][c]) 
                            directionblock[nbr][1]+=1
                    if bckrd<d_parallel[2] and (int(floor(bckrd/d_parallel[1]))==int(floor(nbr/d_parallel[1]))):            
                    
                        if j[1] < d_parallel[1] :
                            
                            for a in range(i[0],i[1]):
                                for c in range(k[0],k[1]):
                                    
                                    if GridMap[a][j[1]][c].ghost:
                                        #print('>>nbr blk, a b c',nbr,'-', GridMap[a][j[1]][c].block, GridMap[a][j[1]][c].x, GridMap[a][j[1]][c].y, GridMap[a][j[1]][c].z,'(',GridMap[a][j[1]][c].indice,GridMap[a][j[1]-1][c].indice,')',GridMap[a][j[1]-1][c].block, GridMap[a][j[1]-1][c].x, GridMap[a][j[1]-1][c].y, GridMap[a][j[1]-1][c].z)
                                        ghostnodes[nbr].append(GridMap[a][j[1]][c])  
                            directionblock[nbr][1]+=1
                    #up down neighbor                    
                    if up>=0 and (int(floor(up/d_parallel[2]))==int(floor(nbr/d_parallel[2]))): 
                        if i[0] > 0 :
                            
                            for b in range(j[0],j[1]):
                                for c in range(k[0],k[1]):                            
                                    if GridMap[i[0]-1][b][c].ghost:
                                        #print('nbr blk, a b c',nbr,'-', GridMap[i[0]-1][b][c].block, GridMap[i[0]-1][b][c].x, GridMap[i[0]-1][b][c].y, GridMap[i[0]-1][b][c].z,'(',GridMap[i[0]-1][b][c].indice,GridMap[i[0]][b][c].indice,')',GridMap[i[0]][b][c].block, GridMap[i[0]][b][c].x, GridMap[i[0]][b][c].y, GridMap[i[0]][b][c].z)
                                        ghostnodes[nbr].append(GridMap[i[0]-1][b][c]) 
                            directionblock[nbr][2]+=1

                    if dwn <=d_parallel[2] and (int(floor(dwn/d_parallel[2]))==int(floor(nbr/d_parallel[2]))):             
                        if i[1] < d_parallel[1] :
                            
                            for b in range(i[0],i[1]):
                                for c in range(k[0],k[1]):
                                    if GridMap[i[1]][b][c].ghost:
                                        #print('>>nbr blk, a b c',nbr,'-', GridMap[i[1]][b][c].block, GridMap[i[1]][b][c].x, GridMap[i[1]][b][c].y, GridMap[i[1]][b][c].z,'(',GridMap[i[1]][b][c].indice,GridMap[i[1]-1][b][c].indice,')',GridMap[i[1]][b][c].block, GridMap[i[1]][b][c].x, GridMap[i[1]][b][c].y, GridMap[i[1]][b][c].z)
                                        ghostnodes[nbr].append(GridMap[i[1]][b][c])                     
                            directionblock[nbr][2]+=1
                    #ghostnodes.append(directionblock)
                    
                    nbr+=1
        
        del  X,Y,Z

        return  ghostnodes, directionblock
    """





