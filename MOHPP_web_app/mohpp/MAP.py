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
import multiprocessing as mp
from mohpp.utilities import sqrt_dist
tabGrid, GridMap, nodeList,Nodelist,ghostzones = [],[],[],[],[]

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
                    tabGrid[i][j][0] =-1
        for i in range(0,depth):
            for j in range(width):
                    tabGrid[i][0][j] =-1
        for i in range(0,height):
            for j in range(0,width):
                    tabGrid[0][i][j] =-1
        #contours
        for i in range(0,depth):
            for j in range(0,height):
                    tabGrid[i][j][width-1] =-1
        for i in range(0,depth):
            for j in range(width):
                    tabGrid[i][height-1][j] =-1
        for i in range(0,height):
            for j in range(0,width):
                    tabGrid[depth-1][i][j] =-1
   
    def accuracy_Test(self):
        
        #tabGrid[0][0][0] = -1
        tabGrid[depth/2-1][height/2-1][width/2 -1] = -1
        #tabGrid[depth/2-1][height/2-1][width/2 -1] = -1
        
    def TinyMAP3D(self):
        #tabGrid[0][0][0] =-1
        
        #obs1
        for i in range(1,3):
            for j in range(5,7):
                for k in range(4,6):
                    tabGrid[i][j][k] =-1     
        
        #d[1] up to down 
        #contours
        for i in range(0,depth):
            for j in range(0,height):
                    tabGrid[i][j][0] =-1
        for i in range(0,depth):
            for j in range(0,height):
                    tabGrid[i][j][width-1] =-1                    
        
        #left to right 
        for i in range(0,depth):
            for j in range(width):
                    tabGrid[i][0][j] =-1
        for i in range(0,depth):
            for j in range(width):
                    tabGrid[i][height-1][j] =-1                    
                  
        for i in range(0,height):
            for j in range(0,width):
                    tabGrid[0][i][j] =-1
        for i in range(0,height):
            for j in range(0,width):
                    tabGrid[depth-1][i][j] =-1
    
    def lesObstacles3d(self):
        
        #obstacle 1
        for k in range(45):
            for i in range(40,80):
                for j in range(0,5):
                    tabGrid[j][i][k] =-1
                    
        for k in range(65):                   
            for i in range(40,80):
                for j in range(55,60):
                    tabGrid[j][i][k] =-1
        for k in range(25):            
            for i in range(75,80):
                for j in range(5,55):
                    tabGrid[j][i][k] =-1

        
        #obstacle 2
        
        for i in range(40,80):
            for j in range(80,85):
                tabGrid[j][i][1] =-1        
        for i in range(40,80):
            for j in range(45,50):
                tabGrid[j][i][1] =-1
        for i in range(80,45):
            for j in range(75,80):
                tabGrid[i][j][1] =-1
    
        #obstacle 3
        for i in range(20,60):
            for j in range(30,35):
                tabGrid[j][i][1] =-1        
        for i in range(20,60):
            for j in range(10,15):
                tabGrid[j][i][1] =-1
        for i in range(5,10):
            for j in range(20,25):
                tabGrid[i][j][1] =-1
        
        #excptected obstacle
        
        for i in range(55,95):
            for j in range(30,35):
                tabGrid[i][j][3] =-2
                
        for i in range(90,95):        
            for j in range(10,35):
                tabGrid[i][j][4] =-2    
        for i in range(50,55):        
            for j in range(10,35):
                tabGrid[i][j][4] =-2  

    def lesObstacles3D(self):
        
        #obs1
        for i in range(0,40):
            for j in range(20,35):
                for k in range(25,35):
                    tabGrid[i][j][k] =-1
        
        #contours
        for i in range (depth):
            for j in range(height):
                tabGrid[i][j][0] = -5
                tabGrid[i][j][width-1] = -5
        
        for i in range (height):
            for j in range (width):
                tabGrid[0][i][j] = -5
                tabGrid[depth-1][i][j] = -5 
                
        for i in range (depth):
            for j in range (width):
                tabGrid[i][0][j] = -5
                tabGrid[i][height-1][j] = -5 
                                         
        
        
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
        #blkx, blky, blkz = 2,2,2
        #length of blocks
        subx, suby, subz = x/blkx,y/blky,z/blkz
        X,Y,Z = [],[],[]
        Nodelist = [[] for _ in range(nbrBlk)]
        print( blkx, blky, blkz ,subx, suby, subz )
        blk_xyz = []
        #time.sleep(5000)
        for i in range(blkx):
            X.append((subx*i,subx*(1+i)))
        #print X
        for i in range(blky):
            Y.append((suby*i,suby*(1+i)))
        #print Y
        for i in range(blkz):    
            Z.append((subz*i,subz*(1+i)))
        nbr,indice = 0,0            
        for db,i in enumerate(X):
            for hb,j in enumerate(Y):
                for wb,k in enumerate(Z):
                    
                    co = 0
                    for a in range(i[0],i[1]):
                        for b in range(j[0],j[1]):
                            for c in range(k[0],k[1]):
                                
                                if i[0] > 0 and GridMap[i[0]][b][c].OBSTACLE != KNOWN_OBSTACLE:
                                    GridMap[i[0]][b][c].ghost = True
                                    
                                if i[1]<x and GridMap[i[1]-1][b][c].OBSTACLE != KNOWN_OBSTACLE:
                                    GridMap[i[1]-1][b][c].ghost = True
                                    
                                if j[0] > 0 and GridMap[a][j[0]][c].OBSTACLE != KNOWN_OBSTACLE:
                                    GridMap[a][j[0]][c].ghost = True
                                    
                                if j[1]<y and GridMap[a][j[1]-1][c].OBSTACLE != KNOWN_OBSTACLE:
                                    GridMap[a][j[1]-1][c].ghost = True
                                    
                                if k[0] > 0 and GridMap[a][b][k[0]].OBSTACLE != KNOWN_OBSTACLE:
                                    GridMap[a][b][k[0]].ghost = True
                                    
                                if k[1]<z and GridMap[a][b][k[1]-1].OBSTACLE != KNOWN_OBSTACLE:    
                                    GridMap[a][b][k[1]-1].ghost = True  
                                    
                                GridMap[a][b][c].block = nbr
                                GridMap[a][b][c].indice = indice
                                indice+=1
                                
                    blk_xyz.append([wb,hb,db])
                    nbr+=1
                    
        blk_d = [wb+1, hb+1, db+1]
        
        d_p = [blk_d[0],blk_d[0]*blk_d[1],blk_d[0]*blk_d[1]*blk_d[2]]
        print blk_d, d_p

        ghostnodes = [[] for _ in range(nbrBlk)]
        neighBlk = [[-1,-1,-1,-1,-1,-1] for _ in range(nbrBlk)]
        nbr = 0 
        d_par = [[0,0,0] for _ in range(nbrBlk)]
        srcObstacles = [[] for _ in range(nbrBlk)]
        
        for ii,i in enumerate(X):
            for jj,j in enumerate(Y):
                for kk,k in enumerate(Z): 
                    co, c0,c2,c1 = 0,0,0,0
                    lft, rgt, frwd, bckrd, up, dwn = nbr-1,nbr+1,nbr-blk_d[0], nbr+blk_d[0], nbr-blk_d[0]*blk_d[1],nbr+blk_d[0]*blk_d[1]
                    
                    c2 = 0
                    for a in range(i[0],i[1]):
                        c1 = 0
                        for b in range(j[0],j[1]):
                            c0 = 0
                            for c in range(k[0],k[1]):
                            
                                #left right
                                if lft>=0 and int(floor(lft/d_p[0]))==int(floor(nbr/d_p[0])):
            
                                    if k[0] > 0 and c == k[0]:
 
                                        if GridMap[a][b][c-1].ghost:
                                            #print GridMap[a][b][c-1].idx 
                                            ghost = deepcopy(GridMap[a][b][c-1])#
                                            ghost.idx = co                                           
                                            ghostnodes[nbr].append(ghost) 
                                            Nodelist[nbr].append(ghost)
                                            co+=1
                                            c0+=1
                                            c1+=1
                                            c2+=1 
                                            neighBlk[nbr][0] =  lft                                                
                        
                                if rgt<d_p[2] and (int(floor(rgt/d_p[0]))==int(floor(nbr/d_p[0]))):
                        
                                    if k[1] < d_p[1] and c +1 ==k[1]:
                         
                                        if GridMap[a][b][k[1]].ghost :
                                            ghost = deepcopy(GridMap[a][b][c+1])
                                            ghost.idx = co                                        
                                            ghostnodes[nbr].append(ghost) 
                                            co+=1
                                            c0+=1
                                            c1+=1
                                            c2+=1                                            
                                            Nodelist[nbr].append(ghost)
                                            neighBlk[nbr][1] = rgt
                                # front back side
                                if frwd>=0 and (int(floor(frwd/d_p[1]))==int(floor(nbr/d_p[1]))):
                      
                                    if j[0] > 0 and b==j[0]:
                             
                                        if GridMap[a][j[0]-1][c].ghost :
                                            
                                            ghost = deepcopy(GridMap[a][b-1][c])
                                            ghost.idx = co                                          
                                            ghostnodes[nbr].append(ghost)
                                            co+=1
                                            c0+=1
                                            c1+=1
                                            c2+=1          
                                            neighBlk[nbr][2] = frwd                                   
                                            #print ' frwd  <<',(a,b,c),(ghost.x,ghost.y,ghost.z),ghost.indice,GridMap[a][b][c].indice,'>>'
                                            Nodelist[nbr].append(ghost)

                                if bckrd<d_p[2] and (int(floor(bckrd/d_p[1]))==int(floor(nbr/d_p[1]))):            
                    
                                    if j[1] < d_p[1] and b+1==j[1]:
                                                                                    
                                        if GridMap[a][j[1]][c].ghost :
                                                                                                                              
                                            ghost = deepcopy(GridMap[a][b+1][c])                                            
                                            ghost.idx = co
                                            ghostnodes[nbr].append(ghost)

                                            co+=1
                                            c0+=1
                                            c1+=1
                                            c2+=1    
                                            neighBlk[nbr][3] = bckrd                                        
                                            #print ' bck  <<',ghost.indice,GridMap[a][b+1][c].indice,'>>', ghost.x, ghost.y, ghost.z,' --',GridMap[a][b+1][c].x,GridMap[a][b+1][c].y,GridMap[a][b+1][c].z
                                            Nodelist[nbr].append(ghost)

                                #up down neighbor                    
                                if up>=0 and (int(floor(up/d_p[2]))==int(floor(nbr/d_p[2]))): 
                                    #print 'testt' ,int(floor(up/d_parallel[2])),int(floor(nbr/d_parallel[2]))
                                    if i[0] > 0 and a ==i[0]:
                            
                                        if GridMap[i[0]-1][b][c].ghost :
                                            ghost = deepcopy(GridMap[a-1][b][c])
                                            ghost.idx = co                                            

                                            ghostnodes[nbr].append(ghost)
                                            co+=1
                                            c0+=1
                                            c1+=1
                                            c2+=1                   
                                            neighBlk[nbr][4] =up                         
                                            Nodelist[nbr].append(ghost)

                                if dwn <=d_p[2] and (int(floor(dwn/d_p[2]))==int(floor(nbr/d_p[2]))):             
                                    
                                    if i[1] < d_p[2] and a+1==i[1]:

                                        if GridMap[i[1]][b][c].ghost :
                                            ghost = deepcopy(GridMap[a+1][b][c])
                                            ghost.idx = co
                                            ghostnodes[nbr].append(ghost)
                                            Nodelist[nbr].append(ghost)
                                            co+=1
                                            c0+=1
                                            c1+=1
                                            c2+=1
                                            neighBlk[nbr][5] =dwn
                                
                                #print'nbr c ',nbr, co
                                GridMap[a][b][c].idx = co
                                co+=1
                                c0+=1
                                c1+=1
                                c2+=1                                
                                if GridMap[a][b][c].OBSTACLE == KNOWN_OBSTACLE:
                                    GridMap[a][b][c].TAG = FORBIDDEN
                                    srcObstacles[nbr].append(GridMap[a][b][c].idx)
                                Nodelist[nbr].append(GridMap[a][b][c]) 
                                
                            d_par[nbr][0] = c0
                        d_par[nbr][1] = c1
                    d_par[nbr][2] = c2            
                        
                    nbr+=1

        for i in range(nbrBlk):
            for nds in ghostnodes[i]:
                nds.sndidx = GridMap[nds.x][nds.y][nds.z].idx
        
        del  X,Y,Z, ghostnodes

        return  neighBlk, d_par, Nodelist,srcObstacles, blk_xyz,blk_d
    
    def BlkNbr(self, x, y, z, nbrBlk):
        blkx, blky, blkz = nbrBlk/(x), nbrBlk/(y/2), nbrBlk/(z/2)       
        #blkx, blky, blkz = nbrBlk/(x), nbrBlk/(y), nbrBlk/(z)#1000
        #blkx, blky, blkz = 2,2,2# for tiny map test 
        #length of blocks
        subx, suby, subz = x/blkx,y/blky,z/blkz
        X,Y,Z = [],[],[]
        Nodelist = [[] for _ in range(nbrBlk)]
        print( blkx, blky, blkz ,subx, suby, subz )
        blk_xyz = []
        #time.sleep(5000)
        for i in range(blkx):
            X.append((subx*i,subx*(1+i)))
        #print X
        for i in range(blky):
            Y.append((suby*i,suby*(1+i)))
        #print Y
        for i in range(blkz):    
            Z.append((subz*i,subz*(1+i)))
        nbr=0            
        for db,i in enumerate(X):
            for hb,j in enumerate(Y):
                for wb,k in enumerate(Z):
                    for a in range(i[0],i[1]):
                        for b in range(j[0],j[1]):
                            for c in range(k[0],k[1]):
                                
                                GridMap[a][b][c].block = nbr
                    nbr+=1        
               
    def processMap_3D(self, seq = 0, nbrBlk=0):       
        global tabGrid, GridMap, nodeList        
        tabGrid = [[[0 for _ in range(width)] for _ in range(height)] for _ in range(depth)]  
        GridMap = [[[0 for _ in range(width)] for _ in range(height)] for _ in range(depth)]
        
        num, nodeList, srcObstacles = 0, [],[]
        print('--> MAP processing...')
    
        '''
        assign node for each cell 
        '''
        for d in range(depth):
            for h in range(height):
                for w in range(width):
                    GridMap[d][h][w] = node3(d, h, w)
        
        self.accuracy_Test()
        self.defineObstacles3D()
      
        if seq ==0:#MSFM or FM
            
            for d in range(depth):
                for h in range(height):
                    for w in range(width):    
                        GridMap[d][h][w].indice = num
                        nodeList.append(GridMap[d][h][w])            
                        if GridMap[d][h][w].OBSTACLE == KNOWN_OBSTACLE or GridMap[d][h][w].OBSTACLE == CONTOUR:
                            GridMap[d][h][w].TAG = FORBIDDEN
                            srcObstacles.append(GridMap[d][h][w].indice)            
                        num+=1
                        
            j = nodeList[srcObstacles[0]]
            
            print 'nodeobs', j.x,j.y,j.z
            for i in nodeList:
                i.anltic  = sqrt_dist((j.x-i.x),(j.y-i.y),(j.z-i.z))
                #print (i.x, i.y, i.z),'-',(j.x,j.y,j.z),'=', sqrt_dist((j.x-i.x),(j.y-i.y),(j.z-i.z))              
            del GridMap, tabGrid
            print'               MAP ready', num, len(nodeList), len(srcObstacles)#, len(srcObstacles)
    
            return [],None, nodeList,[], srcObstacles,0,0,None                                    

        elif seq==1:#VZLHGS
            nbrBlk = 500#8#
            self.BlkNbr(depth, height, width, nbrBlk)
            
            for d in range(depth):
                for h in range(height):
                    for w in range(width):    
            
                        GridMap[d][h][w].indice = num
                        nodeList.append(GridMap[d][h][w])            
                        
                        if GridMap[d][h][w].OBSTACLE == KNOWN_OBSTACLE:
                            GridMap[d][h][w].TAG = FORBIDDEN
                            srcObstacles.append(GridMap[d][h][w].indice)            
                        num+=1
                        
            j = nodeList[srcObstacles[0]]
            print 'nodeobs', j.x,j.y,j.z
            for i in nodeList:
                i.anltic  = sqrt_dist((j.x-i.x),(j.y-i.y),(j.z-i.z))
                #print (i.x, i.y, i.z),'-',(j.x,j.y,j.z),'=', sqrt_dist((j.x-i.x),(j.y-i.y),(j.z-i.z))     
                #i.anltic2 = pow((i.x-j.x),2)*1/100+1/20*pow((i.y-j.y),2)+1/20*pow((i.z-j.z),2) 
                                      
            del GridMap, tabGrid
            return [],None, nodeList,[], srcObstacles,nbrBlk,0,None             
        
        elif seq==2:
            #import multiprocessing as mp
            nbrBlk =500#2*mp.cpu_count()#500#mp.cpu_count()
            ghostzones, d_par,Nodelist,srcObstacles,blkxyz,blk_d = self.segmention(depth, height, width, nbrBlk)
            #print#[Nodelist[0][i].block for i in range(len(Nodelist[0]))]

             
            for d in range(depth):
                for h in range(height):
                    for w in range(width):    
                        nodeList.append(GridMap[d][h][w])
            
            del GridMap, tabGrid
            print('               MAP ready')#, len(srcObstacles))
        
            return ghostzones,d_par,Nodelist, nodeList, srcObstacles,nbrBlk,blkxyz,blk_d
                                
        
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





