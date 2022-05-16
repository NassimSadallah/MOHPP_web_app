'''
Created on Jul 8, 2021

@author: nassim
'''
import time
from PIL import Image
from Node import node, node3
from utilities import height, width, depth, FORBIDDEN, NO_OBSTACLE, KNOWN_OBSTACLE, UNDETECTED_OBSTACLE

tabGrid, GridMap, nodeList = [],[],[]



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
                tabGrid[i][j][0] = -2
                tabGrid[i][j][width-1] = -2
        
        for i in range (height):
            for j in range (width):
                tabGrid[0][i][j] = -2
                tabGrid[depth-1][i][j] = -2 
                
            for i in range (depth):
                for j in range (width):
                    tabGrid[i][0][j] = -2
                    tabGrid[i][height-1][j] = -2 
                                             
        
        
        for i in range(0,35):
            for j in range(45,55):
                for k in range(25,35):
                    tabGrid[i][j][k] =-1#ligne(axe Y), colonnes(axe X), profondeur(axe Z)
        
        for i in range(0,45):
            for j in range(35,55):
                for k in range(45,65):
                    tabGrid[i][j][k] =-1
                    
        for i in range(0,40):
            for j in range(65,75):
                for k in range(25,55):
                    tabGrid[i][j][k] =-1#ligne(axe Y), colonnes(axe X), profondeur(axe Z)

        for i in range(0,45):
            for j in range(85,95):
                for k in range(25,35):
                    tabGrid[i][j][k] =-1#ligne(axe Y), colonnes(axe X), profondeur(axe Z)
        """
        #obs2

        for i in range(0,25):
            for j in range(25,40):
                for k in range(75,95):
                    tabGrid[i][j][k] =-1#ligne(axe Y), colonnes(axe X), profondeur(axe Z)
        
        
        for i in range(0,30):
            for j in range(50,6):
                for k in range(80,90):
                    tabGrid[i][j][k] =-1#ligne(axe Y), colonnes(axe X), profondeur(axe Z)
        
        for i in range(0,45):
            for j in range(65,75):
                for k in range(70,80):
                    tabGrid[i][j][k] =-1#ligne(axe Y), colonnes(axe X), profondeur(axe Z)
        
        for i in range(0,40):
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

        for i in range(0,19):
            for j in range(85,95):
                for k in range(70,85):
                    tabGrid[i][j][k] =-1#ligne(axe height), colonnes(axe X), profondeur(axe Z)
        """
    def defineObstacles3D(self):
     
        for l in range(depth):
            for c in range(height):
                for d in range(width):
                    
                    if tabGrid[l][c][d]==-1:
                        GridMap[l][c][d].OBSTACLE = KNOWN_OBSTACLE
                    
                    elif tabGrid[l][c][d]==-2:
                        GridMap[l][c][d].OBSTACLE = UNDETECTED_OBSTACLE
            
                    else:
                        GridMap[l][c][d].OBSTACLE = NO_OBSTACLE  
    
    def segmention(self, x, y, z,nbrBlk):
        
        blkx, blky, blkz = nbrBlk/(x*20), nbrBlk/(2*y), nbrBlk/(2*z)
        
        #length of blocks
        subx, suby, subz = x/blkx,y/blky,z/blkz
        X,Y,Z = [],[],[]
        nbr = 0
        print( blkx, blky, blkz ,subx, suby, subz )
        #time.sleep(5)
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
                                co+=1
                                #print('abc ',a,b,c,' nbr ', nbr)  
                                GridMap[a][b][c].block = nbr     
                    #print(nbr, co)
                    nbr+=1          
        #print(nbr)
        time.sleep(1000)
        return nbr
    
    def processMap_3D(self, seq = 0, nbrBlk=0):       
        global tabGrid, GridMap, nodeList        
        tabGrid = [[[0 for _ in range(width)] for _ in range(height)] for _ in range(depth)]  
        GridMap = [[[0 for _ in range(width)] for _ in range(height)] for _ in range(depth)]
        
        num, nodeList, srcObstacles = 1, [], []
        print('--> MAP processing...',len(GridMap))
    
        '''
        assign node for each cell 
        '''
        for d in range(depth):
            for h in range(height):
                for w in range(width):
                    GridMap[d][h][w] = node3(d, h, w)
        
        self.MiniMAp3D()
        self.defineObstacles3D()
        
        if seq==1:
            nbrBlk =10000
            self.segmention(depth, height, width, nbrBlk)
        
        for d in range(depth):
            for h in range(height):
                for w in range(width):
                    GridMap[d][h][w].indice = num
                    
                    if GridMap[d][h][w].OBSTACLE == KNOWN_OBSTACLE:
                        GridMap[d][h][w].TAG = FORBIDDEN
                        srcObstacles.append(GridMap[d][h][w].indice)
                    nodeList.append(GridMap[d][h][w])
                    num+=1
        print('               MAP ready', len(nodeList))
        #time.sleep(5)
        return nodeList, srcObstacles,nbrBlk







