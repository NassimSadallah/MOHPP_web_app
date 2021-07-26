'''
Created on Jul 8, 2021

@author: nassim
'''

from PIL import Image
from Node import node
from utilities import FORBIDDEN, NO_OBSTACLE, KNOWN_OBSTACLE, UNDETECTED_OBSTACLE

tabGrid, GridMap = [],[]

'''
Read a binary map with static obstacles
The image should be at the same package as the current script
'''

def readMap(r, r2,binaryMap):

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
define which pixels are refering to the obstacles
''' 
def defineObstacles( r, r2): #on definit les noeuds qui sont des obstacle
    
    for l in range(r2):
        for c in range(r):
            if tabGrid[l][c]==-1:
                GridMap[l][c].OBSTACLE = KNOWN_OBSTACLE
                
            elif tabGrid[l][c]==-2:
                GridMap[l][c].OBSTACLE = UNDETECTED_OBSTACLE
                
            else:
                GridMap[l][c].OBSTACLE = NO_OBSTACLE       
                
'''
    SEQ is a variable which defines if the algorithm is parallelized or not:
        seq=0: one block execution (default value)
        seq=1: use sequenced list of nodes (shorter lenghth to reduce the computation time)
        seq=2: use multi threading (not tested yet)
        seq=3: use multiprocesses on a CPU (not finished)
        seq=4: use multiprocesses on a GPU (not finished)
                                                                
'''
def processMap( r, r2, binaryMap, seq = 0):
    global tabGrid, GridMap
    tabGrid = [[0 for i in range(r)] for j in range(r2)]
    GridMap = [[0 for i in range(r)] for j in range(r2)]
    num, nodeList, srcObstacles = 0, [], []
    print('--> MAP processing...')

    '''
    assign node for each cell 
    '''
    for i in range(r2):
        for j in range(r):
            GridMap[i][j] = node(i, j)
    
    readMap(r,r2,binaryMap)
    
    defineObstacles(r,r2)
    
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
    return nodeList, srcObstacles
        
        
                   
         











