'''
Created on Jul 8, 2021

@author: nassim
'''
from math import floor, cos, sin, pi, sqrt, log10
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from heapq import heappop
from Tkinter import *
import time, os, sys
from multiprocessing import Queue
from pylab import *
from pygame import color


global que
que = Queue()

risks = 0
UNDETECTED_OBSTACLE =-1 
NO_OBSTACLE = 0
KNOWN_OBSTACLE = 1
CONTOUR = -5
FORBIDDEN = -1 
NEW_FORBIDDEN = -2 
INFINI=9999.0
width, height, depth = 100, 100,100# named in nodes as : z,y,x
echelle = 2
FAR, KNOWN, FROZEN, FROZEN_FIX = -1,0,1,2
globalPath = []
#d_ = [width, width*height]
d_ =[width, width*height, width*height*depth] 
d_parallels = [12,12*12,12*12*12]
d_parallel = [10,100,2000]

def coordinatesToIndex(current_coordinates = [],d_ = [], dim = 3):
    
    idx = current_coordinates[0]
    for i in range(1,dim):
        #print current_coordinates[i], d_[i-1]
        idx += current_coordinates[i] * d_[i-1]

    return idx

def indexToCoordinates3D(inde, dim):
    coord = [0 for _ in range(dim)]
    coord[2] = int(floor(inde/d_[1]))
    aux = inde-coord[2]*d_[1]
    coord[1] = int(floor(aux/d_[0]))
    aux -= coord[1]*d_[0]
    coord[0] = aux 
    
    return coord

def indexToCoordinates(idx, d_ = [-1, -1]):
    current_coordinates = [-1, -1]
    current_coordinates[1] = int(floor(idx/d_[0]))
    current_coordinates[0] = (idx-current_coordinates[1]*d_[0])      
    return current_coordinates

def getmaxcost(l):

    m = 0
    for i in l:
        if m < i.cost and i.cost < INFINI:
            m = i.cost  
    return m

def getRisk(l):
    global risks    
    for i in l:
        if i.cost != INFINI and risks< i.risk:      
            risks = i.risk
    return risks 

def getmaxVel(l):
    m = 0
    for i in l:
        if m < i.v:# and not i.v == INFINI:
            m = i.v
        
    return m

def getNorth_East_Down(current, nextStep, down, default_altitude):
    
    dNorth = round((-(nextStep[1])+(current[1])),5)
    dEast =round((nextStep[0])-(current[0]),5)
    dDown = round((default_altitude-down) ,5)   
    
    return dNorth, dEast, dDown

def sqrt_dist(a,b,c = 0):
    
    return round(sqrt(a*a+b*b+c*c),7)

def wavePlot(r= -1, r2= -1,nodes = [],path = [], start = [], goal = []):
    zVelocity = [[0 for _ in range(r)] for _ in range(r2)]
    xvec = np.linspace(0.,r,r)
    yvec = np.linspace(0.,r2,r2)                               
    x,y = np.meshgrid(xvec, yvec)
    
    ind = 0
    cou2 = r2-1
    
    o = getmaxcost(nodes)    
    
    while not cou2==0:
        for c in range(0,r):
            if nodes[ind].cost >o:
                nodes[ind].cost = o+10
            zVelocity[cou2][c] = nodes[ind].cost
            ind+=1 
        cou2-=1                
    
    plt.imshow(zVelocity, cmap='Greys',  interpolation='nearest')
    
    for i in path:
       
        i[1] = r2-i[1]    
        plt.plot(i[0],i[1],'m+')
    
    plt.plot(start.abscice,r2-start.colonne,'^k:')
    plt.plot(goal.abscice,r2-goal.colonne,'go')
    plt.contourf(x, y, zVelocity ,200)                             
    plt.colorbar() 
    
    plt.show()    

def isEmptyList(l):
    
    for i in l:
        if i !=[]:
            #print i
            return False    
    return True

def locateBestIdx(l):
    
    b =(INFINI,INFINI)
    ids = -1
    all = []
    for co,i in enumerate(l):
        if i !=[] and i[0]<b:
            b = i[0]
            ids = co
        if i !=[]:    
            all.append(i[0])
    #return the popped best element in the hosting list

    if ids!=-1:
        best = heappop(l[ids])
        #print 'best and all ', best, l
        return best
        
        
    else:
        print 'NO best element '
        return -1

def UavHeading(heading):
    return heading*pi/180

def degToCartesian(deg):
    return deg*pi/180

def DetectUnexpectedObs(sensType, sensors, UavOrient, curIdx, nodes, extendedObs, safety_margin, sensRange, d_):
    
    isDetected, brake = False, False
    
    cur = [nodes[curIdx].abscice,nodes[curIdx].colonne] # save the x, y coordinates of the current position of the UAV    
    
    if sensType == 'hcsr04':
        sensorsValues=sensors#call the method which read the sensors' values sent from the arduino nano
        sN, sE, sS, sW, sNE, sSE, sSW, sNW = [0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]

        if sensorsValues !=[]: 
            
            if sensorsValues[0] >=.8 and sensorsValues[0] < sensRange:
    
                #transforms the radian values to cartesian one
                sN= [int(round(int(sensorsValues[0])*cos(2*pi-UavOrient+pi/2),0)),-int(round(int(sensorsValues[0])*sin(2*pi-UavOrient+pi/2),0))]   
                obs = [cur[0]+sN[0],cur[1]+sN[1]]
                    
                if ((obs>=[0,0]) and(obs<=[width, height])):
                      
                    curobs = nodes[coordinatesToIndex(obs, d_)]
                
                    if curobs.OBSTACLE != KNOWN_OBSTACLE and curobs.TAG != FORBIDDEN and sN !=[0,0]:            
                        
                        isDetected = True
                        curobs.OBSTACLE = KNOWN_OBSTACLE
                        curobs.TAG = NEW_FORBIDDEN
                        curobs.cost = INFINI
                        extendedObs.append(curobs)
                                    
                        if sensorsValues[0]<=safety_margin:
                            
                            brake = True
            
            if sensorsValues[1] >=.8 and sensorsValues[1] < safety_margin:
                
                sE= [int(round(int(sensorsValues[1])*cos(2*pi-UavOrient),0)),int(round(int(sensorsValues[1])*sin(2*pi-UavOrient),0))]
                
                obs = [cur[0]+sE[0],cur[1]+sE[1]]
                    
                if ((obs>=[0,0]) and(obs<=[width, height])):
                    curobs = nodes[coordinatesToIndex(obs, d_)]
                    
                    if curobs.OBSTACLE != KNOWN_OBSTACLE and curobs.TAG != FORBIDDEN and sE !=[0,0]:            
                        
                        isDetected = True
                        curobs.OBSTACLE = KNOWN_OBSTACLE
                        curobs.TAG = NEW_FORBIDDEN
                        curobs.cost = INFINI
                        extendedObs.append(curobs)
                                    
                        if sensorsValues[1]<=safety_margin:
                            
                            brake = True
            
            if sensorsValues[2] >=.8 and sensorsValues[2] < sensRange:
                
                sS= [int(round(int(sensorsValues[2])*cos(2*pi-UavOrient+pi*3/2),0)),-int(round(int(sensorsValues[2])*sin(2*pi-UavOrient+pi*3/2),0))]          
                obs = [cur[0]+sS[0],cur[1]+sS[1]]
                    
                if ((obs>=[0,0]) and(obs<=[width, height])):
                
                    curobs = nodes[coordinatesToIndex(obs, d_)]
                    
                    if curobs.OBSTACLE != KNOWN_OBSTACLE and curobs.TAG != FORBIDDEN and sS !=[0,0]:            
                        
                        isDetected = True
                        curobs.OBSTACLE = KNOWN_OBSTACLE
                        curobs.TAG = NEW_FORBIDDEN
                        curobs.cost = INFINI
                        extendedObs.append(curobs)
                                    
                        if sensorsValues[2]<=safety_margin:
                      
                            brake = True
            
            if sensorsValues[3] >=.8 and sensorsValues[3] < sensRange:
                
                sW= [int(round(int(sensorsValues[3])*cos(2*pi-UavOrient+pi),0)),int(round(int(sensorsValues[3])*sin(2*pi-UavOrient+pi),0))]            
                obs = [cur[0]+sW[0],cur[1]+sW[1]]
                    
                if ((obs>=[0,0]) and(obs<=[width, height])):
                    
                    curobs = nodes[coordinatesToIndex(obs, d_)]
                    
                    if curobs.OBSTACLE != KNOWN_OBSTACLE and curobs.TAG != FORBIDDEN and sW !=[0,0]:            
                        
                        isDetected = True
                        curobs.OBSTACLE = KNOWN_OBSTACLE
                        curobs.TAG = NEW_FORBIDDEN
                        curobs.cost = INFINI
                        extendedObs.append(curobs)
                                    
                        if sensorsValues[3]<=safety_margin:
                          
                            brake = True
                            
            if sensorsValues[4] >=.8 and sensorsValues[4] < sensRange:
                
                sNE=[int(round(int(sensorsValues[4])*cos(2*pi-UavOrient+pi/4),0)),-int(round(int(sensorsValues[4])*sin(2*pi-UavOrient+pi/4),0))] 
                
                obs = [cur[0]+sNE[0],cur[1]+sNE[1]]
                    
                if ((obs>=[0,0]) and(obs<=[width, height])):            
                    
                    curobs = nodes[coordinatesToIndex(obs, d_)]
                    
                    if curobs.OBSTACLE != KNOWN_OBSTACLE and curobs.TAG != FORBIDDEN and sNE !=[0,0]:            
                        
                        isDetected = True
                        curobs.OBSTACLE = KNOWN_OBSTACLE
                        curobs.TAG = NEW_FORBIDDEN
                        curobs.cost = INFINI
                        extendedObs.append(curobs)
                                    
                        if sensorsValues[4]<=safety_margin:
                            
                            brake = True
                    
            if sensorsValues[5] >=.8 and sensorsValues[5] < sensRange:
                
                sSE=[int(round(int(sensorsValues[5])*cos(2*pi-UavOrient+7*pi/4),0)),-int(round(int(sensorsValues[5])*sin(2*pi-UavOrient+7*pi/4),0))]
                obs = [cur[0]+sSE[0],cur[1]+sSE[1]]
                    
                if ((obs>=[0,0]) and(obs<=[width, height])):            
                
                    curobs = nodes[coordinatesToIndex(obs, d_)]
                    
                    if curobs.OBSTACLE != KNOWN_OBSTACLE and curobs.TAG != FORBIDDEN and sSE !=[0,0]:            
                        
                        isDetected = True
                        curobs.OBSTACLE = KNOWN_OBSTACLE
                        curobs.TAG = NEW_FORBIDDEN
                        curobs.cost = INFINI
                        extendedObs.append(curobs)
                                    
                        if sensorsValues[5]<=safety_margin:
                          
                            brake = True
                
            if sensorsValues[6] >=.8 and sensorsValues[6] < sensRange:
                
                sSW=[int(round(int(sensorsValues[6])*cos(2*pi-UavOrient+5*pi/4),0)),-int(round(int(sensorsValues[6])*sin(2*pi-UavOrient+5*pi/4),0))] 
                obs = [cur[0]+sSW[0],cur[1]+sSW[1]]
                    
                if ((obs>=[0,0]) and(obs<=[width, height])):            
                                
                    curobs = nodes[coordinatesToIndex(obs, d_)]
                    
                    if curobs.OBSTACLE != KNOWN_OBSTACLE and curobs.TAG != FORBIDDEN and sSW !=[0,0]:            
                        
                        isDetected = True
                        curobs.OBSTACLE = KNOWN_OBSTACLE
                        curobs.TAG = NEW_FORBIDDEN
                        curobs.cost = INFINI
                        extendedObs.append(curobs)
                                    
                        if sensorsValues[6]<=safety_margin:
                          
                            brake = True
                        
            if sensorsValues[7] >=.8 and sensorsValues[7] < sensRange:
    
                sNW=[int(round(int(sensorsValues[7])*cos(2*pi-UavOrient+3*pi/4),0)),-int(round(int(sensorsValues[7])*sin(2*pi-UavOrient+3*pi/4),0))]      
                obs = [cur[0]+sNW[0],cur[1]+sNW[1]]
                    
                if ((obs>=[0,0]) and(obs<=[width, height])):             
                    
                    curobs = nodes[coordinatesToIndex(obs, d_)]
                    
                    if curobs.OBSTACLE != KNOWN_OBSTACLE and curobs.TAG != FORBIDDEN and sNW !=[0,0]:            
                        
                        isDetected = True
                        curobs.OBSTACLE = KNOWN_OBSTACLE
                        curobs.TAG = NEW_FORBIDDEN
                        curobs.cost = INFINI
                        extendedObs.append(curobs)
                                    
                        if sensorsValues[7]<=safety_margin:
                           
                            brake = True
                        
            #print  sensorsValues, sN,sE,sS,sW,sNE, sSE, sSW, sNW   
            sensorsValues = []  
        
    elif sensType=='lidar':#if lidar is the integrated sensor

        risks = getRisk(nodes)
        sensorsValues = sensors.items()        
        cartesval = set()
            
        for s in sensorsValues:
            
            if s[1]<4.0 and s[1] >=0.50:
                
                cartesianS = [int(round(int(s[1])*cos(2*pi-UavOrient+degToCartesian(s[0])),0)),-int(round(int(s[1])*sin(2*pi-UavOrient+degToCartesian(s[0])),0))]
                obs = (cur[0]+cartesianS[0],cur[1]+cartesianS[1])
                cartesval.add(obs)
                
        for obs in cartesval:        
                
                curobs = nodes[coordinatesToIndex(obs, d_)]

                       
                
                if (curobs.OBSTACLE != KNOWN_OBSTACLE and curobs.TAG != FORBIDDEN and cartesianS !=[0,0]):            
                    
                    isDetected = True
                    curobs.OBSTACLE = KNOWN_OBSTACLE
                    curobs.TAG = NEW_FORBIDDEN
                    curobs.cost = INFINI
                    extendedObs.append(curobs)
                                
                    if s[1]<=safety_margin:
                       
                        brake = True
                
                    #print obs, curobs.indice

    #extendedObs = list(set(extendedObs))        
    return extendedObs, isDetected, brake#cartesval#
          
def dessiner(noeud,col, grille):
           
            grille.create_rectangle(noeud.abscice*echelle,noeud.colonne*echelle,noeud.abscice*echelle+echelle,
                                        noeud.colonne*echelle+echelle,outline = col,fill=col)

def dessinerPath(x,y,co, grille, fenetre):
 
            
            grille.create_line(x[0]*echelle,x[1]*echelle,y[0]*echelle,y[1]*echelle,width= 2,fill=co)
            fenetre.update()    
   
def drawTK(Nodes, path = globalPath):
    fenetre =Tk()
    fenetre.title("MAP MOPP")
    grille = Canvas(fenetre,width=width*2,height=height*2,bg='grey')
    grille.pack()
    grille.delete(ALL)
    for i in Nodes:
        #ffcccb ffd5d3 ffdddc ffe6e5 ffeeed fff6f6 ffffff
        pourcent = round(i.v,2)*100
        if pourcent ==0:
            dessiner(i, 'black', grille)
        elif pourcent ==1:
            dessiner(i, 'white', grille)
        else: 
            coleur = "grey"+str(int(pourcent))
            dessiner(i, coleur, grille)
        if i.OBSTACLE ==UNDETECTED_OBSTACLE:
            dessiner(i, 'grey', grille)
    if path !=[]:
        x = path[0]
        path.remove(path[0])
        for i in path:
            dessinerPath(x, i, 'magenta', grille, fenetre)
            dessiner(Nodes[coordinatesToIndex([100,74], d_)], 'green', grille)
            dessiner(Nodes[coordinatesToIndex([20,50], d_)], 'red', grille)
            fenetre.update()
            x= i  
    fenetre.mainloop()     

def MAP3DShow(node, sobs,path):
    x,y,z = [],[],[]
    x1,y1,z1 = [],[],[]
    coord = [-1,-1,-1]
    maxcost = getmaxcost(node)
    
    for i in node:
       
        if i.OBSTACLE ==KNOWN_OBSTACLE:
            coord=indexToCoordinates3D(i.indice,3)
    
            x.append(coord[0])
            y.append(coord[1])
            z.append(coord[2])
    
    if len(path) != 0:
        
        for i in path:
            x1.append(i[0])
            y1.append(i[1])
            z1.append(i[2])
    
    
    ax = plt.axes(projection ='3d')
    
    ax.set_xlabel('width')
    ax.set_ylabel('height')
    ax.set_zlabel('deep (z+)');
    
    plt.xlim(0,201)
    plt.ylim(0,201)
    
    s = ax.contour(x,y,z,c='black')
    
    s.set_edgecolors = s.setfacecolors = lambda *args:None
    
    if len(path) !=0:
        
        ax.plot3D(x1, y1, z1,'orange',linewidth = 2)
    
    plt.show()

def Drawcurves(s0,s1,s2):
    #print len(s1)   
    #plt.ylim(0,2)
    from scipy.interpolate import spline
    
    x = np.linspace(0,len(s0),len(s0))
    #ps = spline(x,s0,x)
    #plt.plot(x, ps)
                
    plt.plot(s0,'b',label = 'velocity')
    #plt.plot(s1,'r', label = 'VZLHGS')
    #plt.plot(s2,'g', label = 'BBMSFM')
    
    plt.ylabel('velocity (m/s)')
    plt.xlabel('nodes')#(u'\u03B1 coefficient')
    plt.title("Velocity profil")
    plt.legend(loc='center left', shadow=False, fontsize='large')
    plt.show()
    
    
    return 0

def view3D(x,y,z,c):
    from matplotlib import cbook
    from matplotlib import cm
    from matplotlib.colors import LightSource
    import matplotlib.pyplot as plt
    import numpy as np
    
    
    
    """
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri   
    
    name_color_map = "blues" 
    name_color_map_surface = 'Greens';  # Colormap for the 3D surface only.
    index_x = 0; index_y = 1; index_z = 2; index_c = 3;
    list_name_variables = ['x', 'y', 'z', 'c'];
    fig = plt.figure(); 
    ax = fig.add_subplot(111, projection='3d');
    ax.set_xlabel(list_name_variables[index_x]); ax.set_ylabel(list_name_variables[index_y]);
    ax.set_zlabel(list_name_variables[index_z]);
    plt.title('%s in fcn of %s, %s and %s' % (list_name_variables[index_c], list_name_variables[index_x], list_name_variables[index_y], list_name_variables[index_z]) );
    
    # In this case, we will have 2 color bars: one for the surface and another for 
    # the "scatter plot".
    # For example, we can place the second color bar under or to the left of the figure.
    choice_pos_colorbar = 2;
    
    #The scatter plot.
    img = ax.scatter(x, y, z, c = c, cmap = name_color_map);
    cbar = fig.colorbar(img, shrink=0.5, aspect=5); # Default location is at the 'right' of the figure.
    cbar.ax.get_yaxis().labelpad = 15; cbar.ax.set_ylabel(list_name_variables[index_c], rotation = 270);
    
    # The 3D surface that serves only to connect the points to help visualize 
    # the distances that separates them.
    # The "alpha" is used to have some transparency in the surface.
    surf = ax.plot_trisurf(x, y, z, cmap = name_color_map_surface, linewidth = 0.2, alpha = 0.25);
    
    # The second color bar will be placed at the left of the figure.
    if choice_pos_colorbar == 1: 
        #I am trying here to have the two color bars with the same size even if it 
        #is currently set manually.
        cbaxes = fig.add_axes([1-0.78375-0.1, 0.3025, 0.0393823, 0.385]);  # Case without tigh layout.
        #cbaxes = fig.add_axes([1-0.844805-0.1, 0.25942, 0.0492187, 0.481161]); # Case with tigh layout.
    
        cbar = plt.colorbar(surf, cax = cbaxes, shrink=0.5, aspect=5);
        cbar.ax.get_yaxis().labelpad = 15; cbar.ax.set_ylabel(list_name_variables[index_z], rotation = 90);
    
    # The second color bar will be placed under the figure.
    elif choice_pos_colorbar == 2: 
        cbar = fig.colorbar(surf, shrink=0.75, aspect=20,pad = 0.05, orientation = 'horizontal');
        cbar.ax.get_yaxis().labelpad = 15; cbar.ax.set_xlabel(list_name_variables[index_z], rotation = 0);
    #end
    plt.show();    
    """

def MAP3Dplot(node, blk = None, path = None, seq = None):
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    import numpy as np
    
    
    maxvalue = 0
    x, y, z, co = [],[],[],[]
    if seq==0:
        return 0
    
    elif seq==1:
        anltic = []
        maxvalue=getmaxcost(node)    
        
        print maxvalue
           
        for i in node:
            #if i.OBSTACLE ==KNOWN_OBSTACLE:
                coord=indexToCoordinates3D(i.indice,3)
        
                x.append(coord[0])
                y.append(coord[1])
                z.append(coord[2])
                if i.cost>maxvalue:
                    i.cost = maxvalue+15
                co.append(i.cost)
                anltic.append(i.anltic)
        print len(x), len(co)
        
    elif seq==2:
        for i in range(blk):
            if maxvalue<getmaxcost(node[i]):    
                maxvalue = getmaxcost(node[i])
      
        for j in range(blk):   
            for i in node[j]:
               
                coord=indexToCoordinates3D(i.indice,3)
        
                x.append(coord[0])
                y.append(coord[1])
                z.append(coord[2])
                co.append(i.v)
        
    x1,y1,z1 = [],[],[]
    if len(path) != 0:
        
        for i in path:
            x1.append(i[0])
            y1.append(i[1])
            z1.append(i[2])
    
    print('plotting')
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    color_map = cm.ScalarMappable(cmap=cm.Blues_r)
    clr = np.asarray(z)
    color_map.set_array(clr)
    plt.xlabel('x')
    plt.ylabel('y')

    ax.contourf(x,y,25, c= anltic)#, cmap=plt.hot())#,s = 100)#, color ='blue')  
    #s = ax.scatter(x,y,c=anltic, cmap=plt.hot())
    #s.set_array(co)
    #for p in path:
        
    #    ax.plot3D(X[0],p[1],p[2],'orange',linewidth = 2)
    #ax.plot3D(x1,y1,z1,'orange',linewidth= 2)
    #plt.plot(x1[-2:-1],y1[-2:-1],z1[-2:-1],'go')
    #plt.colorbar(color_map)
    plt.show() 
    print 'end'  
    
    """
    co = []
    for i in node:
        if i.cost<INFINI:
            co.append(i.cost)
  
  
    
    # create grid of x, y, z coordinates
    x, y, z = np.meshgrid(np.arange(100), np.arange(100), np.arange(100))
    
    # calculate cost values
    cost = co
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # plot contour
    contour = ax.contour(x, y, z, cost, cmap='coolwarm')
    fig.colorbar(contour, label='Cost')
    
    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')
    ax.set_title('3D Environment with Cost')
    
    plt.show()
   
    
    
     
    """
    """        
    X, Y = np.meshgrid(np.asarray(x),np.asarray(y))
    Z = np.asarray(z)#np.sinc(np.sqrt(X*X+Y*Y))
    # this is the value to use for the color
    V = np.asarray(c)
    
    # create the figure, add a 3d axis, set the viewing angle
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    
    # here we create the surface plot, but pass V through a colormap
    # to create a different color for each patch
    ax.plot_surface(X, Y, Z,c=V, cmap = plt.hot())        
    plt.show()  
    """

def saveMAP(a, ti, fold, file = "", seq=-1):

    with open(os.path.join(os.path.dirname(os.path.abspath(__file__))+"/../binarymaps/"+fold+"/"+file+".txt"),'a') as mapfile:
        
        mapfile.write("%s-%s\n" %(str(a).strip(),str(ti).strip()))
        mapfile.close()  
        print 'map saved'
    
def loadTimeComptation(fold, fi, seq):
    computation = []
    if seq==0:
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__))+"/../binarymaps/"+fold+"/"+fi+".txt"),'r') as mapfile:
            
            content = mapfile.readlines()
            
            for l in content:
            
                _time =l.split('-')
                computation.append(float(_time[1])) 
 
            mapfile.close()  
            print 'SEQ0: Time charged '
            
    else:
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__))+"/../binarymaps/seqBlk/"+"Seq_MAP.txt"),'r') as mapfile:
            content = mapfile.readlines()
            
            for l in content:
            
                _time =l.split('-')
                computation.append(float(_time[1]))
 
            mapfile.close()  
            print 'SEQ4: Time charged ',computation[0:10]
    return computation        
import pickle  
import bz2   
       
def SavePklMap(PIK, nodes, fold = ""):#, seq = -1):
    
    #if seq==0:#one block exexution of MSFM
    print 'dumping map...'
    f= bz2.BZ2File(os.path.join(os.path.dirname(os.path.abspath(__file__))+"/../binarymaps/"+fold+"/"+PIK),"wb")
    pickle.dump(nodes, f)
    f.close()
    print 'saved!'
    """        
    elif seq==1:
        print 'dumping map...'
        f= bz2.BZ2File(os.path.join(os.path.dirname(os.path.abspath(__file__))+"/../binarymaps/vblhgs/"+PIK),"wb")
        pickle.dump(nodes, f)
        f.close()
        print 'saved!' 
        
    else:
        print 'dumping map...'
        f= bz2.BZ2File(os.path.join(os.path.dirname(os.path.abspath(__file__))+"/../binarymaps/seqBlk/"+PIK),"wb")
        pickle.dump(nodes, f)
        f.close()
        print 'saved!'    
    """
     
def loadPklMAP(PIK,fold):
    print 'loading ', fold, PIK,'...'
    t = time.time()
    f= bz2.BZ2File(os.path.join(os.path.dirname(os.path.abspath(__file__))+"/../binarymaps/"+fold+"/"+PIK),"rb")
    map = pickle.load(f)
    f.close()
    print 'loaded !', time.time()-t
    #print map[0][0],map[20][45].cost 
    return map
    #time.sleep(1000000)

def comparativestudy_CT():
    
    fm = loadTimeComptation("Fm", "Fm", 0)
    #fms= loadTimeComptation("FmStar", "FmStar", 0)
    #msfm= loadTimeComptation("MSFM", "MSFM", 0)
    msfms= loadTimeComptation("MSFMStar", "MSFMStar", 0)
    vzlhgs = loadTimeComptation("vblhgs", "vblhgs", 0) 
    vzlhgss= loadTimeComptation("vblhgsStar", "vblhgsStar", 0)
    bbmsfm= loadTimeComptation("bbMSFM", "bbMSFM", 0)
    bbmsfms = loadTimeComptation("bbMSFMStar", "bbMSFMStar", 0)
    plotResults(fm,msfms,vzlhgs,vzlhgss,bbmsfm,bbmsfms)

def plotResults(*args):

    #plt.ylim(0,160)

    plt.ylabel('time (sec)')
    plt.xlabel(u"\u03B1 (%)")
    linestyleset = ['-','-','-','-',':','-.']
    colorset = ['b','r','g','b','r','g']
    labelset = ['fm','msfms','vblhgs','vblhgss','bbmsfm','bbmsfms']
    for i,arg in enumerate(args):
        plt.plot(arg, colorset[i], linestyle =linestyleset[i], label=labelset[i])
    
    """    
    plt.plot(Rpiapp2[3], 'b', label='Rpi_FM* ')
    plt.plot(Rpiapp2[1], 'r', label='Rpi_MSFM*')
    plt.plot(Rpiapp2[5], 'g', label='Rpi_VBLHGS*')
    
    plt.plot(Nvidiaapp2[5], 'b', linestyle = '-.', label='Nvi_FM*') 
    plt.plot(Nvidiaapp2[3], 'r', linestyle = '-.', label='Nvi_MSFM*')
    plt.plot(Nvidiaapp2[0], 'g', linestyle = '-.', label='Nvi_VBLHGS*')
    """
       
    plt.legend(loc='upper', shadow=False)
    plt.show()   
    
def contourplot(*args):
    return 0   
    
    
    