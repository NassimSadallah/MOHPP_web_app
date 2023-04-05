'''
Main python script : First the binary map is extracted to a grid, 
second, the velocity and the travel time cost are assigned with the Fast Marching 
(in this version only 8/16 neighbors are considered)
Third, the gradient descent is used to compute the optimal path
Finally, we instantiate a vehicle (either within the SITL or on a real drone and command it

Created on Jul 12, 2021

@author: nassim

'''
from cmath import sqrt

"""
import heapq
import time
H = [(21,1),(45,78),(3,5),(22,1),(14,98),(11,54)]
# Covert to a heap
heapq.heapify(H)
print(H)

# Add element
heapq.heappush(H,(18,21))
heapq.heappush(H,(68,0))
#H.sort()
idx = H.index((18,21), )

print('push & idx',idx, H)
H[idx] = (0,0)
print(H)
#heapq._siftdown(H, 0, idx)
for _ in range(5):
    s=heapq.heappop(H)
    print(H)

time.sleep(10000)
"""

import os,time, eel, sys
from mohpp import MAP, CDMap, utilities, OFPSearch, ONPSearch
from UavAndSensors import VehiclesMethods as VeMeth
from mohpp.utilities import  DetectUnexpectedObs, d_,UavHeading,saveMAP,loadTimeComptation, wavePlot, drawTK, globalPath, MAP3DShow, coordinatesToIndex,indexToCoordinates3D,\
    SavePklMap, loadPklMAP, MAP3Dplot, Drawcurves, view3D, comparativestudy_CT
from UavAndSensors.Sensors import Sensors
from math import floor
import numpy as np
import multiprocessing as mp



sitl_connect ='127.0.0.1:14550'
real_connect ='/dev/ttyAMA0' 
UAV = None
#global variable to store the sensed data and return it for the need
sensedData = {}
                                      #z,y,x(3rd dim)
start_coordinates, goal_coordinates = [2,2,2],[88,59,88]#[z,y,x][60,96],[105,38]#colonne/ligne (East/North)
startd_3, goald_3 = [],[]
nextStep, current = [-1.0, -1.0],[-1.0, -1.0]
plannedPath, Nodes, CDM, extendedObs, start_index, goal_index = [],[],[], [], -1,-1
extendedObs, GhostSet = [],[]

my_options = {
    'mode': "None",
    'host': '0.0.0.0',
    'port': 8080,
}

eel.init('webapp')

sensorType, sensor = '', None#'lidar'#'lidar'hcsr04
#comparativestudy_CT()
def reSet():
        
    return 0

alpha =0.0

@eel.expose
def getTVal():
    
    r= dataLecture()
    return r[6]


@eel.expose
def connect(sens):
    
    global UAV, sensorType, sensor
    sensorType = sens
    
    sensor = Sensors(sensorType)  
    if sensor.health=='BAD':
        return 'Bad Signal'
    UAV = VeMeth.UAV().connect_to_vehicle(sitl_connect, 921600)
    location = [UAV.location.global_frame.lat, UAV.location.global_frame.lon] 
    battery = UAV.battery.level
    alt = UAV.location.global_relative_frame.alt
    spd = UAV.airspeed
    head = UAV.heading
    mode = UAV.mode.name
    armable = UAV.is_armable
    eel.spawn(readData_fork)
    
    return location, armable, battery, mode, head, alt, round(spd,2),sensor.health

@eel.expose
def Launch():
    
    global plannedPath, Nodes, goal_index, CDM, start_index, extendedObs,GhostSet, directBlk     
    alpha = 0.2#0.01* int(str(eel.saturation()()))
    #map1= loadTimeComptation(1)
    #map2 = loadTimeComptation(2)
    #map0 = [(1.82*map2[i]+1.53*map1[i]) for i in range(100)]
    #print map0[0:10]
    #Drawcurves(map0, map1, map2)#np.asarray(map0), np.asarray(map), np.asarray(map1))     
    #reads the binary map from the binarymaps package
    #binMap = os.path.join(os.path.dirname(os.path.abspath(__file__))+"/binarymaps/simulation.png")
    
    #defines the obstacles and return the corresponding indexed node list
    #Nodes, srcObs, block = MAP.MAP_2D.processMap(width, height, binMap, seq =1, nbr_blocks=25)
    GhostSet,d_par, nodeslist, Nodes,srcObs, block,blkxyz,blk_d = MAP.MAP_3D().processMap_3D(seq=1, )
    #MAP3DShow(Nodes, srcObs, [])
    #nodenp = np.array(size=(len(Nodes),5))
    #np.reshape(nodenp, (len(Nodes),5))
    #print(len(srcObs), block)#nodeslist[0])
    #time.sleep(10000)
    
    #MAP3Dplot(node=nodeslist,blk=block,path=[],seq= 1)
    #gets the corresponding indices of the cells start and goal
    start_index = coordinatesToIndex(start_coordinates,d_, 3)
    goal_index = coordinatesToIndex(goal_coordinates, d_,3)
    print start_coordinates, goal_coordinates,start_index, goal_index#, nodeslist[start_index].TAG,nodeslist[goal_index].TAG
    from mohpp import FastMarching
    #MAP3DShow(nodeslist, srcObs,[])
    #computes the velocity and the travel time at each node of the map
    seq, Frstapp = 1,True
    
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    
    x = [i for i in range(100)]
    y = [i for i in range(100)]
    z = [i for i in range(100)]
    cost = [i*1.5 for i in range(100)]
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    scatter = ax.scatter(x, y, z, c=cost)
    fig.colorbar(scatter, label='Cost')
    
    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')
    ax.set_title('3D Environment with Cost')
    
    plt.show()
    
    time.sleep(10000)
    """


    
    t=time.time()
    #VelMap=FastMarching.MSfm3D_SeqPar(srcObs,-1,None,-1,nodeslist,d_par,False, block)#loadPklMAP("VBLHGSmap20.0.pkl", "vblhgs")#
    VelMap=FastMarching.VBLHGS3D(srcObs,-1,nodeslist,d_,False, block)#loadPklMAP("VBLHGSmap20.0.pkl", "vblhgs")#
    print'time ', time.time() - t
    t=time.time()
    ertx1, ertx12 = [],[]
    ertx21, ertx22 = [],[]
    cost, anals = [],[]    
    
    for i in VelMap.nodesList:
        coords = utilities.indexToCoordinates3D(i.indice, 3)
        
        #print i.indice,(i.x,i.y,i.z), coords, sqrt(pow(coords[0]-49,2)+pow(coords[1]-49,2)+pow(coords[2]-49,2)),i.anltic, i.full
        ertx1.append(round(abs(i.cost - i.anltic),7))
        ertx12.append(pow(abs(i.cost - i.anltic),2))  
        cost.append(i.cost)
        anals.append(i.anltic)
        
    sums = sum(ertx1)
    s=0
    for i in ertx1:
        s+=i
    print 'sum vs s:', sums, s    
    print cost[117640:117660]
    print anals[117640:117660]
    print ertx1[117640:117660] 
    print (sum(ertx12)/len(ertx12))
    euc = round((sum(ertx12)/len(ertx12)),5)
    l11,l12,l13 = sums/len(ertx1), sqrt(euc), max(ertx1) 
    

    print 'L11',l11,' L12',l12,' L13',l13
    #MAP3Dplot(node=VelMap.nodesList,blk=block,path=[],seq= 1)
    time.sleep(1000000)
    vmap, vmp=CDMap.get_Vel_Cost(nodes = VelMap.nodesList,node=  Nodes ,alpha = alpha, Velocity=1.0, block=block, seq=2)
    #saveMAP(alpha*100, VelMap.tim, "vblhgs", "vblhgs", 0)
    #SavePklMap("VBLHGSmap"+str(int(alpha*100))+""+".pkl", vmap, "vblhgs")
    del VelMap
    print time.time() - t
    t=time.time()
    VelMaps=FastMarching.VZFm3D([start_index],goal_index,vmp,d_,True, block)
    print time.time() - t
    print VelMaps.tim
    plannedPath, velocity_prof = OFPSearch.Gradient_3D(goal_index, start_index, VelMaps.nodesList,None, 3, 1)
    MAP3Dplot(node=VelMaps.nodesList,blk=block,path=plannedPath,seq= 1)
    time.sleep(1000000)
    
    
    _map = loadPklMAP("Fm"+str(5)+".pkl", "Fm") 
    #print _map[start_index].indice,_map[goal_index].indice,_map[start_index].TAG,_map[goal_index].TAG
    #VelMap=FastMarching.VBLHGS3D([start_index],goal_index,nodeslist,d_,False, block)
    #vmap=CDMap.get_Vel_Cost(nodes = VelMap.nodesList, alpha = alpha, Velocity=1)
    VelMap=FastMarching.Fm3D([start_index],goal_index,_map,d_,True, block)
    
    #VelMap=FastMarching.VBLHGS3D([start_index],goal_index,map,d_,True, True, 1, block)
    #VelMaps=FastMarching.VBLHGS3D([start_index],goal_index,map,d_,True, False, 1, block)
    plannedPath, velocity_prof = OFPSearch.Gradient_3D(goal_index, start_index, VelMap.nodesList,None, 3, 1)
    #Drawcurves(velocity_prof, 0, 0)
    #MAP3DShow(VelMap.nodesList, srcObs, plannedPath)
    MAP3Dplot(node=[],blk=block,path=plannedPath,seq= 1)
    
    print('sleeping', VelMap.tim)
    
    time.sleep(10000000)#os._exit(-1)
    if __name__=='__main__':
        
        for i in range(1,101):
 

            from mohpp import FastMarching

            
            alpha  = i/100
            print i
            
            if seq ==0:
                if Frstapp:
                    VelMap=FastMarching.VZFm3D(srcObs,-1,nodeslist,d_,False, 0)
                    vmap=CDMap.get_Vel_Cost(nodes = VelMap.nodesList, alpha = alpha, Velocity=1)
                    saveMAP(i, VelMap.tim, "Fm", "Fm", 0)
                    SavePklMap("Fm"+str(i)+""+".pkl", vmap, "Fm")
                    del VelMap, vmap 
                
                
                else:
                    _map = loadPklMAP("VBLHGSmap"+str(i)+".pkl", "vblhgs")
                    #print _map[start_index].TAG,_map[goal_index].TAG
                    VelMap=FastMarching.Fm3D([start_index],goal_index,_map,d_,True, 0)
                    #VelMap=FastMarching.MSfm3D(srcObs,None,nodeslist,d_,False, 0)
                    vmap = []
                    #vmap=CDMap.get_Vel_Cost(nodes = VelMap.nodesList, alpha = alpha, Velocity=1)
                    saveMAP(i, VelMap.tim, "MSFMStar", "MSFMStar", 0)
                    #SavePklMap("MSFM"+str(i)+""+".pkl", vmap, "MSFM", 0)
                    del VelMap, _map    
                    #time.sleep(10000)  
                          
            elif seq==1:
                _map = loadPklMAP("VBLHGSmap"+str(i)+".pkl", "vblhgs")
                
                VelMap=FastMarching.VBLHGS3D([start_index],goal_index,_map,d_,True, block)
                #VelMap=FastMarching.VBLHGS3D(srcObs,None,nodeslist,d_,True, False, 1, block)
                vmap = []
                #vmap=CDMap.get_Vel_Cost(nodes = VelMap.nodesList, alpha = alpha, Velocity=1)
                #saveMAP(i, VelMap.tim, 1)
                #SavePklMap("VBLHGSStar"+str(i)+""+".pkl", VelMap, vblhgs,1)
                #view3D(x, y, z, c)
                #del VelMap, vmap                   
                
            elif seq==2:
                
                _map = loadPklMAP("map"+str(i)+".pkl", "bbMSFM")
                   
                goal = Nodes[goal_index]
                sblk = Nodes[start_index].block
                gblk = goal.block
                startObs = []
                startObs.append(Nodes[start_index].idx)
                
                #time.sleep(10)
                
                VelMap=FastMarching.MSfm3D_SeqPar(startObs,sblk, goal,gblk, _map, GhostSet, d_par,True, block,blkxyz,blk_d)                
                #state = [False for _ in range(block)]
                #VelMap=FastMarching.MSfm3D_SeqPar([start_index],goal_index, -1,-1, nodeslist, GhostSet, d_par,True, False, state, block)
                #VelMap=FastMarching.MSfm3D_SeqPar(srcObs,-1, -1,-1, nodeslist, GhostSet, d_par,True, False, state, block)
                vmap = []
                #vmap=CDMap.get_Vel_Cost(nodes = VelMap.nodesList, alpha = alpha, Velocity=1, block = block)
                saveMAP(i, VelMap.tim,"bbMSFMStar", "bbMSFMStar", seq=2)
                #SavePklMap("map"+str(i)+""+".pkl", vmap, "bbMSFM")
                del VelMap, vmap, goal, sblk, gblk, startObs
        #loadPklMAP("map1.pkl")
           
        
        """
        goal = Nodes[goal_index]
        sblk = Nodes[start_index].block
        gblk = goal.block
        startObs = []
        startObs.append(Nodes[start_index].idx)
        print startObs,Nodes[start_index].block,Nodes[start_index].idx
        #time.sleep(10)
        
        CDM=FastMarching.MSfm3D_SeqPar(startObs,sblk, goal,gblk, vmap, GhostSet, d_par,True, True, state, block,blkxyz,blk_d)
        """
        """
        queue = mp.Queue()
        
        def child_init(stPts, g, nlist,d_,):
            print("initializing ...")
            global d_parallel, startObs, nodeslist,Nodes, block, goal_index
            startObs = stPts
            nodeslist = nlist
            d_parallel = d_
            goal_index = g
        
        
        pool = mp.Pool(block, initializer = child_init, initargs = (startObs,-1,nodeslist, d_parallel))# for i in range(4)]
        #print("well finished ",time.time() - tic)
        res = [pool.apply_async(FastMarching.ParallelMSFM(), (state[i],queue, True, False,1,1, i)) for i in range(block)]    
        
        #for i in range(block):
        #    pool.apply_async(FastMarching.ParallelMSFM(state[i],queue,startObs[i], -1, [], nodeslist[i], d_parallel, True, False, 1, 1, i))
         
        pool = mp.Pool(processes=block)   
        res = [pool.apply_async(FastMarching.ParallelMSFM(state[i],queue,startObs[i], -1, [], nodeslist[i], d_parallel, True, False, 1, 1, i)) for i in range(block)]
        #pool = [mp.Process(target =FastMarching.ParallelMSFM(startObs[i], -1, Nodes,nodeslist[i], d_parallel, True, False, 1, 1, i)) for i in range(block)]
    
        del startObs
        #del Nodes

        r =0# [queue.get() for _ in pool]  
        print 'r ',r
        time.sleep(10000)      
        g_index = Nodes[goal_index].idx 
        s_index = Nodes[start_index].idx   
        #CDM = [pool.apply_async(FastMarching.ParallelMSFM(queue,[s_index], g_index, [], nodeslist[i], d_parallel, True, False, 1, 1, i)) for i in range(block)]
        #CDM = CDMap.get_Vel_Cost(3, Nodes, srcObs, start_index, [goal_index], alpha, 1.0, d_, seq=1,block=block)
    #wavePlot(d_[0], d_[1]/d_[0], CDM)
        """
        #plannedPath = OFPSearch.Gradient_3D(start_index, goal_index, CDM.nodesList,Nodes, 1)
    #drawTK(Nodes, [start_coordinates,goal_coordinates])
    MAP3DShow(CDM, srcObs, plannedPath)
    return 'Done !'

@eel.expose
def take_off(h):
    
    global plannedPath, Nodes, goal_index, CDM,start_index, extendedObs 
    
    VeMeth.UAV().takeoff(float(str(h)), UAV)
    heading = UavHeading(UAV.heading)
    default_alt = UAV.location.local_frame.down
    
    nextStep = plannedPath[0]
    plannedPath.remove(plannedPath[0])
    nodeIdx = Nodes[utilities.coordinatesToIndex([int(floor(nextStep[0])),int(floor(nextStep[1]))], d_)].indice

    extendedObs, isDetected, brake = DetectUnexpectedObs(sensorType, sensor.getSensorValues(sensorType), UAV.heading, nodeIdx, Nodes, extendedObs, 1.5, 4, d_)

    #Land(UAV)
    
    #lkj = input("ok next test : ")#extendedObs
    globalPath.append(nextStep)
    isReplanning = False   
    
    while nodeIdx !=goal_index:
            
        #isReplanning = False
        if not isReplanning:
            
            current = nextStep
            
            nextStep = plannedPath[0]
            nodeIdx = Nodes[utilities.coordinatesToIndex([int(floor(nextStep[0])),int(floor(nextStep[1]))], d_)].indice
            nodeVel = Nodes[utilities.coordinatesToIndex([int(floor(nextStep[0])),int(floor(nextStep[1]))], d_)].v
            plannedPath.remove(plannedPath[0])
            n, e, d = utilities.getNorth_East_Down(current, nextStep, UAV.location.local_frame.down, default_alt)
            
            print n, e, d#, nodeVel,utilities.sqrt_dist(n, e, d)
            #set the appropriate speed at which the UAV should travel through the point
            UAV.airspeed = nodeVel
            #send command with the north, east, down( -z) distance to move on 
            VeMeth.UAV().send_NED_velocity(n, e,d, UAV)
            globalPath.append(nextStep)
            #pause the script for the corresponding travel time
            time.sleep(utilities.sqrt_dist(n, e, d))
            #nex = input('next')
            
        else:
            #we recompute the global path once we bypassed the dynamic threats
            current = nextStep
            globalPath.append(nextStep)
            plannedPath = OFPSearch.Gradient(nodeIdx, goal_index, CDM, d_)
            
        #globalPath.append(nextStep)
        
        #sensing the surrounding area with the embedded sensors
        extendedObs, isDetected, brake = DetectUnexpectedObs(sensorType, sensor.getSensorValues(sensorType), UAV.heading, nodeIdx, Nodes, extendedObs, 1.5, 4, d_)
        #print isDetected, brake, len(extendedObs)
        #brake = False
        if brake:#if brake is triggered, we must switch to online process
            print 'Into ONPS ...'
            nextStep = ONPSearch.processONPS(nodeIdx, goal_index, heading,extendedObs, isDetected,brake, Nodes, d_, sensor.getSensorValues(sensorType), sensor, UAV, default_alt)
            isReplanning = True
        
        nodeIdx = Nodes[utilities.coordinatesToIndex([int(floor(nextStep[0])),int(floor(nextStep[1]))], d_)].indice
    time.sleep(1)
    #print globalPath
    drawTK(Nodes, globalPath)
    #wavePlot(width, height, Nodes, globalPath, Nodes[start_index], Nodes[goal_index])
    #Land(UAV)
    
    return 'ok'

def readData_fork():
    global sensedData, sensor, sensorType
    sensedData ={}
    sensedData = sensor.getSensorValues(sensorType)
    #print sensedData
first = True
@eel.expose
def dataLecture():
    global sensedData, first, sensor, sensorType

    sensedData ={}
    sensedData = sensor.getSensorValues(sensorType)  
    location = [UAV.location.global_frame.lat, UAV.location.global_frame.lon] 
    battery = UAV.battery.level
    alt = UAV.location.global_relative_frame.alt
    spd = UAV.airspeed
    head = UAV.heading
    mode = UAV.mode.name
    #theta=sensedData#sensor.getSensorValues(sensorType)#, dist = sens.getLidarValues()

    #print theta#, dist
    return location, battery, head, round(alt,1), round(spd,2), mode, sensedData#theta#, dist

@eel.expose
def Land(uav):
    print('Landing mode activated ...')
    VeMeth.UAV.Land(uav)
    
@eel.expose
def test1(h):
    global plannedPath, Nodes, goal_index, CDM,start_index, extendedObs 
    
    VeMeth.UAV().takeoff(float(str(h)), UAV)
    print('alt at ',UAV.location.local_frame.down )
    time.sleep(5)
    Land(UAV) 
    return 'ok'   
@eel.expose
def test2(h):
    global start_coordinates, goal_coordinates
    start_coordinates, goal_coordinates = [100,75],[100,70]#[60,96],[105,38]#colonne/ligne (East/North)
    Launch()
    take_off(h)
    return 'ok'

@eel.expose
def TESTMODE():
    global UAV
    MODE = ['LAND','STABILIZE','GUIDED']
    for mode in MODE:
        VeMeth.UAV()._vehicle_mode(mode, UAV)

@eel.expose
def stopMOHPP():
    return sys.exit(0)

#if __name__=='__main__':
    
eel.start('index.html', my_options, block = True)



"""


'''
connection to vehicle and controlling SITL:('127.0.0.1:14550', 921600), Real: ('/dev/ttyAMA0', baud=921600) 
'''

time.sleep(3)

VeMeth.UAV().takeoff(5.0, UAV)
heading = UavHeading(UAV.heading)
default_alt = UAV.location.local_frame.down
nextStep = plannedPath[0]
plannedPath.remove(plannedPath[0])
nodeIdx = Nodes[utilities.coordinatesToIndex([int(floor(nextStep[0])),int(floor(nextStep[1]))], d_)].indice

isReplanning = False

while nodeIdx !=goal_index:
    
    if not isReplanning:
        
        current = nextStep
        nextStep = plannedPath[0]
        nodeIdx = Nodes[utilities.coordinatesToIndex([int(floor(nextStep[0])),int(floor(nextStep[1]))], d_)].indice
        nodeVel = Nodes[utilities.coordinatesToIndex([int(floor(nextStep[0])),int(floor(nextStep[1]))], d_)].v
        plannedPath.remove(plannedPath[0])
        n, e, d = utilities.getNorth_East_Down(current, nextStep, UAV.location.local_frame.down, default_alt)
        print n, e, d, nodeVel,utilities.sqrt_dist(n, e, d)
        #set the appropriate speed at which the UAV should travel through the point
        UAV.airspeed = nodeVel
        #send command with the north, east, down( -z) distance to move on 
        VeMeth.UAV().send_NED_velocity(n, e,d, UAV)
        #pause the script for the corresponding travel time
        time.sleep(utilities.sqrt_dist(n, e, d))
        
    else:
        #we recompute the global path once we bypassed the dynamic threats
        current = nextStep
        plannedPath = OFPSearch.Gradient(nodeIdx, goal_index, CDM, d_)
        
    globalPath.append(nextStep)
    
    #sensing the surrounding area with the embedded sensors
    extendedObs, isDetected, brake = DetectUnexpectedObs(Server.getTVal(), UAV.heading, nodeIdx, Nodes, extendedObs, 1.5, 4, d_)
    print isDetected, brake, len(extendedObs)
    #brake = False
    if brake:#if brake is triggered, we must switch to online process
        print 'replanning ...'
        nextStep = ONPSearch.processONPS(nodeIdx, goal_index, heading,extendedObs, isDetected, Nodes, d_, Server.getTVal())
        isReplanning = True

    nodeIdx = Nodes[utilities.coordinatesToIndex([int(floor(nextStep[0])),int(floor(nextStep[1]))], d_)].indice
"""
