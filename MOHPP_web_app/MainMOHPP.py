'''
Main python script : First the binary map is extracted to a grid, 
second, the velocity and the travel time cost are assigned with the Fast Marching 
(in this version only 8/16 neighbors are considered)
Third, the gradient descent is used to compute the optimal path
Finally, we instantiate a vehicle (either within the SITL or on a real drone and command it

Created on Jul 12, 2021

@author: nassim

'''
import os,time, eel, sys
from mohpp import MAP, CDMap, utilities, OFPSearch, ONPSearch
from UavAndSensors import VehiclesMethods as VeMeth
from mohpp.utilities import  DetectUnexpectedObs, height, width,depth,d_parallel, d_,UavHeading, wavePlot, drawTK, globalPath, MAP3DShow, coordinatesToIndex,indexToCoordinates3D
from UavAndSensors.Sensors import Sensors
from math import floor
import numpy as np
import multiprocessing as mp
from __builtin__ import True


sitl_connect ='127.0.0.1:14550'
real_connect ='/dev/ttyAMA0' 
UAV = None
#global variable to store the sensed data and return it for the need
sensedData = {}
                                      #z,y,x(3rd dim)
start_coordinates, goal_coordinates = [12,17,17],[87,95,39]#[z,y,x][60,96],[105,38]#colonne/ligne (East/North)
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

    #reads the binary map from the binarymaps package
    binMap = os.path.join(os.path.dirname(os.path.abspath(__file__))+"/binarymaps/simulation.png")
    
    #defines the obstacles and return the corresponding indexed node list
    #Nodes, srcObs, block = MAP.MAP_2D.processMap(width, height, binMap, seq =1, nbr_blocks=25)
    GhostSet,directBlk, nodeslist, Nodes,srcObs, block,dp, hei, wid = MAP.MAP_3D().processMap_3D(seq=1, )
    #MAP3DShow(Nodes, srcObs, [])
    #nodenp = np.array(size=(len(Nodes),5))
    #np.reshape(nodenp, (len(Nodes),5))
    #print(len(srcObs), block)#nodeslist[0])
    #time.sleep(10000)
    
    
    #gets the corresponding indices of the cells start and goal
    start_index = coordinatesToIndex(start_coordinates, 3)
    goal_index = coordinatesToIndex(goal_coordinates, 3)
   
    #MAP3DShow(Nodes, srcObs,[])
    #computes the velocity and the travel time at each node of the map
    if __name__=='__main__':
        
        startObs = [[] for _ in range(block)]

        state = [False for _ in range(block)]
        for i in srcObs:
            #print i, Nodes[i].idx, Nodes[i].block
            startObs[Nodes[i].block].append(Nodes[i].idx)
            #state[Nodes[i].block] = True
            
        #print startObs[0],startObs[50], startObs[450]
        from mohpp import FastMarching
        initime = time.time()
        ite = 0
        
        print d_parallel 
        d_par = [[dp+directBlk[i][0],(dp+directBlk[i][0])*(hei+directBlk[i][1]),\
            (dp+directBlk[i][0])*(hei+directBlk[i][1])*(wid+directBlk[i][2])] for i in range(block)] 
        
        print len(nodeslist[0]),directBlk[0][2],d_par[0], len(nodeslist[1]),d_par[1],len(nodeslist[224]), d_par[224] 
        time.sleep(1000)    
        
        VelMap=FastMarching.MSfm3D_SeqPar(startObs,-1, -1,-1, nodeslist, GhostSet, d_par,True, False, state, block)
        
        vmap = []
        
        vmap=CDMap.get_Vel_Cost(nodes = VelMap.nodesList, alpha = alpha, Velocity=1, block = block)
        
        goal = Nodes[goal_index]
        sblk = Nodes[start_index].block
        gblk = goal.block
        startObs = []
        startObs.append(Nodes[start_index].idx)
        print startObs,Nodes[start_index].block,Nodes[start_index].idx
        
        CDM=FastMarching.MSfm3D_SeqPar(startObs,sblk, goal,gblk, vmap, GhostSet, d_par,True, False, state, block)
       
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
        plannedPath = OFPSearch.Gradient_3D(start_index, goal_index, CDM.nodesList,Nodes, 3)
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
