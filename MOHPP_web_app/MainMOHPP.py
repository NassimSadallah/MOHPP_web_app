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
from mohpp.utilities import  DetectUnexpectedObs, height, width, d_,UavHeading
from UavAndSensors.Sensors import Sensors
from math import floor

sitl_connect ='127.0.0.1:14550'
real_connect ='/dev/ttyAMA0' 
UAV = None
start_coordinates, goal_coordinates = [60,96],[105,35]#EastWest/NorthSouth
nextStep, current = [-1.0, -1.0],[-1.0, -1.0]
plannedPath, Nodes, CDM,globalPath, extendedObs, start_index, goal_index = [],[],[],[], [], -1,-1
extendedObs = []
my_options = {
    'mode': "chrome-app", #or "chrome-app",
    'port': 8000

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
    
    return location, armable, battery, mode, head, alt, round(spd,2),sensor.health

@eel.expose
def Launch():
    
    global plannedPath, Nodes, goal_index, CDM,globalPath, start_index, extendedObs     
    alpha = 0.01* int(str(eel.saturation()()))

    #reads the binary map from the binarymaps package
    binMap = os.path.join(os.path.dirname(os.path.abspath(__file__))+"/binarymaps/Cerist.png")
    
    #defines the obstacles and return the corresponding indexed node list
    Nodes, srcObs, block = MAP.processMap(width, height, binMap, seq =1, nbr_blocks=25)
    
    #gets the corresponding indices of the cells start and goal
    start_index = utilities.coordinatesToIndex(start_coordinates, d_)
    goal_index = utilities.coordinatesToIndex(goal_coordinates, d_)
    
    globalPath = []
    #computes the velocity and the travel time at each node of the map
    CDM = CDMap.get_Vel_Cost(Nodes, srcObs, start_index, [goal_index], alpha, 1.0, d_, seq=1,block=block)
    #wavePlot(d_[0], d_[1]/d_[0], CDM)
    plannedPath = OFPSearch.Gradient(start_index, goal_index, CDM, d_)
    return 'Done !'

@eel.expose
def take_off(h):
    
    global plannedPath, Nodes, goal_index, CDM,globalPath,start_index, extendedObs 
    
    VeMeth.UAV().takeoff(float(str(h)), UAV)
    heading = UavHeading(UAV.heading)
    default_alt = UAV.location.local_frame.down
    
    nextStep = plannedPath[0]
    plannedPath.remove(plannedPath[0])
    nodeIdx = Nodes[utilities.coordinatesToIndex([int(floor(nextStep[0])),int(floor(nextStep[1]))], d_)].indice
    extendedObs, isDetected, brake = DetectUnexpectedObs(sensorType, sensor.getSensorValues(sensorType), UAV.heading, nodeIdx, Nodes, extendedObs, 1.5, 4, d_)
    
    return extendedObs
    
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
            #pause the script for the corresponding travel time
            time.sleep(utilities.sqrt_dist(n, e, d))
            nex = input('next')
            
        else:
            #we recompute the global path once we bypassed the dynamic threats
            current = nextStep
            plannedPath = OFPSearch.Gradient(nodeIdx, goal_index, CDM, d_)
            
        globalPath.append(nextStep)
        
        #sensing the surrounding area with the embedded sensors
        extendedObs, isDetected, brake = DetectUnexpectedObs(sensorType, sensor.getSensorValues(sensorType), UAV.heading, nodeIdx, Nodes, extendedObs, 1.5, 4, d_)
        #print isDetected, brake, len(extendedObs)
        #brake = False
        if brake:#if brake is triggered, we must switch to online process
            print 'Into ONPS ...'
            nextStep = ONPSearch.processONPS(nodeIdx, goal_index, heading,extendedObs, isDetected, Nodes, d_, sensor.getSensorValues(sensorType), sensor, UAV, default_alt)
            isReplanning = True
        
        nodeIdx = Nodes[utilities.coordinatesToIndex([int(floor(nextStep[0])),int(floor(nextStep[1]))], d_)].indice
    time.sleep(1)
    #Land(UAV)
    
    return 'ok'
    
@eel.expose
def dataLecture():
    
    location = [UAV.location.global_frame.lat, UAV.location.global_frame.lon] 
    battery = UAV.battery.level
    alt = UAV.location.global_relative_frame.alt
    spd = UAV.airspeed
    head = UAV.heading
    mode = UAV.mode.name
    theta=sensor.getSensorValues(sensorType)#, dist = sens.getLidarValues()
    print theta
    #print theta#, dist
    return location, battery, head, round(alt,1), round(spd,2), mode, theta#, dist

@eel.expose
def Land(uav):
    print('Landing mode activated ...')
    VeMeth.UAV.Land(uav)
    

@eel.expose
def testLidar(usb):
    from rplidar import RPLidar
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib.animation as animation
    PORT_NAME =usb# '/dev/ttyUSB'
    DMAX = 7000
    IMIN = 0
    IMAX = 100
    def update_line(num, iterator, line):
        scan = next(iterator)
        offsets = np.array([(np.radians(meas[1]), meas[2]) for meas in scan])
        line.set_offsets(offsets)
        intens = np.array([meas[0] for meas in scan])
        #print intens
        #print offsets
        
        line.set_array(intens)
        #time.sleep(2)
        return line,
    def run():
        lidar = RPLidar(PORT_NAME)
        print lidar.get_info()
        fig = plt.figure()
        ax = plt.subplot(111, projection='polar')
        line = ax.scatter([0, 0], [0, 0], s=5, c=[IMIN, IMAX],
                               cmap='Greys_r', lw=4)
        ax.set_rmax(DMAX)
        ax.grid(True)
        iterator = lidar.iter_scans(scan_type='normal')
        ani = animation.FuncAnimation(fig, update_line,
            fargs=(iterator, line), interval=1)
        plt.show()
        lidar.stop()
        lidar.disconnect()
    run()

#testLidar(sens.ser)

@eel.expose
def stopMOHPP():
    return sys.exit(0)

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
