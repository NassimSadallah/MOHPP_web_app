'''
Main python script : First the binary map is extracted to a grid, 
second, the velocity and the travel time cost are assigned with the Fast Marching 
(in this version only 8/16 neighbors are considered)
Third, the gradient descent is used to compute the optimal path
Finally, we instantiate a vehicle (either within the SITL or on a real drone and command it
 
Created on Jul 12, 2021

@author: nassim

'''

from mohpp import MAP, CDMap, utilities, OFPSearch, ONPSearch
from UavAndSensors import VehiclesMethods as VeMeth
from mohpp.utilities import  DetectUnexpectedObs, height, width, d_,UavHeading
from UavAndSensors.Sensors import Sensors
import os,time
from math import floor
import Server


sitl_connect ='127.0.0.1:14550'
real_connect ='/dev/ttyAMA0' 

start_coordinates, goal_coordinates = [151,57],[20,50]
nextStep, current = [-1.0, -1.0],[-1.0, -1.0]
extendedObs = []


#reads the binary map from the binarymaps package
binMap = os.path.join(os.path.dirname(os.path.abspath(__file__))+"/binarymaps/simulation.png")

#defines the obstacles and return the corresponding indexed node list
Nodes, srcObs, block = MAP.processMap(width, height, binMap, seq =1, nbr_blocks=25)

#gets the corresponding indices of the cells start and goal
start_index = utilities.coordinatesToIndex(start_coordinates, d_)
goal_index = utilities.coordinatesToIndex(goal_coordinates, d_)

globalPath = []
#computes the velocity and the travel time at each node of the map
CDM = CDMap.get_Vel_Cost(Nodes, srcObs, start_index, [goal_index], 0.3, 1.0, d_, seq=1,block=block)
#wavePlot(d_[0], d_[1]/d_[0], CDM)
plannedPath = OFPSearch.Gradient(start_index, goal_index, CDM, d_)


'''
connection to vehicle and controlling SITL:('127.0.0.1:14550', 921600), Real: ('/dev/ttyAMA0', baud=921600) 
'''

UAV = VeMeth.UAV().connect_to_vehicle(sitl_connect, 921600)
Server.init()
time.sleep(3)
Server.getTVal()
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
