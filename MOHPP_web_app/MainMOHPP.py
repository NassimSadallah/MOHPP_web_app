'''
Main python script : First the binary map is extracted to a grid, 
second, the velocity and the travel time cost are assigned with the Fast Marching (in this version only 8 neighbors are considered)
Third, the gradient descent is used to compute the optimal path
Finally, we instantiate a vehicule (either within the SITL or on a real drone and command it
 
Created on Jul 12, 2021

@author: nassim

'''
from mohpp import MAP, CDMap, utilities, OFPSearch
from UavAndSensors import commandsUAV, VehiclesMethods as VeMeth
from mohpp.utilities import wavePlot
import os,time
from math import floor



width, height = 200, 150
d_ = [width, width*height]
start_coordinates, goal_coordinates = [151,57],[20,50]
nextStep, current = [-1.0, -1.0],[-1.0, -1.0]

#reads the binary map from the path
binMap = os.path.join(os.path.dirname(os.path.abspath(__file__))+"/binarymaps/simulation.png")

#defines the obstacles and return the corresponding indexed node list
Nodes, srcObs, block = MAP.processMap(width, height, binMap, seq =1, nbr_blocks=25)

#gets the corresponding indices of the cells start and goal
start_index = utilities.coordinatesToIndex(start_coordinates, d_)
goal_index = utilities.coordinatesToIndex(goal_coordinates, d_)


#computes the velocity and the travel time at each node of the map
CDM = CDMap.get_Vel_Cost(Nodes, srcObs, start_index, [goal_index], 0.3, 1.0, d_, seq=1,block=block)
wavePlot(d_[0], d_[1]/d_[0], CDM)
globalPath = OFPSearch.Gradient(start_index, goal_index, CDM, d_)

'''
connection to vehicle and controlling 
'''

UAV = VeMeth.UAV().connect_to_vehicle('127.0.0.1:14550', 921600)
default_alt = VeMeth.UAV().takeoff(5.0, UAV)

nextStep = globalPath[0]
globalPath.remove(globalPath[0])
nodeIdx = Nodes[utilities.coordinatesToIndex([int(floor(nextStep[0])),int(floor(nextStep[1]))], d_)].indice

while nodeIdx !=goal_index:
    
    current = nextStep
    nextStep = globalPath[0]
    nodeVel = Nodes[utilities.coordinatesToIndex([int(floor(nextStep[0])),int(floor(nextStep[1]))], d_)].v
    globalPath.remove(globalPath[0])
    n, e, d = utilities.getNorth_East_Down(current, nextStep, UAV.location.local_frame.down, default_alt)
    UAV.airspeed = nodeVel
    VeMeth.UAV().send_NED_velocity(n, e,d, UAV)
    print n, e, d, nodeVel
    time.sleep(utilities.sqrt_dist(n, e, d))
    nodeIdx = Nodes[utilities.coordinatesToIndex([int(floor(nextStep[0])),int(floor(nextStep[1]))], d_)].indice
