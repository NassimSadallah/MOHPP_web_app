'''
Created on Jul 8, 2021

@author: nassim
'''
from UavAndSensors import VehiclesMethods as VeMethods





def getReadyUAV(connection, baud_rate):

    vehicle = VeMethods.UAV().connect_to_vehicle('127.0.0.1:14550', baud_rate=921600)
    VeMethods.UAV().init_vehicle(vehicle)
    print '          <<  Ready to arm !  >> '
    return vehicle 

def cmdTakeOff(h, UAV):
    VeMethods.UAV().takeoff(h, UAV)
    

def move_InDirections(velx, vely, velz, vehicle):
    VeMethods.UAV().send_NED_velocity(velx, vely, velz, vehicle)
    
