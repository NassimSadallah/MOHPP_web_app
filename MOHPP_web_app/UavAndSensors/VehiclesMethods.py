'''
Created on Jul 8, 2021

@author: nassim
'''
from __future__ import print_function
from dronekit import connect, VehicleMode#, LocationGlobal, LocationGlobalRelative
from pymavlink import mavutil
#import os
import time
from __builtin__ import True

class UAV:
    '''
    Methods to interract with the connected UAV
    '''
    def connect_to_vehicle(self,connection_string,baud_rate):
        print('connection to UAV on ', connection_string)
        vehicle = connect(connection_string,wait_ready=True, baud=baud_rate)
        self.init_vehicle(vehicle)
            
        return vehicle
    
    def init_vehicle(self,vehicle):

        while not vehicle.is_armable:
            print(" Waiting for vehicle initialization...")
            time.sleep(1.5)
        
        vehicle.mode = VehicleMode("GUIDED")
        print(' setting UAV mode to GUIDED...')
        while not vehicle.mode =="GUIDED":
            
            vehicle.mode = VehicleMode("GUIDED")
            time.sleep(1)
        
    def takeoff(self,height,vehicle):
        '''
        Arming the vehicle, height in meters and vehicle is the instantiated object
        '''
        vehicle.armed = True
        while not vehicle.armed:
            print('Waiting for arming...')
            time.sleep(1.5)
        
        '''
        takeoff to the altitude h
        '''
        vehicle.simple_takeoff(height)
        while vehicle.location.global_relative_frame.alt <=height*0.98:
            print('Actual altitude:',vehicle.location.global_relative_frame.alt)
            if not vehicle.armed: 
                print('MOTORS DISARMED...Aborting Mission')
                break
            time.sleep(1.5)    
        return round(vehicle.location.local_frame.down,5)
    
    def Land(self,vehicle):
        vehicle.mode = VehicleMode("LAND")
        
    def send_NED_velocity(self,vel_x, vel_y, vel_z,vehicle):    
        '''
        Move vehicle in direction based on specified velocity vectors (x, y, z).
        Vehicle is the instantiated object.
        '''
        
        msg = vehicle.message_factory.set_position_target_local_ned_encode(
            0,       # time_boot_ms (not used)
            0, 0,    # target system, target component
            mavutil.mavlink.MAV_FRAME_LOCAL_NED, # frame
            0b0000111111000111, # type_mask (only speeds enabled)
            0, 0, 0, # x, y, z positions (not used)
            vel_x, vel_y, vel_z, # x, y, z velocity in m/s
            0, 0, 0, # x, y, z acceleration (not supported yet, ignored in GCS_Mavlink)
            0, 0)    # yaw, yaw_rate (not supported yet, ignored in GCS_Mavlink)
        
        vehicle.send_mavlink(msg)

        
        
