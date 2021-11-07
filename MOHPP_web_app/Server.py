'''
Created on Oct 13, 2021

@author: nassim
'''

import eel

from mohpp import MAP, CDMap, utilities, OFPSearch, ONPSearch
from UavAndSensors import VehiclesMethods as VeMeth
from mohpp.utilities import  DetectUnexpectedObs, height, width, d_,UavHeading
from UavAndSensors.Sensors import Sensors
import os,time
from math import floor


nextStep, current = [-1.0, -1.0],[-1.0, -1.0]
extendedObs = []

lastSensorsValues = [0 for i in range(8)]
sensedArea = 0#Sensors.Sensors()
uav = None
my_options = {
    'mode': "None", #or "chrome-app",
    'port': 8000

}

eel.init('webapp')

@eel.expose
def getTVal():
    while True:
        sens= SensorsValues()
        
        eel.sleep(0.1)
        return sens

@eel.expose
def SensorsValues():
    while True:    
        values =  sensedArea.getSensorsValues()        
        
        eel.sleep(.01)
        return values

@eel.expose
def connect():
    sitl_connect ='127.0.0.1:14550'
    real_connect ='/dev/ttyAMA0' 
    UAV = VeMeth.UAV().connect_to_vehicle(sitl_connect, 921600)
    location = UAV.location.global_frame
    battery = UAV.battery
    alt = UAV.location.local_frame.down
    spd = UAV.airspeed
    lst_hBeat = UAV.last_heartbeat
    head = UAV.heading
    mode = UAV.mode.name
    armable = UAV.is_armable
    return location, armable, battery, mode, lst_hBeat, head, alt, spd

def init():
            
    eel.spawn(SensorsValues)

    eel.start('index.html', my_options, block = False)
#192.168.0.1:8080/miniGCS.html
eel.start('index.html', my_options, block = True)