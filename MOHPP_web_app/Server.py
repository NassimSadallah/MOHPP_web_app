'''
Created on Oct 13, 2021

@author: nassim
'''

import eel, time
from UavAndSensors import Sensors
from threading import Thread

Thread = None
lastSensorsValues = [0 for i in range(8)]
sensedArea = Sensors.Sensors()

my_options = {
    'mode': "None", #or "chrome-app",

}

eel.init('webapp')

@eel.expose
def getTVal():
    while True:
        sens= SensorsValues()
        
        eel.sleep(0.1)
        return sens


@eel.expose    
def send_coordinates(uav):
    return uav.location.global_relative_frame.lon,uav.location.global_relative_frame.lat


#@eel.expose
def SensorsValues():
    while True:    
        values =  sensedArea.getSensorsValues()        
        
        eel.sleep(.01)
        return values

def init():
            
    eel.spawn(SensorsValues)

    eel.start('miniGCS.html', my_options, block = False)
#192.168.0.1:8080/miniGCS.html
