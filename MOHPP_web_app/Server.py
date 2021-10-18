'''
Created on Oct 13, 2021

@author: nassim
'''
import eel, time
from UavAndSensors import Sensors
from threading import Thread

Thread = None

sensedArea = Sensors.Sensors()

my_options = {
    'mode': "chrome", #or "chrome-app",
    'host': '10.0.0.252',
    'port': 8080,
}

eel.init('webapp')



@eel.expose
def getTimess():
    return time.strftime('%c')

@eel.expose    
def send_coordinates(uav):
    return uav.location.global_relative_frame.lon,uav.location.global_relative_frame.lat

@eel.expose
def SensorsValues():
    while True:
        values =  sensedArea.getSensorsValues()
        print values
        return values
        eel.sleep(.2)


eel.start('miniGCS.html', options=my_options, suppress_error=True)
