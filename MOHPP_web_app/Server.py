'''
Created on Oct 13, 2021

@author: nassim
'''
import eel, time

my_options = {
    'mode': "None", #or "chrome-app",
    'host': 'localhost',
    'port': 8080,
}
eel.init('webapp')



@eel.expose
def getTimess():
    return time.strftime('%c')
@eel.expose    
def send_coordinates(uav):
    return uav.location.global_relative_frame.lon,uav.location.global_relative_frame.lat






eel.start('miniGCS.html', options=my_options, suppress_error=True)


