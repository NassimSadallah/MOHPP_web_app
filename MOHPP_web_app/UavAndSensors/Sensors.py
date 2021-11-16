'''
Created on Jul 8, 2021

@author: nassim
'''

import time   
import serial, serial.tools.list_ports
from rplidar import RPLidar



class Sensors(object):
    
    def __init__(self):
        
        self.ser = self.initSensors('lidar')
        lidar = RPLidar(self.ser)
        self.getLidarValues(lidar)

    def initSensors(self, sen):
        print('looking for USB communication ...')

        myports = [tuple(p) for p in list(serial.tools.list_ports.comports())]
        for p in myports:
            
            if 'ttyUSB' in p[0]:
                port = p
                ser = None
                if sen=='hcsr04':
                    ser = serial.Serial(port[0], 9600, timeout=1)
                    ser.flush()
                    
                elif sen=='lidar':
                    ser = port[0]
                
                print('USB serial communication activated on port!', port[0])
                return ser
   
        ser = None
        print('No active USB port is detected')
        return ser
            

    def getSensorsValues(self):
        
        if self.ser != None:
            sensorsValues = []
            while True:
                if self.ser.in_waiting >6:
                    line = self.ser.readline().decode('ISO-8859-1').rstrip()
                    splitline = line.split(',')
          
                    for l in splitline:
                        
                        sensorsValues.append(round(int(l)*0.01,2))  
                        
                    return sensorsValues    
        else:
            exit            

    def getLidarValues(self, lidar):
        
        while True:
            info = lidar.get_info()
            print(info)
        
        health = lidar.get_health()
        print(health)
        
        for i, scan in enumerate(lidar.iter_scans()):
            print('%d: Got %d measurments' % (i, len(scan)))
            #print(lidar.iter_scans())
            if i > 500:
                break
        
        lidar.stop()
        lidar.stop_motor()
        lidar.disconnect()        
                