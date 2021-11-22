'''
Created on Jul 8, 2021

@author: nassim
'''

import time,sys   
import serial, serial.tools.list_ports
from rplidar import RPLidar
import numpy as np
import eel


incr = 0    

class Sensors(object):
    
    def __init__(self, sensor = None):
        
        if sensor =='lidar':
            self.ser = self.initSensors(sensor)
            self.iterator = 0
            try:
                self.lidar = RPLidar(self.ser)
                print(self.lidar.get_health())
                self.iterator=self.lidar.iter_scans(scan_type='normal', max_buf_meas=10000)
                
                #self.getLidarValues(self.lidar)
            except:
            
                sys.exit('No Lidar - PROCESS ABORTED')
            
        elif sensor =='hcsr04':
            self.ser = self.initSensors(sensor)  
            
    def initSensors(self, sen):#sen for sensor type: hs-rc04 or Lidar
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

    '''
    return the sensed distances: value of sensor is either lidar or hcsr04
    '''
    def getSensorValues(self, sensor):
        
        if sensor =='lidar':
            theta={}
            scan = next(self.iterator)

            for meas in scan:#np.array([meas[1] for meas in scan])
                theta[meas[1]] = meas[2]
            
            return theta#, Dist
        
        elif sensor=='hcsr04':
            
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
                exit('No serial connection') 

    """    
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
    """
            

                