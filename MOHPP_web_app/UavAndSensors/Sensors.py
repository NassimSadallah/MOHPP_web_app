'''
Created on Jul 8, 2021

@author: nassim
'''

import time   
import serial, serial.tools.list_ports


class Sensors(object):
    
    def __init__(self):
        
        self.ser = self.initSensors()

    def initSensors(self):
        print('looking for USB communication ...')

        myports = [tuple(p) for p in list(serial.tools.list_ports.comports())]
        for p in myports:
            
            if 'ttyUSB' in p[0]:
                port = p
                ser = serial.Serial(port[0], 9600, timeout=1)
                ser.flush()
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
