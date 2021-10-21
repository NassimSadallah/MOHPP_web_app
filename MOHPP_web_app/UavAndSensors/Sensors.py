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
        try:
            
            myports = [tuple(p) for p in list(serial.tools.list_ports.comports())]
            for p in myports:
               
                if 'ttyUSB' in p[0]:
                    port = p
            ser = serial.Serial(port[0], 9600, timeout=1)
            ser.flush()
            print('USB serial communication activated on port!', port[0])
            return ser
        except:
            self.ser =None
            print('No active USB port is detected')
            

    def getSensorsValues(self):
        sensorsValues = []
        if self.ser != None:
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
    def sensArea(self):
               
        for sensor in self.sensors:
            distanceRound = 0
            i = 5
        
            print "Measurement started for " + sensor['ID'] + ", Ctrl+z to cancle the measurement";
                
            while i>0:
                GPIO.output( sensor['TRIG'], GPIO.LOW);
        
                time.sleep(MEASURE_INTERVAL_TIME); #DELAY
        
                GPIO.output(sensor['TRIG'], GPIO.HIGH);
        
                time.sleep(SENSOR_SETTLE_TIME);
        
                GPIO.output(sensor['TRIG'], GPIO.LOW);
        
                while GPIO.input(sensor['ECHO']) == 0:
                    pulse_start = time.time();
        
                while GPIO.input(sensor['ECHO']) == 1:
                    pulse_end = time.time();
        
                pulse_duration = pulse_end - pulse_start;
        
                distance = pulse_duration * MEASURE_REFERENCE;
                distanceRound += round(distance, 2);
                i-=1
                
            distanceRound = distanceRound/5.0  
            sensor['VALUE'] = int(floor(distanceRound))  
        
        return self.sensors
    




UavOrient = 0#degToRad(0)
    
while True:
    if ser.in_waiting >6:
        line = ser.readline().decode('utf-8').rstrip()
        splitline = line.split(',')

        for l in splitline:
            sensorsValues.append(round(int(l)*0.01,2))  

 
        sN= [int(round(int(sensorsValues[0])*cos(2*pi-UavOrient+pi/2),0)),-int(round(int(sensorsValues[0])*sin(2*pi-UavOrient+pi/2),0))]
        
        sE= [int(round(int(sensorsValues[1])*cos(2*pi-UavOrient),0)),int(round(int(sensorsValues[1])*sin(2*pi-UavOrient),0))]
        
        sS= [int(round(int(sensorsValues[2])*cos(2*pi-UavOrient+pi*3/2),0)),-int(round(int(sensorsValues[2])*sin(2*pi-UavOrient+pi*3/2),0))]
        
        sW= [int(round(int(sensorsValues[3])*cos(2*pi-UavOrient+pi),0)),int(round(int(sensorsValues[3])*sin(2*pi-UavOrient+pi),0))]
        
        sNE=[int(round(int(sensorsValues[4])*cos(2*pi-UavOrient+pi/4),0)),-int(round(int(sensorsValues[4])*sin(2*pi-UavOrient+pi/4),0))] 
        sSE=[int(round(int(sensorsValues[5])*cos(2*pi-UavOrient+7*pi/4),0)),-int(round(int(sensorsValues[5])*sin(2*pi-UavOrient+7*pi/4),0))] 
        sSW=[int(round(int(sensorsValues[6])*cos(2*pi-UavOrient+5*pi/4),0)),-int(round(int(sensorsValues[6])*sin(2*pi-UavOrient+5*pi/4),0))] 
        sNW=[int(round(int(sensorsValues[7])*cos(2*pi-UavOrient+3*pi/4),0)),-int(round(int(sensorsValues[7])*sin(2*pi-UavOrient+3*pi/4),0))]  

        print  sensorsValues[4:8],  sNE, sSE, sSW, sNW       
    
        sensorsValues = []    
    
    time.sleep(1)
        
        
"""        
        

