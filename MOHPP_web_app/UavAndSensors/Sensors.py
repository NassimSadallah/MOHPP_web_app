'''
Created on Jul 8, 2021

@author: nassim
'''

import RPi.GPIO as GPIO
import time
from threading import Thread
from multiprocessing import Process
from math import floor


GPIO.setmode(GPIO.BCM)

GPIO.setwarnings(False) # for disable warnings in terminal

# time for sensor to settle
SENSOR_SETTLE_TIME = 0.00001

MEASURE_INTERVAL_TIME = 0.1 # time delay to measure (min 15miliseconds)                 

# max distance threshold for sensors to react (in cm)
MAX_DISTANCE_THRESHOLD = 500.0

# Speed of sound at sea level = 343 m/s or 34300 cm/s
MEASURE_REFERENCE = 17150


class Sensors(object):
    
    def __init__(self):
        
        self.sensors = []
        '''
        sensors with pin configuration, we define the trig, echo pins 
        and ID for the embedded sensors (8 ultrasonic hc-rs04)  
        '''
        self.sensor1 = {'ID': 's1', 'TRIG': 27, 'ECHO': 12, 'VALUE': -1.0}
        self.sensors.append(self.sensor1) # add to the list
        
        self.sensor2 = {'ID': 's2', 'TRIG': 22, 'ECHO': 16, 'VALUE': -1.0}
        self.sensors.append(self.sensor2) # add to the list
        
        self.sensor3 = {'ID': 's3', 'TRIG': 10, 'ECHO': 20, 'VALUE': -1.0}
        self.sensors.append(self.sensor3) # add to the list
        
        self.sensor4 = {'ID': 's4', 'TRIG': 9, 'ECHO': 21, 'VALUE': -1.0}
        self.sensors.append(self.sensor4) # add to the list
        
        self.sensor5 = {'ID': 's5', 'TRIG': 2, 'ECHO': 24, 'VALUE': -1.0}
        self.sensors.append(self.sensor5) # add to the list
        
        self.sensor6 = {'ID': 's6', 'TRIG': 3, 'ECHO': 25, 'VALUE': -1.0}
        self.sensors.append(self.sensor6) # add to the list
        
        self.sensor7 = {'ID': 's7', 'TRIG': 4, 'ECHO': 8, 'VALUE': -1.0}
        self.sensors.append(self.sensor7) # add to the list
        
        self.sensor8 = {'ID': 's8', 'TRIG': 17, 'ECHO': 7, 'VALUE': -1.0}
        self.sensors.append(self.sensor8) # add to the list
        

        self.initSensors()
        print('initialization of pins done !')


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
    
    
    def initSensors(self):
        
        if len(self.sensors) > 0:
            for sensor in self.sensors:
        
                GPIO.setup( sensor['ECHO'], GPIO.IN )#Sensor's echo pins shoud be in
                
                GPIO.setup( sensor['TRIG'], GPIO.OUT )#Sensor's trig pins should be out
        
        return True
