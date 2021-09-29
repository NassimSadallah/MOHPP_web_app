'''
Created on Jul 8, 2021

@author: nassim
'''
import RPi.GPIO as GPIO
import time
from threading import Thread
from multiprocessing import Process


GPIO.setmode(GPIO.BCM)

GPIO.setwarnings(False) # for disable warnings in terminal

# time for sensor to settle
SENSOR_SETTLE_TIME = 0.00001

MEASURE_INTERVAL_TIME = 0.1 # time delay to measure (min 15miliseconds)                 

# max distance threshold for sensors to react (in cm)
MAX_DISTANCE_THRESHOLD = 5.0

# Speed of sound at sea level = 343 m/s or 34300 cm/s
MEASURE_REFERENCE = 17150


class Sensors(object):
    
    def __init__(self):
        
        self.sensors = []
        '''
        sensors with pin configuration, we define the trig, echo pins 
        and ID for the embedded sensors (8 ultrasonic hc-rs04)  
        '''
        self.sensor1 = {'ID': 'sensor1', 'TRIG': 2, 'ECHO': 24}
        self.sensors.append(self.sensor1) # add to the list
        
        self.sensor2 = {'ID': 'sensor2', 'TRIG': 3, 'ECHO': 25}
        self.sensors.append(self.sensor2) # add to the list
        
        self.sensor3 = {'ID': 'sensor3', 'TRIG': 4, 'ECHO': 8}
        self.sensors.append(self.sensor3) # add to the list
        
        self.sensor4 = {'ID': 'sensor4', 'TRIG': 17, 'ECHO': 7}
        self.sensors.append(self.sensor4) # add to the list
        
        self.sensor5 = {'ID': 'sensor5', 'TRIG': 27, 'ECHO': 12}
        self.sensors.append(self.sensor5) # add to the list
        
        self.sensor6 = {'ID': 'sensor6', 'TRIG': 22, 'ECHO': 16}
        self.sensors.append(self.sensor6) # add to the list
        
        self.sensor7 = {'ID': 'sensor7', 'TRIG': 10, 'ECHO': 20}
        self.sensors.append(self.sensor7) # add to the list
        
        self.sensor8 = {'ID': 'sensor8', 'TRIG': 9, 'ECHO': 21}
        self.sensors.append(self.sensor8) # add to the list

        self.initSensors()


def sensArea(self):

    for sensor in self.sensors:
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
            distanceRound = round(distance, 2);
            i-=1
            
        if(distanceRound < MAX_DISTANCE_THRESHOLD):
            print "Distance of sensor "+ sensor['ID'] + " : ", distanceRound, "cm";
                
        else:
            print "Out of range "+ sensor['ID'] + " : ", distanceRound, "cm";
    
    return 0


def initSensors(self):
    
    if len(self.sensors) > 0:
        for sensor in self.sensors:
    
            GPIO.setup( sensor['ECHO'], GPIO.IN )#Sensor's echo pins shoud be in
            
            GPIO.setup( sensor['TRIG'], GPIO.OUT )#Sensor's trig pins should be out
    
    
    
    
    
    
    
    
    
    
    