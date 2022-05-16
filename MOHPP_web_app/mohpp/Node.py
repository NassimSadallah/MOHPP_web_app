'''
Created on Jul 8, 2021

@author: nassim
'''
from utilities import FAR, INFINI
class node():# class node
    def __init__(self,x,y):
        
        self.colonne = x # height
        self.abscice = y #width
        self.long = None
        self.lat = None
        self.OBSTACLE = 0
        
        self.cost = INFINI
        self.type =FAR
        self.v = 1.0
        self.block = 0 # the containing block
        self.indice = -1
        self.cros = False
        self.TAG = None
        self.risk = 0
        self.full = 0
        self.hCost = 0
        self.G = 0
        self.F = 0
        self.H = 0
        self.parent = self
        self.visited = False
        
class node3():
    def __init__(self, x, y, z):
        self.x = x#depth
        self.y = y#width
        self.z = z#height
        #longitude and latitude coordinates
        self.long = None
        self.lat = None
        self.TAG = None
        self.OBSTACLE = 0
        self.cost = INFINI
        self.type =FAR
        self.v = 1.0
        self.block = 0 # the containing block
        self.indice = -1
        self.cros = False
        self.safety_index = 0
        self.full = 0
        self.hCost = 0
        
        self.G = 0
        self.F = 0
        self.H = 0
        self.parent = self
        self.visited = False
