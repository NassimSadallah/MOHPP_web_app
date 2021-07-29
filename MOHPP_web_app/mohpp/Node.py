'''
Created on Jul 8, 2021

@author: nassim
'''
class node():# class node
    def __init__(self,x,y):
        
        self.colonne = x # height
        self.abscice = y #width
        self.long = None
        self.lat = None
        self.OBSTACLE = 0
        
        self.cost = 99999999999.0
        self.type ='F'
        self.v = 1.0
        self.block = 0 # the containing block
        self.indice = -1
        self.cros = False
        self.TAG = None
        self.full = 0
        self.hCost = 0
        self.G = 0
        self.F = 0
        self.H = 0
        self.parent = self
        self.visited = False
        
    