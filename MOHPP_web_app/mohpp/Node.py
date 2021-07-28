'''
Created on Jul 8, 2021

@author: nassim
'''
class node():# class node
    def __init__(self,x,y):
        
        self.colonne = x
        self.abscice = y # elle va s'incrementer sur l'axe x (les case des colonnes )
        self.long = None
        self.lat = None
        self.OBSTACLE = 0
        
        self.cost = 99999999999.0
        self.type ='F'
        self.v = 1.0
        self.block = 0 # the hosting block
        self.indice = -1
        self.cros = False
        self.TAG = None
        
        self.G = 0
        self.F = 0
        self.H = 0
        self.parent = self
        self.visited = False
        
    