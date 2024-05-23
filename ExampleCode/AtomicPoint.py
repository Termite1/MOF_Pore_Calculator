# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 19:10:37 2023

@author: samda
"""

import math
import VanDerWaalsRadius as vdwrf

class AtomicPoint:
    
    #AtomicPoint Constructor
    
    """
    Constructor for the AtomicPoint class
    
    @param x X float coordinate of atomic point
    @param y Y float coordinate of atomic point
    @param z Z float coordinate of atomic point
    @param element String atomic symbol of element, all caps
    """
    def __init__(self, x, y, z, element, **kwargs):
        self.x = x
        self.y = y
        self.z = z
        self.element = element
        self.VanDerWaalsRadius = vdwrf.findRadius(element)
        
        if 'bondType' in kwargs:
            self.bondType = kwargs['bondType']
            if self.bondType == '':
                self.bondType = 'H_'
        else: 
            self.bondType = 'H_'
            
        if 'vdwr' in kwargs:
            self.VanDerWaalsRadius = kwargs['vdwr']
    

    """
    Calculates the distance between two given AtomicPoints
    
    @param self One AtomicPoint point of calculation
    @param point2 Second AtomicPoint of calculation
    @param includeVanDerWaalRadius Boolean indicateing whether or not the Van 
           Der Waal radius of the ppoints should be taken into consideration
    @return float Distance between the two atomic points
    """
    def distBetween(self, point2, includeVanDerWaalRadius):
        if(includeVanDerWaalRadius == True):
            temp_x = self.x - point2.x 
            temp_y = self.y - point2.y
            temp_z = self.z - point2.z
            return (math.hypot(temp_x,temp_y,temp_z) - (self.VanDerWaalsRadius + point2.VanDerWaalsRadius))
        elif(includeVanDerWaalRadius == False):
            temp_x = self.x - point2.x
            temp_y = self.y - point2.y
            temp_z = self.z - point2.z
            return math.hypot(temp_x,temp_y,temp_z)
        
    """
    Calculates the distance between two given AtomicPoints, but only in the 
    (x,y) plane.
    
    @param self One AtomicPoint point of calculation
    @param point2 Second AtomicPoint of calculation
    @param includeVanDerWaalRadius Boolean indicateing whether or not the Van 
           Der Waal radius of the ppoints should be taken into consideration
    @return float Distance between the two atomic points
    """
    def distBetween2D(self, point2, includeVanDerWaalRadius):
        if(includeVanDerWaalRadius == True):
            temp_x = self.x - point2.x 
            temp_y = self.y - point2.y
            return (math.hypot(temp_x,temp_y) - (self.VanDerWaalsRadius + point2.VanDerWaalsRadius))
        elif(includeVanDerWaalRadius == False):
            temp_x = self.x - point2.x
            temp_y = self.y - point2.y
            return math.hypot(temp_x,temp_y)
        
    """
    Prints coordinates of an AtomicPoint.
    
    @param self The AtomicPoint
    """    
    def printCoords(self):
        print(f'({self.x:.3f},{self.y:.3f},{self.z:.3f}) with vdwr = {self.VanDerWaalsRadius:.3f}')