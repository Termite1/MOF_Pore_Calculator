# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 13:58:11 2024

Checks if there is a line of sight between two points, taking into account obscuring spheres

@author: samda
"""
from AtomicPoint import AtomicPoint
import numpy as np
import math

def getMass(a):
    """
    Returns average mass of atom given the element 

    Parameters
    ----------
    a : AtomicPoint
        AtomicPoint elemental mass is being returned for.

    Returns
    -------
    float
        Average mass of given element.

    """
    
    massDict = {"H": 1.0078,
                "C": 12.011,
                "O": 15.999,
                "ZN": 65.380}
    
    return massDict[a.element]


def centerOfMass(pointList, useMass):
    '''
    Function calculates the (x,y,z) center of mass given a list of AtmicPoints

    Parameters
    ----------
    pointList : List of AtomicPoints
        Structure center of mass is being found for.
    useMass : boolean
        If true, uses mass to calculate center of mass. If false, finds 
        geometric center of mass.

    Returns
    -------
    AtomicPoint location of center of mass (with no vdwr)
    '''
    #If useMass is true
    if useMass:
        #Keeps track of totals
        total = {'x': 0, 'y': 0, 'z': 0}
        massTotal = 0
        
        #adds x,y,y of each atom in pointList to totals. Accounts for mass
        #track mass as well
        for p in pointList:
            temp_mass = getMass(p)
            
            total['x'] += p.x * temp_mass
            total['y'] += p.y * temp_mass
            total['z'] += p.z * temp_mass
        
            massTotal += temp_mass
        
        #Devides by total mass to get center of mass
        total['x'] = total['x']/massTotal
        total['y'] = total['y']/massTotal
        total['z'] = total['z']/massTotal
        
        return AtomicPoint(total['x'], total['y'], total['z'], "")
        
    #if use mass is false
    else:
        #Keeps track of totals
        total = {'x': 0, 'y': 0, 'z': 0}
        
        #adds x,y,y of each atom in pointList to totals
        for p in pointList:
            total['x'] += p.x
            total['y'] += p.y
            total['z'] += p.z
        
        #Devides by number of atoms in list to get average
        total['x'] = total['x']/len(pointList)
        total['y'] = total['y']/len(pointList)
        total['z'] = total['z']/len(pointList)
        
        return AtomicPoint(total['x'], total['y'], total['z'], "")
    
    
def magnitude(v):
    '''
    Calculates magnitude of given vector.

    Parameters
    ----------
    v : numpy vector
        vector.

    Returns
    -------
    float
        vector magnitude
    '''
    return math.sqrt(sum(pow(element, 2) for element in v))
    
    
def lineOfSightCheck(a, b, pointList, **kwargs):
    '''
    Function checks if there is an unobstructed line of sight between points a 
    and b, or if that line of sight is obscured by atoms in the point list

    Parameters
    ----------
    a : AtomicPoint
        Point 1.
    b : AtomicPoint
        Point 2.
    pointList : List of AtomicPoints
        Molecular structure being inspected.
    **kwargs
        ignoreIndex : int
            pointList index to ignore (only takes 1 for now)

    Returns
    -------
    Boolean.
        True if unobstructed line of sight, false otherwise
        Would not be hard to include obstructing atoms
    '''
    
    AB = np.array([b.x-a.x, b.y-a.y, b.z-a.z]) #Vector line between point a and point b
    mag_AB = magnitude(AB) #Vector magnitude
    
    #repeats for every AtomicPoint in the pointList
    for p in pointList:
        #Checks to see if indexes should be ignronered, and if so, skips listed indexes
        if 'ignoreIndex' in kwargs: 
            if pointList.index(p) == kwargs['ignoreIndex']: continue
        
        AP = np.array([p.x-a.x, p.y-a.y, p.z-a.z]) #Vector line between point a and point b
        
        cp = np.cross(AB, AP) #Cross product of AB and AP
        
        #Shortest distance between p and AB
        CD = magnitude(cp) / mag_AB
        
        #Makes shure the distance is larger than the Van Der Waals radius.
        #if it isn't the line passes throgh the atom, and the atom obscures the line of sight
        #between the two points
        if CD < p.VanDerWaalsRadius:
            return False
        
    #If line not obstructed, return True
    return True
        
    
    