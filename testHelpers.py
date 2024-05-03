# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 09:06:13 2024

Helper functions for testing growCircle

@author: samda
"""
from AtomicPoint import AtomicPoint
import matplotlib.pyplot as plt


def nearestAtom(a, pointList, excluded):
    '''Function finds the nearest AtomicPoint to point a, ignoring any atoms
    given in the excluded list. VanDerWaals radius is acounted for.
    
    @peram a AtomicPoint
    @peram pointList of AtomicPoints to search
    @peram excluded list of AtomicPoints to ignore
    @return Nearest AtomicPoint to a ignoring excluded atoms
    '''
    
    #add **kwargs for distBetween and distBetween2D options?
    
    nearest_dist = 100000
    nearest_atom = AtomicPoint(0, 0, 0, "")

    #Checks all AtomicPoints in the given list
    for p in pointList:
        exclusion_check = False #Can't break outer for loop from inner for loop
        #Checks if the AtomicPoint is excluded and if so skips it
        for e in excluded:
            if p == e: 
                exclusion_check = True
                break #Don't need to continue loop afer match is found
        if exclusion_check: continue
        
        if a.distBetween2D(p, True) < nearest_dist:
            nearest_atom = p
            nearest_dist = a.distBetween2D(p, True)
            
    return nearest_atom


def plotAtom(p, ax, **kwargs):
    ''' 
    Function helps automate creation and ploting of circle patchs when graphing 
    atoms. Circles are ploted onto sublot ax. Atoms are assumed to have been 
    projected into 2D using planerize and planerizePoint from lcfProcessing.
    
    @peram p AtomicPoint to be graphed
    @peram ax plt subplot circles are being attached to 
    **kwargs 
        color : change color of the circle, default red
        radius : change radius ofr circle, default p.VanDerWaalsRadius
    '''
    c_color = 'r'
    radius = p.VanDerWaalsRadius
    
    if 'color' in kwargs:
        c_color = kwargs['color']
    
    if 'radius' in kwargs:
        radius = kwargs['radius']
    
    circle = plt.Circle((p.x, p.y), radius, color=c_color)
    ax.add_patch(circle)
    
    
def atomCollide(a, pointList, **kwargs):
    '''
    Function checks if given AtomicPoint intersects with any other atoms in the
    pointList. Returns True if it does so, false otherwise.
    
    @peram a AtomicPoint being checked
    @peram pointList List of AtomicPoints that p is being chekd against
    **kwargs
        count : Boolean, if true count the number intersecting atoms
    @return True if p itersects with any of the atoms in the pointList, returns 
    false otherwise
        If 'count' = True, then answer will be returned as [Boolean, int count of atoms intersecting]
    '''
    if 'count' in kwargs:
        if kwargs['count'] == True:
            iCount = 0
            intersecting = False
            
            for p in pointList:
                if a.distBetween(p, True) < 0: 
                    intersecting = True
                    iCount += 1
            
            if intersecting: return [True, iCount]
            else: return [False, 0]
    
    for p in pointList:
        if a.distBetween(p, True) < 0: return True
    return False
    

