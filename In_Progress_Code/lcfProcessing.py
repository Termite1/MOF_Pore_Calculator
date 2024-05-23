# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 13:42:29 2023

@author: samda
"""

from AtomicPoint import AtomicPoint
import numpy as np
import math


def normalize(v):
    '''
    Function normailzes a given vector and returns it.
    
    @peram v Vector in list form
    @peram Unit vector in direction of v in the form of a list
    '''
    magnitude = 0
    #calculates magnitude of vector from components, runs for all components
    for component in v:
        magnitude += (float(component)**2)
    final_magnitude = math.sqrt(magnitude)
    #Devides each component of the vector by the vector's magnitude
    normalized_vector = [x/final_magnitude for x in v]
    return normalized_vector
    

def dot_product(a, b):
    '''
    Function calculates the dot product of two vectors (lists) or equal length
    
    @peram a Vector 1
    @peram b Vector 2
    @return Dot product as list
    '''
    total = 0
    for i in range (len(a)):
        total += (a[i] * b[i])
    return total
    

def planarize(a, b, c):
    '''
    Function takes 3 points with x, y, y components and finds the 2d plane that intersects them all.
    Returns the principle values to convert (x, y, z) points into points in the new plane. Updated
    function does not use matric calculations.
    
    Based on this: https://www.quora.com/What-is-the-fastest-way-to-find-the-equation-of-a-plane-given-three-points
    
    @param a Point 1 with (x, y, z) components
    @param b Point 2 with (x, y, z) components
    @param c Point 3 with (x, y, z) components
    @return [Plane Origin as AtomicPoint, Surface Normal Unit Vector, Axis Unit Vector 1, Othoganal Axis Unit Vector 2, 
                 Constant for cartesian planar equation with regards to surface normal k]
    '''
    ab = [b.x-a.x,b.y-a.y,b.z-a.z]  #vector from a --> b
    ac = [c.x-a.x,c.y-a.y,c.z-a.z]  #vector from a --> c
    
    n_ab = normalize(ab) #normaized ab vector
    n_ac = normalize(ac) #normaized ac vector
    
    #Cross product vector between ab and ac
    surface_normal = [n_ab[1]*n_ac[2] - n_ab[2]*n_ac[1], -1*(n_ab[0]*n_ac[2] - n_ab[2]*n_ac[0]), n_ab[0]*n_ac[1] - n_ab[1]*n_ac[0]]
    n_surface_normal = normalize(surface_normal)
    
    #Second planar axis along with ab, cross product of surface normal and ab
    in_plane_normal = [n_surface_normal[1]*n_ab[2] - n_surface_normal[2]*n_ab[1], -1*(n_surface_normal[0]*n_ab[2] - n_surface_normal[2]*n_ab[0]), n_surface_normal[0]*n_ab[1] - n_surface_normal[1]*n_ab[0]]
    n_in_plane_normal = normalize(in_plane_normal)
    
    #Planar constant --> u1 ∗ x + u2 ∗ y + u3 ∗ z + K = 0
    k = -1 * (a.x*surface_normal[0] + a.y*surface_normal[1] + a.z*surface_normal[2])
    
    #Declares planar origin = point a
    origin = AtomicPoint(a.x, a.y, a.z, "")
    
    return [origin, n_surface_normal, n_ab, n_in_plane_normal, k]
    

def planarizePoint(a, pv):
    '''
    Function takes a (x,y,z) point and a set of values produced by planarize, and projects that point
    onto the given plane, returning its (u,w) position. 

    Code from here: https://stackoverflow.com/questions/23472048/projecting-3d-points-to-2d-plane

    @param a Point with (x,y,z) values to be projected onto given plane
    @param pv Principle planar values calculated by planarize. [Plane Origin as AtomicPoint, Surface Normal Unit Vector, Axis Unit Vector 1, 
                                                                Othoganal Axis Unit Vector 2, Constant for cartesian planar equation with regards to surface normal k]
    @return AtomicPoint with (u,w) position as projected onto the given plane. s will be set to the out of plane seperation value (not sure if sign will be correct) (x,y,z) --> (u,v,s)
    '''
    origin = pv[0]
    surface_normal = pv[1]
    axis1 = pv[2]
    axis2 = pv[3]
    
    origin_vector = [a.x - origin.x, a.y - origin.y, a.z - origin.z]
    
    s = dot_product(surface_normal, origin_vector)
    u = dot_product(axis1, origin_vector)
    v = dot_product(axis2, origin_vector)

    solution = AtomicPoint(u, v, s, "", vdwr=a.VanDerWaalsRadius)
    return solution     


def nearestAtom(a, pointList, excluded, **kwargs):
    '''Function finds the nearest AtomicPoint to point a, ignoring any atoms
    given in the excluded list. VanDerWaals radius is acounted for. Default is
    distBetween calcualted in 3D space, for 2D add kwargs dimension='2D'
    
    @peram a AtomicPoint
    @peram pointList of AtomicPoints to search
    @peram excluded list of AtomicPoints to ignore
    **kwargs
        dimension : if = 2D uses distBetween2D as opposed to distBetween. 
                    if = 3D, uses distBetween as normal
    @return Nearest AtomicPoint to a ignoring excluded atoms
    '''
    
    #add **kwargs for distBetween and distBetween2D options?
    
    nearest_dist = 100000
    nearest_atom = AtomicPoint(0, 0, 0, "")

    if 'dimension' in kwargs:
        if kwargs['dimension'] == '2D':
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
        
        elif kwargs['dimension'] == '3D':
            #Checks all AtomicPoints in the given list
            for p in pointList:
                exclusion_check = False #Can't break outer for loop from inner for loop
                #Checks if the AtomicPoint is excluded and if so skips it
                for e in excluded:
                    if p == e: 
                        exclusion_check = True
                        break #Don't need to continue loop afer match is found
                if exclusion_check: continue
                
                if a.distBetween(p, True) < nearest_dist:
                    nearest_atom = p
                    nearest_dist = a.distBetween(p, True)
            return nearest_atom
        
        #Error input handling
        else:
            print("Error: nearestAtom function received incorrect kwargs argument")
            return nearest_atom
    
    #Checks all AtomicPoints in the given list
    for p in pointList:
        exclusion_check = False #Can't break outer for loop from inner for loop
        #Checks if the AtomicPoint is excluded and if so skips it
        for e in excluded:
            if p == e: 
                exclusion_check = True
                break #Don't need to continue loop afer match is found
        if exclusion_check: continue
        
        if a.distBetween(p, True) < nearest_dist:
            nearest_atom = p
            nearest_dist = a.distBetween(p, True)
            
    return nearest_atom


def atomCollide(a, pointList, **kwargs):
    '''
    Function checks if given AtomicPoint intersects with any other atoms in the
    pointList. Returns True if it does so, false otherwise.
    
    
    Currently only has distBetween2D, will eventually want to add kwargs condition
    for 3D distBetween
    
    
    @peram a AtomicPoint being checked
    @peram pointList List of AtomicPoints that p is being chekd against
    **kwargs
        count : Boolean, if true count the number intersecting atoms
        record : Boolean, if true record the intersecting atoms in a list. count
                 must be true for record to record all intersecting atoms
        threshold : Float, threshold under which the distance between the spheres
                has to be for it to count as a collision. -0.00001 is the default
    @return True if p itersects with any of the atoms in the pointList, returns 
    false otherwise
        If 'count' = True, then answer will be returned as [Boolean, int count of atoms intersecting]
        If 'record' = True, then list of atoms recorded will be appended to the end on an existing 
            answer list
    '''
    recordList = []
    recording = False
    if 'record' in kwargs:
        if kwargs['record'] == True:
            recording = True
    
    threshold = -0.00001
    if 'threshold' in kwargs:
        threshold = kwargs['threshold']
    
    if 'count' in kwargs:
        if kwargs['count']:
            iCount = 0
            intersecting = False
            
            for p in pointList:
                if a.distBetween2D(p, True) < threshold: 
                    intersecting = True
                    
                    
                    if recording:
                        recordList.append(p)
                    iCount = iCount + 1
            
            if intersecting: 
                if recording: return [True, iCount, recordList]
                else: return [True, iCount]
            else: 
                return [False, 0]
    
    for p in pointList:
        if a.distBetween2D(p, True) < threshold: return True
    return False