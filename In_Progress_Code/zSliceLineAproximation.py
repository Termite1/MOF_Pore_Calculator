# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 12:05:05 2023

This set of functions approximates the central line traveling through a pore
by find the largest circle that can fit in the pore at each slize of the z axis
at regular intervals (0.1 Angstroms?). The line is approximated by defining a 
path through the center of each calculated circle.

@author: samda
"""

from AtomicPoint import AtomicPoint
import math
import LargestCircleFinder as lcf
import numpy as np
import PoreCalculatorV0Functions as pcf
import lcfProcessing as proc


'''
Calculates the radius of the circle in the z plane created by the z plane at a 
certain value cutting thorugh an atom.

@param p Atomic Point
@param z (x,y) plane location along the z axis
@return AtomicPoint with VanDerWaalsRadius calculated for the z slice
'''
def zSlice(p, z):
    #Checks if z falls in the radius of the point, accounts for negative z
    if (abs(p.z) + p.VanDerWaalsRadius < abs(z) or abs(p.z) - p.VanDerWaalsRadius > abs(z)):
        newPoint = AtomicPoint(p.x, p.y, p.z, p.element)
        newPoint.VanDerWaalsRadius = 0
    else:
        a = abs(p.z - z)
        c = p.VanDerWaalsRadius
        b = math.sqrt(c**2 - a**2)
    
        newPoint = AtomicPoint(p.x, p.y, p.z, p.element)
        newPoint.VanDerWaalsRadius = b
    
    return newPoint


'''
Given a list of points and a z value returns a list of points and their radius
at that slice of the (x,y) plane.

@param pointList List of atomic points
@param z level being calculated for
@return List of points that fall in that (x,y) plane with their adjusted radius
'''
def sliceList(pointList, z):
    i = 0 #number of points for the adjusted list
    
    #Counts the number of eligible points in the list
    for p in pointList:
        if (abs(z) > abs(p.z) - p.VanDerWaalsRadius and abs(z) < abs(p.z) + p.VanDerWaalsRadius):
            i += 1

    #Creates new list
    newList = [None] * i
    i = 0
    
    #Generated the adjusted list
    for p in pointList:
        if (abs(z) > abs(p.z) - p.VanDerWaalsRadius and abs(z) < abs(p.z) + p.VanDerWaalsRadius):
            newPoint = zSlice(p, z)
            newList[i] = newPoint
            i += 1
    
    return newList


def approxCircle(pointList, calcPoint, n):
    """
    Function finds the largest circle that can find in the location defined by
    finding the n nearest atoms to the calcPoint and then calculating the 
    largest surrounded circle that can fit in that area. Requires all the 
    AtomicPoints in the pointList to have the same z to work in a 2D slice
    of the xy axes. z of calcPoint needs to be the same as the z of the slice.

    Parameters
    ----------
    pointList : List of AtomicPoints
        List of AtomicPoints describing the position of the atoms.
    calcPoint : AtomicPoint
        AtomicPoint that the calculation grows from.
    n: int
        Number of atoms nearest to the calcPoint to calculate from.

    Returns
    -------
    List with [0] the AtomicPoint details of the calculated circle, and [1] a 
    list of the 3 AtomicPoints that define the circle about the calcPoint
    """
    #Get n nearest atoms to the calcPoint
    limitedList = pcf.listLimiter2D(calcPoint, pointList, n, setZ=calcPoint.z)
    
    largestCircle = AtomicPoint(0,0,calcPoint.z,"")
    surroundingAtoms = [None, None, None]
    
    #Checks every combination of the given atoms
    for a in range(0, len(limitedList)):
        for b in range(1, len(limitedList)):
            if b == a: continue
            for c in range(2, len(limitedList)):
                if c == a or c == b: continue
                
                newCircle = lcf.threeAtomCircumcenter(limitedList[a], limitedList[b], limitedList[c])
                #Makes sure the circle is surrounded
                if lcf.simpleCircleSurrounded(limitedList[a], limitedList[b], limitedList[c], newCircle) == False:
                    continue
                
                #make sure circle created doesn't collide with any atoms
                #for atom in limitedList:
                #    if newCircle.distBetween2D(atom, True) < 0: continue
                if proc.atomCollide(newCircle, limitedList): continue
                
                #keeps the circle if its radius is larger than the current
                #largest circle
                if newCircle.VanDerWaalsRadius > largestCircle.VanDerWaalsRadius:
                    largestCircle = newCircle
                    surroundingAtoms[0] = limitedList[a]
                    surroundingAtoms[1] = limitedList[b]
                    surroundingAtoms[2] = limitedList[c]
    
    return [largestCircle, surroundingAtoms]


def pathGen(pointList, **kwargs):
    '''
    Given a list of points defining a pore generates an approximation of the line
    the travels through the center of the pore.

    This function could run into errors if there are U-shaped protrusions cutting one side to a dead end
    and the other to the pore exit or a smaller pore pathway. Could try going from the other side at 
    that point.

    @param pointList List of Atomic Points
    
    kwargs 
        smallest_r : if present, records and returns smallest radius encountered during path
    
    @return List of AtomicPoints approximating line traveling through the center of the pore.
    '''
    
    stepSize = 0.3 #Size in Angstroms of the distance between steps
    n = 40 #Number of nearest atoms to grab when performing approxCircle on slices with anomoulus answers
    #Should take ~1 sec to calculate
    bounds = lcf.boundsFinder(pointList)
    
    zMax = bounds[4]
    zMin = bounds[5]
    
    #List of z axis slices to run calculations on
    zList = np.arange(zMin, zMax, stepSize)
    
    numPathPoints = len(zList)
    pathList = [None] * numPathPoints 
    i = 0 #index for pathList and zList
    
    smallest_r = 10000 #records smallest radius encountered during path
    
    #starting location for calc point
    calcPoint = AtomicPoint((bounds[0]/2 + bounds[1]/2)*2, (bounds[2]/2 + bounds[3]/2)*2, zList[i], '')
    
    #calculates at each z value in the zList
    for z in zList:
        
        calcPoint.z = z
        
        pointSet = sliceList(pointList, z)
        if len(pointSet) < 3:
            #print(f"Not enough atoms at z = {z} to form a circle. Only {len(pointSet)} atoms.")
            pathList[i] = calcPoint
        #calcuates the location of the cicle in the slice
        newPoint = lcf.growCircle(pointSet, calcPoint, bounds)
        newPoint[1].VanDerWaalsRadius = newPoint[0]
        #checks if circle is valid
        if newPoint[0] > 0:
            #Checks if start of calculation is within the bounds of the solution,
            #hopefully will stop anomolous answers from occuring
            check = True
            
            if newPoint[1].distBetween2D(calcPoint, True) > 0:
                check = False
                '''
                #Use approxCircle when anomoulus answer detected
                result = approxCircle(pointSet, calcPoint, n)
                #Makes new calcPoint
                pathList[i] = result[0]
                calcPoint = AtomicPoint(result[0].x, result[0].y, 0, "")
                '''
            #-0.1 is a large collision dectection threshold, but for this approximation
            #only the most egregious errors need to be caught and corrected, i.e. full out intersections
            #as to minor collisions due to programatic issues with the growCircle function
            elif proc.atomCollide(newPoint[1], pointSet, threshold=-0.1):
                check = False
            #Intended to make sure circle is acutally touching the atoms the bound it, not necessary at this time
            #elif lcf.threeAtomRadiusCheck(newPoint[2][0], newPoint[2][1], newPoint[2][2], newPoint[1], circle=True) == False:
                #check = False
            
            '''
            else:
                pathList[i] = newPoint[1]            
                #makes the new calcPoint the center of the calculated circle 
                calcPoint = AtomicPoint(newPoint[1].x, newPoint[1].y, 0, "")
                calcPoint.VanDerWaalsRadius = 0
            '''
            
            if check == False:    
                #Use approxCircle when anomoulus answer detected
                result = approxCircle(pointSet, AtomicPoint(0,0,z,""), n)
                #Makes new calcPoint
                pathList[i] = result[0]
                layerAns = result[0]
                
                if 'smallest_r' in kwargs:
                    if result[0].VanDerWaalsRadius < smallest_r: smallest_r = result[0].VanDerWaalsRadius
                
                calcPoint = AtomicPoint(result[0].x, result[0].y, 0, "")
            else:
                pathList[i] = newPoint[1]            
                
                if 'smallest_r' in kwargs:
                    if newPoint[0] < smallest_r: smallest_r = newPoint[0]
                
                #makes the new calcPoint the center of the calculated circle 
                calcPoint = AtomicPoint(newPoint[1].x, newPoint[1].y, 0, "")
            
        else:
            '''
            #Just shifts the calcPoint up in the z axis, but really if this happens you have an error
            print(f'Collision error at z = {z}, radius = {newPoint[0]}, newPoint ({newPoint[1].x}, {newPoint[1].y}, {newPoint[1].z}), calcPoint ({calcPoint.x},{calcPoint.y},{calcPoint.z})')
            pathList[i] = calcPoint
            '''
            #Use approxCircle when anomoulus answer detected
            result = approxCircle(pointSet, AtomicPoint(0,0,z,""), n)
            #Makes new calcPoint
            
            pathList[i] = result[0]
            layerAns = result[0]
            
            if 'smallest_r' in kwargs:
                if result[0].VanDerWaalsRadius < smallest_r: smallest_r = result[0].VanDerWaalsRadius
            
            calcPoint = AtomicPoint(result[0].x, result[0].y, 0, "")
        
        i += 1
    
    if 'smallest_r' in kwargs: return [pathList, smallest_r]
    return pathList