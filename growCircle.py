# -*- coding: utf-8 -*-
"""
Created on Sun Jan 21 16:07:56 2024

@author: samda
"""

import LargestCircleFinder as lcf
import lcfProcessing as proc
from AtomicPoint import AtomicPoint
import PoreCalculatorV0Functions as pcf


def approxCircle(pointList, calcPoint, n):
    """
    Function linds the largest circle that can fin in the location defined by
    finding the n nearest atoms to the calcPoint and then calculating the 
    largest surrounded circle that can fit in that area. Requires all the 
    AtomicPoints in the pointList to have the same z to work in a 2D slice
    of the xy axes. z of calcPoint needs to be the same as the z of the slice.

    Seems to be working properly now?

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
    for a in range(0, len(limitedList) - 1):
        for b in range(1, len(limitedList) - 1):
            if b == a: continue
            for c in range(2, len(limitedList) - 1):
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
                

def growCircle(pointList, calcPoint, bounds):
    """
    Given a list of AtomicPoints, finds the largest circle that can fit in the 
    locations defined, as calculated from the calcPoint. Intended to be used to 
    find the largest circle that can fit through the pore of a Metal Organic 
    Framework, with atoms described by an input pdb file.  

    @param pointList List of AtomicPoints describing the position of the atoms
    @param calcPoint AtomicPoint that the calculation grows from
    @param bounds List of floats holdng the x,y bounds of the pointArray [Xmax, Xmin, Ymax, Ymin]
    @return List containing [0] the float radius of the circle grown, [1] the AtomicPoint details of the calcPoint, 
            and [2] a list of the 3 AtomicPoints that define the circle about the calcPoint
    """
    #initializes the displacement vector [x,y]
    displacement_vector = [0, 0]
    maxCycles = 100 #Maximum number of calculation cicles, beyond this something is probably wrong
    
    #Find nearest atom to calcPoint
    ap1 = proc.nearestAtom(calcPoint, pointList, [], dimension='2D')
    ap2 = None
    ap3 = None
    #Set circle radius to distance between point and nearest atom
    calcPoint.VanDerWaalsRadius = calcPoint.distBetween2D(ap1, True)
    
    #Checks if calcPoint is inside the ap1 vdwr and cancels calculation if so
    if calcPoint.VanDerWaalsRadius < 0:
        return [-3, AtomicPoint(0, 0, 0, "")] 
    
    #List of step-sizes, could probably make something more dynamic, but this should cover my basis for the moment
    scaleList = [1, 0.1, 0.01, 0.001, 0.0001, 0.00001]
    
    #Calculate displacement_vector (inverse of directional vector from calcPoint to ap1)
    displacement_vector[0] = calcPoint.x - ap1.x
    displacement_vector[1] = calcPoint.y - ap1.y  
    
    #normalize displacement_vector (magnitude set to 1)
    displacement_vector = proc.normalize(displacement_vector)
    
    #what I am currently using as logic for the while loop
    calculating = True
    stepIndex = 0
    calcCycles = 0
    
    #Second Circle Snap calculation --> 1 point growth
    #Could probably abstract this to its own function
    while calculating == True:
        
        #makes sure function doesn't get caught in infinite loops and just cancles the program if things go on for too long
        if calcCycles > maxCycles: return [-4, AtomicPoint(0, 0, 0, "")] 
        calcCycles = calcCycles + 1
        
        #Find x and y displacement
        x_disp = displacement_vector[0] * scaleList[stepIndex]
        y_disp = displacement_vector[1] * scaleList[stepIndex]
        
        #Calculate new position with respect to the calcPoint
        new_x = calcPoint.x + x_disp
        new_y = calcPoint.y + y_disp

        #out of bounds checking
        if(new_x > bounds[0] or new_x < bounds[1] or new_y > bounds[2] or new_y < bounds[3]): 
            return [-1, AtomicPoint(0, 0, 0, '')] #error return
        
        #Make new location into an AtomicPoint
        new_point = AtomicPoint(new_x, new_y, 0, "")
        new_point.VanDerWaalsRadius = new_point.distBetween2D(ap1, True)
        
        #Checks if new point collides with any atoms and if so how many
        collideResult = proc.atomCollide(new_point, pointList, count=True, record=True)
        
        if collideResult[0] == False: #No collisions, update calcPoint as new_point
            calcPoint = new_point
        elif collideResult[0] == True and collideResult[1] == 1: #Only collided with 1 atom, move to next step
            calculating = False
            ap2 = collideResult[2][0]
            
            #calculates location of new point
            #might want to walk through edge cases of twoAtomCircle
            calcPoint = lcf.twoAtomCircle(ap1, ap2, calcPoint)
            calcPoint.VanDerWaalsRadius = ap1.distBetween2D(calcPoint, True)
        else: #More than 1 atom collided with, reduces step size then repeates
            #reset, use smaller step size incriment
            if stepIndex < 5:
                stepIndex = stepIndex + 1
                continue
            else:
                #Should maybe handle this diffeently, but unlikely to come up due
                #to small step sizes avalible
                return [-2, AtomicPoint(0, 0, 0, '')] #error return

            
    #continue to three atom circle from here
    #Reset values
    calculating = True
    stepIndex = 0
    calcCycles = 0
    
    while calculating == True:
        #In this case, ap1 and ap2 are the circles the algorithm is 'growing' from
        #The vector direction of growth is perpendicular to the line created from
        #the points of the calculating circle touching the atomic points
        
        #makes sure function doesn't get caught in infinite loops and just cancels the program if things go on for too long
        #May need to add a better method of error handling though
        if calcCycles > maxCycles: return [-5, AtomicPoint(0, 0, 0, "")] 
        calcCycles = calcCycles + 1
        
        #Calculate directional vectors from ap1 --> calcPoint and ap2 --> calcPoint
        v1 = [calcPoint.x - ap1.x, calcPoint.y - ap1.y]
        v2 = [calcPoint.x - ap2.x, calcPoint.y - ap2.y]
        
        #normalize the directional vectors
        v1 = proc.normalize(v1)
        v2 = proc.normalize(v2)
        
        #average the directional vectors
        displacement_vector = [0.5 * (v1[0] + v2[0]), 0.5 * (v1[0] + v2[1])]
        displacement_vector = proc.normalize(displacement_vector)
        
        #Find x and y displacement
        x_disp = displacement_vector[0] * scaleList[stepIndex]
        y_disp = displacement_vector[1] * scaleList[stepIndex]
        
        #Calculate new position with respect to the calcPoint
        new_x = calcPoint.x + x_disp
        new_y = calcPoint.y + y_disp
    
        #out of bounds checking
        if(new_x > bounds[0] or new_x < bounds[1] or new_y > bounds[2] or new_y < bounds[3]): 
            return [-1, AtomicPoint(0, 0, 0, '')] #error return
        
        #Make new location into an AtomicPoint
        new_point = AtomicPoint(new_x, new_y, 0, "")
        new_point.VanDerWaalsRadius = new_point.distBetween2D(ap1, True)
        
        #Checks if new point collides with any atoms and if so how many
        collideResult = proc.atomCollide(new_point, pointList, count=True, record=True)
        
        #No collisions, update calcPoint as new_point
        if collideResult[0] == False: calcPoint = new_point
        elif collideResult[0] == True and collideResult[1] > 1:
            #reset, use smaller step size incriment
            if stepIndex < 5:
                stepIndex = stepIndex + 1
                continue
            else:
                #Do things here
                #Should maybe handle this diffeently, but unlikely to come up due
                #to small step sizes avalible
                #First find the nearest 3 points. If they sourround the atom
                #three atom circle and return that as the answer
                return [-2, AtomicPoint(0, 0, 0, '')] #error return
        elif collideResult[0] == True and collideResult[1] == 1:
            #A possible answer has been found
            calculating = False
            
            #Set third point the the most recent one collided with
            ap3 = collideResult[2][0]
            
            #Generate the new circle
            new_point = lcf.threeAtomCircumcenter(ap1, ap2, ap3)
            
            #Check if the new circle is sourrounded by the three points
            surround = lcf.circleSurrounded(ap1, ap2, ap3, new_point)
            
            #if the circle is surrounded, return it and any related information
            if surround[0]:
                #Technically don't have to return radius, that is imbeded in the vdwr
                return [new_point.VanDerWaalsRadius, new_point, [ap1, ap2, ap3]]
            
            #if the circle isn't surrounded, grab the outermost of the surroudning points,
            #and update the location of the calcPoint
            else:
                ap1 = surround[1]
                ap2 = surround[2]
                calcPoint = new_point
                
                calculating = True
                stepIndex = 0