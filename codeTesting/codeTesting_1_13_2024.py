# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 21:33:06 2024

Testing of zSlice code for frames other than 111 and 121, getting collision error

@author: samda
"""

import zSliceLineAproximation as zsla
import pdbPoreReader
import LargestCircleFinder as lcf
import numpy as np
from AtomicPoint import AtomicPoint
import matplotlib.pyplot as plt
import sys
import math

#Possible create a file with some functions I commonly use to varify things graphing wise to expidite things
#Issue is with how tight circles are resolved with the end forcing to the 3 nearest even in bad
#situations

#Remove quadrant check for current circleCalc function until heuristic can be properly defined
#Get Grimm the atoms the were bad for his vs my code
#fix growCircle. Do step by step walkthrough of that function for the bad frame

file = '05 3x3x3 middle pore frame 112.pdb'
pdbList = pdbPoreReader.atomReaderPDB(file)


stepSize = 0.03 #Size in Angstroms of the distance between steps
bounds = lcf.boundsFinder(pdbList)

zMax = bounds[4]
zMin = bounds[5]

#List of z axis slices to run calculations on
zList = np.arange(zMin, zMax, stepSize)

numPathPoints = len(zList)
pathList = [None] * numPathPoints 
i = 0 #index for pathList and zList

#starting location for calc point
calcPoint = AtomicPoint((bounds[0]/2 + bounds[1]/2)*2, (bounds[2]/2 + bounds[3]/2)*2, zList[i], '')


j = 0


#calculates at each z value in the zList
for z in zList:
    
    calcPoint.z = z
    
    pointSet = zsla.sliceList(pdbList, z)
    if len(pointSet) < 3:
        print(f"Not enough atoms at z = {z} to form a circle. Only {len(pointSet)} atoms.")
        pathList[i] = calcPoint
    #calcuates the location of the cicle in the slice
    
    
    if j == 31:
        
        calcPoint.printCoords()
        print(z)
        
        #Find the atomic point nearest to the calc point, accounting for vdwr
            #Initializes the points
        ap1 = pointSet[0] #should initialze with different values
        ap2 = pointSet[0] #could do AtomicPoint(abs(Xmax)+abs(Xmin), abs(Ymax)+abs(Ymin), 0, "")
        ap3 = pointSet[0]
        current_radius = calcPoint.distBetween2D(ap1, True)
        maxCycles = 100
        
        #Adjustment for zSlice
        originZ = calcPoint.z
        
        
        #First circle snap calculation --> 0 point growth
        for p in pointSet:
            if(calcPoint.distBetween2D(p, True) < calcPoint.distBetween2D(ap1, True)):
                ap1 = p
                current_radius = calcPoint.distBetween2D(ap1, True)
        
        
        
        #Checks if calcPoint is inside the ap1 vdwr and cancels calculation if so
        #Should probably add something that instead find the nearest open area and
        #continues the claulcation from there
        if calcPoint.distBetween2D(ap1, False) < ap1.VanDerWaalsRadius: 
            print(-3)
            #return [-3, AtomicPoint(0, 0, 0, "")] 
            sys.exit()
        
        #List of step-sizes, could probably make soemthing more dynamic, but this should cover my basis for the moment
        scaleList = [1, 0.1, 0.01, 0.001, 0.0001, 0.00001]
        
        #Used to calculate angle from atomic point to calc point, vdwr not needed
        trig_x = ap1.x - calcPoint.x
        trig_y = ap1.y - calcPoint.y
        theta = math.acos(abs(trig_x) / math.hypot(trig_x, trig_y))
        
        #what I am currently using as logic for the while loop
        calculating = True
        stepIndex = 0
        calcCycles = 0
        
        
        k = 0
        
        #Second Circle Snap calculation --> 1 point growth
        #Could probably abstract this to its own function
        while calculating == True:
            
            #makes sure function doesn't get caught in infinite loops and just cancles the program if things go on for too long
            if calcCycles > maxCycles: 
                print(-4)
                #return [-4, AtomicPoint(0, 0, 0, "")] 
                sys.exit()
            calcCycles = calcCycles + 1
            
            #X and Y displacments for how far to move c
            disp_x = scaleList[stepIndex] * math.cos(theta)
            disp_y = scaleList[stepIndex] * math.sin(theta)        
            #Grow away from that nearest point
            
            new_x = calcPoint.x + (-lcf.sign(trig_x) * disp_x)
            new_y = calcPoint.y + (-lcf.sign(trig_y) * disp_y)
        
            #out of bounds checking
            if(new_x > bounds[0] or new_x < bounds[1] or new_y > bounds[2] or new_y < bounds[3]): 
                
                #print(f'1.1: new_x = {new_x}, new_y = {new_y}')
                print(-1)
                #return [-1, AtomicPoint(0, 0, 0, '')] #error return
                sys.exit()
        
            new_point = AtomicPoint(new_x, new_y, 0, "")
            new_radius = new_point.distBetween2D(ap1, True)
        
            eligible_points_counter = 0
            
            #Checks every point in the  given atomic array, and sees if any other the 
            #first contact point are at the edge or inside the bounds of the new radius
            for p in pointSet:
                if(new_point.distBetween2D(p, True) <= new_radius and ap1 != p):
                    eligible_points_counter = eligible_points_counter + 1
                    ap2 = p
                    #Track points
            
            #Only proceeds if only 1 new point is encountered. Otherwise decriments step 
            #size and checks again.
            #If only 1 new point is found (ignoring the contact point), snaps the circle
            #to includee the new contact point, the proceeds to next step of calculation
            if(eligible_points_counter > 1):
                #reset, use smaller incriment
                if stepIndex < 5:
                    stepIndex = stepIndex + 1
                    continue
                else:
                    #print("Incrimenting Error.")
                    print(-2)
                    #return [-2, AtomicPoint(0, 0, 0, '')] #error return
                    sys.exit()
            elif(eligible_points_counter == 1):
                #perform circle snap calculation for the two points
                calculating = False
                
                #calculates location of new point, function above
                calcPoint = lcf.twoAtomCircle(ap1, ap2, calcPoint)
                current_radius = ap1.distBetween2D(calcPoint, True)
                
            else:
                calcPoint = new_point
                current_radius = new_radius #Technicaly don't need to track this here
        
         
        
        #Third Circle Snap calculation --> 2 (or 3) point growth
        
        calculating = True
        stepIndex = 0
        calcCycles = 0
        
        while calculating == True:
            #In this case, ap1 and ap2 are the circles the algorithm is 'growing' from
            #The vector direction of growth is perpendicular to the line created from
            #the points of the calculating circle touching the atomic points
            
            bounds = lcf.boundsFinder(pointSet)

            fig, ax = plt.subplots()
            for p in pointSet:
                circle = plt.Circle((p.x, p.y), p.VanDerWaalsRadius, color='r')
                ax.add_patch(circle)
            
            circle2 = plt.Circle((calcPoint.x, calcPoint.y), current_radius, color='blue')
            ax.add_patch(circle2)
            
            circleAp1 = plt.Circle((ap1.x, ap1.y), ap1.VanDerWaalsRadius, color='k')
            ax.add_patch(circleAp1)
            
            circleAp2 = plt.Circle((ap2.x, ap2.y), ap2.VanDerWaalsRadius, color='purple')
            ax.add_patch(circleAp2)
            
            circleAp3 = plt.Circle((ap3.x, ap3.y), ap3.VanDerWaalsRadius, color='orange')
            ax.add_patch(circleAp3)
            
            ax.set_xlim(bounds[1], bounds[0])
            ax.set_ylim(bounds[3], bounds[2])

            plt.show()
            
            print("Code check")
            
            #makes sure function doesn't get caught in infinite loops and just cancles the program if things go on for too long
            if calcCycles > maxCycles: 
            
                '''
                Circle Surround seems to be the failing test, or at least certain circle configurations can't seem to grow properly'
                
                if eligible_points_counter == 1:
                    return
                '''
                print(-5)
                #return [-5, AtomicPoint(0, 0, 0, "")] 
                sys.exit()
        
            
            calcCycles = calcCycles + 1
            
            '''
            #Might be worth abstracting vector calculation to a funtion or using the methods
            #python provides to simplify, should test this section.
            '''
            
            
            
            #Vector from ap1 toward calcPoint
            vp1cX = calcPoint.x - ap1.x
            vp1cY = calcPoint.y - ap1.y
            #Unit Vectors
            uvp1cX = vp1cX / math.sqrt(vp1cX**2 + vp1cY**2)
            uvp1cY = vp1cY / math.sqrt(vp1cX**2 + vp1cY**2)
            
            #Vector from ap2 toward calcPoint
            vp2cX = calcPoint.x - ap2.x
            vp2cY = calcPoint.y - ap2.y
            #Unit Vectors
            uvp2cX = vp2cX / math.sqrt(vp2cX**2 + vp2cY**2)
            uvp2cY = vp2cY / math.sqrt(vp2cX**2 + vp2cY**2)
            
            #Add vectors to get average vector direction of travel
            #Don't need to average the vectors if I normalize afterwards anyways, right?
            vX = (uvp1cX + uvp2cX) / 2
            vY = (uvp1cY + uvp2cY) / 2
            
            #Convert to unit vectors again
            uvX = vX / math.sqrt(vX**2 + vY**2)
            uvY = vY / math.sqrt(vX**2 + vY**2)
            
            #New displacement in x and y
            disp_x = scaleList[stepIndex] * uvX
            disp_y = scaleList[stepIndex] * uvY
            
            #Get new location of calculation
            new_x = calcPoint.x + disp_x
            new_y = calcPoint.y + disp_y
            
            #out of bounds checking
            if(new_x > bounds[0] or new_x < bounds[1] or new_y > bounds[2] or new_y < bounds[3]): 
                print(-1)
                #return [-1, AtomicPoint(0, 0, 0, '')] #error return
                sys.exit()
            
            new_point = AtomicPoint(new_x, new_y, 0, "")
            new_radius = new_point.distBetween2D(ap1, True)
        
            eligible_points_counter = 0
            
            #Checks every point in the given atomic array, and sees if any other 
            #point at or inside the bounds of the new radius other than ap1 and ap2
            for p in pointSet:
                if(new_point.distBetween2D(p, True) <= new_radius and ap1 != p and ap2 != p):
                    eligible_points_counter += 1                
                    ap3 = p
            
            #Only proceeds if only 1 new point is encountered. Otherwise decriments step 
            #size and checks again.
            #If only 1 new point is found (ignoring the contact point), snaps the circle
            #to includee the new contact point
            if(eligible_points_counter > 1):
                #reset, use smaller incriment
                if stepIndex < 5:
                    stepIndex = stepIndex + 1
                    continue
                else:
                    #if it gets to the state where there are still too many eligible points at
                    #0.001 step incriments, just find the nearest point and calculate the circle from there.
                    for p in pointSet:
                        if(new_point.distBetween2D(p, True) <= new_point.distBetween2D(ap3, True) and ap1 != p and ap2 != p):
                            eligible_points_counter += 1                
                            ap3 = p
                    
                    calcPoint = lcf.threeCircleCircumcenter(ap1, ap2, ap3)
                    current_radius = ap1.distBetween2D(calcPoint, True)
                    calcPoint.z = originZ 
                    
                    
                    bounds = lcf.boundsFinder(pointSet)

                    fig, ax = plt.subplots()
                    for p in pointSet:
                        circle = plt.Circle((p.x, p.y), p.VanDerWaalsRadius, color='r')
                        ax.add_patch(circle)
                    
                    circle2 = plt.Circle((calcPoint.x, calcPoint.y), current_radius, color='blue')
                    ax.add_patch(circle2)
                    
                    circleAp1 = plt.Circle((ap1.x, ap1.y), ap1.VanDerWaalsRadius, color='k')
                    ax.add_patch(circleAp1)
                    
                    circleAp3 = plt.Circle((ap3.x, ap3.y), ap3.VanDerWaalsRadius, color='orange')
                    ax.add_patch(circleAp3)
                    
                    circleAp2 = plt.Circle((ap2.x, ap2.y), ap2.VanDerWaalsRadius, color='purple')
                    ax.add_patch(circleAp2)
                    
                    ax.set_xlim(bounds[1], bounds[0])
                    ax.set_ylim(bounds[3], bounds[2])

                    plt.show()
                    
                    print("Code check")
                    
                    
                    print("Valid end!")
                    #return [current_radius, calcPoint, [ap1, ap2, ap3]]
                    sys.exit()
                
            elif(eligible_points_counter == 1 and ap1 != ap2 and ap1 != ap3 and ap2 != ap3):
                #perform circle snap calculation for the three points
                calculating = False
                
                #Now need to check if the points bound the calcPoint on all side/limit growth
                #If they do thats the final circle. If they don't and all three point's fall 
                #within a <180 degree arc, find the outermost circles of the three circles 
                #and repeat this loop with those points as the new ap1 and ap2
                surround = lcf.circleSurrounded(ap1, ap2, ap3, calcPoint)
                
                if surround[0]:
                    #Current return is in the form of a list with radius and calcPount, but could just 
                    #make vdwr the current radius if I wanted to 
                    #calculates location of new point, function above
                    calcPoint = lcf.threeCircleCircumcenter(ap1, ap2, ap3)
                    current_radius = ap1.distBetween2D(calcPoint, True)
                    
                    calcPoint.z = originZ 
                    
                    #return [current_radius, calcPoint, [ap1, ap2, ap3]]
                    print("Valid end 2!")
                    sys.exit()
                #if circle is not surounded, continue to grow the circle, the 2 outermost of the 
                #3 defining atomic points servining as the new ap1 and ap2
                else:                
                    ap1 = surround[1]
                    ap2 = surround[2]
                    
                    calculating = True
                    stepIndex = 0
                
            else:
                calcPoint = new_point
                current_radius = new_radius #Technicaly don't need to track this here
                            
        sys.exit() #end catch
        
    else:
        newPoint = lcf.growCircle(pointSet, calcPoint, bounds)
        newPoint[1].VanDerWaalsRadius = newPoint[0]
    #checks if circle is valid
    
    
    if j == 33:
        '''
        bounds = lcf.boundsFinder(pointSet)

        fig, ax = plt.subplots()
        for p in pointSet:
            circle = plt.Circle((p.x, p.y), p.VanDerWaalsRadius, color='r')
            ax.add_patch(circle)
        
        circle2 = plt.Circle((calcPoint.x, calcPoint.y), 1, color='blue')
        ax.add_patch(circle2)
        
        ax.set_xlim(bounds[1], bounds[0])
        ax.set_ylim(bounds[3], bounds[2])

        plt.show()
        '''
        break
    
    
    elif newPoint[0] > 0:
        pathList[i] = newPoint[1]            
        #makes the new calcPoint the center of the calculated circle 
        calcPoint = newPoint[1]
        calcPoint.VanDerWaalsRadius = 0
    else:
        #Just shifts the calcPoint up in the z axis, but really if this happens you have an error
        #print(f'Collision error at z = {z}, radius = {newPoint[0]}, newPoint ({newPoint[1].x}, {newPoint[1].y}, {newPoint[1].z}), calcPoint ({calcPoint.x},{calcPoint.y},{calcPoint.z})')
        
        
        bounds = lcf.boundsFinder(pointSet)

        fig, ax = plt.subplots()
        for p in pointSet:
            circle = plt.Circle((p.x, p.y), p.VanDerWaalsRadius, color='r')
            ax.add_patch(circle)
        
        circle2 = plt.Circle((calcPoint.x, calcPoint.y), 1, color='blue')
        ax.add_patch(circle2)
        
        ax.set_xlim(bounds[1], bounds[0])
        ax.set_ylim(bounds[3], bounds[2])

        plt.show()
        
        print(j)
        
        break
        
    
        
        pathList[i] = calcPoint
    
    i += 1
    
    j += 1