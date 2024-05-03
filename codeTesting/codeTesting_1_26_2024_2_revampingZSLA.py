# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 14:45:57 2024

zsla test

@author: samda
"""
from AtomicPoint import AtomicPoint
import pdbPoreReader
import zSliceLineAproximation as zsla
import LargestCircleFinder as lcf
import numpy as np
import growCircle as gc

#Environemnt set up
file = '05 3x3x3 middle pore frame 111.pdb'
pdbList = pdbPoreReader.atomReaderPDB(file)

#zsla.pathGen(pdbList)


#PathGen Function
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

#calculates at each z value in the zList
for z in zList:
    
    calcPoint.z = z
    
    pointSet = zsla.sliceList(pdbList, z)
    if len(pointSet) < 3:
        print(f"Not enough atoms at z = {z} to form a circle. Only {len(pointSet)} atoms.")
        pathList[i] = calcPoint
    #calcuates the location of the cicle in the slice
    
    #Getting div by 0 errors compaired to lcf.circleGrow
    newPoint = gc.growCircle(pointSet, calcPoint, bounds)
    
    
    newPoint[1].VanDerWaalsRadius = newPoint[0]
    #checks if circle is valid
    if newPoint[0] > 0:
        pathList[i] = newPoint[1]            
        #makes the new calcPoint the center of the calculated circle 
        calcPoint = newPoint[1]
        calcPoint.VanDerWaalsRadius = 0 #Setting results vdwr to 0 as well
    else:
        #Just shifts the calcPoint up in the z axis, but really if this happens you have an error
        print(f'Collision error at z = {z}, radius = {newPoint[0]}, newPoint ({newPoint[1].x}, {newPoint[1].y}, {newPoint[1].z}), calcPoint ({calcPoint.x},{calcPoint.y},{calcPoint.z})')
        pathList[i] = calcPoint
    
    i += 1


#Validation
for point in pathList:
    point.printCoords()