# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 16:48:15 2024

Walking through the growCircle function to see what the issues are

@author: samda
"""
import LargestCircleFinder as lcf
import lcfProcessing as proc
from AtomicPoint import AtomicPoint
import pdbPoreReader
import zSliceLineAproximation as zsla
import testHelpers as th
import matplotlib.pyplot as plt
import time
import growCircle as gc
import PoreCalculatorV0Functions as pcf
import numpy as np


#Check if previous point is in bounds of followup circle, if not, use approxCircle, 
#otherwise use growCricle.
#Check by using a.distBetween2D(b, False), if dist < a.vdwr, b inside circle a


#Environemnt set up
file = '05 3x3x3 middle pore frame 111.pdb'
pdbList = pdbPoreReader.atomReaderPDB(file)

#z = 30.586000000000034
#inputPoint = AtomicPoint(-0.4385308775314408,-0.4077724205622452,30.586000000000034, "")

#pointList = zsla.sliceList(pdbList, z)
#bounds = lcf.boundsFinder(pointList)

n=30

start = time.time()

#code testing
#result = gc.approxCircle(pointSet, inputPoint, n)


stepSize = 0.1 #Size in Angstroms of the distance between steps
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
    
    
    print(f'{z}')
    
    
    calcPoint.z = z
    
    pointSet = zsla.sliceList(pdbList, z)
    if len(pointSet) < 3:
        print(f"Not enough atoms at z = {z} to form a circle. Only {len(pointSet)} atoms.")
        pathList[i] = calcPoint
    #calcuates the location of the cicle in the slice
    limitedList = pcf.listLimiter2D(calcPoint, pointSet, n, setZ=calcPoint.z)
    newPoint = gc.approxCircle(limitedList, calcPoint, n)
    
    pathList[i] = newPoint[0]            
    #makes the new calcPoint the center of the calculated circle 
    calcPoint = newPoint[0]
    calcPoint.VanDerWaalsRadius = 0
        
    i += 1
    #print(i)
    #getting math domain error on cycle 9 (8 last number before error)

end = time.time()
print(f'Function took {end - start} seconds to run.')


#limitedAtoms = pcf.listLimiter2D(inputPoint, pointSet, n)
#result[0].printCoords()

'''
fig, ax = plt.subplots()

ax.set_xlim(bounds[1], bounds[0])
ax.set_ylim(bounds[3], bounds[2])

for p in pointSet:
    th.plotAtom(p, ax)
    
for q in limitedAtoms:
    th.plotAtom(q, ax, color='orange')
    
for v in result[1]:
    th.plotAtom(v, ax, color='purple')
    
th.plotAtom(result[0], ax, color='blue')

plt.show()
'''