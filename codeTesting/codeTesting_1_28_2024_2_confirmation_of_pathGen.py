# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 12:31:46 2024

pathGen testing after approxCircle fix

@author: samda
"""
import zSliceLineAproximation as zsla
import pdbPoreReader as pdb
import LargestCircleFinder as lcf
import lcfProcessing as proc
import numpy as np
from AtomicPoint import AtomicPoint
import matplotlib.pyplot as plt
import sys
import math
import testHelpers as th
import growCircle as gc
import time
#for the animation
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter


file = 'MOF-5 middle pore frames 111-121.pdb'
#pdbList = pdbPoreReader.atomReaderPDB(file)

frames = pdb.multiframePDBReader(file)

#print(len(frames)) - 11 - [0-10]

pdbList = frames[9]

#IT WORKS!!! Do time test tomorrow and try to get the updates loaded into pathGen
#generate animation as well?
#Once pathGen works, move onto consecutive pore calc


start = time.time()


stepSize = 0.1 #Size in Angstroms of the distance between steps
n = 40 #Number of nearest atoms to grab when performing approxCircle on slices with anomoulus answers
#Should take ~1 sec to calculate
bounds = lcf.boundsFinder(pdbList)

zMax = bounds[4]
zMin = bounds[5]


#List of z axis slices to run calculations on
zList = np.arange(zMin, zMax, stepSize)

numPathPoints = len(zList)
pathList = [None] * numPathPoints 


#Just for the animation creation
zPointsList=[]


i = 0 #index for pathList and zList

#starting location for calc point
calcPoint = AtomicPoint((bounds[0]/2 + bounds[1]/2)*2, (bounds[2]/2 + bounds[3]/2)*2, zList[i], '')

#calculates at each z value in the zList
for z in zList:
    
    calcPoint.z = z
    
    pointSet = zsla.sliceList(pdbList, z)
    
    #Just for the animation creation
    zPointsList.append(pointSet)
    
    if len(pointSet) < 3:
        #print(f"Not enough atoms at z = {z} to form a circle. Only {len(pointSet)} atoms.")
        pathList[i] = calcPoint
    #calcuates the location of the cicle in the slice
    newPoint = lcf.growCircle(pointSet, calcPoint, bounds)
    newPoint[1].VanDerWaalsRadius = newPoint[0]
    
    layerAns = newPoint[1]
    #checks if circle is valid
    if newPoint[0] > 0:
        #Checks if start of calculation is within the bounds of the solution,
        #hopefully will stop anomolous answers from occuring
        check = True
        
        if newPoint[1].distBetween2D(calcPoint, True) > 0:
            check = False
        #-0.1 is a large collision dectection threshold, but for this approximation
        #only the most egregious errors need to be caught and corrected, i.e. full out intersections
        #as to minor collisions due to programatic issues with the growCircle function
        elif proc.atomCollide(newPoint[1], pointSet, threshold=-0.1):
            check = False
        #Intended to make sure circle is acutally touching the atoms the bound it, not necessary at this time
        #elif lcf.threeAtomRadiusCheck(newPoint[2][0], newPoint[2][1], newPoint[2][2], newPoint[1], circle=True) == False:
            #check = False
        
        if check == False:    
            #Use approxCircle when anomoulus answer detected
            result = zsla.approxCircle(pointSet, AtomicPoint(0,0,z,""), n)
            #Makes new calcPoint
            pathList[i] = result[0]
            layerAns = result[0]
            
            calcPoint = AtomicPoint(result[0].x, result[0].y, 0, "")
        else:
            pathList[i] = newPoint[1]            
            #makes the new calcPoint the center of the calculated circle 
            calcPoint = AtomicPoint(newPoint[1].x, newPoint[1].y, 0, "")
    else:
        #Use approxCircle when anomoulus answer detected
        result = zsla.approxCircle(pointSet, AtomicPoint(0,0,z,""), n)
        #Makes new calcPoint
        
        #Try printing each of the atoms in the point set then set up testing 
        #approx circle in a new file
        
        
        #result[0].printCoords()
        
        pathList[i] = result[0]
        layerAns = result[0]
        
        calcPoint = AtomicPoint(result[0].x, result[0].y, 0, "")
    
    
    i += 1
    
    
    
        
    

    #for s in result[1]:
    #    th.plotAtom(s, ax, color="purple")
        
    


    plt.show()
    
    #layerAns.printCoords()

end = time.time()


print(end-start)
print(len(pathList))

#Make animation
#Works!

#Add small central black circle for central path vizualization
'''
fig, ax = plt.subplots(1,1)
fig.set_size_inches(5,5)
fig.tight_layout(rect=[0, 0.03, 1, 0.95])

def animate(i):
    ax.clear()
    
    zPoints = zPointsList[i]
    resultCircle = pathList[i]

    th.plotAtom(resultCircle, ax, color="blue")

    for zp in zPoints:
        th.plotAtom(zp, ax)

    ax.set_xlim(bounds[1], bounds[0])
    ax.set_ylim(bounds[3], bounds[2])

ani = FuncAnimation(fig, animate, frames=len(pathList), repeat=False)

plt.close()

ani.save("zSlice-Frame-113-walkthrough.gif", dpi=300, writer=PillowWriter(fps=10))
'''

#for layer in pathList:
#    layer.printCoords()

'''
fig, ax = plt.subplots()

for p in pointList:
    th.plotAtom(p, ax)
    
th.plotAtom(result[0], ax, color="blue")

for s in result[1]:
    th.plotAtom(s, ax, color="purple")
    
ax.set_xlim(bounds[1], bounds[0])
ax.set_ylim(bounds[3], bounds[2])


plt.show()
'''