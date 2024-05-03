# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 15:24:24 2024

Generating histogram of largest angle when threeAtom circle returns non-answers

@author: samda
"""
import LargestCircleFinder as lcf
import lcfProcessing as proc
import pdbPoreReader as pdb
from AtomicPoint import AtomicPoint
import PoreCalculatorV0Functions as pcf


file = '05 3x3x3 middle pore frame 121.pdb'
pdbList = pdb.atomReaderPDB(file)
bounds = lcf.boundsFinder(pdbList)

z = (bounds[5] + bounds[4]) / 2

center = AtomicPoint(0, 0, z, "")
pointList = pcf.listLimiter(center, pdbList, 100)

badSolutions = []

for a in range(0, len(pointList) - 1):
    for b in range(1, len(pointList) - 1):
        if b == a: continue
        for c in range(2, len(pointList) - 1):
            if c == a or c == b: continue

            #Grimm code
            circle = lcf.threeAtomCircumcenter(pointList[a], pointList[b], pointList[c])
            
            #check if distance to atpms from circle center is the same
            dist1 = circle.distBetween(pointList[a], True)
            dist2 = circle.distBetween(pointList[b], True)
            dist3 = circle.distBetween(pointList[c], True)
            
            check1 = False
            check2 = False
            check3 = False
            
            if abs(dist1 - dist2) > 0.001: check1 = True    
            if abs(dist2 - dist3) > 0.001: check2 = True
            if abs(dist1 - dist3) > 0.001: check3 = True
            
            #Answer incorrect, find and return the largest angle as well as the surrounding atoms
            if check1 or check2 or check3:
                side1 = pointList[a].distBetween(pointList[b], False)
                side2 = pointList[b].distBetween(pointList[c], False)
                side3 = pointList[a].distBetween(pointList[c], False)
                
                angleList = []
                angleList.append(lcf.cosLawAngle(side1, side2, side3))
                angleList.append(lcf.cosLawAngle(side2, side3, side1))
                angleList.append(lcf.cosLawAngle(side3, side1, side2))
                
                #append solution to list
                badSolutions.append(max(angleList)) #, [pointList[a], pointList[b], pointList[c]] add to get atoms

#Write the bad solutions down
with open('outputAngles.txt', 'a') as f:
    for line in badSolutions:
        f.write(str(line) + '\n')
                