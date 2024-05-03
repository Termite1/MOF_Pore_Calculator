# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 21:36:05 2024

@author: samda
"""

import pdbPoreReader
import lcfProcessing as proc
from AtomicPoint import AtomicPoint
import zSliceLineAproximation as zsla
import AtomicGraph as ag
import TPAAGraphFunctions as tpaa
import LargestCircleFinder as lcf
import PoreCalculatorV0Functions as pcf

file = '05 3x3x3 middle pore frame 121.pdb'

pdbList = pdbPoreReader.atomReaderPDB(file)

prebounds = lcf.boundsFinder(pdbList)
#calculates center of bounded area, generates list of pore points from there
#could probably come up with something more eligant or just use circle calc. Circle calc won't work with occluded z axis though
center = AtomicPoint((prebounds[0] + prebounds[1]) / 2, (prebounds[2] + prebounds[3]) / 2, (prebounds[4] + prebounds[5]) / 2, "")

#generate limited list from pore file
pointList = pcf.listLimiter(center, pdbList, 100)

#generates sliceList
sliceList = zsla.pathGen(pdbList)
for p in sliceList:
    p.VanDerWaalsRadius = p.VanDerWaalsRadius * 0.2
    
#Creates graph for checking what caps the calculated circles touch
myGraph = ag.AtomicGraph()
tpaa.genTPAACapGraph(myGraph, pdbList)


a = AtomicPoint(-1.863, -2.096, 33.412, 'H')
b = AtomicPoint(3.938, -0.983, 35.325, 'C')
c = AtomicPoint(1.780, 3.137, 34.804, 'C')


solution = pcf.circleCalculation(a, b, c, pdbList, sliceList, myGraph)
print(solution[0])
