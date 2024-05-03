# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 09:06:00 2024

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
import testHelpers as th
import growCircle as gc

#Environment Setup - Start

file = '05 3x3x3 middle pore frame 112.pdb'
pdbList = pdbPoreReader.atomReaderPDB(file)

z = 30.586000000000034
inputPoint = AtomicPoint(-0.4385308775314408,-0.4077724205622452,30.586000000000034, "")

pointSet = zsla.sliceList(pdbList, z)
bounds = lcf.boundsFinder(pointSet)

#ap1 = th.nearestAtom(inputPoint, pointSet, [])
#currentRadius = inputPoint.distBetween2D(ap1, True)

#ap1.printCoords()

#Environment Setup - End


#Do full walkthrough of function --> look at and plot circles

result = gc.growCircle(pointSet, inputPoint, bounds)
print(result[0])

#Plot Creation - Start

#result = AtomicPoint(-3.4188238400801056,-2.9120287052730207,32.145228655530836,"")
#p1 = AtomicPoint(-1.187,-4.751,32.03, "")
#p2 = AtomicPoint(-5.19,-0.249,31.768, "")


fig, ax = plt.subplots()

for p in pointSet:
    th.plotAtom(p, ax)

#th.plotAtom(result, ax, color="blue", radius=1.9985843668253662)
#th.plotAtom(p1, ax, color="orange", radius=0.8955305689925394)
#th.plotAtom(p2, ax, color="purple", radius=1.221833049152002)

'''
th.plotAtom(ap1, ax, color="orange")
#th.plotAtom(point2, ax, color="purple")
#th.plotAtom(point3, ax, color="black")
#th.plotAtom(theoryCenter, ax, radius=currentRadius, color="blue")
th.plotAtom(inputPoint, ax, color="blue", radius=currentRadius)


x1 = np.array([inputPoint.x, -3 * math.cos(theta)])
y1 = np.array([inputPoint.y, -3 * math.sin(theta)])
plt.plot(x1, y1, color = "orange")

x2 = np.array([inputPoint.x, ap1.x])
y2 = np.array([inputPoint.y, ap1.y])
plt.plot(x2, y2, color = "purple")
'''

ax.set_xlim(bounds[1], bounds[0])
ax.set_ylim(bounds[3], bounds[2])

#ax.set_xlim(-7.5, 0)
#ax.set_ylim(-5, 0)

plt.show()