# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 10:28:02 2024

Testing/walkthrough of growCircle function

And math for circleApprox

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


#Environment setup
file = '05 3x3x3 middle pore frame 112.pdb'
pdbList = pdbPoreReader.atomReaderPDB(file)

#z = 30.586000000000034
#calcPoint = AtomicPoint(-0.4385308775314408,-0.4077724205622452,30.586000000000034, "")
z = 29.896000000000008

#calcPoint = AtomicPoint(-0.4385308775314408,-0.4077724205622452,29.896000000000008,"")

pointList = zsla.sliceList(pdbList, z)
bounds = lcf.boundsFinder(pointList)

'''
#Math for circleApprox
When z =  29.896000000000008
A=-2.2047839195375136
B=-6.011207097735543
C=-0.6870376605527184
D=1.1441527440879364
E=0.0
F=0.0
r1=0.9051403206133399


#print((A*B + C*D + E*F - r1)**2) #squared
print(133.68438561230516 - 158.69700168126502)

#print((A**2 + C**2 + E**2 -1) * (B**2 + D**2 + F**2 - r1**2)) #squared

print() #squared

#All devided by
print(A**2 + C**2 + E**2 - 1)

#show moleculeas
a = AtomicPoint(5.028,1.462,29.896000000000008,"", vdwr=0.9051403206133399)
b = AtomicPoint(3.589,6.814,29.896000000000008, "", vdwr = 1.4094818196770078)
c = AtomicPoint(2.959,7.207,29.896000000000008,"", vdwr = 0.29047375096559236)
'''



'''

atomList = [a,b,c]

fig, ax = plt.subplots()

for p in pointList:
    th.plotAtom(p, ax)
    
th.plotAtom(a, ax, color="purple")
th.plotAtom(b, ax, color="blue")
th.plotAtom(c, ax, color="orange")

ax.set_xlim(bounds[1], bounds[0])
ax.set_ylim(bounds[3], bounds[2])
  
plt.show()
'''