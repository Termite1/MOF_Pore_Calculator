# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 09:47:36 2024

Testing approxCircle on frame with few atoms

@author: samda
"""
from AtomicPoint import AtomicPoint
import zSliceLineAproximation as zsla
import PoreCalculatorV0Functions as pcf
import LargestCircleFinder as lcf
import lcfProcessing as proc
import matplotlib.pyplot as plt
import testHelpers as th
import pdbPoreReader as pdb


file = 'MOF-5 middle pore frames 111-121.pdb'
#pdbList = pdbPoreReader.atomReaderPDB(file)

frames = pdb.multiframePDBReader(file)

pdbList = frames[2]

z = 31.806000000000118
calcPoint = AtomicPoint(-0.6083174926580314,0.8871399979064402,0, "")
calcPoint.z = z

pointSet = zsla.sliceList(pdbList, z)
bounds = lcf.boundsFinder(pointSet)
n=40

result_1=zsla.approxCircle(pointSet, calcPoint, n)


fig, ax = plt.subplots()

for p in pointSet:
    th.plotAtom(p, ax)

th.plotAtom(result_1[0], ax, color="blue")

for atom in result_1[1]:
    th.plotAtom(atom, ax, color='purple')

ax.set_xlim(bounds[1], bounds[0])
ax.set_ylim(bounds[3], bounds[2])

#ax.set_xlim(-7.5, 0)
#ax.set_ylim(-5, 0)

plt.show()



result_2 = lcf.growCircle(pointSet, calcPoint, bounds)
result_2[1].VanDerWaalsRadius=result_2[0]
result_2[1].printCoords()
for atom in result_2[2]:
    atom.printCoords()




fig, ax = plt.subplots()

for p in pointSet:
    th.plotAtom(p, ax)

th.plotAtom(result_2[1], ax, color="blue")

for atom in result_2[2]:
    th.plotAtom(atom, ax, color='purple')

ax.set_xlim(bounds[1], bounds[0])
ax.set_ylim(bounds[3], bounds[2])

#ax.set_xlim(-7.5, 0)
#ax.set_ylim(-5, 0)

plt.show()