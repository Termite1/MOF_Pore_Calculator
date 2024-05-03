# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 22:01:12 2024

@author: samda
"""

#Make new file, just run CircleCalculator with atoms that give known solution.
#See what error return is, then use that to find the bug(s)

#Also open old files and give grimm example of overlapping circles where the model breaks

from AtomicPoint import AtomicPoint
import lcfProcessing as proc
import matplotlib.pyplot as plt
import LargestCircleFinder as lcf

a = AtomicPoint(3, 3, 3, "C")
b = AtomicPoint(3, -3, 3, "C")
c = AtomicPoint(-3, 0, -3, "C")
pointList = [a,b,c]
center = lcf.threeAtomCircumcenter(a, b, c)

pvals = proc.planarize(a, b, c)

newList = []

for p in pointList:
    temp_p = proc.planarizePoint(p, pvals)
    newList.append(temp_p)

p_center = proc.planarizePoint(center, pvals)

print(lcf.simpleCircleSurrounded(newList[0], newList[1], newList[2], p_center))

bounds = lcf.boundsFinder(newList)


fig, ax = plt.subplots()

for p in newList:
    p.printCoords()
    circle = plt.Circle((p.x, p.y), p.VanDerWaalsRadius, color='r')
    ax.add_patch(circle)

p_center.printCoords()
center_circle = plt.Circle((p_center.x, p_center.y), p_center.VanDerWaalsRadius, color='blue')
ax.add_patch(center_circle)

ax.set_xlim(bounds[1], bounds[0])
ax.set_ylim(bounds[3], bounds[2])

plt.show()