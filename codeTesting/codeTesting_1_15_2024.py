# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 15:27:58 2023

@author: samda
"""

import LargestCircleFinder as lcf
from AtomicPoint import AtomicPoint 
import pdbPoreReader as pdb
import PoreCalculatorV0Functions as pcf
import lcfProcessing as proc
import matplotlib.pyplot as plt
import time

#Run only for circles that are sourrounded?


file = '05 3x3x3 middle pore frame 121.pdb'
pdbList = pdb.atomReaderPDB(file)
bounds = lcf.boundsFinder(pdbList)

z = (bounds[5] + bounds[4]) / 2

center = AtomicPoint(0, 0, z, "")
pointList = pcf.listLimiter(center, pdbList, 100)


start1 = time.time()


for a in range(0, len(pointList) - 1):
    for b in range(1, len(pointList) - 1):
        if b == a: continue
        for c in range(2, len(pointList) - 1):
            if c == a or c == b: continue
            
            #My code
            plane = proc.planarize(pointList[a], pointList[b], pointList[c])
            pa = proc.planarizePoint(pointList[a], plane, False)
            pb = proc.planarizePoint(pointList[b], plane, False)
            pc = proc.planarizePoint(pointList[c], plane, False)


            new_point_list = [pa, pb, pc]
            bounds = lcf.boundsFinder(new_point_list)


            #projected Center is still in 2D plane
            pcenter = lcf.threeCircleCircumcenter(pa, pb, pc)
            r_1 = pcenter.distBetween2D(pa, True)

            #circle center in (x,y,z) coords
            center_1 = proc.planarizePoint(pcenter, plane, True)
            center_1.VanDerWaalsRadius = r_1

end1 = time.time()
start2 = time.time()

for a in range(0, len(pointList) - 1):
    for b in range(1, len(pointList) - 1):
        if b == a: continue
        for c in range(2, len(pointList) - 1):
            if c == a or c == b: continue

            #Grimm code
            circle_2 = lcf.threeAtomCircumcenter(pointList[a], pointList[b], pointList[c])             

#end2 = time.time()

#print(f'My code took {end1 - start1} sec to run for 100*99*98 entries')
#print(f"Grimm's math took {end2 - start2} sec to run for 100*99*98 entries")



            #Check if circles are the same
            if abs(circle_2.x - center_1.x) > 0.001 or abs(circle_2.y - center_1.y) > 0.001 or abs(circle_2.z - center_1.z) > 0.001:
                dist11 = center_1.distBetween(pointList[a], True)
                dist12 = center_1.distBetween(pointList[b], True)
                dist13 = center_1.distBetween(pointList[c], True)
                
                dist21 = circle_2.distBetween(pointList[a], True)
                dist22 = circle_2.distBetween(pointList[b], True)
                dist23 = circle_2.distBetween(pointList[c], True)
                
                check11 = False
                check12 = False
                check13 = False
                check21 = False
                check22 = False
                check23 = False
                
                if abs(dist11 - dist12) < 0.001: check11 = True    
                if abs(dist12 - dist13) < 0.001: check12 = True
                if abs(dist11 - dist13) < 0.001: check13 = True
                if abs(dist21 - dist22) < 0.001: check21 = True
                if abs(dist22 - dist23) < 0.001: check22 = True
                if abs(dist21 - dist23) < 0.001: check23 = True

                valid_1 = True
                valid_2 = True
                
                if check11 == False or check12 == False or check13 == False: valid_1 = False
                if check21 == False or check22 == False or check23 == False: valid_2 = False
                
                #if valid_1 == valid_2: continue
                if valid_1 == True and valid_2 == True: continue 
                
                print("Error. My code found:")
                center_1.printCoords()
                center_1.VanDerWaalsRadius = 0

                print("")

                print(center_1.distBetween(pointList[a], True))
                print(center_1.distBetween(pointList[b], True))
                print(center_1.distBetween(pointList[c], True))

                print("")

                print("Grimm's math found:")
                circle_2.printCoords()            
                circle_2.VanDerWaalsRadius = 0

                print("")

                print(circle_2.distBetween(pointList[a], True))
                print(circle_2.distBetween(pointList[b], True))
                print(circle_2.distBetween(pointList[c], True))

                print("")

                print("Defining Atoms:")
                pointList[a].printCoords()
                pointList[b].printCoords()
                pointList[c].printCoords()
                print("") 
                print('')
                



'''
a = 0
b = 1
c = 2


atom1 = AtomicPoint(3.025,3.768,35.044, '', vdwr=1.7)
atom2 = AtomicPoint(4.474,-0.757,37.688, '', vdwr=1.7)
atom3 = AtomicPoint(3.86,0.873,36.345, '', vdwr=1.1)

pointList = [atom1, atom2, atom3]


#My code
plane = proc.planarizeV2(pointList[a], pointList[b], pointList[c])
pa = proc.planarizePointV2(pointList[a], plane, False)
pb = proc.planarizePointV2(pointList[b], plane, False)
pc = proc.planarizePointV2(pointList[c], plane, False)


new_point_list = [pa, pb, pc]
bounds = lcf.boundsFinder(new_point_list)


#projected Center is still in 2D plane
pcenter = lcf.threeCircleCircumcenter(pa, pb, pc)
r_1 = pcenter.distBetween2D(pa, True)

#circle center in (x,y,z) coords
center_1 = proc.planarizePointV2(pcenter, plane, True)
center_1.VanDerWaalsRadius = r_1

#Grimm code
circle_2 = lcf.threeAtomCircumcenter(pointList[a], pointList[b], pointList[c])      


#Check if circles are the same
if abs(circle_2.x - center_1.x) > 0.001 or abs(circle_2.y - center_1.y) > 0.001 or abs(circle_2.z - center_1.z) > 0.001:
    print("Error. My code found:")
    center_1.printCoords()
    temp1 = center_1.VanDerWaalsRadius
    center_1.VanDerWaalsRadius = 0

    print("")

    print(center_1.distBetween(pointList[a], True))
    print(center_1.distBetween(pointList[b], True))
    print(center_1.distBetween(pointList[c], True))

    print("")

    print("Grimm's math found:")
    circle_2.printCoords()      
    temp2 = circle_2.VanDerWaalsRadius      
    circle_2.VanDerWaalsRadius = 0

    print("")

    print(circle_2.distBetween(pointList[a], True))
    print(circle_2.distBetween(pointList[b], True))
    print(circle_2.distBetween(pointList[c], True))

    print("")

    print("Defining Atoms:")
    pointList[a].printCoords()
    pointList[b].printCoords()
    pointList[c].printCoords()
    print("") 
    
    center_1.VanDerWaalsRadius = temp1
    circle_2.VanDerWaalsRadius = temp2
    

print('Planarized Points:')
pa.printCoords()
pb.printCoords()
pc.printCoords()

print("")
pcircle = proc.planarizePointV2(center_1, plane, False)
pcircle.printCoords()
pcircle2 = proc.planarizePointV2(circle_2, plane, False)
pcircle2.printCoords()


fig, ax = plt.subplots()           
   
ax.set_ylim(30,40)
ax.set_xlim(12,22)


pltCircle1 = plt.Circle((pa.x, pa.y), pa.VanDerWaalsRadius)          
pltCircle2 = plt.Circle((pb.x, pb.y), pb.VanDerWaalsRadius)
pltCircle3 = plt.Circle((pc.x, pc.y), pc.VanDerWaalsRadius)

pltCircleCircle = plt.Circle((pcircle.x, pcircle.y), pcircle.VanDerWaalsRadius, color="r")
pltCircleCircle2 = plt.Circle((pcircle2.x, pcircle2.y), pcircle2.VanDerWaalsRadius, color="r")


ax.add_patch(pltCircleCircle)


ax.add_patch(pltCircle1)
ax.add_patch(pltCircle2)
ax.add_patch(pltCircle3)

ax.plot

'''
                
                


