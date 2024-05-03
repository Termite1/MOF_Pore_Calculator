# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 11:17:50 2024

@author: samda
"""

import PoreCalculatorV0Functions as pcf
import pdbPoreReader as pdb
from AtomicPoint import AtomicPoint
import MonteCarloRandomSamplingV0 as mcrs
import consecutiveFramesPoreCalculator as cfpc
import zSliceLineAproximation as zsla
import AtomicGraph as ag
import TPAAGraphFunctions as tpaa
import cylinderGraph as cg
import LargestCircleFinder as lcf

prev_sol = AtomicPoint( -0.658 ,  0.173 , 35.106, "")

file = "MOF-5 middle pore frames 111-161.pdb"
file2 = "allAns_111-161.pdb"

all_frames = pdb.multiframePDBReader(file)
all_Ans = pdb.atomReaderPDB(file2)

#Check frame 135, not getting answer I should. Same for 117
#Think about stats concerning surounding atoms. Might be worth looking at, see if 
#there are any obvious patterns or limitations. Try view to center of mass next?
current_frame = all_frames[16]

#get nearest 50 atoms
frame_atoms = cfpc.atomFinderNum(prev_sol, 50, current_frame)

#approx central path
sliceList = zsla.pathGen(current_frame)
#Reducing vdwr size
for p in sliceList:
    p.VanDerWaalsRadius = p.VanDerWaalsRadius * 0.4


#Creates graph for checking what caps the calculated circles touch
myGraph = ag.AtomicGraph()
tpaa.genTPAACapGraph(myGraph, current_frame)

#This is producing duplicates (x6) - fix
sol = pcf.poreCalcL(frame_atoms, current_frame, sliceList, myGraph, allTopAns=True)

#real function
#sol2 = mcrs.MonteCarloSet(50, 1, current_frame, center=prev_sol, numPoints=50, percentVDWR=0.4, allTopAns=34)

for ans in sol:
    ans[1].printCoords()
    print(ans[3][0])
    print(ans[3][1])
    print(ans[3][2])
    print("")

bounds = lcf.boundsFinder(current_frame)

'''
for ans in sol:
    new_node = cg.cylinderNode(ans[1], ans[2][0], ans[2][1], ans[2][2])
    new_node.ap.printCoords()
'''

myCGraph = cg.cylinderGraph()
myCGraph.genGraph(current_frame, sol, bounds)

print("")

print(len(myCGraph.myVerticies))
print(len(myCGraph.myAdjList))


'''
print("")
cluster_list = pcf.clusterPointList(sol, 34, "solution_lowMaxR")
print(len(cluster_list))
print(cluster_list[0][0])


#Sorts the list in decending order by radius
cluster_list = pcf.sortList(cluster_list, "solution_topRadius")
print("")
for item in cluster_list:
    print(item[0])
'''

#print("")
#print(sol2[0])