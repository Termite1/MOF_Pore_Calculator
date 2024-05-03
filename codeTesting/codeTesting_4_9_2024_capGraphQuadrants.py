# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 13:55:50 2024

@author: samda
"""
import pdbPoreReader as pdb
import AtomicGraph as ag
import TPAAGraphFunctions as tpaa

file = "MOF-5 middle pore frames 111-161.pdb"
all_frames = pdb.multiframePDBReader(file)


frame = all_frames[23]

myGraph = ag.AtomicGraph()
tpaa.genTPAACapGraph(myGraph, frame)

'''
for v in myGraph.myVerticies:
    v.ap.printCoords()
    print(v.quadrant)
    print(v.index)
    print("")
'''