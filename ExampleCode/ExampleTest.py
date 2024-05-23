# -*- coding: utf-8 -*-
"""
Created on Wed May 22 14:23:21 2024

Example application of the code to frame 114 of the MOF-5 TPAA capped simulation
to find the largest solid sphere that can pass through the pore this frame.

@author: samda
"""
from AtomicPoint import AtomicPoint
import pdbPoreReader as pdb
import ExampleCode as ec

file = 'test_frame_114.pdb'
n = 50 #Nearest number of atoms to the center to be selected for further calculation
z = 35 #z partition height
center = AtomicPoint(-0.180, 0.777, 35.103, "") #location of solution to previous frame

pointList = pdb.atomReaderPDB(file) #Turns pdb file into a list of AtomicPoints

#Finds the largest solid sphere that can pass through the pore
solution = ec.poreCalculator(pointList, center, n, z) 

#vdwr is the radius (Van der waals radius for atoms)
solution[1].printCoords()
print("")
pdb.pdbWrite([solution[1]],1)