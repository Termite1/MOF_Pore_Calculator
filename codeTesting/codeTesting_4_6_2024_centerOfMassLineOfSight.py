# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 12:15:32 2024

for testing centerOfMassLineOfSight code

@author: samda
"""
import centerOfMassLineOfSight as mlos
from AtomicPoint import AtomicPoint
import pdbPoreReader as pdb


file = "MOF-5 middle pore frames 111-161.pdb"

frames = pdb.multiframePDBReader(file)

current_frame = frames[0]

res1 = mlos.centerOfMass(current_frame, True)
res2 = mlos.centerOfMass(current_frame, False)

res1.printCoords()
res2.printCoords()