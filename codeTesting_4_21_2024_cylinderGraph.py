# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 14:18:01 2024

@author: samda
"""

import pdbPoreReader as pdb
import zSliceLineAproximation as zsla

file = "MOF-5 middle pore frames 111-161.pdb"
file2 = "allAns_111-161.pdb"

all_frames = pdb.multiframePDBReader(file)
all_Ans = pdb.atomReaderPDB(file2)

count = -1

for frame in all_frames:
    count += 1
    if count == 0: continue
    
    path = zsla.pathGen(frame, smallest_r=True)
    print(path[1])

    #sol = mcrs.MonteCarloSet(50, 1, frame, center=all_Ans[count - 1], numPoints=50, percentVDWR=0.4, topAns=20)