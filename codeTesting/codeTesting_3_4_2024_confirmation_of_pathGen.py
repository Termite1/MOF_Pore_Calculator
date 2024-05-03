# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 15:17:08 2024

@author: samda
"""
import zSliceLineAproximation as zsla
import pdbPoreReader as pdb
import time


file = 'MOF-5 middle pore frames 141-151.pdb'
#pdbList = pdbPoreReader.atomReaderPDB(file)

all_frames = pdb.multiframePDBReader(file)

#print(len(frames)) - 11 - [0-10]

current_frame = all_frames[3]

'''
totalStart = time.time()

for frame in all_frames:
    #IT WORKS!!! Do time test tomorrow and try to get the updates loaded into pathGen
    #generate animation as well?
    #Once pathGen works, move onto consecutive pore calc

    start = time.time()

    pathList = zsla.pathGen(frame)

    end = time.time()
    
    print(end - start)
    
totalEnd = time.time()
print("")
print(totalEnd - totalStart)
'''

pathList = zsla.pathGen(current_frame)

pdb.pdbWrite(pathList, 269)