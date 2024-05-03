# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 15:07:57 2024

Created to verify answers aquired over consecutive frames were the same as 
answers produced from indevidual frames.

@author: samda
"""

import pdbPoreReader as pdb
import consecutiveFramesPoreCalculator as cfpc
import MonteCarloRandomSamplingV0 as mcrs
import time
from AtomicPoint import AtomicPoint


file = "MOF-5 middle pore frames 111-161.pdb"

frames = pdb.multiframePDBReader(file)

#mcrs_solutions = []

'''
i = 0

for frame in frames:
    result = mcrs.MonteCarloSet(40, 15, frame)
    mcrs_solutions.append(result)
    
    i += 1
    print(i)
'''

start = time.time()
#Consec uses a pop function which modifies the frame list by removing the first entry
consec_solutions = cfpc.multiframePoreCalculator(frames)
end = time.time()

#Get calculated sphere location, result[1] for each


#("Monte Carlo Random Sampling Solutions")
#for i in mcrs_solutions:
#    print(i[0])
    
#print("")
#for i in mcrs_solutions:
#    print(f"{i[1].x:.3f}, {i[1].y:.3f}, {i[1].z:.3f}")

'''
print("")
for i in mcrs_solutions:
    for j in i[2]:
        j.printCoords()
    print("")
'''


print("")
print("Consecutive Frames Solutions")
for i in consec_solutions:
    print(f"{i[0]:.3f}")



print("")
for i in consec_solutions:
    pdb.pdbWrite([i[1]], 269)
    #print(f"{i[1].x:.3f}, {i[1].y:.3f}, {i[1].z:.3f} --- {j}")


'''
j = 111

print("")
for i in consec_solutions:
    print(j)
    j += 1
    pdb.pdbWrite(i[2], 270)
    print("")
'''

'''
print("i1")
for i in consec_solutions:
    print(f"{i[3][0]:.3f}")

print("")
print("i2")
for i in consec_solutions:
    print(f"{i[3][1]:.3f}")

print("")
print("i3")
for i in consec_solutions:
    print(f"{i[3][2]:.3f}")
'''

    
print("")
print(f"Took consec frames {end - start} sec to run.")