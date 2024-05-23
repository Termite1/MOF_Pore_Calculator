# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 22:04:37 2024

File contains code that calculates the largest solid sphere that can fit through 
the given pore over multiple consecutive frames.

@author: samda
"""
import pdbPoreReader as pdb
import consecutiveFramesPoreCalculator as cfpc
import time


file = 'MOF-5 middle pore frames 111-121.pdb'
       
tStart = time.perf_counter()

frames = pdb.multiframePDBReader(file)
solutions = cfpc.multiframePoreCalculator(frames)

for i in solutions:
    print(i[0])

tEnd = time.perf_counter()
time_diff = tEnd - tStart


print(f'The program took {time_diff:.3f} seconds to run')