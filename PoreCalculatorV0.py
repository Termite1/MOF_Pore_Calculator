# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 12:42:34 2023

@author: samda
"""

import PoreCalculatorV0Functions as pcf
import MonteCarloRandomSamplingV0 as mcrs
import time
from AtomicPoint import AtomicPoint
import pdbPoreReader as pdb
 
file = '05 3x3x3 middle pore frame 112.pdb'
       
tStart = time.perf_counter()

#result = pcf.poreCalcF(file)
result = mcrs.MonteCarloSet(100, 1, None, fileName=file)

tEnd = time.perf_counter()
time_diff = tEnd - tStart

print(f'Largest circle is located at ({result[1].x:.3f},{result[1].y:.3f},{result[1].z:.3f})')
print(f'It has radius = {result[0]}')
print(f'The points it is bounded by are ({result[2][0].x:.3f},{result[2][0].y:.3f},{result[2][0].z:.3f}), ({result[2][1].x:.3f},{result[2][1].y:.3f},{result[2][1].z:.3f}), and ({result[2][2].x:.3f},{result[2][2].y:.3f},{result[2][2].z:.3f})')
print(result[2][0].element)
print(result[2][1].element)
print(result[2][2].element)
print(f'The program took {time_diff:.3f} seconds to run')


pdb.pdbWrite(result[2], 280)
