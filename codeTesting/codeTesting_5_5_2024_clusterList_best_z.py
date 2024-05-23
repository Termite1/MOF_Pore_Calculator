# -*- coding: utf-8 -*-
"""
Created on Sun May  5 10:24:25 2024

@author: samda
"""

import pdbPoreReader as pdb
import consecutiveFramesPoreCalculator as cfpc
import csv
import MonteCarloRandomSamplingV0 as mcrs

file = "MOF-5 middle pore frames 111-161.pdb"
all_frames = pdb.multiframePDBReader(file)

file2 = "allAns_111-161.pdb"
all_Ans = pdb.atomReaderPDB(file2)

suspected_solutions_r = [2.993,2.867,2.486,2.819,2.895,3.051,3.289,3.548,2.983,3.069,2.761,3.02,3.355,3.442
,2.578,2.681,1.936,2.355,2.513,3.019,2.68,3.332,3.047,3.052,2.876,2.94,2.558,2.455,2.688,2.593,2.652,2.511
,1.938,1.926,2.193,2.593,3.226,3.047,2.193,1.644,2.422,2.518,2.326,2.454,2.664,2.22,2.334,2.219,2.128,2.66,2.577]

all_sol = []

with open('clusterList_z_optimization_4.csv', 'w+', newline='') as csvfile:
    #fieldnames = ['z', 'diff_avg', 'diff_max', 'diff_max_index', 'num_0', 'num<=0.2, num>0.2']

    writer = csv.writer(csvfile)

    #Writes file header
    #writer.writerow(fieldnames)

    #runs from 32.5 to 37
    for z_init in range(130,148):
        z = z_init/4 
        
        #writer.writerow(["{:.2f}".format(z)])
        
        #sol_list = [] #cfpc.multiframePoreCalculator(all_frames,z_split=z)
        
        sol_list = cfpc.multiframePoreCalculator(all_frames, z_split=z)
        sol_list = [[z]] + sol_list
        
        #all_sol.append(sol_list)
        
        
        #print("")
        for i in sol_list:
            #pdb.pdbWrite([i[1]], 269)
            writer.writerow([i[0]])
        
        writer.writerow([""])
        
        '''
        discrepancy_list = []
        diff_tot = 0
        
        num_0 = 0 #number with no diff
        num_l2 = 0 #number with diff <= 0.2
        num_g2 = 0 #number with diff > 0.2
        
        
        for i in range(len(sol_list)):
            diff = abs(sol_list[i][0] - suspected_solutions_r[i])
            discrepancy_list.append(diff)
            diff_tot += diff
            
            if diff < 0.001: num_0 += 1
            elif diff > 0.2: num_g2 += 1
            else: num_l2 += 1
            
        diff_avg = diff_tot / len(sol_list)
        diff_max = max(discrepancy_list)
        diff_max_index = discrepancy_list.index(diff_max)
        '''
    
        #row = [z, diff_avg, diff_max, num_0, num_l2, num_g2]
    '''
    #runs through each row
    for i in range(len(all_sol[0])):
        this_line = [] #used to construct write row statement
        
        #runs through each column
        for j in range(len(all_sol)):
            temp = all_sol[j][i][0]
            this_line.append(temp)
        
        writer.writerow(this_line)
        '''
    