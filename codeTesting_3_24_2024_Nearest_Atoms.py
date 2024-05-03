# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 12:08:23 2024

@author: samda
"""

from AtomicPoint import AtomicPoint
import pdbPoreReader as pdb
import consecutiveFramesPoreCalculator as cfpc
import MonteCarloRandomSamplingV0 as mcrs
import csv
import centerOfMassLineOfSight as mlos
import lcfProcessing as proc
import math
import numpy as np
import LargestCircleFinder as lcf
import AtomicGraph as ag
import TPAAGraphFunctions as tpaa
import MOF5_Index_Identifiers as m5ii

prev_sol = AtomicPoint(0.214  , 0.489 , 36.107, "")

file = "MOF-5 middle pore frames 111-161.pdb"
file2 = "allAns_111-161.pdb"

all_frames = pdb.multiframePDBReader(file)
all_Ans = pdb.atomReaderPDB(file2)

#Frame 161
#current_frame = all_frames[50]

#nearest_aList = cfpc.atomFinderNum(prev_sol, 40, current_frame)

#pdb.pdbWrite(nearest_aList, 1)

#sol = mcrs.MonteCarloSet(40, 1, current_frame, center=prev_sol, numPoints=40, percentVDWR=0.3)

#print(sol[0])

#pdb.pdbWrite([sol[1]], 41)

'''
sol = mcrs.MonteCarloSet(50, 1, current_frame, center=prev_sol, numPoints=50, percentVDWR=0.4, topAns=20)

for i in sol:
    print(f'{i[0]:.3f}')


print("")
all_sol = []
for i in sol:
    all_sol.append(i[1])

pdb.pdbWrite(all_sol, 1)
'''

'''
#for having all the centers of mass and geometric centers
cm_list = []
gc_list = []
for frame in all_frames:
    center_of_mass = mlos.centerOfMass(frame, True)
    geometric_center = mlos.centerOfMass(frame, False)
    
    cm_list.append(center_of_mass)
    gc_list.append(geometric_center)
    
pdb.pdbWrite(cm_list, 111)

print("")
pdb.pdbWrite(gc_list, 111)  
'''


frame_base = 111
count = -1


#Writes topAns analysis file
with open('topAns_sa_index_analysis.csv', 'w+', newline='') as csvfile:
    #fieldnames = ['Frame Count', 'r', 'x', 'y', 'z', 'dist_cm_xy', 'dist_cm_z', 'dist_gc_xy', 'dist_gc_z', '', 'x_1', 'y_1', 'z_1', 'atom_type_1', 'quadrant_1', 'angle_1', 'see_cm_1', 'dist_cm_xy_1', 'dist_cm_z_1', 'see_gc_1', 'dist_gc_xy_1', 'dist_gc_z_1', '', 'x_2', 'y_2', 'z_2', 'atom_type_2', 'quadrant_2', 'angle_2', 'see_cm_2', 'dist_cm_xy_2', 'dist_cm_z_2', 'see_gc_2', 'dist_gc_xy_2', 'dist_gc_z_2', '', 'x_3', 'y_3', 'z_3', 'atom_type_3', 'quadrant_3', 'angle_3', 'see_cm_3', 'dist_cm_xy_3', 'dist_cm_z_3', 'see_gc_3', 'dist_gc_xy_3', 'dist_gc_z_3']
    fieldnames = ['Frame Count', 'r', 'x', 'y', 'z', '', 'angle_1', 'index_1', 'location_1', '', 'angle_2', 'index_2', 'location_2', '', 'angle_3', 'index_3', 'location_3']
    
    writer = csv.writer(csvfile)

    #Writes file header
    writer.writerows([fieldnames])
    
    #Generates graph for the frame - only need to do once, indexes are the same even if position isn't
    #myGraph = ag.AtomicGraph()
    #tpaa.genTPAACapGraph(myGraph, all_frames[0])
    
    for frame in all_frames:
        count += 1
        if count == 0: continue
        
        sol = mcrs.MonteCarloSet(50, 1, frame, center=all_Ans[count - 1], numPoints=50, percentVDWR=0.4, topAns=20)
        
        #center_of_mass = mlos.centerOfMass(frame, True)
        #geometric_center = mlos.centerOfMass(frame, False)
        
        
        #Repeates for each atom (i.e. answer) in sol
        for a in sol:
            #Skips invalid answers
            if a[1].z == 0: continue
            
            #surface_normal = np.array(proc.planarize(a[2][0], a[2][1], a[2][2])[1])
            #xy_normal = np.array([0,0,1])
            
            #theta = math.asin(mlos.magnitude(np.cross(surface_normal, xy_normal))/(mlos.magnitude(surface_normal) * mlos.magnitude(xy_normal)))
            
            '''
            #Sphere distances to center of mass and geometric center
            dist_cm_xy = center_of_mass.distBetween2D(a[1], False)
            dist_cm_z = center_of_mass.z - a[1].z
            
            dist_gc_xy = geometric_center.distBetween2D(a[1], False)
            dist_gc_z = geometric_center.z - a[1].z
            
            #[Frame Count, r, x, y, z]
            sphere = [frame_base + count, a[0], a[1].x, a[1].y, a[1].z, dist_cm_xy, dist_cm_z, dist_gc_xy, dist_gc_z]
            '''

            sphere = [frame_base + count, a[0], a[1].x, a[1].y, a[1].z]
            
            #Get sides of triangle bounding atoms create round sphere
            side_a = a[2][0].distBetween(a[2][1], False)
            side_b = a[2][0].distBetween(a[2][2], False)
            side_c = a[2][1].distBetween(a[2][2], False)
            
            #angles atoms make in triangle
            angle_0 = lcf.cosLawAngle(side_a, side_b, side_c)
            angle_2 = lcf.cosLawAngle(side_b, side_c, side_a)
            angle_1 = lcf.cosLawAngle(side_a, side_c, side_b)
            
            '''
            #initializes quadrants
            quadrant_0 = None
            quadrant_1 = None
            quadrant_2 = None
            
            #Matches indexes to quadrants. Runs through every atom in the molecular structure
            for vertex in myGraph.myVerticies:
                if vertex.index == a[3][0]: quadrant_0 = vertex.quadrant
                if vertex.index == a[3][1]: quadrant_1 = vertex.quadrant
                if vertex.index == a[3][2]: quadrant_2 = vertex.quadrant
            
            #See center of mass or geometric center, and distances
            see_cm_0 = mlos.lineOfSightCheck(a[2][0], center_of_mass, frame)
            dist_cm_xy_0 = center_of_mass.distBetween2D(a[2][0], False)
            dist_cm_z_0 = center_of_mass.z - a[2][0].z
            
            see_gc_0 = mlos.lineOfSightCheck(a[2][0], geometric_center, frame)
            dist_gc_xy_0 = geometric_center.distBetween2D(a[2][0], False)
            dist_gc_z_0 = geometric_center.z - a[2][0].z
            
            
            see_cm_1 = mlos.lineOfSightCheck(a[2][1], center_of_mass, frame)
            dist_cm_xy_1 = center_of_mass.distBetween2D(a[2][1], False)
            dist_cm_z_1 = center_of_mass.z - a[2][1].z
            
            see_gc_1 = mlos.lineOfSightCheck(a[2][1], geometric_center, frame)
            dist_gc_xy_1 = geometric_center.distBetween2D(a[2][1], False)
            dist_gc_z_1 = geometric_center.z - a[2][1].z
            
            
            see_cm_2 = mlos.lineOfSightCheck(a[2][2], center_of_mass, frame)
            dist_cm_xy_2 = center_of_mass.distBetween2D(a[2][2], False)
            dist_cm_z_2 = center_of_mass.z - a[2][2].z
            
            see_gc_2 = mlos.lineOfSightCheck(a[2][2], geometric_center, frame)
            dist_gc_xy_2 = geometric_center.distBetween2D(a[2][2], False)
            dist_gc_z_2 = geometric_center.z - a[2][2].z
            
            data = []
            
            
            #Go about constructing data line, do so in order of decreasing angle
            #[x, y, z, atom_type, angle, see center of mass, see geometric center] - need to do quadrants soon
            if angle_0 > angle_1 and angle_0 > angle_2:
                data = ['', a[2][0].x, a[2][0].y, a[2][0].z, a[2][0].bondType, quadrant_0, angle_0, see_cm_0, dist_cm_xy_0, dist_cm_z_0, see_gc_0, dist_gc_xy_0, dist_gc_z_0]
                if angle_1 > angle_2:
                    data_1 = ['', a[2][1].x, a[2][1].y, a[2][1].z, a[2][1].bondType, quadrant_1, angle_1, see_cm_1, dist_cm_xy_1, dist_cm_z_1, see_gc_1, dist_gc_xy_1, dist_gc_z_1]
                    data_2 = ['', a[2][2].x, a[2][2].y, a[2][2].z, a[2][2].bondType, quadrant_2, angle_2, see_cm_2, dist_cm_xy_2, dist_cm_z_2, see_gc_2, dist_gc_xy_2, dist_gc_z_2]
                else:
                    data_1 = ['', a[2][2].x, a[2][2].y, a[2][2].z, a[2][2].bondType, quadrant_2, angle_2, see_cm_2, dist_cm_xy_2, dist_cm_z_2, see_gc_2, dist_gc_xy_2, dist_gc_z_2]
                    data_2 = ['', a[2][1].x, a[2][1].y, a[2][1].z, a[2][1].bondType, quadrant_1, angle_1, see_cm_1, dist_cm_xy_1, dist_cm_z_1, see_gc_1, dist_gc_xy_1, dist_gc_z_1]
            elif angle_1 > angle_0 and angle_1 > angle_2:
                data = ['', a[2][1].x, a[2][1].y, a[2][1].z, a[2][1].bondType, quadrant_1, angle_1, see_cm_1, dist_cm_xy_1, dist_cm_z_1, see_gc_1, dist_gc_xy_1, dist_gc_z_1]
                if angle_0 > angle_2:
                    data_1 = ['', a[2][0].x, a[2][0].y, a[2][0].z, a[2][0].bondType, quadrant_0, angle_0, see_cm_0, dist_cm_xy_0, dist_cm_z_0, see_gc_0, dist_gc_xy_0, dist_gc_z_0]
                    data_2 = ['', a[2][2].x, a[2][2].y, a[2][2].z, a[2][2].bondType, quadrant_2, angle_2, see_cm_2, dist_cm_xy_2, dist_cm_z_2, see_gc_2, dist_gc_xy_2, dist_gc_z_2]
                else:
                    data_1 = ['', a[2][2].x, a[2][2].y, a[2][2].z, a[2][2].bondType, quadrant_2, angle_2, see_cm_2, dist_cm_xy_2, dist_cm_z_2, see_gc_2, dist_gc_xy_2, dist_gc_z_2]
                    data_2 = ['', a[2][0].x, a[2][0].y, a[2][0].z, a[2][0].bondType, quadrant_0, angle_0, see_cm_0, dist_cm_xy_0, dist_cm_z_0, see_gc_0, dist_gc_xy_0, dist_gc_z_0]
            else:
                data = ['', a[2][2].x, a[2][2].y, a[2][2].z, a[2][2].bondType, quadrant_2, angle_2, see_cm_2, dist_cm_xy_2, dist_cm_z_2, see_gc_2, dist_gc_xy_2, dist_gc_z_2]
                if angle_0 > angle_1:
                    data_1 = ['', a[2][0].x, a[2][0].y, a[2][0].z, a[2][0].bondType, quadrant_0, angle_0, see_cm_0, dist_cm_xy_0, dist_cm_z_0, see_gc_0, dist_gc_xy_0, dist_gc_z_0]
                    data_2 = ['', a[2][1].x, a[2][1].y, a[2][1].z, a[2][1].bondType, quadrant_1, angle_1, see_cm_1, dist_cm_xy_1, dist_cm_z_1, see_gc_1, dist_gc_xy_1, dist_gc_z_1]
                else:
                    data_1 = ['', a[2][1].x, a[2][1].y, a[2][1].z, a[2][1].bondType, quadrant_1, angle_1, see_cm_1, dist_cm_xy_1, dist_cm_z_1, see_gc_1, dist_gc_xy_1, dist_gc_z_1]
                    data_2 = ['', a[2][0].x, a[2][0].y, a[2][0].z, a[2][0].bondType, quadrant_0, angle_0, see_cm_0, dist_cm_xy_0, dist_cm_z_0, see_gc_0, dist_gc_xy_0, dist_gc_z_0]
            
            #construct data line
            full_line = sphere + data + data_1 + data_2
            '''
            
            #indexes a[3][0-2]
            
            location_0 = m5ii.identifyMOF5Index((frame_base + count), a[3][0])
            location_1 = m5ii.identifyMOF5Index((frame_base + count), a[3][1])
            location_2 = m5ii.identifyMOF5Index((frame_base + count), a[3][2])
            
            data = []
            
            
            #Go about constructing data line, do so in order of decreasing angle
            #[x, y, z, atom_type, angle, see center of mass, see geometric center] - need to do quadrants soon
            if angle_0 > angle_1 and angle_0 > angle_2:
                data = ['', angle_0, a[3][0], location_0]
                if angle_1 > angle_2:
                    data_1 = ['', angle_1, a[3][1], location_1]
                    data_2 = ['', angle_2, a[3][2], location_2]
                else:
                    data_1 = ['', angle_2, a[3][2], location_2]
                    data_2 = ['', angle_1, a[3][1], location_1]
            elif angle_1 > angle_0 and angle_1 > angle_2:
                data = ['', angle_1, a[3][1], location_1]
                if angle_0 > angle_2:
                    data_1 = ['', angle_0, a[3][0], location_0]
                    data_2 = ['', angle_2, a[3][2], location_2]
                else:
                    data_1 = ['', angle_2, a[3][2], location_2]
                    data_2 = ['', angle_0, a[3][0], location_0]
            else:
                data = ['', angle_2, a[3][2], location_2]
                if angle_0 > angle_1:
                    data_1 = ['', angle_0, a[3][0], location_0]
                    data_2 = ['', angle_1, a[3][1], location_1]
                else:
                    data_1 = ['', angle_1, a[3][1], location_1]
                    data_2 = ['', angle_0, a[3][0], location_0]
            
            #construct data line
            full_line = sphere + data + data_1 + data_2
            
            #write it
            writer.writerows([full_line])
            
        #write spacer
        writer.writerows([['']])
