# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 17:55:45 2023

Set of functions that find the largest solid sphere that can pass through the 
given pore, finding the results from a number of consecutive frames produced by 
a molecular dynamics simulation.

@author: samda
"""

from AtomicPoint import AtomicPoint
import zSliceLineAproximation as zsla
import PoreCalculatorV0Functions as pcf
import AtomicGraph as ag
import TPAAGraphFunctions as tpaa
import LargestCircleFinder as lcf
import MonteCarloRandomSamplingV0 as mcrs


def atomFinderNum(p, n, atom_list):
    '''Function finds the n closest atoms to point p from the given list of atoms.
    Atoms in atom_list and point p should be in the form of AtomicPoints.
    
    @peram p AtomicPoint
    @peram n Number of atoms to find and return
    @peram atom_list List of atoms to find atoms from
    @return List of n nearest atoms to point p
    '''
    atoms_dist = []
    center = AtomicPoint(p.x, p.y, p.z, "")
    
    #Repeats for each atom in the given list of atom
    for atom in atom_list:
        distance = center.distBetween(atom, True)
        atoms_dist.append([distance, atom])
    
    atoms_dist.sort()
    #grabs nearest n atoms
    solution_list = atoms_dist[:n]
    
    solutions = []
    #Grabs the atoms from the sortted list and returns them
    for pair in solution_list:
        solutions.append(pair[1])
    
    return solutions 


def atomFinderDist(p, dist, atom_list):
    '''Function finds the atoms within dist angstroms of sphere p from the 
    given list of atoms. Atoms in atom_list and sphere p should be in the form 
    of AtomicPoints, with p.VanDerWaalsRadius equal to the sphere's radius.
    
    @peram p AtomicPoint
    @peram dist Distance within to find atoms
    @peram atom_list List of atoms to find atoms from
    @return List atoms within dist angstroms of sphere p
    '''
    solutions = []

    #Repeats for each atom in the given list of atom
    for atom in atom_list:
        #checks if atoms is within dist angstroms of p, if so adds it to the list
        if p.distBetween(atom, True) < dist:
            solutions.append(atom)
            
    return solutions
    

def multiframePoreCalculator(frame_list):
    '''Function finds the largest solid sphere that can pass through each of 
    the given frames of the pore. Input is produced the from the positional 
    frame data from a .pdb file. The multiframePDBReader function from the 
    pdbPoreReader file will produce the list from a .pdb file.
    
    @param frame_list List of list of AtomicPoints produced with multiframePDBReader 
    @return List of solid sphere coordinates and radii in the form of an AtomicPoint. 
        The largest solid sphere that can pass through pore of each given frame. [[AtomicPoint, [List of 3 bounding atoms]], ...]
    '''

    solutions = []
    n = 50 #number of atoms to get for successive frame calculations
    vp = 0.4 #Percent of zSlice VanDerWaals radius used when checking if point in the pore are valid
    
    frame_1 = frame_list.pop(0)
    solution_1 = [0,0,0]
    
    #Makes sure no 0 or -1 answers are passed
    while solution_1[0] <= 0:
        solution_1 = mcrs.MonteCarloSet(100, 1, frame_1, percentVDWR=0.3)
    solutions.append(solution_1)
    
    prev_solution = solution_1[1]
    prev_frame = frame_1
    
    #Generates graph for the frame - only need to do once, indexes are the same even if position isn't
    myGraph = ag.AtomicGraph()
    tpaa.genTPAACapGraph(myGraph, frame_1)
    
    #For C_R1 checking. Grabs indexes of r1 carbons
    r1list = []
    for v in myGraph.myVerticies:
        if v.ap.bondType == "C_R1":
            r1list.append(v.index)
    
    #Finds solution for each given frame
    for frame in frame_list:
        #Gets limited list of atoms for frame calculation. Alternative function is atomFinderDist
        frame_atoms = atomFinderNum(prev_solution, n, frame)
    
        sliceList = zsla.pathGen(frame) #generates sliceList

        for p in sliceList: #Adjusts vdwr appropriatly 
            p.VanDerWaalsRadius = p.VanDerWaalsRadius * vp
        
        #Calculates the solution
        frame_solution_set = pcf.poreCalcL(frame_atoms, frame, sliceList, myGraph, zbExclude=r1list, allTopAns=True)
        #frame_solution = frame_solution_set[1] #Gets atomic point of solution circle center
        #frame_solution.VanDerWaalsRadius = frame_solution_set[0] #set's solution sphere vdwr to calculated radius


        #Splits list into two based on some z value and returns section
        #with the lower maximum radius
        cut_list = pcf.clusterPointList(frame_solution_set, 34, "solution_lowMaxR")
            
        #Sorts the list in acending order by radius
        cut_list = pcf.sortList(cut_list, "solution_topRadius")
            
        #returns entry with max radius of that section
        frame_solution_set = cut_list[-1]


        #Runs Monte Carlo if there is no answer. Repeats until there is an answer
        while frame_solution_set[0] <= 0:
            frame_solution_set = mcrs.MonteCarloSet(100, 1, frame, percentVDWR=0.4)

        #Calculate position change of surrounding atoms between frames
        saPosChange = pcf.distChange(prev_frame, frame, frame_solution_set[3])
        #replcae indexes with atoms positional change
        frame_solution_set[3] = saPosChange

        ''' Need to add index tracking to verticies
        #Have option to calculate full change in postion of each cap between frames.
        for i in range(1,4):
            for v in myGraph.myVerticies:
                #Only gets a particular quadrant at a time
                if v.quadrant == i:
        '''

        solutions.append(frame_solution_set)
        prev_solution = frame_solution_set[1]
        prev_frame = frame
        
    return solutions