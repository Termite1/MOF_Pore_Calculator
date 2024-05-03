# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 20:00:41 2023

@author: samda
"""

from AtomicPoint import AtomicPoint
import PoreCalculatorV0Functions as pcf
import LargestCircleFinder as lcf
import random
import pdbPoreReader
import zSliceLineAproximation as zsla
import AtomicGraph as ag
import TPAAGraphFunctions as tpaa


'''
Redundant function due to having random.sample as a function. Leaving in for 
now just in case a file make use of it else where. Remove during next major
refactor.

@param n Number of points to sample
@param pointList List of AtomicPoints to sample from
@return List [] size n of Atomic Points
'''
def randomSamplingSet (n, pointList):
    samplelist = random.sample(pointList, n)
    return samplelist

'''
Function calculates the largest solid sphere that can fit through the pore described
by the pointList. Samples random sets of points to calculate from.

@param ss Sample size of random sampling
@param n Number of iteration of sampling
@param pointList List of AtomicPoints to sample from
@return List with [float largest circle Radius, largest Circle Center (x,y,z), 
                   and a list of the three points bounding the largest circle] 
                  of the minimum sphere generated over the iterations
'''
def MonteCarloSet (ss, n, pointList, **kwargs):
    #Number of points that will actually be used to calculate the largest sphere that
    #can fit through the pore, can change number in future.
    if 'numPoints' in kwargs:
        N = kwargs['numPoints']
    else: N = 100
    
    #Percent of zSlice VanDerWaals radius used when checking if point in the pore are valid
    if 'percentVDWR' in kwargs:
        vp = kwargs['percentVDWR']
    else: vp = 0.2
    
    #read points from file
    if 'fileName' in kwargs:
        pdbList = pdbPoreReader.atomReaderPDB(kwargs['fileName'])
    else:
        pdbList = pointList
    
    prebounds = lcf.boundsFinder(pdbList)
    #calculates center of bounded area, generates list of pore points from there
    #could probably come up with something more eligant or just use circle calc. Circle calc won't work with occluded z axis though
    center = AtomicPoint((prebounds[0] + prebounds[1]) / 2, (prebounds[2] + prebounds[3]) / 2, (prebounds[4] + prebounds[5]) / 2, "")
    
    #Allows manual setting of calculation center
    if 'center' in kwargs:
        center = kwargs['center']
    
    #generate limited list from pore file
    limitedList = pcf.listLimiter(center, pdbList, N)
    
    #generates sliceList
    sliceList = zsla.pathGen(pdbList)
    for p in sliceList:
        p.VanDerWaalsRadius = p.VanDerWaalsRadius * vp
    
    solutionsList = [None] * n
    i = 0 #index
    invalid = -1
    
    #Creates graph for checking what caps the calculated circles touch
    myGraph = ag.AtomicGraph()
    tpaa.genTPAACapGraph(myGraph, pdbList)
    
    #For C_R1 checking. Grabs indexes of r1 carbons
    r1list = []
    for v in myGraph.myVerticies:
        if v.ap.bondType == "C_R1":
            r1list.append(v.index)
    
    #Runs n times
    while i < n:
        sampleList = random.sample(limitedList, ss)
        
        if 'topAns' in kwargs:
            solutionsList[i] = pcf.poreCalcL(sampleList, pdbList, sliceList, myGraph, zbExclude=r1list, topAns=kwargs['topAns']) #iterates and calculates over the list
            return solutionsList[i] #topAns only deisnged to return answers from 1 frame at the moment
        if 'allTopAns' in kwargs:
            #answers = pcf.poreCalcL(sampleList, pdbList, sliceList, myGraph, zbExclude=r1list, allTopAns=True)
            answers = pcf.poreCalcL(sampleList, pdbList, sliceList, myGraph, allTopAns=True)
            
            
            #Splits list into two based on some z value and returns section
            #with the lower maximum radius
            cut_list = pcf.clusterPointList(answers, kwargs["allTopAns"], "solution_lowMaxR")
            
            #Sorts the list in acending order by radius
            cut_list = pcf.sortList(cut_list, "solution_topRadius")
            
            #returns entry with max radius of that section
            solutionsList[i] = cut_list[-1]
            
        else: solutionsList[i] = pcf.poreCalcL(sampleList, pdbList, sliceList, myGraph, zbExclude=r1list) #iterates and calculates over the list
        
        #print(f"Radius = {solutionsList[i][0]:.3f} at ({solutionsList[i][1].x:.3f}, {solutionsList[i][1].y:.3f}, {solutionsList[i][1].z:.3f})")
        if solutionsList[i][0] <= 0:
            solutionsList[i][0] = invalid #makes sure there are no values equal to or below 0.
            #Radius that are are set invalidly high
        i += 1 #increments index
    
    radiusList = [None] * n
    i = 0
    
    for a in solutionsList:
        radiusList[i] = solutionsList[i][0]
        i += 1
    
    maxrad = max(radiusList) #gets smallest solution calculated
    solution = None
    
    for a in solutionsList:
        if maxrad == a[0]:
            solution = a
    
    return solution
