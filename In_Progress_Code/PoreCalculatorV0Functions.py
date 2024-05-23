# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 11:15:03 2023

This file contains all the main functions needed to find hte largest solid 
sphere that can pass through a given pore described by a pdb file. Major function
include:
    The CircleCalculation function that calculates a sphere based on 3
    given AtomicPoints and performs a number of checks to confirm that the point
    generated is valid.
    
    poreCalcF that generates the result from a pdb file
    
    poreCalcL that generates the result from a number of given perameters

@author: samda
"""

import LargestCircleFinder as lcf
import lcfProcessing as proc
from AtomicPoint import AtomicPoint
import heapq
import pdbPoreReader
import zSliceLineAproximation as zsla
import AtomicGraph as ag
import TPAAGraphFunctions as tpaa


def circleCalculation(a, b, c, pointList, sliceList, capGraph, **kwargs):
    '''
    Function takes 3 atoms, finds the sphere they define, then checks if that sphere
    is valid in the context. If the sphere is valid, its center and radius are returned. 
    Otherwise an invalide result will be returned (radius < 0)

    @param a Atomic Point a
    @param b Atomic Point b
    @param c Atomic Point c
    @param pointList List of points in the system the circle is being calculated from
    @param sliceList list of points approximating line traveling through the center of the pore
    @param capGraph Graph of the AtomicPoints comprising the caps on the pore
    @peram **kwargs
        exclude_C_R1 List of ints
            Used for C_R1 exclusion. Solutions that are below 2 or more atoms 
            within the given list of indexs are considered invalid
    @return List with [float Circle Radius, Circle Center (x,y,z), a list of the 
                       three atoms bounding the sphere, and a list of the bounding atom's indexes from the pointList]
    '''
    center = lcf.threeAtomCircumcenter(a, b, c)
    
    r = center.VanDerWaalsRadius
    
    #Check if candidate sphere radius is large enough. Haven't seen an answer yet with a vdwr < 1.5, keeping
    #limit at 1 for now
    if r <= 1: return [-7, center]
    
    #Trying a euclid distance check. So far no answers have been more than 2.75
    #xy dist from the center of the pore (0,0). Flat limiting at 3 for now, not sure if totally applicable,
    #could see missing outlier
    if center.distBetween2D(AtomicPoint(0, 0, 0, ""), False) > 3:
        return [-8, center]
    
    #Checks if candidate sphere intersects other atoms, returns invalid result if so
    for p in pointList:
        if p.distBetween(center, True) < 0 and p != a and p != b and p != c:
            return [-3, center]
    
    #Takes three given atoms and finds the values for the defined plane
    plane  = proc.planarize(a, b, c)
    
    #Converts atomic and candidate sphere coordinates to coordinates of new plane
    pa = proc.planarizePoint(a, plane)
    pb = proc.planarizePoint(b, plane)
    pc = proc.planarizePoint(c, plane)
    pcenter = proc.planarizePoint(center, plane)

    #Checks if the candidate is surrounded, returns invalid if not
    if lcf.simpleCircleSurrounded(pa, pb, pc, pcenter) == False:
        return [-1, center]
    

    #checks if sphere intersects with the line traveling through the center of the pore at least once
    pathCheck = False
    
    #Runs for every slice in z made with the pathGen function from zsla
    for z_slice in sliceList:
        #adjusts for z location
        tp = zsla.zSlice(center, z_slice.z)
        
        #Only runs if center circle cuts through relevant z plane
        if tp.VanDerWaalsRadius != 0:
            #checks if the circle from center intersects with the sliceList circle
            if z_slice.distBetween2D(tp, True) < 0: 
                pathCheck = True
            
    if pathCheck == False:
        return [-4, center]
    
    
    #Grabs list index of surounding atoms
    i1 = pointList.index(a)
    i2 = pointList.index(b)
    i3 = pointList.index(c)

    #Temporarily removing quadratn check until it the logic for it is fixed

    '''
    #checks what quadrants each of the atoms defining the circle belong to. 
    #Spheres are only valid if they are made from atoms from three different
    #caps or caps that are arranged diaganol to each other.
    atomList = [a, b, c]
    quadrantList = [0, 0, 0]
    i = 0
    
    quadrantCheck = False
    #repeats 3 times, 1 for each atom
   
    #Could possibly simplify by doing the graph first each frame then setting quadrents right there
    #reduces repeats of the next step! --> Do this!
    
    #Don't want quadrants to be apart of atomic points
    for p in atomList:
        for v in capGraph.myVerticies:
            #records the quadrants of the atomic points
            if v.ap == p:
                quadrantList[i] = v.quadrant

    q_bools = [False] * 4
    for atom in quadrantList:
        q_bools[atom - 1] = True
    
    #Checks if the sphere is valid
    #Sphere is valid if atoms come from 3 different caps --> 3 diff caps automatically have 1 diagonal pair
    #sphere is valid if atoms come from diagonal caps
    if q_bools[0] and q_bools[3]: quadrantCheck = True
    elif q_bools[1] and q_bools[2]: quadrantCheck = True
    
    if quadrantCheck == False:
        [-5, center]
    '''
    
    rc = 0
    #Check solution where answers are below a C_R1 carbon (or other atoms as indicated by the given indexes)
    if "exclude_C_R1" in kwargs:
        rc = 0 #R1 count
        r1list = kwargs["exclude_C_R1"]
        for carbon in r1list:
            #make temp point after z check
            temp_point = AtomicPoint(center.x, center.y, 0, "") 
            #only counts when solution is below caps
            if pointList[carbon].z > center.z:
                if temp_point.distBetween2D(pointList[carbon], True) < r:
                    rc += 1
    #Only exclude answer if center below 2 or more C_R1 carbons            
    if rc >= 2:
        return [-6, center]
    
    
    return [r, center, [a,b,c], [i1, i2, i3]]


'''
Function takes a list of AtomicPoints and returns a new list of length size containing 
the size closest points to AtomicPoint a

@param a Atomic Point
@param pointList List of points to be checked
@param size number of points to be returned
@return List of size Atomic Points nearest to point a from the given pointList
'''
def listLimiter(a, pointList, size):
    distList = [None] * len(pointList)
    distDict = {}
    limitedList = [None] * size
    
    i = 0; #list index
    
    #For loop runs through each itme in list, calulates distance to given point
    #records that point in a new list, and records the distance:AtomicPoint pair 
    #with the distance as the key in a dictionary
    for p in pointList:
        d = a.distBetween(p, True)
        
        #update dictionary
        add = {d : p}
        distDict.update(add)
        
        distList[i] = d
        i = i + 1
        
    #Finds the produces a size length list of the smallest distances in the distance list    
    valList = heapq.nsmallest(size, distList)
    
    i = 0
    #Fills limited list with the AtomicPoints associate with the size smallest distances.
    for v in valList:
        limitedList[i] = distDict.pop(v)
        i = i + 1
        
    return limitedList


def listLimiter2D(a, pointList, n, **kwargs):
    """
    Function takes a list of AtomicPoints and returns a new list of length n 
    containing the n closest points to AtomicPoint a. Only accounts for x and y
    axes, ignoring z. If kwargs has setZ, can set a z for the produced atoms in 
    the list.

    Parameters
    ----------
    a : AtomicPoint
        Point distance calculations are performed from.
    pointList : List of AtomicPoints
        List of AtomicPoints to be checked.
    n : int
        Number of nearest atoms to be returned.
    **kwargs : 
        setZ : Float
            Sets z value of returned atoms in list to the setZ value.

    Returns
    -------
    limitedList : List of AtomicPoints
        List of n nearest atoms to AtomicPoint a.

    """
    distList = [None] * len(pointList)
    distDict = {}
    limitedList = [None] * n
    
    #Tracks if setting z or not
    if 'setZ' in kwargs:
        zSet = True
    else: zSet = False
    
    #makes sure n isn't larger than the given list
    if n >= len(pointList):
        if zSet:
            for j in range(0, len(pointList)):
                tempAtom = pointList[j]
                tempAtom.z = kwargs['setZ']
                distList[j] = tempAtom
            return distList
        else: return pointList
    
    i = 0; #list index
    
    #For loop runs through each itme in list, calulates distance to given point
    #records that point in a new list, and records the distance:AtomicPoint pair 
    #with the distance as the key in a dictionary
    for p in pointList:
        d = a.distBetween2D(p, True)
        
        #update dictionary
        add = {d : p}
        distDict.update(add)
        
        distList[i] = d
        i = i + 1
        
    #Finds the produces a size length list of the smallest distances in the distance list    
    valList = heapq.nsmallest(n, distList)
    
    i = 0
    #Fills limited list with the AtomicPoints associate with the size smallest distances.
    for v in valList:
        tempAtom = distDict.pop(v)
        
        if zSet: tempAtom.z = kwargs['setZ']
        
        limitedList[i] = tempAtom
        i = i + 1
        
    return limitedList


def insertPoint(p, l):
    '''
    Function updates list by inserting AtomicPoints into the list and keeping them
    organized by van der Waals radius in decending order. Doens't increase size of
    list, so items may be lost if they are too small.

    Parameters
    ----------
    p : AtomicPoint
        AtomicPoint being inserted into the list.
    l : List of AtomicPoints
        List of AtomicPoints new point is being inserted into.

    Returns
    -------
    None. The list is just updated.

    '''
    a = p
    
    #Run through each item in the given list
    for i in range(len(l)):
        b = l[i]
        #If the selected AP's vdwr is larger the one being looked at, the points
        #swapped. 
        
        #Helps ensure answers aren't duplicated in the list
        if f"{a[0]:.3f}" == f"{b[0]:.3f}": break
        
        #performs check
        if a[0] - b[0] > 0.001:
            l[i] = a
            a = b


def sortList(pointList, method):
    '''
    Sorts given list of AtomicPoints by some method. General methods will be
    based on insertion sort for now. Code derived from https://www.geeksforgeeks.org/sorting-algorithms-in-python/#

    Parameters
    ----------
    pointList : List of AtomicPoints
        List of AtomicPoints to sort.
    method : String
        Describes sorting method.
            "topRadius" : Sorts from largest to smallest VanDerWaalsRadius 
            "solution_topRadius" : Sorts solution list from largest to smallest
                calculated radius
        
    Returns
    -------
    Sorted list of AtomicPoints or list of solutions.
    '''
    
    if method == "topRadius":        
        
        # Outer loop to traverse on len(pointList)  
        for i in range (1, len(pointList)):
            a = pointList[i]
            
            # Move elements of list1[0 to i-1], 
            # which are greater to one position
            # ahead of their current position 
            j = i - 1
            
            while j >= 0 and a.VanDerWaalsRadius < pointList[j].VanDerWaalsRadius:
                pointList[j + 1] = pointList[j]
                j -= 1
                
            pointList[j + 1] = a
            
        return pointList
    
    if method == "solution_topRadius":        
        
        # Outer loop to traverse on len(pointList)  
        for i in range (1, len(pointList)):
            a = pointList[i]
            
            # Move elements of list1[0 to i-1], 
            # which are greater to one position
            # ahead of their current position 
            j = i - 1
            
            while j >= 0 and a[0] < pointList[j][0]:
                pointList[j + 1] = pointList[j]
                j -= 1
                
            pointList[j + 1] = a
            
        return pointList
    

#may want to return to 10 largest spheres to check?
def poreCalcL(pointList, originList, sliceList, myGraph, **kwargs):
    '''Function tests every permutation of 3 points from the limited point list for valid maximum sphere results

    @param pointList to iterate over
    @param originList list of points defining original system
    @param sliceList list of points approximating line traveling through the center of the pore
    @peram **kwargs
        zbExclude list of ints. Used in circleCalc, answers below atoms with given indexes considered invalid
        topAns int. Returns top int answers
        allTopAns None. Returns all top answers if present in kwargs. topAns
            has priority
    @return List with [float largest circle Radius, largest Circle Center (x,y,z), 
                       a list of the three points bounding the largest circle,
                       and the indexes of the bounding points from the originList]
    '''
    
    lcc = AtomicPoint(0, 0, 0, '') #largest circle center
    lcr = 0                        #largest circle radius
    lcbp = [AtomicPoint(0, 0, 0, '')] * 3 #largest circle bounding points
    bpi = [0, 0, 0]
    
    #For returning top topAns solution.
    topAns = []
    if 'topAns' in kwargs:
        topAns = [None] * kwargs['topAns']
        #Fills list with dummy van der waals radii
        for i in range(len(topAns)):
            topAns[i] = [0, AtomicPoint(0, 0, 0, '')]
    
    #do multithreading instead here eventually
    for a in range(0, len(pointList)):
        for b in range(a+1, len(pointList)):
            if b == a: continue
            for c in range(b+1, len(pointList)):
                val = None
                
                if 'zbExclude' in kwargs:
                    val = circleCalculation(pointList[a], pointList[b], pointList[c], originList, sliceList, myGraph, exclude_C_R1=kwargs['zbExclude'])
                else:
                    val = circleCalculation(pointList[a], pointList[b], pointList[c], originList, sliceList, myGraph)
                
                if 'topAns' in kwargs:
                    insertPoint(val, topAns)
                elif val[0] > lcr:
                    #redundant, could just save val
                    lcr = val[0] #largest circle radius
                    lcc = val[1] #largest circle center
                    lcbp = val[2] #largest circle bounding points
                    bpi = val[3] #bounding points indexes
                    
                if 'allTopAns' in kwargs:
                    if val[0] > 0:
                        topAns.append(val)

    if 'topAns' in kwargs: return topAns
    elif 'allTopAns' in kwargs: 
        #Check to make sure there is an answer. If not, dummy value returned 
        #just in case. Inside list since topAns is a list format
        if len(topAns) == 0: return [[lcr, lcc, lcbp, bpi]]
        return topAns
    else: return [lcr, lcc, lcbp, bpi]


def distChange(f1, f2, indexList):
    '''
    Function return the euclidiean change in position for a list of atoms 
    between one frame and the next, provided the atoms have the same index in 
    each frame.

    Parameters
    ----------
    f1 : List of AtomicPoints
        First frame.
    f2 : List of AtomicPoints
        Second frame.
    indexList : List of ints
        Indicies of atoms being examined.

    Returns
    -------
    ???: List of floats
        Eucidean change in position for the related atom.

    '''
    
    #Holds list of position changes
    posChangeList = []
    
    #Repeats for each index in the index list
    for i in indexList:
        a = f1[i] #Grabs atom's postion in first frame
        b = f2[i] #Grabs atom's postion in second frame
        
        #Calculates distance between them. Ignored van der Waal's radius
        d = a.distBetween(b, False)
        #Adds distance to list
        posChangeList.append(d)
    
    return posChangeList


def clusterPointList(pointList, z, c):
    '''
    Function takes a list of AtomicPoints or list of solutions and splits it 
    based on some z value. It then returns one of the list based on some 
    operation deterimned by string perameter c.
    
    In theory it wouldn't be difficult to account for splitting based on values
    other than z.

    Parameters
    ----------
    pointList : List of AtomicPoints or list of solutions
        List of AtomicPoints to be split.
    z : float
        z value to split list at.
    c : String
        Operation used to determine how to split list
            "lowMaxR": Returns the list with lowest maximum radius
            "solution_lowMaxR": Returns the solution list with lowest maximum 
                radius

    Returns
    -------
    List of AtomicPoints or list of solutions
    '''
    #Initializes lists
    topList = []
    botList = []
    
    if c == "lowMaxR":
        for i in pointList:
            if i.z > z: topList.append(i)
            else: botList.append(i)
    if c == "solution_lowMaxR":
        for i in pointList:
            if i[1].z > z: topList.append(i)
            else: botList.append(i)
    
    #Checks if both lists have entries. If one doesn't, the other is returned
    #by default    
    if len(topList) == 0: return botList
    if len(botList) == 0: return topList
    
    #Returns the list with lowest maximum radius
    if c == "lowMaxR":
        #top and bottom max radius for seperating based on operation "lowMaxR"
        trMax = 0
        brMax = 0
        
        #Gets max top radius - Probably won't work, fix with something similar to what was done below
        for a in topList:
            if a.VanDerWaalsRadius > trMax: trMax = a.VanDerWaalsRadius
        
        #Gets max bottom radius
        for a in botList:
            if a.VanDerWaalsRadius > brMax: trMax = a.VanDerWaalsRadius
    
        #Returns list with lower max radius
        if trMax > brMax: return botList
        else: return topList
        
    #Returns the solution list with lowest maximum radius
    if c == "solution_lowMaxR":
        #top and bottom max radius for seperating based on operation "lowMaxR"
        trMax = sortList(topList, "solution_topRadius")[-1][0]
        brMax = sortList(botList, "solution_topRadius")[-1][0]
        
        #Returns list with lower max radius
        if trMax > brMax: return botList
        else: return topList