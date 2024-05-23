# -*- coding: utf-8 -*-
"""
Created on Wed May 22 13:13:19 2024

Example of key code used to calculate that largest solid sphere that can pass 
through a MOF's pores

@author: samda
"""
import math
from matplotlib import path
from AtomicPoint import AtomicPoint
import heapq


def poreCalculator(pointList, center, n, z):
    '''
    Function calculates the largest solid sphere that can fit through the pore described
    by the pointList.

    Parameters
    ----------
    pointList : List of AtomicPoints
        Pore being explored
    center : AtomicPoint
        Location from which nearest AtomicPoints are found.
    n : int
        Number of AtomicPoints to select for calculation
    z : float
        z Height partition used for clusterPointList function.

    Returns
    -------
    solution_set : [float largest circle Radius, largest Circle Center (x,y,z), 
                       a list of the three points bounding the largest circle,
                       and the indexes of the bounding points from the originList]
        Information about best solution candidate sphere.
    '''
    #generate limited list from pore file
    limitedList = listLimiter(center, pointList, n)
    
    #calculate all answers
    allAns = poreCalcL(limitedList, pointList)
    
    #select best answer
    cut_list = clusterPointList(allAns, z)
    cut_list = sortList(cut_list)
    solution_set = cut_list[-1]
    
    return solution_set


def listLimiter(a, pointList, size):
    '''
    Function takes a list of AtomicPoints and returns a new list of length size containing 
    the size closest points to AtomicPoint a

    @param a Atomic Point
    @param pointList List of points to be checked
    @param size number of points to be returned
    @return List of size Atomic Points nearest to point a from the given pointList
    '''
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


def boundsFinder(pointArray):
    '''
    Function takes a pointArray and produces the maxmum and minimum x, y, z values
    in an array in the form of [Xmax, Xmin, Ymax, Ymin, Zmax, Zmin]

    @param pointArray Array of points with x, y, z components
    @return Array of floats [Xmax, Xmin, Ymax, Ymin, Zmax, Zmin]
    '''
    
    #Bounding Checker Portion of the Function. This could be performed outside
    #of the circleGrow function and added as an additional input
    Xmax = pointArray[0].x 
    Xmin = pointArray[1].x 
    Ymax = pointArray[0].y 
    Ymin = pointArray[1].y 
    Zmax = pointArray[0].z 
    Zmin = pointArray[1].z 
    
    for a in pointArray: 
        if a.x > Xmax: 
            Xmax = a.x 
        elif a.x < Xmin: 
            Xmin = a.x 
        if a.y > Ymax: 
            Ymax = a.y 
        elif a.y < Ymin: 
            Ymin = a.y 
        if a.z > Zmax: 
            Zmax = a.z 
        elif a.z < Zmin: 
            Zmin = a.z 
    
    return [Xmax, Xmin, Ymax, Ymin, Zmax, Zmin]


def poreCalcL(pointList, originList):
    '''Function tests every permutation of 3 points from the limited point list for valid candidate sphere results

    @param pointList to iterate over
    @param originList list of points defining original system
    @return List with [float largest circle Radius, largest Circle Center (x,y,z), 
                       a list of the three points bounding the largest circle,
                       and the indexes of the bounding points from the originList]
    '''
    lcc = AtomicPoint(0, 0, 0, '') #largest circle center
    lcr = 0                        #largest circle radius
    lcbp = [AtomicPoint(0, 0, 0, '')] * 3 #largest circle bounding points
    bpi = [0, 0, 0]
    
    #For returning top solutions
    topAns = []
        
    #runs through every 3 atom combination in the point list
    for a in range(0, len(pointList)):
        for b in range(a+1, len(pointList)):
            if b == a: continue
            for c in range(b+1, len(pointList)):
                #Calculates if the 3 atoms produce a condidate sphere
                val = circleCalculation(pointList[a], pointList[b], pointList[c], originList)
                
                #If candidate sphere valid, adds it to possible answers list
                if val[0] > 0:
                    topAns.append(val)
                    
    #Checks to make sure there is an answer. If not, dummy value returned 
    #just in case. Inside list since topAns is a list format
    if len(topAns) == 0: return [[lcr, lcc, lcbp, bpi]]
    return topAns


def circleCalculation(a, b, c, pointList):
    '''
    Function takes 3 atoms, finds the sphere they define, then checks if that sphere
    is valid in the context. If the sphere is valid, its center and radius are returned. 
    Otherwise an invalide result will be returned (radius < 0)

    @param a Atomic Point a
    @param b Atomic Point b
    @param c Atomic Point c
    @param pointList List of points in the system the circle is being calculated from
    @return List with [float Circle Radius, Circle Center (x,y,z), a list of the 
                       three atoms bounding the sphere, and a list of the bounding atom's indexes from the pointList]
    '''
    center = threeAtomCircumcenter(a, b, c)
    
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
    plane = planarize(a, b, c)
    
    #Converts atomic and candidate sphere coordinates to coordinates of new plane
    pa = planarizePoint(a, plane)
    pb = planarizePoint(b, plane)
    pc = planarizePoint(c, plane)
    pcenter = planarizePoint(center, plane)

    #Checks if the candidate is surrounded, returns invalid if not
    if simpleCircleSurrounded(pa, pb, pc, pcenter) == False:
        return [-1, center]
    
    
    #Grabs list index of surounding atoms
    i1 = pointList.index(a)
    i2 = pointList.index(b)
    i3 = pointList.index(c)
    
    return [r, center, [a,b,c], [i1, i2, i3]]


def clusterPointList(pointList, z):
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

    Returns
    -------
    List of solutions (see circle calculation return)
    '''
    #Initializes lists
    topList = []
    botList = []
    
    #runs through every possible solution
    for i in pointList:
        if i[1].z > z: topList.append(i)
        else: botList.append(i)
    
    #Checks if both lists have entries. If one doesn't, the other is returned
    #by default    
    if len(topList) == 0: return botList
    if len(botList) == 0: return topList
    
    #Returns the solution list with lowest maximum radius

    #top and bottom max radius for seperating based on operation "solution_topRadius"
    trMax = sortList(topList)[-1][0]
    brMax = sortList(botList)[-1][0]
        
    #Returns list with lower max radius
    if trMax > brMax: return botList
    else: return topList


def sortList(pointList):
    '''
    Sorts given list of AtomicPoints by some method. General methods will be
    based on insertion sort for now. Code derived from https://www.geeksforgeeks.org/sorting-algorithms-in-python/#

    Parameters
    ----------
    pointList : List of AtomicPoints
        List of AtomicPoints to sort.
        
    Returns
    -------
    Sorted list of AtomicPoints or list of solutions.
    '''
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


def threeAtomCircumcenter(p1, p2, p3):
    """
    Function for calculating the center of a circle described by how it touches three non-linear atoms.
    
    Explinations behind equations used in three_atoms_touching.pdf
    
    Mathmactical formula for the coordinates of a circle center based on https://math.stackexchange.com/questions/3100828/calculate-the-circle-that-touches-three-other-circles
    
    @param p1 AtomicPoint 1
    @param p2 AtomicPoint 2
    @param p3 AtomicPoint 3
    @return AtomicPoint describing the center of the circle defined by the 3 atoms
    """
    r1 = p1.VanDerWaalsRadius
    r2 = p2.VanDerWaalsRadius
    r3 = p3.VanDerWaalsRadius
    
    a = (p2.y - p1.y) * (p3.z - p1.z) - (p2.z - p1.z) * (p3.y - p1.y)
    b = (p2.z - p1.z) * (p3.x - p1.x) - (p2.x - p1.x) * (p3.z - p1.z)
    c = (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x)
    d =  a * p1.x + b * p1.y + c * p1.z
    
    l = p2.x - p1.x + (a/c) * p1.z - (a/c) * p2.z
    m = p2.y - p1.y + (b/c) * p1.z - (b/c) * p2.z
    n = r1 - r2
    p_ = r1**2 - r2**2 + p2.x**2 - p1.x**2 + p2.y**2 - p1.y**2 + p2.z**2 - p1.z**2
    p = 0.5*p_ + (d/c) * p1.z - (d/c) * p2.z

    q = p3.x - p1.x + (a/c) * p1.z - (a/c) * p3.z
    s = p3.y - p1.y + (b/c) * p1.z - (b/c) * p3.z
    t = r1 - r3
    w_ = r1**2 - r3**2 + p3.x**2 - p1.x**2 + p3.y**2 - p1.y**2 + p3.z**2 - p1.z**2
    w = 0.5*w_ + (d/c) * p1.z - (d/c) * p3.z
    
    A = (s*n - m*t) / (s*l - m*q)
    B_ = (s*p - m*w) / (s*l - m*q)
    B = B_ - p1.x
    C = (q*n - l*t) / (q*m - l*s)
    D_ = (q*p - l*w) / (q*m - l*s)
    D = D_ - p1.y
    E = (-b/c) * C - (a/c) * A
    F = (d/c) - (b/c) * D_ - (a/c) * B_ - p1.z
    
    if ((A*B + C*D + E*F - r1)**2 - (A**2 + C**2 + E**2 -1) * (B**2 + D**2 + F**2 - r1**2)) < 0:
        #print("Error: Negative value calculated in uncertainty square root.")
        #Error return
        return AtomicPoint(0, 0, 0, '', vdwr=-1)
    
    #Getting a negative of a sqrt
    value = -(A*B + C*D + E*F - r1) / (A**2 + C**2 + E**2 - 1)
    uncertainty = math.sqrt((A*B + C*D + E*F - r1)**2 - (A**2 + C**2 + E**2 -1) * (B**2 + D**2 + F**2 - r1**2)) / (A**2 + C**2 + E**2 - 1)

    r_ans1 = value + uncertainty
    r_ans2 = value - uncertainty

    x_1 = A*r_ans1 + B_
    x_2 = A*r_ans2 + B_

    y_1 = C*r_ans1 + D_
    y_2 = C*r_ans2 + D_

    z_1 = (d/c) - (b/c) * y_1 - (a/c) * x_1
    z_2 = (d/c) - (b/c) * y_2 - (a/c) * x_2
    
    #Solution 1
    sol_1 = AtomicPoint(x_1, y_1, z_1, "") 
    sol_1.VanDerWaalsRadius = abs(r_ans1)
    #Checks if the circle touches the surface of each atom, eliminate answers 
    #from linear or near-linear combination of atoms
    rCheck1 = threeAtomRadiusCheck(p1, p2, p3, sol_1) 
    
    #solution 2
    sol_2 = AtomicPoint(x_2, y_2, z_2, "") 
    sol_2.VanDerWaalsRadius = abs(r_ans2)
    #Checks if the circle touches the surface of each atom, eliminate answers 
    #from linear or near-linear combination of atoms
    rCheck2 = threeAtomRadiusCheck(p1, p2, p3, sol_2)


    #Ensure that the surrounded circle is returned rather than the surrounding circle
    if rCheck1 and rCheck2: 
        if abs(r_ans1) < abs(r_ans2):
            return sol_1
        else:
            return sol_2
    elif rCheck1 and rCheck2 == False:
        return sol_1
    elif rCheck1 == False and rCheck2:
        return sol_2
    else:
        if abs(r_ans1) < abs(r_ans2):
            return sol_1
        else:
            return sol_2
        
        
def threeAtomRadiusCheck(p1, p2, p3, sphere, **kwargs):
    '''Function takes three atoms and a sphere and checks if the sphere touches
    the surface of each of the atoms.
    
    @param p1 AtomicPoint 1
    @param p2 AtomicPoint 2
    @param p3 AtomicPoint 3
    @peram sphere sphere being checked described by an AtomicPoint
    @return True if sphere touches the surface of each of the atoms, false otherwise
    '''
    #Only runs if 2D o
    if 'circle' in kwargs:
        if kwargs['circle'] == True:
            if abs(sphere.distBetween2D(p1, True)) < 0.01:
                if abs(sphere.distBetween2D(p2, True)) < 0.01:
                    if abs(sphere.distBetween2D(p3, True)) < 0.01:
                        return True
    
    else:
        if abs(sphere.distBetween(p1, True)) < 0.001:
            if abs(sphere.distBetween(p2, True)) < 0.001:
                if abs(sphere.distBetween(p3, True)) < 0.001:
                    return True
    return False


def simpleCircleSurrounded(p1, p2, p3, c):
    """
    Function checks if a circle c is surrounded by points p1, p2, p3. A surrounded circle 
    would be unable to grow. Returns True if c is surrounded, False otherwise.

    @param p1 AtomicPoint of bounding circle 1
    @param p2 AtomicPoint of bounding circle 2
    @param p3 AtomicPoint of bounding circle 3
    @param c AtomicPoint of the center of circle c
    @return Boolean. True is surrounded, False otherwise
    """
    
    p = path.Path([(p1.x,p1.y),(p2.x,p2.y),(p3.x,p3.y)])
    point = [c.x,c.y]
    
    if p.contains_point(point):
        return True
    
    return False


def planarize(a, b, c):
    '''
    Function takes 3 points with x, y, y components and finds the 2d plane that intersects them all.
    Returns the principle values to convert (x, y, z) points into points in the new plane. Updated
    function does not use matric calculations.
    
    Based on this: https://www.quora.com/What-is-the-fastest-way-to-find-the-equation-of-a-plane-given-three-points
    
    @param a Point 1 with (x, y, z) components
    @param b Point 2 with (x, y, z) components
    @param c Point 3 with (x, y, z) components
    @return [Plane Origin as AtomicPoint, Surface Normal Unit Vector, Axis Unit Vector 1, Othoganal Axis Unit Vector 2, 
                 Constant for cartesian planar equation with regards to surface normal k]
    '''
    ab = [b.x-a.x,b.y-a.y,b.z-a.z]  #vector from a --> b
    ac = [c.x-a.x,c.y-a.y,c.z-a.z]  #vector from a --> c
    
    n_ab = normalize(ab) #normaized ab vector
    n_ac = normalize(ac) #normaized ac vector
    
    #Cross product vector between ab and ac
    surface_normal = [n_ab[1]*n_ac[2] - n_ab[2]*n_ac[1], -1*(n_ab[0]*n_ac[2] - n_ab[2]*n_ac[0]), n_ab[0]*n_ac[1] - n_ab[1]*n_ac[0]]
    n_surface_normal = normalize(surface_normal)
    
    #Second planar axis along with ab, cross product of surface normal and ab
    in_plane_normal = [n_surface_normal[1]*n_ab[2] - n_surface_normal[2]*n_ab[1], -1*(n_surface_normal[0]*n_ab[2] - n_surface_normal[2]*n_ab[0]), n_surface_normal[0]*n_ab[1] - n_surface_normal[1]*n_ab[0]]
    n_in_plane_normal = normalize(in_plane_normal)
    
    #Planar constant --> u1 ∗ x + u2 ∗ y + u3 ∗ z + K = 0
    k = -1 * (a.x*surface_normal[0] + a.y*surface_normal[1] + a.z*surface_normal[2])
    
    #Declares planar origin = point a
    origin = AtomicPoint(a.x, a.y, a.z, "")
    
    return [origin, n_surface_normal, n_ab, n_in_plane_normal, k]


def planarizePoint(a, pv):
    '''
    Function takes a (x,y,z) point and a set of values produced by planarize, and projects that point
    onto the given plane, returning its (u,w) position. 

    Code from here: https://stackoverflow.com/questions/23472048/projecting-3d-points-to-2d-plane

    @param a Point with (x,y,z) values to be projected onto given plane
    @param pv Principle planar values calculated by planarize. [Plane Origin as AtomicPoint, Surface Normal Unit Vector, Axis Unit Vector 1, 
                                                                Othoganal Axis Unit Vector 2, Constant for cartesian planar equation with regards to surface normal k]
    @return AtomicPoint with (u,w) position as projected onto the given plane. s will be set to the out of plane seperation value (not sure if sign will be correct) (x,y,z) --> (u,v,s)
    '''
    origin = pv[0]
    surface_normal = pv[1]
    axis1 = pv[2]
    axis2 = pv[3]
    
    origin_vector = [a.x - origin.x, a.y - origin.y, a.z - origin.z]
    
    s = dot_product(surface_normal, origin_vector)
    u = dot_product(axis1, origin_vector)
    v = dot_product(axis2, origin_vector)

    solution = AtomicPoint(u, v, s, "", vdwr=a.VanDerWaalsRadius)
    return solution   


def normalize(v):
    '''
    Function normailzes a given vector and returns it.
    
    @peram v Vector in list form
    @peram Unit vector in direction of v in the form of a list
    '''
    magnitude = 0
    #calculates magnitude of vector from components, runs for all components
    for component in v:
        magnitude += (float(component)**2)
    final_magnitude = math.sqrt(magnitude)
    #Devides each component of the vector by the vector's magnitude
    normalized_vector = [x/final_magnitude for x in v]
    return normalized_vector


def dot_product(a, b):
    '''
    Function calculates the dot product of two vectors (lists) or equal length
    
    @peram a Vector 1
    @peram b Vector 2
    @return Dot product as list
    '''
    total = 0
    for i in range (len(a)):
        total += (a[i] * b[i])
    return total