# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 18:21:26 2023

The purpose of program is to find the largest circle that can fit in a given
perimeter defined by an array of (x,y) coordinate values. 

@author: samda
"""

import math
from matplotlib import path

from AtomicPoint import AtomicPoint


def sign(n):
    """
    Returns the sign of a given number, or 0 if it is equal to zero

    @param n A number
    @return 1 if n > 0, -1 if n < 0, or 0 if n = 0
    """
    
    if n == 0: return 0
    elif n > 0: return 1
    elif n < 0: return -1
    

def cosLawAngle(a, b, c):
    """
    Calcuates the angle opening onto side c of trignle with sides a, b, c using the
    Cosine Law of triangles

    @param a Float side of a triangle
    @param b Float side of a triangle
    @param c Float side of a triangle
    @return Float theta representing the angle opening onto side c
    """
    
    return math.acos((c**2 - a**2 - b**2)/(-2 * a * b))


def twoAtomCircle(p1, p2, v):
    """
    Takes two atomic points with (x,y) coordinates and a vector from one point toward
    the direction the center of circle c is moving in (v). Acounts for Van der Waal
    radii. Returns the new location of circle center c. Equation documentation is in
    the Samuel Darer OneNote on the pore algorithm.

    @param p1 Atomic point 1, the point c was originally growing away from
    @param p2 Atomic point 2, the point circle c just collided with
    @param v Atomic point of c's last location, a vector in the direction c was traveling
    @return AtomicPoint Point c with its new location
    """
    
    #make p1 the center of my new coordiate system
    tp1 = AtomicPoint(0, 0, 0, p1.element)
    tp1.VanDerWaalsRadius = p1.VanDerWaalsRadius
    
    tp2 = AtomicPoint(p2.x - p1.x, p2.y - p1.y, 0, p2.element)
    tp2.VanDerWaalsRadius = p2.VanDerWaalsRadius
    
    tv = AtomicPoint(v.x - p1.x, v.y - p1.y, 0, "")
    
    #Added to increase readability
    vdwr1 = tp1.VanDerWaalsRadius
    vdwr2 = tp2.VanDerWaalsRadius
    
    h = tp1.distBetween2D(tp2, False)
    
    o_theta = math.acos(tp2.x / math.sqrt( tp2.x**2 + tp2.y**2 ) )
    theta = math.acos(tv.x / math.sqrt( tv.x**2 + tv.y**2 ) )
    
    r = (1/2) * (vdwr1**2 - vdwr2**2 + h**2 - 2 * vdwr1 * h * math.cos(theta - o_theta)) / (vdwr2 - vdwr1 + h * math.cos(theta - o_theta))
    
    a = r + vdwr1
    
    #Find new point, and account for the fact that the calculation was centered on p1
    c = AtomicPoint(a * math.cos(sign(tv.x) * theta) + p1.x, a * math.sin(sign(tv.y) * theta) + p1.y, 0, "")
    
    return c



def threeCircleCircumcenter(p1, p2, p3):
    """
    Given 3 AtomicPoints each defining the center of a circle, finds the location
    representing the center of the circle the points define

    Mathmactical formula for the coordinates of a circle center was aquired from https://math.stackexchange.com/questions/3100828/calculate-the-circle-that-touches-three-other-circles

    @param a AtomicPoint 1
    @param b AtomicPoint 2
    @param c AtomicPoint 3
    @return AtomicPoint describing the center of the circle defined by the 3 points. z value of returned point = p1.z
    """
    
    vdwr1 = p1.VanDerWaalsRadius
    vdwr2 = p2.VanDerWaalsRadius
    vdwr3 = p3.VanDerWaalsRadius
    
    Ka = -(vdwr1**2) + (vdwr2**2) + (p1.x**2) - (p2.x**2) + (p1.y**2) - (p2.y**2)
    Kb = -(vdwr1**2) + (vdwr3**2) + (p1.x**2) - (p3.x**2) + (p1.y**2) - (p3.y**2)
    
    D = (p1.x * (p2.y - p3.y)) + (p2.x * (p3.y - p1.y)) + (p3.x * (p1.y - p2.y))
    A0 = ((Ka * (p1.y - p3.y)) + (Kb * (p2.y - p1.y))) / (2 * D)
    B0 = -1 * ((Ka * (p1.x - p3.x)) + (Kb * (p2.x - p1.x))) / (2 * D)
    A1 = -1 * ((vdwr1 * (p2.y - p3.y)) + (vdwr2 * (p3.y - p1.y)) + (vdwr3 * (p1.y - p2.y))) / D
    B1 = ((vdwr1 * (p2.x - p3.x)) + (vdwr2 * (p3.x - p1.x)) + (vdwr3 * (p1.x - p2.x))) / D

    C0 = ((A0 - p1.x)**2) + ((B0 - p1.y)**2) - (vdwr1**2)
    C1 = (A1 * (A0 - p1.x)) + (B1 * (B0 - p1.y)) - vdwr1
    C2 = (A1**2) + (B1**2) - 1
    
    r = ((-1 * C1) - math.sqrt((C1**2) - (C0 * C2))) / C2
    
    new_point = AtomicPoint(A0 + (A1 * r), B0 + (B1 * r), p1.z, "")
    
    return new_point


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
    

def threeAtomCircumcenter(p1, p2, p3):
    """
    Function for calculating the center of a circle described by how it touches three non-linear atoms.
    Prof. Grimm's version of the equation, explinations in paper/file on three atom circumcenter
    
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
    rCheck1 = threeAtomRadiusCheck(p1, p2, p3, sol_1) #Checks if the circle touches the surface of each atom
    
    #solution 2
    sol_2 = AtomicPoint(x_2, y_2, z_2, "") 
    sol_2.VanDerWaalsRadius = abs(r_ans2)
    rCheck2 = threeAtomRadiusCheck(p1, p2, p3, sol_2) #Checks if the circle touches the surface of each atom

    '''
    print("Grimm's code answers 1 and 2")
    sol_1.printCoords()
    sol_2.printCoords()
    print('')
    '''
    
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
    

def circleSurrounded(p1, p2, p3, c):
    """
    Function checks if a circle c is surrounded by points p1, p2, p3. A surrounded circle 
    would be unable to grow. Returns True if c is surrounded, False otherwise.

    @param p1 AtomicPoint of bounding circle 1
    @param p2 AtomicPoint of bounding circle 2
    @param p3 AtomicPoint of bounding circle 3
    @param c AtomicPoint of the center of circle c
    @return List. [0] is a Boolean, True if circle c is surrounded by p1, p2, and p3, 
            otherwise returns False. If [0] is False, [1] is one edge of the 
            arc of three points about c and [2] is the other edge 
    """
    
    p = path.Path([(p1.x,p1.y),(p2.x,p2.y),(p3.x,p3.y)])
    point = [c.x,c.y]
    
    if p.contains_point(point):
        return [True, p1, p2]
    
    #Triangle 1 - p1 excluded
    a = p2.distBetween2D(c, False)
    b = p3.distBetween2D(c, False)
    cp = p2.distBetween2D(p3, False)
    p1_theta = cosLawAngle(a, b, cp)
    
    #Triangle 2 - p2 excluded
    a = p1.distBetween2D(c, False)
    b = p3.distBetween2D(c, False)
    cp = p1.distBetween2D(p3, False)
    p2_theta = cosLawAngle(a, b, cp)
    
    #Triangle 3 - p3 excluded
    a = p1.distBetween2D(c, False)
    b = p2.distBetween2D(c, False)
    cp = p1.distBetween2D(p2, False)
    p3_theta = cosLawAngle(a, b, cp)
    
    #Find the largest angle. Return the atoms the the angle isn't attached to
    n = max(p1_theta, p2_theta, p3_theta)
    if (n == p1_theta): return [False, p2, p3]
    elif (n == p2_theta): return [False, p1, p3]
    elif (n == p3_theta): return [False, p1, p2]


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
    
    #print ('Xmax = %f, Xmin = %f, Ymax = %f, Ymin = %f' % (Xmax,Xmin,Ymax,Ymin)) 
    return [Xmax, Xmin, Ymax, Ymin, Zmax, Zmin]


def inBounds(a, bounds):
    '''
    Function checks if an AtomicPoint falls withing the given (x,y,z) bounds

    @param a AtomicPoint
    @param bounds
    @return Boolean. Returns True if point is in bounds, False otherwise.
    '''
    
    if(a.x > bounds[0] or a.x < bounds[1] or a.y > bounds[2] or a.y < bounds[3] or a.z > bounds[4] or a.z < bounds[5]): 
        return False
    else: return True


def growCircle (pointArray, calcPoint, bounds):
    """
    Given a list of AtomicPoints, finds the largest circle that can fit in the 
    locations defined, as calculated from the calcPoint. Intended to be used to 
    find the largest circle that can fit through the pore of a Metal Organic 
    Framework, with atoms described by an input pdb file.  

    @param pointArray List of AtomicPoints describing the position of the atoms
    @param calcPoint AtomicPoint that the calculation grows from
    @param bounds Array of floats holdng the x,y bounds of the pointArray [Xmax, Xmin, Ymax, Ymin]
    @return List containing [0] the float radius of the circle grown, [1] the AtomicPoint details of the calcPoint, 
            and [2] a list of the 3 AtomicPoints that define the circle about the calcPoint
    """
    
    #Find the atomic point nearest to the calc point, accounting for vdwr
        #Initializes the points
    ap1 = pointArray[0] #should initialze with different values
    ap2 = pointArray[0] #could do AtomicPoint(abs(Xmax)+abs(Xmin), abs(Ymax)+abs(Ymin), 0, "")
    ap3 = pointArray[0]
    current_radius = calcPoint.distBetween2D(ap1, True)
    maxCycles = 100
    
    #Adjustment for zSlice
    originZ = calcPoint.z
    
    
    #First circle snap calculation --> 0 point growth
    for p in pointArray:
        if(calcPoint.distBetween2D(p, True) < calcPoint.distBetween2D(ap1, True)):
            ap1 = p
            current_radius = calcPoint.distBetween2D(ap1, True)
    
    #Checks if calcPoint is inside the ap1 vdwr and cancels calculation if so
    #Should probably add something that instead find the nearest open area and
    #continues the claulcation from there
    if calcPoint.distBetween2D(ap1, False) < ap1.VanDerWaalsRadius: return [-3, AtomicPoint(0, 0, 0, "")] 
    
    #List of step-sizes, could probably make soemthing more dynamic, but this should cover my basis for the moment
    scaleList = [1, 0.1, 0.01, 0.001, 0.0001, 0.00001]
    
    #Used to calculate angle from atomic point to calc point, vdwr not needed
    trig_x = ap1.x - calcPoint.x
    trig_y = ap1.y - calcPoint.y
    theta = math.acos(abs(trig_x) / math.hypot(trig_x, trig_y))
    
    #what I am currently using as logic for the while loop
    calculating = True
    stepIndex = 0
    calcCycles = 0
    
    #Second Circle Snap calculation --> 1 point growth
    #Could probably abstract this to its own function
    while calculating == True:
        
        #makes sure function doesn't get caught in infinite loops and just cancles the program if things go on for too long
        if calcCycles > maxCycles: return [-4, AtomicPoint(0, 0, 0, "")] 
        calcCycles = calcCycles + 1
        
        #X and Y displacments for how far to move c
        disp_x = scaleList[stepIndex] * math.cos(theta)
        disp_y = scaleList[stepIndex] * math.sin(theta)        
        #Grow away from that nearest point
        
        new_x = calcPoint.x + (-sign(trig_x) * disp_x)
        new_y = calcPoint.y + (-sign(trig_y) * disp_y)
    
        #out of bounds checking
        if(new_x > bounds[0] or new_x < bounds[1] or new_y > bounds[2] or new_y < bounds[3]): 
            
            #print(f'1.1: new_x = {new_x}, new_y = {new_y}')
            
            return [-1, AtomicPoint(0, 0, 0, '')] #error return
    
        new_point = AtomicPoint(new_x, new_y, 0, "")
        new_radius = new_point.distBetween2D(ap1, True)
    
        eligible_points_counter = 0
        
        #Checks every point in the  given atomic array, and sees if any other the 
        #first contact point are at the edge or inside the bounds of the new radius
        for p in pointArray:
            if(new_point.distBetween2D(p, True) <= new_radius and ap1 != p):
                eligible_points_counter = eligible_points_counter + 1
                ap2 = p
                #Track points
        
        #Only proceeds if only 1 new point is encountered. Otherwise decriments step 
        #size and checks again.
        #If only 1 new point is found (ignoring the contact point), snaps the circle
        #to includee the new contact point, the proceeds to next step of calculation
        if(eligible_points_counter > 1):
            #reset, use smaller incriment
            if stepIndex < 5:
                stepIndex = stepIndex + 1
                continue
            else:
                #print("Incrimenting Error.")
                return [-2, AtomicPoint(0, 0, 0, '')] #error return
        elif(eligible_points_counter == 1):
            #perform circle snap calculation for the two points
            calculating = False
            
            #calculates location of new point, function above
            calcPoint = twoAtomCircle(ap1, ap2, calcPoint)
            current_radius = ap1.distBetween2D(calcPoint, True)
            
        else:
            calcPoint = new_point
            current_radius = new_radius #Technicaly don't need to track this here
        
    #Third Circle Snap calculation --> 2 (or 3) point growth
    
    calculating = True
    stepIndex = 0
    calcCycles = 0
    
    while calculating == True:
        #In this case, ap1 and ap2 are the circles the algorithm is 'growing' from
        #The vector direction of growth is perpendicular to the line created from
        #the points of the calculating circle touching the atomic points
        
        #makes sure function doesn't get caught in infinite loops and just cancles the program if things go on for too long
        if calcCycles > maxCycles: 
            
            '''
            Circle Surround seems to be the failing test, or at least certain circle configurations can't seem to grow properly'
            
            if eligible_points_counter == 1:
                return
            '''
            
            return [-5, AtomicPoint(0, 0, 0, "")] 
        calcCycles = calcCycles + 1
        
        '''
        #Might be worth abstracting vector calculation to a funtion or using the methods
        #python provides to simplify, should test this section.
        '''
        
        #Vector from ap1 toward calcPoint
        vp1cX = calcPoint.x - ap1.x
        vp1cY = calcPoint.y - ap1.y
        #Unit Vectors
        uvp1cX = vp1cX / math.sqrt(vp1cX**2 + vp1cY**2)
        uvp1cY = vp1cY / math.sqrt(vp1cX**2 + vp1cY**2)
        
        #Vector from ap2 toward calcPoint
        vp2cX = calcPoint.x - ap2.x
        vp2cY = calcPoint.y - ap2.y
        #Unit Vectors
        uvp2cX = vp2cX / math.sqrt(vp2cX**2 + vp2cY**2)
        uvp2cY = vp2cY / math.sqrt(vp2cX**2 + vp2cY**2)
        
        #Add vectors to get average vector direction of travel
        #Don't need to average the vectors if I normalize afterwards anyways, right?
        vX = (uvp1cX + uvp2cX) / 2
        vY = (uvp1cY + uvp2cY) / 2
        
        #Convert to unit vectors again
        uvX = vX / math.sqrt(vX**2 + vY**2)
        uvY = vY / math.sqrt(vX**2 + vY**2)
        
        #New displacement in x and y
        disp_x = scaleList[stepIndex] * uvX
        disp_y = scaleList[stepIndex] * uvY
        
        #Get new location of calculation
        new_x = calcPoint.x + disp_x
        new_y = calcPoint.y + disp_y
        
        #out of bounds checking
        if(new_x > bounds[0] or new_x < bounds[1] or new_y > bounds[2] or new_y < bounds[3]): 
            return [-1, AtomicPoint(0, 0, 0, '')] #error return
        
        new_point = AtomicPoint(new_x, new_y, 0, "")
        new_radius = new_point.distBetween2D(ap1, True)
    
        eligible_points_counter = 0
        
        #Checks every point in the given atomic array, and sees if any other 
        #point at or inside the bounds of the new radius other than ap1 and ap2
        for p in pointArray:
            if(new_point.distBetween2D(p, True) <= new_radius and ap1 != p and ap2 != p):
                eligible_points_counter += 1                
                ap3 = p
        
        #Only proceeds if only 1 new point is encountered. Otherwise decriments step 
        #size and checks again.
        #If only 1 new point is found (ignoring the contact point), snaps the circle
        #to includee the new contact point
        if(eligible_points_counter > 1):
            #reset, use smaller incriment
            if stepIndex < 5:
                stepIndex = stepIndex + 1
                continue
            else:
                #if it gets to the state where there are still too many eligible points at
                #0.001 step incriments, calc a different way
                return [-6, AtomicPoint(0, 0, 0, '')] #error return
            
        elif(eligible_points_counter == 1 and ap1 != ap2 and ap1 != ap3 and ap2 != ap3):
            #perform circle snap calculation for the three points
            calculating = False
            
            #Now need to check if the points bound the calcPoint on all side/limit growth
            #If they do thats the final circle. If they don't and all three point's fall 
            #within a <180 degree arc, find the outermost circles of the three circles 
            #and repeat this loop with those points as the new ap1 and ap2
            surround = circleSurrounded(ap1, ap2, ap3, calcPoint)
            
            if surround[0]:
                #Current return is in the form of a list with radius and calcPount, but could just 
                #make vdwr the current radius if I wanted to 
                #calculates location of new point, function above
                calcPoint = threeCircleCircumcenter(ap1, ap2, ap3)
                current_radius = ap1.distBetween2D(calcPoint, True)
                
                calcPoint.z = originZ 
                
                return [current_radius, calcPoint, [ap1, ap2, ap3]]
            #if circle is not surounded, continue to grow the circle, the 2 outermost of the 
            #3 defining atomic points servining as the new ap1 and ap2
            else:                
                ap1 = surround[1]
                ap2 = surround[2]
                
                calculating = True
                stepIndex = 0
            
        else:
            calcPoint = new_point
            current_radius = new_radius #Technicaly don't need to track this here
            
    #for error handling, the if statement above is the real return statement.
    return