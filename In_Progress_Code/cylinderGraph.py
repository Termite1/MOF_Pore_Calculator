# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 22:28:01 2024

Code used to form and search cylinder graph network for determining the largest
sphere that can pass through a pore

@author: samda
"""

from AtomicPoint import AtomicPoint
import lcfProcessing as proc
import LargestCircleFinder as lcf
import numpy as np
import math
import random


class cylinderNode:
    def __init__(self, ap, a, b, c, **kwargs):
        '''
        Initilizes cylinder node. Used to perform cylinder graph search.

        Parameters
        ----------
        ap : AtomicPoint
            AtomicPoint describing the sphere that defines this cylinder
        a : AtomicPoint
            One of the AtomicPoints surrounding sphere ap
        b : AtomicPoint
            One of the AtomicPoints surrounding sphere ap
        c : AtomicPoint
            One of the AtomicPoints surrounding sphere ap
        **kwargs
            no_vector : Any value
                if present, ignores vector calculation

        Returns
        -------
        None.

        '''
        #AtomicPoint
        self.ap = ap
        self.r = ap.VanDerWaalsRadius #double stating just so value is easier to work with
        
        self.surrounding = [a,b,c]
        
        #This is direction of surface normal discribing given sphere (as a unit vector)
        if 'no_vector' in kwargs:
            self.vector = None
        elif 'vector' in kwargs:
            self.vector = kwargs['vector']
        else: self.vector = proc.normalize(proc.planarize(a, b, c)[1])
        
        #Need to compute these
        self.heightP = None #Positive height, vector in +z direction
        self.heightN = None ##Positive height, vector in -z direction, could probably make both arbitrary
       
        #If True indicates that this cylinder connects to the top or bottom of the pore
        self.startConnection = False
        self.endConnection = False
        
       
        
class cylinderEdge:
    def __init__(self, a, b, r):
        '''
        Represent an edge between two cylinder, or the cylinder that can 
        connect two cylinder nodes.

        Parameters
        ----------
        a : cylinderNode
            First cylinderNode.
        b : cylinderNode
            second cylinderNode.
        r : float
            Radius of connecting cylinder.

        Returns
        -------
        None.
        '''
        self.nodes = [a, b]
        self.r = r
        

class cylinderGraph:
    def __init__(self):
        '''
        Used to perform cylinder search. Cylinder derived from calculated spheres
        are nodes while cylinders concting nodes are edges.

        Returns
        -------
        None.

        '''
        self.myVerticies = []
        self.myAdjList = []
        
        #Used for exploring network
        self.accessable = False
        self.reached = False
        
        #Will need to create/initilize these
        self.start = None #Something representing the top of the pore
        self.end = None #something representing the bottom of the pore
        
    
    def genGraph(self, pointList, nodeList, bounds):
        '''
        Generates cylinder graph network for exploration of pore with best-
        first search.

        Parameters
        ----------
        pointList : List of AtomicPoints
            Atomic environemnt being explored.
        nodeList : List of AtomicPoints
            Sphere nodes that will be turned into cylindrical nodes.
        bounds : List of floats [Xmax, Xmin, Ymax, Ymin, Zmax, Zmin]
            Bounds of the pore in the positive and negative directions of each 
            axis, based on the atoms the furthest in each of those directions.

        Returns
        -------
        None.

        '''
        #Generates cylinder node list from the AtomicPoints in the nodeList
        for i in range(len(nodeList)):
            new_node = cylinderNode(nodeList[i][1], nodeList[i][2][0], nodeList[i][2][1], nodeList[i][2][2])
            self.myVerticies.append(new_node)
            cylinderHeight(pointList, self.myVerticies[i], bounds)
            
            print('')
            new_node.ap.printCoords()
            print(new_node.heightN + new_node.heightP)
            print('')
        
        #Run through every pair of nodes
        for j in range(len(self.myVerticies)):
            for k in range(j+1, len(self.myVerticies)):
                #Checks if there is viable connection
                temp_r = cylinderConnection(pointList, self.myVerticies[j], self.myVerticies[k], True)
                
                if temp_r == 0: continue
                #if so, creates edge and adds it to the adjacency list    
                temp_edge = cylinderEdge(self.myVerticies[j], self.myVerticies[k], temp_r)
                self.myAdjList.append(temp_edge)
        

def exploreCylinderGraph(cylinderGraph):
    '''
    Function uses best-first search to explore given cylinder network. Returns
    the radius of the smallest cylinder passed through (the radius of the 
    largest sphere that could actually trravel through the network.)

    Parameters
    ----------
    cylinderGraph : cylinderGraph
        The graph being explored.

    Returns
    -------
    float
        Returns the radius of the smallest cylinder passed through (the radius 
        of the largest sphere that could actually trravel through the network.)

    '''
    reached_nodes = []
    #accessable_nodes = []
    
    best_r = 10000 #Tracks current passthrough radius
    
    end_reached = False
    
    #Creates a starting node
    dummy_point = AtomicPoint(0, 0, 0, '')
    startNode = cylinderNode(dummy_point, None, None, None, no_vector=True)
    startNode.reached = True
    
    #Adds edges connecting nodes with a start connection to the start
    for node in cylinderGraph.myVerticies:
        if node.startConnection == True:
            cylinderGraph.myAdjList.append(cylinderEdge(startNode, node, node.r))
            
            #not sure I need
            #accessable_nodes.append(node)
    
    #Adds start node to reached list
    reached_nodes.append(startNode)
    
    #Continues until end is reached
    while end_reached == False:
        #Check all reached nodes for an end connection. If present, end function
        for rn in reached_nodes:
            if rn.endConnection == True:
                return best_r
        
        #tracks current best edge and r of that edge
        tracking_r = 0
        tracking_edge = None
        #Run through edges that have only 1 reached node, keep edge with largest r that is still
        #less than or equal to best_r
            #See if I can limit r to 3 decimal floats to equals works better? The might just be for print
        for node in reached_nodes:
            for edge in cylinderGraph.myAdjList:
                if node in edge.nodes:
                    #skips edges where neither nodes have already been reached
                    if edge.nodes[0].reached == False and edge.nodes[1].reached == False: continue
                    #skips edges where both nodes have already been reached
                    if edge.nodes[0].reached == True and edge.nodes[1].reached == True: continue
                    
                    #tracks largest r edge whith r less than or equal to current best_r
                    if edge.r <= best_r and edge.r > tracking_r:
                        tracking_r = edge.r
                        tracking_edge = edge
        
        #grabs newly reached node
        new_node = None
        if tracking_edge.nodes[0].reached == False: new_node = tracking_edge.nodes[0]
        else: new_node = tracking_edge.nodes[1]
        
        new_node.reached = True
        reached_nodes.append(new_node)
        
        #Adds all adjacent nodes not in accessible list to accessible list (don't actully need to do this?)
        '''
        for edge in cylinderGraph.myAdjList:
            if new_node in edge.nodes:
                other_node = None
                if edge.nodes[0] == new_node: other_node = edge.nodes[1]
                else: other_node = edge.nodes[0]
                
                if other_node.accessable
        '''


def cylinderHeight(pointList, cn, bounds):
    '''
    Takes a given cylinderNode, determines its positive and negative height,
    and sets those values.

    Parameters
    ----------
    pointList : List of AtomicPoints
        Envronment cylinder is present in.
    cn : cylinderNode
        Cylinder height is being calculated for.
    bounds : List of floats [Xmax, Xmin, Ymax, Ymin, Zmax, Zmin]
        Bounds of the pore in the positive and negative directions of each 
        axis, based on the atoms the furthest in each of those directions.

    Returns
    -------
    None
    '''
    #Need to impliment catch for if cn.vector[2] == 0
    
    #Draw line of sight to above and bellow the pore
    topZ = bounds[4] + 2
    botZ = bounds[5] - 2
    
    #Find value cylinder directional vector needs to be multiplied by to reach
    #the set z value
    dist_to_top = (topZ - cn.ap.z) / cn.vector[2]
    dist_to_bot = (botZ - cn.ap.z) / cn.vector[2]

    #Find the x and y position of the point at that z value
    topX = cn.vector[0] * dist_to_top
    botX = cn.vector[0] * dist_to_bot
    
    topY = cn.vector[1] * dist_to_top
    botY = cn.vector[1] * dist_to_bot

    #Generates top and bottom point location as AtomicPoints
    top = AtomicPoint(topX, topY, topZ, '')
    bot = AtomicPoint(botX, botY, botZ, '')

    #Finds the nearest colliding atom in the top or bottom direction, if there is one
    
    print('top')
    
    top_pbound = quickCylinderLineOfSight(pointList, cn.ap, top, cn.surrounding, cn.r - 0.1)
    
    print("")
    print('bot')
    
    bot_pbound = quickCylinderLineOfSight(pointList, cn.ap, bot, cn.surrounding, cn.r - 0.1)    
    print('')


    #If line of sight unobscured, fuction sets cylinder height and notes that it connects to the top or bottom of the pore
    #If line of sight is obscured, calculates distance to obscuring atom and sets cylinder height
    if top_pbound == -1:
        cn.heightP = abs(dist_to_top)
        cn.startConnection = True
    else:
        cn.heightP = top_pbound
    if bot_pbound == -1:
        cn.heightN = abs(dist_to_bot)
        cn.endConnection = True
    else:
        cn.heightN = bot_pbound


def cylinderConnection(pointList, a, b, reportRadius, **kwargs):
    '''
    Function checks if two cylinders can be connected by a third cylinder 
    through a given molecular structure. Line of sight is drawn as a cylinder 
    whose radius can be returned if reportRadius is True. If reportRadius is 
    False, then if a cylinder can connect the two given cylinder nodes True is 
    returned. If the nodes have no line of sight, False is returned.

    Parameters
    ----------
    pointList : List of AtomicPoints
        List of AtomicPoints describing the molecular pore being investigated
    a : cylinderNode
        First cylinderNode.
    b : cylinderNode
        Second cylinderNode.
    reportRadius : boolean
        If True, then if there is a cylinder that can connect the two given
        cylinderNodes, that cylinder's radius is returned. If False, then if
        a cylinder can connect the two given cylinder nodes, only True is 
        returned. No matter what, if no cylinder can connect the nodes then
        False is returned.
    **kwargs
        cycles : int
            Number of generation genetic algorithm goes through, default 10.
            Cycles + 1 generations created.
        population : int
            Number of elements to have in population, default 24. 
            (Value must be divisble by 4).

    Returns
    -------
    If reportRadius == True:
        float (radius of line of sight sphere) or 0 (no line of sight)
    If reportRadius == False
        True (there is line of sight) or 0 (no line of sight)
    '''
    
    #Try using a genetic algorithm? Still need to account for angle difference
    #3 terms, h1 (dist along cylinder 1), h2 (dist along cylinder 2), w2 (shift along cylinder w)
    #cylinder 1 is the cylinder with smaller r
    
    major_mutation_rate = 0.1 #Chance for big shift in position
    major_mutation_varience = 1 #Amount of possible varience. 1 means value is randomly regenerated
    minor_mutation_rate = 0.4 #Chance for small shifts in position
    minor_mutation_varience = 0.05 #Amount of possible varience in either direction. 
                        #0.05 means up to 5% varience of total from current position in either direction
    
    #Size of child population, parent population half child population
    population = 24
    if 'population' in kwargs:
        population = kwargs['population']
    
    #number of generation to run the genetic algorithm through
    cycles = 10
    current_cycle = 0
    if 'cycles' in kwargs:
        cycles = kwargs['cycles']
    
    #Estabilishes maximum radius of connecting cylinder, equal to the radous of the smaller cylinder
    max_r = 0
    lc = None #Large cylinder
    sc = None #Small cylinder
    if a.r > b.r: 
        max_r = b.r
        lc = a
        sc = b
    else: 
        max_r = a.r
        lc = b
        sc = a
    
    L1 = sc.heightN + sc.heightP #Length 1, length along smaller cylinder
    L2 = lc.heightN + lc.heightP #Length 2, length along larger cylinder
    total_W2 = 2 * (lc.r - sc.r) #Width 2, absolute amount of variation possible across larger cylinder's body
                                 #values are going to be half it, but range from -W2/2 to W2/2.
    
    #generate initial parent population
    children = []
    
    i = 0 #index tracking addition of parents to parent list
    while i < population:
        #Generates random values for intial parents
        new_L1 = random.random() * L1
        new_L2 = random.random() * L2
        new_W2 = (random.random() * total_W2) - (total_W2/2)
        
        #[L1, L2, W2, calculated r]
        new_connection = [new_L1, new_L2, new_W2, 1000]
        children.append(new_connection)
        i += 1
    
    #Genetic algorithm, runs until cycle # of iteration have occured
    while current_cycle < cycles:
        #Calculates largest open connecting cylinder that can be created along the path described
        for cylinder in children:
            #Skips calculation for values that have already been calculated
            if cylinder[3] != 1000: continue 
            
            #Generate AtomicPoint locations on the fly
            s_point = cylinderPoint(sc, cylinder[0], max_r) #point on smaller cylinder
            il_point = cylinderPoint(lc, cylinder[1], max_r) #intial point on larger cylinder
            
            #Vector line from s_point toward il_point 
            iAB = np.array([il_point.x-s_point.x, il_point.y-s_point.y, il_point.z-s_point.z]) 
            
            #regnerates point on larger cylinder accounting for W2
            l_point = cylinderPoint(lc, cylinder[1], max_r, v=iAB, w=cylinder[2]) 
            
            #Calculates number of atoms obscuring generated connection 
            cylinder[3] = cylindricalLineOfSightV2(pointList, s_point, l_point, max_r)
            if cylinder[3] == 0: #If unobscured connection found, ends function
                                     #and returns it immediatly
                if reportRadius: return max_r
                else: return True
            '''
            #Ends calculation if cylinder exists and isn't reporting radius
            #if reportRadius == False and cylinder[3] > 0: return True
            '''
                
        #Sorts children list by number of obscuring atoms in acending order
        children.sort(key = lambda x: x[3])
        
        #Selects parent population (half of children population)
        new_pop = int(population / 2)
        parents = children[:new_pop]
        
        #Saves top two results to be added to subsequent generation with no changes
        first = parents[0]
        second = parents[1]
        
        #shuffles parent list
        random.shuffle(parents)
        
        children = []
        
        #Generate children (crossover between pairs of parents) and cause mutations
        j = 0 #index for children generation
        while j < new_pop:
            #Gets parents
            p1 = parents[j]
            p2 = parents[j+1]
            
            #Need to fill out crossover function
            new_children = crossover(p1, p2, major_mutation_rate, major_mutation_varience, minor_mutation_rate, minor_mutation_varience, L1, L2, total_W2)
            
            #Adds newly created children to children list
            for child in new_children:
                children.append(child)
            
            j += 2
        
        #Adds previous best two results to children list
        children.append(first)
        children.append(second)
        
        current_cycle += 1
    
    
    #Performs a last check on the resulting children
    #Calculates largest open connecting cylinder that can be created along the path described
    for cylinder in children:
        #Skips calculation for values that have already been calculated
        if cylinder[3] != 1000: continue 
        
        #Generate AtomicPoint locations on the fly
        s_point = cylinderPoint(sc, cylinder[0], max_r) #point on smaller cylinder
        il_point = cylinderPoint(lc, cylinder[1], max_r) #intial point on larger cylinder
        
        #Vector line from s_point toward il_point 
        iAB = np.array([il_point.x-s_point.x, il_point.y-s_point.y, il_point.z-s_point.z]) 
        
        #regnerates point on larger cylinder accounting for W2
        l_point = cylinderPoint(lc, cylinder[1], max_r, v=iAB, w=cylinder[2]) 
        
        #Calculates radius for generated connection 
        cylinder[3] = cylindricalLineOfSightV2(pointList, s_point, l_point, max_r)
        if cylinder[3] == 0: #If unobscured connection found, ends fucntion
                                 #and returns it immediatly
            if reportRadius: return max_r
            else: return True
    '''
    #Sorts children list by calculated radius in decending order
    children.sort(reverse=True, key = lambda x: x[3])
    
    #Checks if there is at least 1 valid solution (selects the best valid solution)
    if children[0][3] > 0:
        if reportRadius: return children[0][3]
        else: return True
    else: return 0
    '''
    return 0
    
    
def cylindricalLineOfSight(pointList, a, b, r):
    '''
    Function calculates if there is a cylindrical line of sight between two 
    points. If there are spheres (atoms) in the way, the function will return
    the radius of the largest cylinder that stays within the given radius.
    
    #Note for future, see if I can adapt this to account for angle varience.
    #May not need. Think about having option to return location of line-of-sight 
    #cylinder though

    Three atom circumcenter doesn't work the wway it would need to for this right now.

    Parameters
    ----------
    pointList : List of AtomicPoints
        Environment contianing spheres (atoms) the line of sight is being drawn
        through.
    a : AtomicPoint
        First point line of sight is being drawn from.
    b : AtomicPoint
        Second point line of sight is being drawn from.
    r : float
        Radius of cylindrical line of sight being drawn.

    Returns
    -------
    float
        Radius of largest cylindrical line of sight that can be drawn in the given
        radius. 0 if no line of sight can be drawn.

    '''
    
    AB = np.array([b.x-a.x, b.y-a.y, b.z-a.z]) #Vector line between point a and point b
    mag_AB = magnitude(AB) #Vector magnitude
    
    #List of AtomicPoints the cylinder colides with
    collision_list = []
    
    #repeats for every AtomicPoint in the pointList
    for p in pointList:
        #Vector line between point a and point b
        AP = np.array([p.x-a.x, p.y-a.y, p.z-a.z]) 
        
        #Cross product of AB and AP
        cp = np.cross(AB, AP) 
        
        #Shortest distance between p and AB
        CD = magnitude(cp) / mag_AB
        
        #Makes shure the distance is larger than the Van Der Waals radius.
        #if it isn't the line passes throgh the atom, and the atom obscures the line of sight
        #between the two points
        if CD < (p.VanDerWaalsRadius + r):
            collision_list.append(p)
        
    #If line not obstructed, return r
    if len(collision_list) == 0:
        return r
    
    #project all the point collided with into a 2D plane with surface normal AB
    #Generate plane axes and origin
    surface_normal = proc.normalize(AB)
    axis1 = proc.normalize(genOrthoginalVector(AB))
    axis2 = proc.normalize(np.cross(AB, axis1))
    
    #Origin point halfway between a and b
    origin = AtomicPoint(a.x + AB[0], a.y + AB[1], a.z + AB[2], '')
    
    plane_values = [origin, surface_normal, axis1, axis2]
    
    #Project all collided atoms onto generated 2D plane
    projected_list = []
    for point in collision_list:
        projected_list.append(proc.planarizePoint(point, plane_values))
    
    #maybe have some way to eliminate answers to reduce calculation times?
    #check if cylinder completly obscured?
    
    
        #check for collision and throw those out
        #Make sure this calculates the internal 3 atom circle? - test with large circle and
        #two internal smaller circles
    cylinder_circle = AtomicPoint(0, 0, 0, '', vdwr=r)
    
    best_r = 0
    
    #Method requires at least 2 collided atoms, checks for that
    if len(projected_list > 1):
        #Calculate the three atom circumcenter of every pair of collided atom with
        #the circle projection of the surrounding cylinder
        for i in range(len(projected_list)):
            for j in range(i+1, len(projected_list)): #Ensures there is no double checking
                ans = lcf.threeAtomCircumcenter(projected_list[i], projected_list[j], cylinder_circle)
                
                intersection = False
                #Checks if circle intersects other atoms, discards result if so
                for projected_point in pointList:
                    if projected_point.distBetween2D(ans, True) < 0 and projected_point != projected_list[i] and projected_point != projected_list[j]:
                        intersection = True
                        break
                
                #If there is an intersection, skips to next circle
                if intersection: continue
                
                #Keep best radius found
                if ans.VanDerWaalsRadius > best_r: best_r = ans.VanDerWaalsRadius
        
    #Calculate the two atom circumcenter for every collided atom with and the 
    #circle projection of the surrounding cylinder
    dummy_origin = AtomicPoint(0, 0, 0, '')
    for i in range(len(projected_list)):
        d = dummy_origin.distBetween2D(projected_list[i], True)
        
        #I have graphic of this calculation
        ans_radius = (d - projected_list[i].VanDerWaalsRadius + r) / 2
        
        #Calculates unit direction vector toward calculated circle from r1
        r1_u_v = lcf.normalize([-1 * projected_list[i].x, -1 * projected_list[i].y])
        
        #generate atomic point location --> (-r1 unit vector * r1 + radius) + r1(x,y)
        #Check this works
        ans = AtomicPoint(projected_list[i].x + (r1_u_v[0] * (projected_list[i].VanDerWaalsRadius + ans_radius)), projected_list[i].y + (r1_u_v[1] * (projected_list[i].VanDerWaalsRadius + ans_radius)), 0, '', vdwr=ans_radius)

        intersection = False
        #Checks if circle intersects other atoms, discards result if so
        for projected_point in pointList:
            if projected_point.distBetween2D(ans, True) < 0 and projected_point != projected_list[i]:
                intersection = True
                break
        
        #If there is an intersection, skips to next circle
        if intersection: continue
        
        #Keep best radius found
        if ans.VanDerWaalsRadius > best_r: best_r = ans.VanDerWaalsRadius
    
    #could modify this to return the cylinder as well, but that is a later concern
    return best_r


def cylindricalLineOfSightV2(pointList, a, b, r):
    '''
    Function calculates if there is a cylindrical line of sight between two 
    points. If there are spheres (atoms) in the way, the function will count
    the number in the way.
    
    Parameters
    ----------
    pointList : List of AtomicPoints
        Environment contianing spheres (atoms) the line of sight is being drawn
        through.
    a : AtomicPoint
        First point line of sight is being drawn from.
    b : AtomicPoint
        Second point line of sight is being drawn from.
    r : float
        Radius of cylindrical line of sight being drawn.

    Returns
    -------
    int
        Number of atoms between the two points intersecting the given radius.

    '''
    
    AB = np.array([b.x-a.x, b.y-a.y, b.z-a.z]) #Vector line between point a and point b
    mag_AB = magnitude(AB) #Vector magnitude
    
    #List of AtomicPoints the cylinder colides with
    num_collision = 0
    
    #repeats for every AtomicPoint in the pointList
    for p in pointList:
        #Vector line between point a and point b
        AP = np.array([p.x-a.x, p.y-a.y, p.z-a.z]) 
        
        #Cross product of AB and AP
        cp = np.cross(AB, AP) 
        
        #Shortest distance between p and AB
        CD = magnitude(cp) / mag_AB
        
        #ocount collisions
        if CD < (p.VanDerWaalsRadius + r):
            num_collision += 1
        
    return num_collision


def quickCylinderLineOfSight(pointList, a, b, s, r):
    '''
    Function calculates if there is a cylindrical line of sight between two 
    points. If there are spheres (atoms) in the way, the function will return
    the distance to the first atom collided with on the path from point a to 
    point b.

    Parameters
    ----------
    pointList : List of AtomicPoints
        Environment lign of sight is being calculated through.
    a : AtomicPoint
        Point line of sight is being drawn from.
    b : AtomicPoint
        Point line of sight is being drawn to.
    s : List of AtomicPoints
        AtomicPoints to ignore during calculation
    r : float
        Radius of cylindrical line of sight being drawn.

    Returns
    -------
    -1 if line of sight unobscured, the float distance to the first atom 
    collided with if the line of sight is obscured.
    '''
    AB = np.array([b.x-a.x,b.y-a.y,b.z-a.z]) #Vector line being drawn from a to b
    d = proc.normalize(AB) #Unit directional vector along line being drawn
    m_AB = magnitude(AB) #magtnitude of AB
    
    smallest_ans = 10000 #Dummy answer
    
    #Checks against every atom in the list
    for c in pointList:
        if c in s: continue #Skips ignored atoms
        #Line from a to point c
        AC = np.array([c.x-a.x,c.y-a.y,c.z-a.z])
        m_AC = magnitude(AC)
        
        BC = np.array([c.x-b.x,c.y-b.y,c.z-b.z])
        m_BC = magnitude(BC)
        
        t = np.dot(d, AC)        

        p_x = a.x + t * d[0]
        p_y = a.y + t * d[1]
        p_z = a.z + t * d[2]
        
        #Point on line AB closest to c
        p = AtomicPoint(p_x, p_y, p_z, '')

        #Line from c to p
        CP = np.array([p.x-c.x,p.y-c.y,p.z-c.z])
        
        #magnitude of CP, or distance between c and p
        m_CP = magnitude(CP)
        
        #helps check if point is in direction of AB (instead of -AB) or for helping to find
        #height
        #Length along AB line of sight till location CP would lead to is reached
        l = math.sqrt(m_AC**2 - m_CP**2)
        
        #If point doesn't intersect line of sight, next point is moved to
        if r <= m_CP - c.VanDerWaalsRadius: continue
        elif m_BC > m_AB and l > c.VanDerWaalsRadius: continue
        #If point center is inside cylindrical line of sight
        elif r >= m_CP:
            h = l - c.VanDerWaalsRadius #Distance until line of sight collides with point
        #When point radius intersects line of sight
            
        else:
            j = m_CP-r #Distance inside intersected point but outside cylinder
            k = math.sqrt(c.VanDerWaalsRadius**2 - j**2) #Distance in -AB direction to where cylinder first touches this point
        
            h = l - k #Distance until line of sight collides with point
        
        
        c.printCoords()
        print(h)
        
        
        if h < smallest_ans: smallest_ans = h #Updates smallest ans of new answer smaller
        
    #return True if the value didn't change, no points in the way.
    if smallest_ans == 10000: return -1
    else: return smallest_ans #Otherwise return distance to first collided point
        

def cylinderPoint(c, l, r, **kwargs):
    '''
    Coverts cylinder value (length along cylinder) into an
    (x,y,z) AtomicPoint position.

    Parameters
    ----------
    c : cylinderNode
        cylinder AtomicPoint is being created from.
    l : float
        Length along cylinder.
    r : float
        Radius of resulting point.
    **kwargs
        w : float
            Distance across cylinder. Requires v.
        v : numpy vector or float list [x,y,z]
            Vector describing line between the point on the smaller cylinder 
            and the original point on the larger cylinder. Used when 
            regenerating a point for the large cylinder, now acounting for
            movement across cylinder. Requires w.

    Returns
    -------
    AtomicPoint.
        (x,y,z) position described by given cylinder and the length along it.
    '''
    
    #Base point at bottom of cylinder
    bp_x = c.ap.x - (c.vector[0]*c.heightN)
    bp_y = c.ap.y - (c.vector[1]*c.heightN)
    bp_z = c.ap.z - (c.vector[2]*c.heightN)
    
    #generates new coordinate values
    new_x = (c.vector[0] * l) + bp_x
    new_y = (c.vector[1] * l) + bp_y
    new_z = (c.vector[2] * l) + bp_z
    
    #Used when regenerating point
    if 'w' in kwargs and 'v' in kwargs:
        #creates v unit directional vector
        dir_v = proc.normalize(kwargs['v'])
        
        #Generates vector defining movement along W2 axis
        W2_axis = np.cross(c.vector, dir_v)
        
        W2 = kwargs['w'] #Grabs W2 value
        
        #Generates updated values accounting for w
        u_x = new_x + (W2_axis[0] * W2)
        u_y = new_y + (W2_axis[1] * W2)
        u_z = new_z + (W2_axis[2] * W2)
        
        #generates and returns point
        u_point = AtomicPoint(u_x, u_y, u_z, '', vdwr=r)
        return u_point
    else:
        #regular generation of  point
        new_point = AtomicPoint(new_x, new_y, new_z, '', vdwr=r)
        return new_point


def crossover(p1, p2, MMr, MMv, mmr, mmv, L1, L2, W2_total):
    '''
    Function performs crossover, child generation, and child mutation for
    cylinderConnection genetic algorithm.

    Parameters
    ----------
    p1 : [L1 (float), L2 (float), W2 (float), N/A]
        First parent.
    p2 : [L1 (float), L2 (float), W2 (float), N/A]
        Second parent.
    MMr : float
        Major mutation rate. Chance of major mutation out of 1.
    MMv : float
        Percent varience from current value. If 1, value randomly regenerated.
    mmr : float
        Minor muation rate. Chance of minor mutation. Performed alongside
        major muation check (total mutation chance = MMr + mmr)
    mmv : float
        Percent varience from current value in either direction based on total
        value. Will wrap around at ends to prevent values from accumulating 
        near bounds.
    L1 : float
        Total value for L1
    L2 : float
        Total value for L2
    W2_total : float
        Max value for W2. Should be total span (Ex. if -0.5 to 0.5, put 1)

    Returns
    -------
    [L1 (float), L2 (float), W2 (float), 0] x 4
        Creates and returns list of 4 children
    '''
    #fill this out
    
    #Generates indicies to swap
    swap_vals = random.sample(range(1,3), 2)
    
    child_list = []
    max_vals = [L1, L2, W2_total/2]
    
    #Generates 4 children
    for i in range(4):
        child_list.append([0,0,0,0])
    
    n = 0 # index for children
    #Populates the children with values, perfroms swaps in pairs
    for j in swap_vals:
        #Runs through indexes 0-2
        for k in range(3):
            if k == j: #For swap
                child_list[n][k] = p2[k]
                child_list[n+1][k] = p1[k]
            else: #for non-swap
                child_list[n][k] = p1[k]
                child_list[n+1][k] = p2[k]
        n += 2 #incriments by 2 to next pair
    
    #Now to check for mutations in each field of each child
    for child in child_list:
        for m in range(3):
            chance = random.random() #Generates percentage
            
            if chance <= MMr:
                #Perform major mutaion
                child[m] = majorMutation(m, child[m], MMv, max_vals[m])
            elif chance <= MMr + mmr:
                #perform minor mutation
                child[m] = minorMutation(m, child[m], MMv, max_vals[m])
    
    return child_list


def majorMutation(i, pv, MMv, max_val):
    '''
    Function performs major mutation operation on given value for the crossover
    genetic algorithm function.

    Parameters
    ----------
    i : int
        Index so the value being mutated is known.
    pv : float
        Previous value of field
    MMv : float 
        Major mutation varience. If = to 1, value regnerated (value currently 1)
    max_val : float
        Maximum value of field being mutated.

    Returns
    -------
    float
        New mutated value of field.
    '''
    #Checks if operating on W2 field
    full_max = max_val
    if i == 2: full_max = 2 * max_val
    
    new_val = 0 #value that will be returned
    
    #if major mutation varience is 1
    if MMv == 1:
        #regenerates
        new_val = random.random() * full_max
        #Accounts for W2 offset
        if i == 2: new_val = new_val - max_val

    return new_val


def minorMutation(i, pv, mmv, max_val):
    '''
    Function performs minor mutation operation on given value for the crossover
    genetic algorithm function.

    Parameters
    ----------
    i : int
        Index so the value being mutated is known.
    pv : float
        Previous value of field
    mmv : float 
        Minor mutation varience. Varience in either direction as a percentage 
        of the maximum value, so mmv shouldn't be larger than 0.5.
    max_val : float
        Maximum value of field being mutated.

    Returns
    -------
    float
        New mutated value of field.
    '''
    #Checks if operating on W2 field
    full_max = max_val
    if i == 2: full_max = 2 * max_val
    
    new_rand = random.random()
    new_val = pv + ((new_rand * 2 * mmv * full_max) - (mmv * full_max))
    
    #Next check for over or underflow and wrap around if neccissary 
    #for W2
    if i == 2: 
        if new_val > max_val: new_val -= full_max
        elif new_val < -max_val: new_val += full_max
    else: #for the others
        if new_val > max_val: new_val -= max_val
        elif new_val < 0: new_val += max_val
    
    return new_val


def magnitude(v):
    '''
    Calculates magnitude of given vector.

    Parameters
    ----------
    v : numpy vector
        vector.

    Returns
    -------
    float
        vector magnitude
    '''
    return math.sqrt(sum(pow(element, 2) for element in v))


def genOrthoginalVector(v):
    '''
    Given a vector, generates an arbitray second vector that is orthoginal to it

    Parameters
    ----------
    v : numpy array [x, y, z]
        Vector used to generate orthoginal vector.

    Returns
    -------
    numpy array [x, y, z]
        Vector orthoginal to vector v.
    '''
    
    #Orthoginal vector's 1st element set to v's second element
    ux = v[1]
    
    #Orthoginal vector's 2nd element set to v's second element
    uy = v[0]
    
    #Solving the dot product, orthoginal vector's third element equal to the 
    #product of v's first 2 eleemtns time -2 divided by c's third element
    uz = (-2*v[0]*v[1])/v[2]
    
    #Orthoginal vector
    u = np.array([ux,uy,uz])
    return u
    