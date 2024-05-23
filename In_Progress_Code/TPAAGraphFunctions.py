# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 16:06:36 2023

@author: samda
"""

import AtomicGraph as ag
from AtomicPoint import AtomicPoint
from collections import deque


def quadrantCheck(ap):
    """Function takes an atomic point and returns which (x,y) grid quadrant it resides in.
    (x+,y+) = 1, (x-,y+) = 2, (x+,y-) = 3, (x-,y-) = 4. If x and/or y = 0, returns 0
    
    @param ap AtomicPoint
    @return 0, 1, 2, 3, or 4, see conditions above
    """
    
    if ap.x > 0:
        if ap.y > 0: return 1
        elif ap.y < 0: return 3
        elif ap.y == 0: 
            print(f"Qudrant Check Error: y = 0. x = {ap.x}")
            return 0
    elif ap.x < 0:
        if ap.y > 0: return 2
        elif ap.y < 0: return 4
        elif ap.y == 0: 
            print(f"Qudrant Check Error: y = 0. x = {ap.x}")
            return 0
    elif ap.x == 0:
        print(f"Qudrant Check Error: x = 0. y = {ap.y}")
        return 0


def axisCheck(al):
    """Function takes a set of 4 atomic points from a linker and returns which pair of grid (x,y) quadrants it resides in.
    (x+,y+-) = 5, (x+-,y+) = 6, (x-,y+-) = 7, (x+-,y-) = 8. If != 2 quadrants are found, returns 0 and something funcky is going on.
    
    @param al List of 4 AtomicPoints
    @return 0, 5, 6, 7, or 8. see conditions above
    """
    #Each quadrant is associated with quadrant - 1 as its related index to align with the list starting index of 0
    quadList = [False] * 4
    #Goes through each AtomicPoint in the al list
    for ap in al:
        #Finds the quadrant the atom belongs too then set the relvant index to True
        quadList[quadrantCheck(ap) - 1] = True
    
    #Used for checking if there are 2 True in the list
    checkTot = 0
    for check in quadList:
        if check == True:
            checkTot += 1
            
    if checkTot != 2:
        #Something funky is going on if this happens, linker exists in 1 or 3 or 
        #more quadrants, which shouldn't be possible for MOF-5 or similar MOFs 
        #with a pore centered at the origins and linkers oriented along axes
        print(f"Linker Check Error: {checkTot} quadrants detected, only 2 should be associated with a linker")
        return 0
    
    #q 1 and 3
    if quadList[0] and quadList[2]: return 5
    #q 1 and 2
    elif quadList[0] and quadList[1]: return 6
    #q 2 and 4
    elif quadList[1] and quadList[3]: return 7
    #q 3 and 4
    elif quadList[2] and quadList[3]: return 8
        
    print(f"Linker Check Error: Diagonal quadrants detected")
    return 0


def updateQuadrant(myGraph, quad, start):
    """
    Function updates the quadrant of a descrete portion of an AtomicGraph
    connected by edges so the verticies all belong to that quadrant.
    
    @peram myGraph AtomicGraph being operated on that verticies and edges exist in
    @peram quad Qudrant that vertex quadrants will be changed to
    @peram start is an AtomicVertex that belongs to the descrete graph to be updated the the function will start from
    """
    q = deque()
    start.predecessor = True #Just as a way of indicating that this is the beggining node
    q.append(start)
    qLen = 1
    
    #Runs as long as there are items in the queue
    while qLen != 0:
        tempV = q.popleft()
        qLen -= 1    
        
        #Updates vertex's quadrant
        tempV.quadrant = quad
        
        for e in myGraph.myAdjList:
            if e == None: continue #For the end of the list to get it done with quickly. Adj list should end with None values
            
            if myGraph.edgeContains(tempV, e):
                #Grabs paired vertex from edge
                newV = myGraph.pairReturn(tempV, e)
                #Makes sure vertex hasn't been encountered before
                if newV.predecessor == None:
                    #Adds vertex to the queue
                    newV.predecessor == tempV
                    q.append(newV)
                    qLen += 1
    

def genTPAACapGraph(myGraph, pointList):
    """Function generates graph of caps on a MOF with Triphenylacetic acid-based caps.
    Assumes rectangular MOF structure with pore centered on the origin and linkers
    aligned with either the x or y axis.
    
    @param pointList List of AtomicPoints graph(s) are being generated from.
    """
    
    #Length used to calculate what atoms are bonded
    bondDist = 2.6 #Angstroms
    #tLength = len(myGraph.myAdjList)
    
    #Index for vertex and then the edge list
    i = 0
    #Secondary index for the start list, whih grabs the C_3 carbon atoms
    j = 0
    
    myGraph.myVerticies = [None] * len(pointList)
    myGraph.numVerticies = 0
    startList = [None] * 4
    
    #Convert atomic points to graph verticies
    for p in pointList:
        new_vertex = ag.AtomicVertex(p)
        new_vertex.index = pointList.index(p)
        if p.bondType == "C_3":        
            new_vertex.quadrant = quadrantCheck(p)  
            startList[j] = new_vertex
            j += 1
        
        myGraph.myVerticies[i] = new_vertex
        myGraph.numVerticies += 1
        i += 1
    
    q = deque() #Initializes FIFO queue
    qLen = 0
    
    #add each vertex from start list to the queue
    for s in startList:
        q.append(s)
        qLen += 1

    #Runs as long as there are items in the queue
    while qLen != 0:
        tempV = q.popleft()
        qLen -= 1
        
        #checks all the avalible vertices
        for v in myGraph.myVerticies:
            #sees if any are within the listed bond distance
            #if tempV.ap.distBetween(v.ap, False) <= bondDist:
            if tempV.ap.distBetween(v.ap, True) <= 0:
                
                #If a zinc is reached, that is too far. At metal cluster. Skip to next item in queue.
                if v.ap.element == "ZN":
                    continue
                
                #If a vertex has not been encountered yet, it's quadrant is marked and the vertex added to the queue
                if v.quadrant == 0:
                    v.quadrant = tempV.quadrant
                    q.append(v)
                
                #checks if edge between the verticies exists, creates one if there isn't.
                if myGraph.edgeExist(tempV, v) == False:
                    myGraph.makeEdge(tempV, v)
    
    
def genBDCLinkerGraph(myGraph):
    """Function generates graph of BDC (terephthalic acid) based or related 
    linkers on a MOF with dicarboxylic acid connections to a metal center, just 
    zinc for now. Assumes rectangular MOF structure with pore centered on 
    the origin and linkers aligned with either the x or y axis.
    
    Currently performed after creating the cap graph. Assumes graph has already
    created from the atoms pdb frame
    
    
    Currrently untested
    
    
    @param myGraph graph of atoms describing the pore
    """
    #Could probably simplify code
    
    #Length used to calculate what atoms are bonded
    bondDist = 2.4 #Angstroms
    
    q = deque() #Initializes FIFO queue --> keeps track of graph starting locations --> don't need queue here now, could use list and for loop
    qLen = 0
    
    q2 = deque() #Initializes 2nd FIFO queue --> used for graph generation
    q2Len = 0
    
    #Queue all the 'O_R' atoms, they are the oxygens bonded to the zinc 
    #metal clusters. They are used as the base for linker graph generation
    for v in myGraph.myVerticies:
        if v.ap.bondType == 'O_R':
            q.append(v)
            qLen += 1
    
    #Runs as long as there are items in the queue
    while qLen != 0:
        tempV = q.popleft()
        qLen -= 1
        
        #Function only operates on vertices that are not yet part of a graph
        if tempV.quadrant != 0:
            continue
        
        #Vertex added to second queue for graph creation
        q2.append(tempV)
        q2Len += 1
        
        #Tracks O_R oxygens added to graph, if 4 encountered that means graph represent a linker, less than 4 means extraneous atoms
        ox_tracker = []
        
        #Now create graph
        while q2Len != 0:
            temp2V = q2.popleft()
            q2Len -= 1
            
            #9 is currently going to be the quadrant dummy value for and in calculation graph
            #9 as a quadrant label should not exist outside of this function running
            temp2V.quadrant = 9 
            
            #Tracks 'O_R' oxygens encounters, see note above where ox_tracker in initialized
            if temp2V.bondType == 'O_R':
                ox_tracker.append(temp2V)
            
            #checks all the avalible vertices
            for v in myGraph.myVerticies:
                #sees if any are within the listed bond distance
                if tempV.ap.distBetween(v.ap, False) <= bondDist:
                    #If a zinc is reached, that is too far. At metal cluster. Skip to next item in queue.
                    if v.ap.element == "ZN":
                        continue
                    
                    #Function only operates on vertices that are not yet part of the graph yet
                    if v.quadrant != 0:
                        continue
                    
                    #Since above if passed, the v hasn't yet been evaluated by the function yet
                    #Won't have edge to that point yet either, would have edge if point has been operated on
                    #previously
                    q2.append(v)
                    q2Len += 1
                    
                    #checks if edge between the verticies exists, creates one if there isn't.
                    if myGraph.edgeExist(tempV, v) == False:
                        myGraph.makeEdge(tempV, v)
        
        
        oxLen = len(ox_tracker)
        new_quadrant = -1
        #if 4 'O_R" oxygens encountered that means graph represent a linker, less than 4 means extraneous atoms
        if oxLen == 4:
            #Based on geometric locations of O_R atoms, determines quadrant of linker atoms
            new_quadrant = axisCheck(ox_tracker)
        
        #Updates new discrete graph to either have a linker quadrant or an invalid quadrant label of -1 for
        #incomplete portions of linkers that may exist at the edges of the simulation
        updateQuadrant(myGraph, new_quadrant, tempV)