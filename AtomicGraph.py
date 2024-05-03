# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 11:53:10 2023

Code for creating a graph ot atomic points.

@author: samda
"""

import AtomicPoint


def expandList(l):
    """Function takes a list and doubles it in size, adding the items from the old list ot the front of the new list

    @param l List to be expanded
    """
    new_length = len(l) * 2
    new_list = [None] * new_length
    i = 0
    
    #Adds all the items from the old list to the new list
    for item in l:
        new_list[i] = l
        i += 1
    
    return new_list


class AtomicVertex:
    
    def __init__(self, a):
        """Constructor for AtomicVertex class.
        
        @param a AtomicPoint
        """
        #AtomicPoint
        self.ap = a
        self.distance = 10000 #Should make larger if I somehow have really large graphs. Values should be ~infinity
        self.predecessor = None
        self.quadrant = 0
        self.index = -1
    
    
class AtomicGraph:
    
    def __init__(self):
        """Constructor for AtomicGraph class.
            
        @param pointList List of AtomicPoints that will be turned into the graph 
        """
        self.myVerticies = [None]
        self.numVerticies = 0
        self.myAdjList = [None] * 4
        self.numEdges = 0
    
    
    def isFull(self, l):
        """Function checks if a given list is full (no None values present). Returns True if so, False otherwise.
        
        @param
        @return True if list has no None values, False otherwise.
        """
        isFull = False
        
        #Checks every item in the given list
        for item in l:
            if l == None:
                isFull = True
        
        return isFull
    
    
    def makeEdge(self, v1, v2):
        """Function makes an edge between two given
        verticies and records them in the myAdjList
        
        @param v1 AtomicVertex1
        @param v2 AtomicVertex2
        """
        new_edge = (v1, v2)
        
        #Expands myAdjList if it is full
        if self.isFull(self.myAdjList):
            expandList(self.myAdjList)
        
        #Adds new_edge to list at the first empty index encountered
        for item in self.myAdjList:
            if item == None:
                self.numEdges += 1
                item = new_edge
                break
  
  
    def edgeExist(self, v1, v2):
        """Function checks if there is already an edge between v1 and v2.
        Returns True if there is, False otherwise.
        
        @param v1 AtomicVertex1
        @param v2 Atomicvertex2
        @return True if an edge between the verticies exists, False otherwise
        """
        #run through every edge in the list, checks if both vertices are present in the edge
        for e in self.myAdjList:
            if e == None:
                #continue
                return False #Items added sequentially, so this should be fine
                             #no need to go through all the None values
            elif (v1 == e[0] or v1 == e[1]) and (v2 == e[0] or v2 == e[1]):
                return True
        return False
    
    
    def edgeContains(self, v, edge):
        """Function checks if the given edge contains vertex v.
        Returns True if it does, False otherwise.
        
        @param v AtomicVertex
        @param edge as a (v1, v2) tuple
        @return True if given edge contains vertex v, False otherwise
        """
        if edge[0] == v or edge[1] == v: return True
        else: return False
        
        
    def pairReturn(self, v, edge):
        """Function returns counterpart to input vector forming given edge.
        
        @param v AtomicVertex
        @param edge as a (v1, v2) tuple
        @return v counterpart
        """
        if self.edgeContains(v, edge) == False:
            print("Error: Edge does not contain given vertex!")
            return None #Return if edge somehow does not contain v
        
        if edge[0] == v: return edge[1]
        elif edge[1] == v: return edge[0]
        
        
    '''
    def genFullGraph(self, pointList):
        """
        ???
        """
        
        #Length used to calculte what atoms are bonded
        bondDist = 2 #Angstroms
        tLength = len(self.myAdjList)
        
        i = 0
        
        self.myVerticies = [None] * len(pointList)
        
        #Convert atomic points to graph verticies
        for p in pointList:
            self.myVerticies[i] = AtomicVertex(p)
            self.numVerticies += 1
            i += 1
        
        #Generates adjacency list
        for v in self.myVerticies:
            for v2 in self.myVerticies:
                #check distance and create tuple pair of verticies
                if v.ap.distBetween(v2.ap, False) < bondDist:
                    #Grows list when full by doubleing the length
                    if i == tLength:
                        expandList(self.myAdjList)
                    #check if edge already exists
                    #point where ut wiuld be useful to have hash function for edge definitions?
                    #...
                    
                    self.myAdjList[i] = (v, v2)
                    i += 0
    '''
    

    
                
                
                
            
                
            
                