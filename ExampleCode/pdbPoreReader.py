# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 20:52:14 2023

This program reads a pdb file and turns it into array of (x,y) points
representing the perimeter of the pore.

example pdb file is arranged as such
First Line: ??? (Skip?)

ATOM # atom_type [amino_acid?] chain_name sequence_number X Y Z ??? ??? Element

Each field is seperated (tab delimited), so should be able to grab the
info I want easily

Problem: Only atom locations are shown, not bonds. Either the file will
need to be interprited for bond locations, or bonds will be effectivly
neglected. This may be fine, hydrogen atom locations result in very
spiky structure.

Finding the perimeter atoms is really going to be the challenging part
since there is no clear indication of what atoms are connected to what.

Should probably bake in a bounds function that causes the program to only
read points between some set of z values, probably an overload

@author: samda
"""
from AtomicPoint import AtomicPoint
    

def atomReaderPDB (pdb_input):
    """
    Function reads a pdb file and turns it into a list of AtomicPoints

    @param: The function takes the name of a pdb file to process
    @return: Function returns a list of Point2D contianing all the (x,y) 
             coordinates of the atoms in the pdb file
    """
    
    with open(pdb_input) as f:
        lines = f.readlines()
        num_entries = len(lines)
    
        atomicArray = [None] * (num_entries - 2)
    
        #For loop skips the first line then reads to the second to last line
        for a in range(1, num_entries - 1):
            temp_line = lines[a].split() 
            
            #5 is x coord, 6 y coord, 7 z coord
            atomicArray[a-1] = AtomicPoint(float(temp_line[5]), #x coord
                                           float(temp_line[6]), #y coord
                                           float(temp_line[7]), #z coord
                                           temp_line[10],       #element
                                           bondType=temp_line[2]) #bond type      

        return atomicArray  
    

def multiframePDBReader(pdb_input):
    """
    Function reads a pdb file with multiple frames and turns it into a list of lists of AtomicPoints

    @param pdb_input The function takes the name of a pdb file to process
    @return Function returns a list of lists of AtomicPoints contianing all the coordinates of the atoms in the pdb file.
        Each list represents a frame of the pdf file
    """  

    with open(pdb_input) as f:
        frames_list = []
        
        allLines = f.readlines()
        allLines.pop(0) #first line never contains data
    
        #Holds data on current frame
        currentAtomicList = []
    
        for line in allLines:
            temp_line = line.split()
            #checks if current line is a split between frames. Indicated by "END" characters, which is only 1 item.
            #if so, saves current frame, the starts new list for new frame
            if len(temp_line) < 11:
                frames_list.append(currentAtomicList)
                currentAtomicList = []
                continue
            #Checks if the end of the pdb file has been reached. 
            #Some files have additional lines afte the last "END", this catches that
            #Does not save a frame if there is no "END" for it
            elif line == allLines[-1]:
                return frames_list 
            
            #5 is x coord, 6 y coord, 7 z coord, 10 element, 2 varient of element by bond type
            currentAtomicList.append(AtomicPoint(float(temp_line[5]), #x coord
                                           float(temp_line[6]), #y coord
                                           float(temp_line[7]), #z coord
                                           temp_line[10],       #element
                                           bondType=temp_line[2])) #bond type
        
        #Returns list of lists of "END" is on the last line of the file
        return frames_list 


def pdbWrite(pointList, n):
    """
    Function takes a list of points and prints them out in a format suitble to be 
    copied directly into a pdb file.

    @param pointList List of AtomicPoints that will be processed and printed
    @param n Integer initial index of points printed.
    """
    
    #incrementor for index
    i = 0
    
    #repeated for every point in the point list
    for p in pointList:
        j = i + n
        
        #Argument 1
        ATOM = "ATOM"
        
        #Argument 2
        index = str(j)
        
        #Space 1
        l1 = 7 - len(index)
        s1 = " " * l1
        
        #Argument 3 -- Having errors
        bondType = p.bondType
        
        #Space 2 & 3
        l2 = len(bondType)
        if l2 == 2:
            s2 = " " * 2
            s3 = " " * 6
        elif l2 == 3:
            s2 = " " * 2
            s3 = " " * 5
        elif l2 == 4:
            s2 = " "
            s3 = " " * 5
        #don't know if there are bond types longer than 4
        
        #Argument 4
        a4 = "X"
        
        #Space 4
        s4 = " " * 3
        
        #Argument 5
        a5 = "1"
        
        #Argument 6
        X = "%.3f" % p.x
        
        #Space 5
        l5 = 12 - len(X)
        s5 = " " * l5
        
        #Argument 7
        Y = "%.3f" % p.y
        
        #Space 6
        l6 = 8 - len(Y)
        s6 = " " * l6
        
        #Argument 8
        Z = "%.3f" % p.z
        
        #Space 7
        l7 = 8 - len(Z)
        s7 = " " * l7
        
        #Space 8
        s8 = " " * 2
        
        #Argument 9
        a9 = "1.00"
        
        #Space 9
        s9 = " " * 2
        
        #Argument 10
        a10 = "0.00"
        
        #Argument 11
        if p.element == '':
            element = "H"
        else:
            element = p.element
            
        #Space 10
        l10 = 12 - len(element)
        s10 = " " * l10
        
        #Assemble Line
        line = (ATOM, s1, index, s2, bondType, s3, a4, s4, a5, s5, X, s6, Y, s7, Z, s8, a9, s9, a10, s10, element)
        pdbLine = ''.join(line)
        
        #Print Line (Hopefully newline for each)
        print(pdbLine)
        
        #increment index
        i += 1
        
        