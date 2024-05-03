# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 14:48:31 2024

File for keeping track of MOF-5 Cap index placement

@author: samda
"""


def identifyMOF5Index(frame, index):
    '''
    If given an index from a properly bounded MOF-5 frame (-11 --> 11 in x,y, > 28 or 30? in z)
    returns what part of the molecular structure the index atom is a part of.
    
    Parameters
    ----------
    frame : int
        Frame index is from. Function only definitly works from frames 111-161 right now
    index : int
        Index of atom being described.

    Returns
    -------
    String description of atom location in MOF structure.

    '''

    if frame <= 121:
            
        #Frames 111-121

        Cap1_Stalk = ['Cap1_Stalk', 3, 27, 138, 151]
        Cap1_Leaf_A = ['Cap1_Leaf_A', 47, 113, 117, 118, 119, 120, 258, 259, 260, 261, 262]
        Cap1_Leaf_B = ['Cap1_Leaf_B', 46, 111, 112, 114, 115, 116, 253, 254, 255, 256, 257]
        Cap1_Leaf_C = ['Cap1_Leaf_C', 48, 110, 121, 122, 123, 124, 263, 264, 265, 266, 267]
        M_Cluster_1 = ['M_Cluster_1', 5, 9, 19, 126, 130, 134, 141, 148, 160, 168, 177, 181, 185, 191]

        Cap2_Stalk = ['Cap2_Stalk', 1, 25, 140, 152]
        Cap2_Leaf_A = ['Cap2_Leaf_A', 41, 83, 87, 88, 89, 90, 229, 230, 231, 232]
        Cap2_Leaf_B = ['Cap2_Leaf_B', 40, 81, 82, 84, 85, 86, 224, 225, 226, 227, 228]
        Cap2_Leaf_C = ['Cap2_Leaf_C', 42, 80, 91, 92, 93, 94, 233, 234, 235, 236, 237]
        M_Cluster_2 = ['M_Cluster_2', 6, 8, 17, 127, 131, 133, 143, 145, 158, 167, 180, 182, 187, 189]

        Cap3_Stalk = ['Cap3_Stalk', 0, 24, 139, 150]
        Cap3_Leaf_A = ['Cap3_Leaf_A', 39, 65, 76, 77, 78, 79, 219, 220, 221, 222, 223]
        Cap3_Leaf_B = ['Cap3_Leaf_B', 38, 68, 72, 73, 74, 75, 214, 215, 216, 217, 218]
        Cap3_Leaf_C = ['Cap3_Leaf_C', 37, 66, 67, 69, 70, 71, 209, 210, 211, 212, 213]
        M_Cluster_3 = ['M_Cluster_3', 4, 11, 16, 125, 129, 136, 142, 146, 157, 165, 178, 183, 188, 190]

        Cap4_Stalk = ['Cap4_Stalk', 2, 26, 137, 149]
        Cap4_Leaf_A = ['Cap4_Leaf_A', 45, 95, 106, 107, 108, 109, 248, 249, 250, 251, 252]
        Cap4_Leaf_B = ['Cap4_Leaf_B', 44, 98, 102, 103, 104, 105, 243, 244, 245, 246, 247]
        Cap4_Leaf_C = ['Cap4_Leaf_C', 43, 96, 97, 99, 100, 101, 238, 239, 240, 241, 242]
        M_Cluster_4 = ['M_Cluster_4', 7, 10, 18, 28, 128, 132, 135, 144, 147, 159, 166, 179, 184, 186, 192]

        Linker_5 = ['Linker_5', 20, 23, 33, 36, 57, 60, 62, 63, 169, 172, 174, 175, 201, 204, 206, 207]
        Linker_6 = ['Linker_6', 13, 14, 30, 31, 50, 51, 53, 56, 154, 155, 161, 164, 194, 195, 197, 200]
        Linker_7 = ['Linker_7', 21, 22, 34, 35, 58, 59, 61, 64, 170, 171, 173, 176, 202, 203, 205, 208]
        Linker_8 = ['Linker_8', 12, 15, 29, 32, 49, 52, 54, 55, 153, 156, 162, 163, 193, 196, 198, 199]
        
        #All lists for easy processing
        list_of_lists = [Cap1_Stalk, Cap1_Leaf_A, Cap1_Leaf_B, Cap1_Leaf_C, M_Cluster_1, Cap2_Stalk, Cap2_Leaf_A, Cap2_Leaf_B, Cap2_Leaf_C, M_Cluster_2, Cap3_Stalk, Cap3_Leaf_A, Cap3_Leaf_B, Cap3_Leaf_C, M_Cluster_3, Cap4_Stalk, Cap4_Leaf_A, Cap4_Leaf_B, Cap4_Leaf_C, M_Cluster_4, Linker_5, Linker_6, Linker_7, Linker_8]

        for item in list_of_lists:
            if index in item: return item[0]

    elif frame >= 122:

        #Frames 122-161

        Cap1_Stalk = ['Cap1_Stalk', 3, 27, 137, 150]
        Cap1_Leaf_A = ['Cap1_Leaf_A', 46, 112, 116, 117, 118, 119, 258, 259, 260, 261, 262]
        Cap1_Leaf_B = ['Cap1_Leaf_B', 45, 110, 111, 113, 114, 115, 253, 254, 255, 256, 257]
        Cap1_Leaf_C = ['Cap1_Leaf_C', 47, 109, 120, 121, 122, 123, 263, 264, 265, 266, 267]
        M_Cluster_1 = ['M_Cluster_1', 5, 9, 19, 125, 129, 133, 140, 147, 159, 167, 176, 180, 184, 190]

        Cap2_Stalk = ['Cap2_Stalk', 1, 25, 139, 151]
        Cap2_Leaf_A = ['Cap2_Leaf_A', 40, 82, 86, 87, 88, 89, 228, 229, 230, 231, 232]
        Cap2_Leaf_B = ['Cap2_Leaf_B', 39, 80, 81, 83, 84, 85, 223, 224, 225, 226, 227]
        Cap2_Leaf_C = ['Cap2_Leaf_C', 41, 79, 90, 91, 92, 93, 233, 234, 235, 236, 237]
        M_Cluster_2 = ['M_Cluster_2', 6, 8, 17, 126, 130, 132, 142, 144, 157, 166, 179, 181, 186, 188]

        Cap3_Stalk = ['Cap3_Stalk', 0, 24, 138, 149]
        Cap3_Leaf_A = ['Cap3_Leaf_A', 38, 64, 75, 76, 77, 78, 218, 219, 220, 221, 222]
        Cap3_Leaf_B = ['Cap3_Leaf_B', 37, 67, 71, 72, 73, 74, 213, 214, 215, 216, 217]
        Cap3_Leaf_C = ['Cap3_Leaf_C', 36, 65, 66, 68, 69, 70, 208, 209, 210, 211, 212]
        M_Cluster_3 = ['M_Cluster_3', 4, 11, 16, 124, 128, 135, 141, 145, 156, 164, 177, 182, 187, 189]

        Cap4_Stalk = ['Cap4_Stalk', 2, 26, 136, 148]
        Cap4_Leaf_A = ['Cap4_Leaf_A', 44, 94, 105, 106, 107, 108, 248, 249, 250, 251, 252]
        Cap4_Leaf_B = ['Cap4_Leaf_B', 43, 97, 101, 102, 103, 104, 243, 244, 245, 246, 247]
        Cap4_Leaf_C = ['Cap4_Leaf_C', 42, 95, 96, 98, 99, 100, 238, 239, 240, 241, 242]
        M_Cluster_4 = ['M_Cluster_4', 7, 10, 18, 127, 131, 134, 143, 146, 158, 165, 178, 183, 185, 191]

        Linker_5 = ['Linker_5', 20, 23, 32, 35, 56, 59, 61, 62, 168, 171, 173, 174, 200, 203, 205, 206]
        Linker_6 = ['Linker_6', 13, 14, 29, 30, 49, 50, 52, 55, 153, 154, 160, 163, 193, 194, 196, 199]
        Linker_7 = ['Linker_7', 21, 22, 33, 34, 57, 58, 60, 63, 169, 170, 172, 175, 201, 202, 204, 207]
        Linker_8 = ['Linker_8', 12, 15, 28, 31, 48, 51, 53, 54, 152, 155, 161, 162, 192, 195, 197, 198]
        
        #All lists for easy processing
        list_of_lists = [Cap1_Stalk, Cap1_Leaf_A, Cap1_Leaf_B, Cap1_Leaf_C, M_Cluster_1, Cap2_Stalk, Cap2_Leaf_A, Cap2_Leaf_B, Cap2_Leaf_C, M_Cluster_2, Cap3_Stalk, Cap3_Leaf_A, Cap3_Leaf_B, Cap3_Leaf_C, M_Cluster_3, Cap4_Stalk, Cap4_Leaf_A, Cap4_Leaf_B, Cap4_Leaf_C, M_Cluster_4, Linker_5, Linker_6, Linker_7, Linker_8]

        for item in list_of_lists:
            if index in item: return item[0]