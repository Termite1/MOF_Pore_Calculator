# -*- coding: utf-8 -*-
"""
Created on Fri May  3 11:33:52 2024

@author: samda
"""

import LargestCircleFinder as lcf
import pdbPoreReader as pdb

file = "MOF-5 middle pore frames 111-161.pdb"

all_frames = pdb.multiframePDBReader(file)

#Base frame 111. Frame 136
current_frame = all_frames[25]

#Atom indexes for collision error
a1 = current_frame[199]
b1 = current_frame[254]
c1 = current_frame[103]

#Atom indexes for sourround error (small)
a2 = current_frame[120]
b2 = current_frame[114]
c2 = current_frame[110]

#Answer 1
a3 = current_frame[199]
b3 = current_frame[56]
c3 = current_frame[246]

#Answer 2
a4 = current_frame[199]
b4 = current_frame[254]
c4 = current_frame[246]

#Answer 3
a5 = current_frame[254]
b5 = current_frame[59]
c5 = current_frame[246]



s1 = lcf.threeAtomCircumcenter(a1,b1,c1)
s2 = lcf.threeAtomCircumcenter(a2,b2,c2)
s3 = lcf.threeAtomCircumcenter(a3,b3,c3)
s4 = lcf.threeAtomCircumcenter(a4,b4,c4)
s5 = lcf.threeAtomCircumcenter(a5,b5,c5)


write_list = [s1, s2, s3, s4, s5]

pdb.pdbWrite(write_list, 1)

