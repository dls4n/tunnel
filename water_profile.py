#!/usr/bin/env python
#
#
# water_profile: identification of water molecules inside protein cavity determined by CAVER Analyst 2.0
#    Copyright (C) 2025  David Stokes 
#           stokes@nyu.edu
#           Dept. of Biochemistry and Molecular Pharmacology
#           NYU School of Medicine
#           550 First Ave
#           New York, NY 10016   USA
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#


# usage: water_profile.py tunnel.pdb molecule.pdb

#
#
#

from __future__ import print_function  # to enable python2.7
import platform
import warnings
import sys
import os
import math
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO

warnings.filterwarnings('ignore')    # ignore all warnings
#warnings.filterwarnings('ignore', message="Invalid or missing occupancy")

print("water_profile running under python vers: %s" % platform.python_version())
if len(sys.argv) < 2:
    print ('''         water_profile.py tunnel.pdb molecule.pdb
   This program will determine locations of water molecules lying inside a CAVER cavity''')

    try: tunnelpdb = raw_input('tunnel profile PDB file: ')   # this is for python2
    except NameError: tunnelpdb = input('tunnel profile PDB file: ')   # this is for python3
else:
    tunnelpdb = sys.argv[1]

if len(sys.argv) < 3:
    try: moleculepdb = raw_input('molecule PDB file: ')   # this is for python2
    except NameError: moleculepdb = input('molecule PDB file: ')   # this is for python3
else:
    moleculepdb = sys.argv[2]

# Process caver tunnel PDB file so it can be parsed by Bio.PDB
tunnelatomspdb = 'tunnel_conv.pdb'
f = open(tunnelpdb,'r')
o = open(tunnelatomspdb,'w')   # temporary file with only ATOM records

indat = f.readlines()
resnew = 0   # base residue number - to be incremented
for l in indat:
    if 'ATOM' in l:
        resnum = l[22:26]
        if resnew == 0: resnew = int(resnum)
        else: resnew += 1   # increment residue
        oline = l.replace(resnum,"%s" % resnew)
        o.write(oline)

f.close()
o.close()

#dummy = input('finished conversion: ')

# create parser to read PDB files
parser = PDBParser()

# read tunnel atoms from PDB
tstruct = parser.get_structure('caver', tunnelatomspdb)
tmodel = tstruct[0]
tchain = tmodel["T"]    # caver uses chain T

#dummy = input('finished reading tunnel file')

# create dictionary with tunnel atom coordinates
tunnel_dict = {}
running_dist = 0.0    # running distance along the tunnel profile starting at first point
last_coord = (0.0, 0.0, 0.0)   #coord of last point along profile
first_point = True
for residue in tchain:
    resnum = residue.get_id()[1]    # residue number (int)
    for atom in residue:  # there should only be one 'H' atom
        if atom.get_name() == 'H':
            coord = atom.get_coord()     # tuple
            occup = atom.get_bfactor()    # radius of tunnel
            if first_point:  # this is the first point
                dist = 0.0
                last_coord = coord
                first_point = False
            else:
                dist += math.dist(coord,last_coord)
                last_coord = coord
            tunnel_dict[resnum] = [resnum, coord, dist, occup]
        else:
            print ('Encountered unexpected ATOM in tunnel PDB, e.g., not "H"')
            sys.exit(0)

#dummy = input('finished processing tunnel file')

# read molecule from PDB
mstruct = parser.get_structure('molecule', moleculepdb)
mmodel = mstruct[0]

#dummy = input('finished reading molecule PDB')

water_dict = {}
for chain in mmodel:         # loop over entire molecule looking for waters
    for residue in chain:
        if residue.get_resname() == 'HOH':
            # look for closest tunnel atom for dist calculation
            ydist = 1000.   # to find closest tunnel atom
            inside = False    # to signal whether inside the tunnel cavity
      # inside tunnel depends on being within one of the spheres defining the tunnel
      # this may not be the closest sphere
            for tatom in tunnel_dict:
                water_coord = residue['O'].get_coord()
                dist = math.dist(water_coord, tunnel_dict[tatom][1])
                if dist <= tunnel_dict[tatom][3]: inside = True  # water inside tunnel
                if dist < ydist:   # this loop to record closest tunnel atom
                    ydist = dist   # distance from tunnel profile
                    xdist = tunnel_dict[tatom][2]  # coordinate along tunnel profile
                    tres = tunnel_dict[tatom][0]  # index of tunnel atom along profile
                    hohchain = residue.get_full_id()[2]  # chain for water molecule
                    hohres = residue.get_id()[1]    # resnum for water molecule
            if inside:    # record coordinates
                water_dict[hohres] = [xdist,ydist,hohres,hohchain,tres,water_coord]

# write results
r = open('tunnel_waters.dat','w')
print ("    x            y         hoh #  chain    tunnel #")
r.write("tunnel: %s\n" % tunnelpdb)
r.write("model: %s\n" % moleculepdb)
r.write ("    x            y         hoh #  chain    tunnel #\n")
for key,value in water_dict.items():
#    print (value)
    print ("%10.3f %10.3f %10d %5s %10d" % (value[0],value[1],value[2],value[3],value[4]))
    r.write("%10.3f %10.3f %10d %5s %10d\n" % (value[0],value[1],value[2],value[3],value[4]))
r.close()
# write tunnel profile distances
t = open('tunnel_profile.dat','w')
t.write("    res #      dist\n")
for key,value in tunnel_dict.items():
    t.write("%10d %10.3f\n" % (value[0],value[2]))
t.close()

# calculate interwater distances
water_distances = []
for key1,water1 in water_dict.items():
    for key2,water2 in water_dict.items():
        if water1[2] != water2[2] and water1[3] == water2[3]:  # same chain only
            dist = math.dist(water1[5],water2[5])
            if dist < 4: water_distances.append( [water1[2],water2[2],dist,water1[3]] )

# write out inter-water distances
dout = open('tunnel_water_dist.dat','w')
dout.write ("    w1        w2       dist      chain\n")
for d in water_distances:
#    print ("%10d %10d %10.3f %10s" % (d[0],d[1],d[2],d[3]) )
    dout.write("%10d %10d %10.3f %10s\n" % (d[0],d[1],d[2],d[3]) )
dout.close()


#os.system('rm -f %s' % tunnelatomspdb)    # clean up temporary tunnel PDB
            
            





