#!/usr/bin/python2.7

"""
author : Oanh Vu
This script runs on pymol and color each atom based on new b-factors stored on a text file
File contains the spectrum values of each atom in order should be stored in the {ligand}_b.txt file
"""
from pymol import cmd
from pymol import stored
import numpy as np

def newB_mol (ligand, mi, ma):
  '''
    color the atom based on the statistic of the each atom
    
    ligand: name of the molecule object, file contains the spectrum values of each atom in order should be stored in the {ligand}_b.txt file
    
    mi, ma are min and max values for coloring
    '''

  #cmd.load( ligand + ".sdf")
  stored.count = 0
  cmd.iterate("all", "stored.count += 1")
  inFile = open(ligand+"_b.txt", 'r')
  stored.newB = []
  count = 0
  for line in inFile.readlines(): 
    stored.newB.append( float(line) )
    count = count+1
  inFile.close()
  #mi = min(stored.newB)
  #ma = max(stored.newB)
  for atom in range(1, count+1):
    #print atom
    cmd.alter("%s and id %s" %(ligand, atom), "b=%s"%stored.newB.pop(0))
  #cmd.save("%s_newB.sdf"%ligand , "%s"% ligand)
  #cmd.delete(ligand)
  #inFile.close()
  cmd.spectrum("b", "rainbow", ligand, mi, ma)
  #cmd.ramp_new("b", mol, [min, max], "rainbow")

cmd.extend("newB_mol", newB_mol)
