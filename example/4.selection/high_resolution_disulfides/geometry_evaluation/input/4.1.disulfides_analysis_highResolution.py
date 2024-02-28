#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Karen J. Gonzalez
"""

import os
import pandas as pd
from pyrosetta import *
from rosetta import *
from rosetta.core.scoring import *
import argparse
import dill

######################## Functions ########################

def expandSequence(sequence): #Sequence can be a list or a string
    seqExtended = []
    for i in sequence:
        if i == "R":
            seqExtended.append("ARG")
        elif i == "H":
            seqExtended.append("HIS")    
        elif i == "K":
            seqExtended.append("LYS")            
        elif i == "D":
            seqExtended.append("ASP")            
        elif i == "E":
            seqExtended.append("GLU")            
        elif i == "S":
            seqExtended.append("SER")            
        elif i == "T":
            seqExtended.append("THR")            
        elif i == "N":
            seqExtended.append("ASN")    
        elif i == "Q":
            seqExtended.append("GLN")    
        elif i == "C":
            seqExtended.append("CYS")    
        elif i == "G":
            seqExtended.append("GLY")            
        elif i == "P":
            seqExtended.append("PRO")            
        elif i == "A":
            seqExtended.append("ALA")            
        elif i == "V":
            seqExtended.append("VAL")            
        elif i == "I":
            seqExtended.append("ILE")            
        elif i == "L":
            seqExtended.append("LEU")    
        elif i == "M":
            seqExtended.append("MET")    
        elif i == "F":
            seqExtended.append("PHE")    
        elif i == "Y":
            seqExtended.append("TYR")
        elif i == "W":
            seqExtended.append("TRP") 
    return seqExtended

def contractSeq(sequence): #sequence has to be a list with each amino acid as an independent element
    seqExtended = []
    for i in sequence:
        if i == "ARG":
            seqExtended.append("R")
        elif i == "HIS":
            seqExtended.append("H")    
        elif i == "LYS":
            seqExtended.append("K")            
        elif i == "ASP":
            seqExtended.append("D")            
        elif i == "GLU":
            seqExtended.append("E")            
        elif i == "SER":
            seqExtended.append("S")            
        elif i == "THR":
            seqExtended.append("T")            
        elif i == "ASN":
            seqExtended.append("N")    
        elif i == "GLN":
            seqExtended.append("Q")    
        elif i == "CYS":
            seqExtended.append("C")    
        elif i == "GLY":
            seqExtended.append("G")            
        elif i == "PRO":
            seqExtended.append("P")            
        elif i == "ALA":
            seqExtended.append("A")            
        elif i == "VAL":
            seqExtended.append("V")            
        elif i == "ILE":
            seqExtended.append("I")            
        elif i == "LEU":
            seqExtended.append("L")    
        elif i == "MET":
            seqExtended.append("M")    
        elif i == "PHE":
            seqExtended.append("F")    
        elif i == "TYR":
            seqExtended.append("Y")
        elif i == "TRP":
            seqExtended.append("W") 
    return seqExtended 

def find_disulfides(pose):
    disulfides = []
    disulf_vector = utility.vector1_std_pair_unsigned_long_unsigned_long_t() # Create an empty vector for disulfide bonds
    core.conformation.disulfide_bonds(pose.conformation(), disulf_vector) # Detect disulfide bonds in starting structure
    for res in disulf_vector:
        disulfides.append([res[0], res[1]]) #the disulfides are in rosetta numbering
    return disulfides

################################################################################
### Command line arguments 
parser = argparse.ArgumentParser(description='Geometry analysis of high-resolution dataset')
parser.add_argument('pdbs_directory', help = "Directory hosting the entire crystallized high-resolution dataset (relaxed)")
args = parser.parse_args()
### Initiallize pyrosetta
init()  
scorefxn = get_fa_scorefxn()

### Gather high-resolution PDBs
pdbs_dir = args.pdbs_directory
pdbs = [os.path.join(pdbs_dir, i) for i in os.listdir(pdbs_dir) if i.endswith(".pdb")]

### Find disulfides, calculate geometry and energy scores
columns = ["SS_rosetta","fa13", "dihedral_x1","dihedral_x2" ,"dihedral_x3","angle_cb1_s1_s2", "dist_cb1_ca2", "description"]            
data = pd.DataFrame(index= list(range(len(pdbs))), columns=columns)
for i in range(len(pdbs)):
    pdb = pdbs[i]
    pdb_pose =  pose_from_pdb(pdb)
    disulf = find_disulfides(pdb_pose)
    for dslf in disulf:  
        data.loc[i,"description"] = pdb           
        res1 = pdb_pose.pdb_info().pose2pdb(dslf[0])
        res2 = pdb_pose.pdb_info().pose2pdb(dslf[1])        
        data.loc[i,"SS_rosetta"]= str(dslf[0])+"_"+str(dslf[1])
        scorefxn(pdb_pose)
        data.loc[i,"fa13"] = pdb_pose.energies().residue_total_energies(dslf[0])[dslf_fa13]
        ##### calculating bond geometry
        residue = pdb_pose.residue(dslf[0])
        partner = pdb_pose.residue(dslf[1])
        ## atom coordinates
        s1_xyz = residue.xyz("SG")
        s2_xyz = partner.xyz("SG")
        cb1_xyz = residue.xyz("CB")    
        cb2_xyz = partner.xyz("CB")
        ca1_xyz = residue.xyz("CA")
        ca2_xyz = partner.xyz("CA")
        n1_xyz = residue.xyz("N") 
        ###### calculating dihedral angles (torsion angles)
        # dihedral1  = N-Ca-Cb-Sg
        # dihedral2 = Ca-Cb-Sg-Sg
        # dihedral3 = Cb-Sg-Sg-Cb      
        dihedral_x1 = rosetta.numeric.dihedral(n1_xyz, ca1_xyz,cb1_xyz, s1_xyz)
        dihedral_x2 = rosetta.numeric.dihedral(ca1_xyz, cb1_xyz,s1_xyz, s2_xyz)
        dihedral_x3 = rosetta.numeric.dihedral(cb1_xyz, s1_xyz,s2_xyz, cb2_xyz)         
        data.loc[i,"dihedral_x1"] =dihedral_x1
        data.loc[i,"dihedral_x2"] =dihedral_x2
        data.loc[i,"dihedral_x3"] =dihedral_x3              
        # angle_cb1_s1-s2"
        angle_cb1_s1_s2 = rosetta.numeric.angle_degrees(cb1_xyz,s1_xyz,s2_xyz)
        data.loc[i,"angle_cb1_s1_s2"] =angle_cb1_s1_s2       
        # distance cb1_ca2"
        dist_cb1_ca2 = (cb1_xyz - ca2_xyz).norm()
        data.loc[i,"dist_cb1_ca2"] = dist_cb1_ca2
        
 
# Save data
data.dropna(axis=0, how="any", inplace=True)
with open('geometry_high_resolution.pkl', 'wb') as file:
    dill.dump({'data_high_resolution': data}, file)



        
        
            
            
            
    









