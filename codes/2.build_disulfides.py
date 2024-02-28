#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Karen J. Gonzalez
"""
import os
import pandas as pd
import math
import argparse
import subprocess
from pyrosetta import *
from rosetta import *


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

def read_alignment(alignment):  #alignment in fasta format. Make sure the sequence alignment corresponds with the structural alignment
    with open(alignment, "r") as fopen:
        sequences = {}
        for line in fopen:
            line = line.strip()
            if line.startswith(">"):
                seq_name=line[1:]
                sequences[seq_name] = []
            else:
                for aa in line:
                    sequences[seq_name].append(aa)
    sequences = pd.DataFrame.from_dict(sequences)
    return sequences

def add_data_according_to_alignment(protein, alignment_dataframe, data_values, column_name="PDB_#" ):
    new = []
    counter=0
    for index, row in alignment_dataframe.iterrows():
        if alignment_dataframe.loc[index,protein] == "-":
            new.append("-") 
        else:
            new.append(data_values[counter])
            counter+=1                
    new_column = protein + "_"+column_name       
    alignment_dataframe[new_column] = new
    return alignment_dataframe

def sequence_from_pdb(pdb, chain="all", contracted=True):
    positions = []
    seq = []
    chains =[] 
    fopen = open(pdb,"r")
    old = ''
    for line in fopen:
        if line.startswith("ATOM"):
            new = int(line[22:26])
            if new != old:
                ch = line[21] #chain in the pdb
                if chain =="all":
                    chains.append(ch)
                else:
                    if not ch in chain:
                        continue
                    else:
                        chains.append(ch)
                aa = line[17:20]
                seq.append(aa)
                positions.append(new)
                old = new
    if contracted:
        seq = contractSeq(seq)
    return seq, positions, chains

def RMSD(vector1,vector2):
    euclidean_distance = ((vector1[0]-vector2[0])**2) + ((vector1[1]-vector2[1])**2) + ((vector1[2]-vector2[2])**2)
    rmsd = math.sqrt((euclidean_distance/2))
    return rmsd

def ca_coordinates(pdb, chain="A"):  
    data = []
    rosetta_position = 0
    with open(pdb,'r') as fopen:
        for line in fopen:
            if line.startswith("ATOM") and line[13:15] == "CA":
                if line[21] in chain:
                    position = int(line[22:26])
                    rosetta_position+=1
                    sequence = contractSeq([line[17:20]])  
                    data.append([float(line[30:38]), float(line[38:46]), float(line[46:54]), sequence[0], position, rosetta_position])     
        fopen.close()
    coordinates = pd.DataFrame(data, columns=["x","y","z","seq","PDB_#","rosetta_#"])
    coordinates.index = coordinates["rosetta_#"].tolist()
    return coordinates

def round_half_up(n, decimals=0):
    multiplier = 10 ** decimals
    return math.floor(n*multiplier + 0.5 )/multiplier

        
def potential_dislf(pose, ref_resi, distance_cutoff, min_loop):   
    potential_dislf= [] 
    ref_resi_name = pose.residue(ref_resi).name()
    if ref_resi_name != "GLY" and ref_resi_name != "PRO" and ref_resi_name != "CYS:disulfide":
        r1 = pose.residue(ref_resi).xyz("CB")
        total_res =  pose.total_residue()
        for i in range(1,total_res+1):
            resi_type = pose.residue(i).name() 
            if resi_type != "GLY" and resi_type != "PRO" and resi_type != "CYS:disulfide":   
                if ref_resi+min_loop < i or i< ref_resi-min_loop: #check residues that are at least 5 residues apart from the reference residue
                        r2 = pose.residue(i).xyz("CB")
                        distance = r1.distance(r2)
                        if distance <= distance_cutoff:         
                            potential_dislf.append(i)
    return potential_dislf
                                         
################################################
root = os.getcwd()
### Command line arguments 
parser = argparse.ArgumentParser(description='Building potential disulfide bonds')
parser.add_argument('arg_file', help = "File containing all input arguments")
parser.add_argument('pdb_model', help = "PDB model (output from relaxation)")

args = parser.parse_args()
# input arguments from arg_file
arg_dict = {}
with open(args.arg_file, "r") as fopen:
    for line in fopen:
        if line and line[0].isalpha():
            line = line.split("=")
            arg_dict[line[0].strip()] = line[1].strip()

### Obtain ca coordinates
pdb_post = arg_dict["pdb_post"]
monomer_chain_post = arg_dict["post_monomer_ch"]
coordinates_post = ca_coordinates(pdb_post, chain= monomer_chain_post)

pdb_pre_ref = arg_dict["pdb_pre_ref"]
monomer_chain_pre = arg_dict["pre_monomer_ch"]
coordinates_pre = ca_coordinates(pdb_pre_ref, chain= monomer_chain_pre)
# Only consider residue positions that are shared between pre- and postfusion structures
alignment = read_alignment(arg_dict["alignment"])
alignment = alignment.rename(columns={arg_dict["pre_name_alignment"]: "prefusion", arg_dict["post_name_alignment"]: "postfusion"})
alignment = add_data_according_to_alignment("prefusion", alignment, coordinates_pre.loc[:,"rosetta_#"].tolist(), column_name="rosetta_#")
alignment = add_data_according_to_alignment("postfusion", alignment, coordinates_post.loc[:,"rosetta_#"].tolist(), column_name="rosetta_#")
                                                   
idx_shared = []
for index, row in alignment.iterrows():
    if alignment.loc[index,"prefusion"] != "-" and alignment.loc[index,"postfusion"] != "-":
        idx_shared.append(index)
                                                    
# Trim coordinates to shared positions only
coordinates_pre_shared = coordinates_pre.loc[alignment.loc[idx_shared,"prefusion_rosetta_#"].tolist(), :]
coordinates_post_shared = coordinates_post.loc[alignment.loc[idx_shared,"postfusion_rosetta_#"].tolist(),:].set_index(pd.Index(coordinates_pre_shared.index.tolist()), inplace=False) # we need to have the two dataframes with the same indexes. The indexes of prefusion were chosen here  

### Calculate per residue root mean square deviation 
distance = []
for index, row in coordinates_pre_shared.iterrows():
    data1 = coordinates_pre_shared.loc[index, ["x","y","z"]].tolist()
    data2 = coordinates_post_shared.loc[index, ["x","y","z"]].tolist()
    rmsd = RMSD(data1,data2)
    distance.append(rmsd)
    
### Find positions that move more than ~10 Angs
mobile_pre = []
mobile_post = []

for i in range(len(distance)):
    if round_half_up(distance[i]) >= 10 :
        mobile_pre.append(coordinates_pre_shared.iloc[i,-1])
        mobile_post.append(coordinates_post_shared.iloc[i,-1])

### Add flanking regions
'''We have also included regions flanking the mobile areas where the secondary structure differs between the pre- and postfusion structures.'''
''' we search until finding a block of X (custom) amino acids matching their secondary strucures. '''

# Identify secondary structure of input PDBs
init() #initialize pyrosetta
pose_pre_ref = pose_from_pdb(pdb_pre_ref)
pose_post = pose_from_pdb(pdb_post)
dssp = pyrosetta.rosetta.protocols.moves.DsspMover()
dssp.apply(pose_pre_ref)    
dssp.apply(pose_post)
ss_pre = pose_pre_ref.secstruct() #secondary structure prediction for the prefusion state
ss_post = pose_post.secstruct() #secondary structure prediction for the postfusion state

# Add secondary structure information to aligned sequences
alignment = add_data_according_to_alignment("prefusion", alignment, ss_pre, column_name="ss")
alignment = add_data_according_to_alignment("postfusion", alignment, ss_post, column_name="ss")

# Restrict the analysis to shared regions
alignment_shared = alignment.loc[idx_shared,:]

# Find flanking regions 
flanking_pre = []
flanking_post  = []
resi_flanking = int(arg_dict["flanking"])
for i in range(len(mobile_pre)):
    if i+1 < len(mobile_pre):
        if mobile_pre[i] +1 != mobile_pre[i+1]  : #find where the mobile area is discontinuos
            #analyze upstream neighbors
            neighbor_up = mobile_pre[i] +1     
            if neighbor_up in alignment_shared['prefusion_rosetta_#'].tolist():
                idx_neighbor_up = alignment_shared[alignment_shared['prefusion_rosetta_#'] == neighbor_up].index.tolist()[0]
                for idx in range(idx_neighbor_up, alignment_shared.index[-1]+1):
                    if alignment_shared.loc[idx:idx+resi_flanking, "prefusion_ss"].tolist() == alignment_shared.loc[idx:idx+resi_flanking, "postfusion_ss"].tolist():
                        flanking_pre += alignment_shared.loc[idx_neighbor_up:idx+resi_flanking, 'prefusion_rosetta_#'].tolist()
                        flanking_post += alignment_shared.loc[idx_neighbor_up:idx+resi_flanking, 'postfusion_rosetta_#'].tolist()                                                  
                        break
                              
# Update mobile regions with flanking residues                                 
mobile_pre += flanking_pre ; mobile_post += flanking_post        

# Include the fusion peptide residues 
fusion_peptide = arg_dict["fusion_peptide"].split("-")
fp_idx_start = coordinates_pre[coordinates_pre["PDB_#"] == int(fusion_peptide[0])].index.tolist()[0]
fp_idx_end = coordinates_pre[coordinates_pre["PDB_#"] == int(fusion_peptide[1])].index.tolist()[0]
fusion_peptide = coordinates_pre.loc[fp_idx_start:fp_idx_end,"rosetta_#"].tolist()
mobile_pre+= fusion_peptide

mobile_pre = list(set(mobile_pre)); mobile_pre.sort()
mobile_post = list(set(mobile_post)); mobile_post.sort()

### Make disulfide bond and output pdb for relaxation
out_path = os.path.join(root,"disulfide_pdbs")
subprocess.call(["mkdir",out_path]) 
# Load the default ref2015 score function
scorefxn = get_fa_scorefxn()
sfxn_disulfonly = ScoreFunction()
sfxn_disulfonly.set_weight(core.scoring.ScoreType.dslf_fa13, 1.25)  # Full-atom disulfide energy term


model = args.pdb_model
pose_model = pose_from_pdb(model)
#let's analyze only the monomeric structure
pose_model_monomer = pose_model.split_by_chain(1)


for ref_resi in mobile_pre:
    # Find potential disulfides based on Cb-Cb distance
    disulf = potential_dislf(pose_model_monomer, ref_resi, 6, 8) 
    for res2 in disulf:
        pose_out = pose_model_monomer.clone()
        #Build disulfides
        disulfidize_mover = protocols.denovo_design.DisulfidizeMover() 
        disulfidize_mover.set_match_rt_limit(2)
        disulfidize_mover.set_max_disulf_score(3.5)
        disulfidize_mover.make_disulfide(pose_out, ref_resi, res2, True, scorefxn)
        #check disulfide match_rt
        disulf_potential = core.scoring.disulfides.DisulfideMatchingPotential()
        match_rt = disulfidize_mover.check_disulfide_match_rt(pose_out, ref_resi, res2,disulf_potential, False)
        #check disulfide score
        ss_score = disulfidize_mover.check_disulfide_score(pose_out, ref_resi, res2, sfxn_disulfonly, scorefxn)
  
        if match_rt and ss_score:
            #save pdb in output folder
            if os.path.dirname(model) != '':
                name = model.split("/")[-1]
            else:
                name = model
            out_name = name[:-4]+"_"+str(ref_resi)+"_"+str(res2)+".pdb"
            out_path2 = os.path.join(out_path, out_name)   
            pose_out.dump_pdb(out_path2)        

     





































