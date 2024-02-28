#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Karen J. Gonzalez
"""
import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyrosetta import *
from rosetta import *
from rosetta.core.scoring import *
import subprocess
from rosetta.core.scoring.sasa import SasaCalc
import argparse
from Bio.PDB import PDBParser, Superimposer
import dill
import json
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['Arial']

######################## Functions ########################

def expandSequence(sequence): # Sequence can be a list or a string
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

def contractSeq(sequence): # Sequence has to be a list with each amino acid as an independent element
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
        disulfides.append([res[0], res[1]]) #Disulfides are in Rosetta numbering
    return disulfides

def rms_CA(pdb, ref_pdb, model= 0): 
    # Load the reference and designed PDBs 
    p = PDBParser()
    reference_pdb = p.get_structure("reference", ref_pdb)
    designed_pdb = p.get_structure("designed", pdb)   
    # Extract first chain from reference structures
    ref_monomer= next(reference_pdb[model].get_chains(), None) 
    # Get the list of CA atoms from both reference and designed PDBs
    ca_atoms_ref = [atom for atom in ref_monomer.get_atoms() if atom.name == 'CA']
    ca_atoms_design = [atom for atom in designed_pdb.get_atoms() if atom.name == 'CA']
    # Make sure that both structures have the same number of CA atoms
    assert len(ca_atoms_ref) == len(ca_atoms_design), "Different number of CA atoms"
    # Align the structures
    imposer = Superimposer()
    imposer.set_atoms(ca_atoms_ref, ca_atoms_design)    
    # Apply the rotation/translation matrix to the second structure
    imposer.apply(designed_pdb.get_atoms())
    return imposer.rms

def dslf_geometry(dslf): # dslf is a list containing only the two residues forming the disulfide (Rosetta numbering)
    ##### Calculating bond geometry
    residue = pdb_pose.residue(int(dslf[0]))
    partner = pdb_pose.residue(int(dslf[1]))
    ## Atom coordinates
    s1_xyz = residue.xyz("SG")
    s2_xyz = partner.xyz("SG")
    cb1_xyz = residue.xyz("CB")    
    cb2_xyz = partner.xyz("CB")
    ca1_xyz = residue.xyz("CA")
    ca2_xyz = partner.xyz("CA")
    n1_xyz = residue.xyz("N") 
    ###### Calculating dihedral angles (torsion angles)
    # dihedral1  = N-Ca-Cb-Sg
    # dihedral2 = Ca-Cb-Sg-Sg
    # dihedral3 = Cb-Sg-Sg-Cb 
    dihedral_x1 = rosetta.numeric.dihedral(n1_xyz, ca1_xyz,cb1_xyz, s1_xyz)
    dihedral_x2 = rosetta.numeric.dihedral(ca1_xyz, cb1_xyz,s1_xyz, s2_xyz)
    dihedral_x3 = rosetta.numeric.dihedral(cb1_xyz, s1_xyz,s2_xyz, cb2_xyz)                    
    angle_cb1_s1_s2 = rosetta.numeric.angle_degrees(cb1_xyz,s1_xyz,s2_xyz)
    dist_cb1_ca2 = (cb1_xyz - ca2_xyz).norm() 
    return dihedral_x1, dihedral_x2, dihedral_x3, angle_cb1_s1_s2, dist_cb1_ca2

def disulf_relabeling(labels):
    numbering = []
    for idx in labels:
        idx = idx.split("_")
        # Reorganize the labeling so that the smaller residue number is at the beginning 
        if int(idx[0]) < int(idx[1]):
            numb = idx[0]+"_"+idx[1]
        else:
            numb = idx[1]+"_"+ idx[0]                     
        numbering.append(numb)
    return numbering    

def remove_outliers(df, column):
    Q1 = df[column].quantile(0.25)
    Q3 = df[column].quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR
    return df[(df[column] >= lower_bound) & (df[column] <= upper_bound)]

def number_bins(data):
    data.sort()
    n = len(data)
    range_bins = max(data) - min(data)
    number_intervals =  math.sqrt(n) 
    width =  range_bins /number_intervals
    bins = []
    limit = min(data)
    for i in range(round(number_intervals)):
        bins.append(limit)
        limit += width
    bins.append(limit)
    return bins

def allowed_values(value, bin_edges, counts):
    for i in range(len(bin_edges) - 1):
        if counts[i+1] > 0:
            if bin_edges[i] <= value < bin_edges[i+1]:
                return True
    return False            

def customized_histogram(disulfide_data, score_term,  high_resol_data, bins, title, xticks = np.arange(-0.7,1, step=0.1), xlim_right=0, xlim_left=0, figure_size=(2,2) ):    
    high_resol_data = high_resol_data.loc[:,score_term].tolist()
    counts, edges = np.histogram(high_resol_data, bins=bins) 

    data_x = disulfide_data.loc[:, score_term].tolist()
    data_y = disulfide_data.loc[:, "SS_rosetta"]
      
    # Plot designed disulfides scores
    fig, ax1 = plt.subplots(figsize=figure_size)
    # Make the background first
    min_value = min(counts)
    max_value = max(counts)
    normalized = [(i-min_value)/(max_value-min_value) for i in counts]     
    for i in range(0, len(edges)-1):    
        x = (edges[i] + edges[i+1]) / 2
        width = abs(edges[i] - edges[i+1])
        transparency = normalized[i]  # The transparency will correspond to the frequency (height) of each histogram bar    
        ax1.bar(x, len(data_y), width=width,alpha=transparency, color="teal")
    # Scatter plot
    ax1.scatter(data_x, data_y, s=15, c="black")
    
    # X axes
    plt.xticks(xticks, rotation=90, fontsize=5)
    plt.xlim(right=xlim_right, left = xlim_left)
    plt.xlabel(score_term, fontsize=8, labelpad=5)
    # Y axes    
    ax1.tick_params(axis='y', labelsize=5,right= False, left= True, labelright= False, labelleft=True)
    y_lim = len(ax1.get_yticks())
    plt.ylim(bottom=0, top= y_lim-0.5)    
    ax1.xaxis.grid(color='gray', linestyle='dashed', linewidth=0.3) #grid features
    ax1.yaxis.grid(color='gray', linestyle='dashed', linewidth=0.3) #grid features
    ax1.set_axisbelow(True)

    plt.title(title , fontsize=12 )
    plt.tight_layout(pad=0.1)
    plt.savefig(score_term+"_"+ title+ '.pdf', format='pdf')
    plt.close()
    
def RMSD(vector1,vector2):
    euclidean_distance = ((vector1[0]-vector2[0])**2) + ((vector1[1]-vector2[1])**2) + ((vector1[2]-vector2[2])**2)
    rmsd = math.sqrt((euclidean_distance/2))
    return rmsd

def ca_coordinates(pdb, ch="A"):  
    data = []
    rosetta_position = 0
    with open(pdb,'r') as fopen:
        for line in fopen:
            if line.startswith("ATOM") and line[13:15] == "CA":
                if line[21] in ch:
                    position = int(line[22:26])
                    rosetta_position+=1
                    sequence = contractSeq([line[17:20]])  
                    data.append([float(line[30:38]), float(line[38:46]), float(line[46:54]), sequence[0], position, rosetta_position])     
        fopen.close()
    coordinates = pd.DataFrame(data, columns=["x","y","z","seq","PDB_#","rosetta_#"])
    coordinates.index = coordinates["rosetta_#"].tolist()
    return coordinates    
    
def read_alignment(alignment):  # Alignment in fasta format. Make sure the sequence alignment corresponds with the structural alignment
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

def round_half_up(n, decimals=0):
    multiplier = 10 ** decimals
    return math.floor(n*multiplier + 0.5 )/multiplier

### filters 
def per_resi_packing_score(pose, repetitions):
    packing = protocols.pose_metric_calculators.PackstatCalculator(repetitions,True)
    packing_per_res = json.loads(packing.get("residue_packstat", pose))  
    return packing_per_res

def CYS_neighbors(ref_resi, resi2, pose, total_res, cutoff, format_res = "rosetta", ch="A"):
    distance_cutoff = cutoff
    if format_res != "rosetta":
        ref_resi = pose.pdb_info().pdb2pose(ch,resi)
    r1 = pose.residue(ref_resi)
    r1_xyz = r1.xyz(r1.nbr_atom())
    for i in range(1,total_res+1):
        if i != resi2 and i != ref_resi:
            r2 = pose.residue(i)
            resi_type = r2.name() 
            if resi_type == "CYS:disulfide":   
                r2_xyz = r2.xyz(r2.nbr_atom())
                distance = r1_xyz.distance(r2_xyz)
                if distance <= distance_cutoff:         
                    return True            
    return False

def neighborhood(ref_resi, pose, total_res, cutoff, format_res = "rosetta", ch="A"):
    distance_cutoff = cutoff
    if format_res != "rosetta":
        ref_resi = pose.pdb_info().pdb2pose(ch,resi)
    r1 = pose.residue(ref_resi)
    r1_xyz = r1.xyz(r1.nbr_atom())
    neighbors = []
    for i in range(1,total_res+1):
        if i != ref_resi:
            r2 = pose.residue(i)
            r2_xyz = r2.xyz(r2.nbr_atom())
            distance = r1_xyz.distance(r2_xyz)
            if distance <= distance_cutoff:         
                neighbors.append(i)
    return neighbors

def get_sc_sasa(pose, positions, probe_size=2.2):
    sasa = SasaCalc()
    sasa.set_probe_radius(probe_size)
    sasa.calculate(pose)
    sasas = []    
    for res in positions:
        sasas.append(sasa.get_residue_sasa_sc()[res])
    return sasas

def hydrophobic(pose, residue_index):
    hydrophobic_aa = ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP']
    residue_name = pose.residue(residue_index).name3()
    if residue_name in hydrophobic_aa:
        return True
    else:
        return False

def scores_in_neighborhood(list_all_res_scores, list_neighbors): # Only Rosetta numbering
    total = 0
    for i in list_neighbors:
        total += list_all_res_scores[i-1]
    return total
            
################################################################################ 
### Command line arguments 
parser_arg = argparse.ArgumentParser(description='Selection of non-native disulfide bonds')
parser_arg.add_argument('arg_file', help = "File containing all input arguments")
args = parser_arg.parse_args()
home = os.getcwd() 
# Input arguments from arg_file
arg_dict = {}
with open(args.arg_file, "r") as fopen:
    for line in fopen:
        if line and line[0].isalpha():
            line = line.split("=")
            arg_dict[line[0].strip()] = line[1].strip()

## Discard designs where backbone deviated largely from starting structure 
ref_dir = arg_dict["starting_pdbs_dir"]            
designs_dir = arg_dict["designs_dir"]
designs = [i for i in os.listdir(designs_dir) if i.endswith("pdb")]        
rms_data = []
for pdb in designs:    
    ref_pdb = pdb.split("_")
    ref_pdb = "_".join(ref_pdb[:-3])+".pdb"
    ref_pdb = os.path.join(ref_dir,ref_pdb )
    design_pdb = os.path.join(designs_dir,pdb )
    rms = rms_CA(design_pdb, ref_pdb)
    rms_data.append([pdb, rms])
rms_data_df = pd.DataFrame(rms_data, columns=["pdb", "rms"])
                          
filtered_rmsd = rms_data_df[rms_data_df["rms"]<= 0.5]

### Initiallize pyrosetta
init()  
scorefxn = get_fa_scorefxn()

#### Find best scoring disulfides across all structures 
columns = ["SS_rosetta","fa13", "dihedral_x1","dihedral_x2" ,"dihedral_x3","angle_cb1_s1_s2", "dist_cb1_ca2", "description"]                      
disulfides = []
for i in filtered_rmsd.loc[:,"pdb"]:
    ss = i.split("_")
    ss = "_".join(ss[-3:-1])
    disulfides.append(ss)

disulfides = list(set(disulfides))

best_scoring = pd.DataFrame(index = disulfides, columns=columns)
best_scoring.fillna(0.1, inplace=True) # 0.1 as the energy threshold to acept a disulfide

for index, row in filtered_rmsd.iterrows():
    pdb = filtered_rmsd.loc[index,"pdb"]
    pdb_path = os.path.join(designs_dir,pdb)
    pdb_pose =  pose_from_pdb(pdb_path)
    scorefxn(pdb_pose)
    
    dslf = pdb.split("_")
    dslf = [dslf[-3],dslf[-2]]
    dslf_energy = pdb_pose.energies().residue_total_energies(int(dslf[0]))[dslf_fa13]
    dslf_str = "_".join(dslf)

    if dslf_energy < best_scoring.loc[dslf_str, "fa13"]:
        dihedral_x1, dihedral_x2, dihedral_x3, angle_cb1_s1_s2, dist_cb1_ca2 = dslf_geometry(dslf)            
        best_scoring.loc[dslf_str, "SS_rosetta"] = dslf_str
        best_scoring.loc[dslf_str, "fa13"] = dslf_energy           
        best_scoring.loc[dslf_str, "dihedral_x1"] = dihedral_x1
        best_scoring.loc[dslf_str, "dihedral_x2"] = dihedral_x2
        best_scoring.loc[dslf_str, "dihedral_x3"] = dihedral_x3
        best_scoring.loc[dslf_str, "angle_cb1_s1_s2"] = angle_cb1_s1_s2
        best_scoring.loc[dslf_str, "dist_cb1_ca2"] = dist_cb1_ca2
        best_scoring.loc[dslf_str, "description"] = pdb
   
best_scoring = best_scoring[best_scoring["fa13"] < 0.1] # Let's remove the disulfides that did not pass the energy threshold

### Remove duplicated disulfides
remove=[]
for disulf in best_scoring.index:
    resi = disulf.split("_")
    reverse = resi[1]+"_"+resi[0]
    if disulf not in remove and reverse not in remove:
        if reverse in best_scoring.index:
            sc1=best_scoring.loc[disulf,"fa13"] 
            sc2=best_scoring.loc[reverse,"fa13"] 
            if sc1 > sc2:
                remove.append(disulf)
            else:
                remove.append(reverse)
                
best_scoring.drop(index=remove, inplace=True)  

###### Native disulfides
native_pdb = arg_dict["starting_unrelaxed_pdb_dir"] 
native_pose = pose_from_pdb(native_pdb)
native_dslf = find_disulfides(native_pose)

relaxed_native_pdbs = [i for i in os.listdir(ref_dir) if i.endswith(".pdb")]
total_indexes = list(range(len(relaxed_native_pdbs)*len(native_dslf)))
native_data = pd.DataFrame(index=total_indexes, columns=columns[:-1])

idx = -1
for pdb in relaxed_native_pdbs:
    pdb_path = os.path.join(ref_dir,pdb)
    pdb_pose =  pose_from_pdb(pdb_path)
    scorefxn(pdb_pose)
    for dslf in native_dslf:
        idx +=1
        dslf_energy = pdb_pose.energies().residue_total_energies(dslf[0])[dslf_fa13]
        dslf_str = str(dslf[0])+"_"+str(dslf[1])
        dihedral_x1, dihedral_x2, dihedral_x3, angle_cb1_s1_s2, dist_cb1_ca2 = dslf_geometry(dslf)            
        native_data.loc[idx, "SS_rosetta"] = dslf_str
        native_data.loc[idx, "fa13"] = dslf_energy           
        native_data.loc[idx, "dihedral_x1"] = dihedral_x1
        native_data.loc[idx, "dihedral_x2"] = dihedral_x2
        native_data.loc[idx, "dihedral_x3"] = dihedral_x3
        native_data.loc[idx, "angle_cb1_s1_s2"] = angle_cb1_s1_s2
        native_data.loc[idx, "dist_cb1_ca2"] = dist_cb1_ca2
        native_data.loc[idx, "description"] = pdb    

# Group data according to their disulfide bond
numeric_columns = native_data.columns.difference(["SS_rosetta", "description"])
native_data[numeric_columns] = native_data[numeric_columns].apply(pd.to_numeric)
grouped = native_data.iloc[:,:-1].groupby('SS_rosetta')

# Remove outlier values
cleaned_groups = []
for group_name, group_data in grouped:
    for column in numeric_columns:
        group_data = remove_outliers(group_data, column)
    cleaned_groups.append(group_data)
cleaned_data = pd.concat(cleaned_groups)

# Now group by 'SS_rosetta' again and calculate the mean of the cleaned data 
average_native_data = cleaned_data.groupby('SS_rosetta').mean()
    
# Reorganize the disulfide labeling so that the smaller residue number is at the beginning                                                         
average_native_data["SS_rosetta"] = disulf_relabeling(average_native_data.index)
average_native_data.index = average_native_data.loc[:,"SS_rosetta"].tolist()     
best_scoring["SS_rosetta"] = disulf_relabeling(best_scoring.index)
best_scoring.index = best_scoring.loc[:,"SS_rosetta"].tolist()  
# Sort the disulfides
average_native_data = average_native_data.sort_values(by='SS_rosetta', key=lambda x: x.str.split('_').str[0].astype(int))
best_scoring = best_scoring.sort_values(by='SS_rosetta', key=lambda x: x.str.split('_').str[0].astype(int))

### Since this calculation takes long, we can save the data so we don't have to re-run it
with open('designed_disulfides_data.pkl', 'wb') as file:
    dill.dump({'average_native_data': average_native_data, 'best_scoring': best_scoring}, file)


############################### Filtering ###############################
''' If you saved the disulfide data and want to perform just the filtering process, the saved data can be loaded with the code below '''
#with open('designed_disulfides_data.pkl', 'rb') as file:
#    data = dill.load(file)
#average_native_data = data['average_native_data']
#best_scoring = data['best_scoring']

########## Comparison of geometric parameters with high-resolution structures
# Import data from high-resolution disulfides 
geometry_highRes = arg_dict["geometry_highResolution"] 
with open(geometry_highRes, 'rb') as file:
    data = dill.load(file)
data_high_resolution = data['data_high_resolution']

# Based on high-resolution data, define the limits within which the disulfide geometry should fall
limits = {}
for key, item in data_high_resolution.loc[:,"fa13":"dist_cb1_ca2"].iteritems():
    bins = number_bins(data_high_resolution.loc[:,key].tolist())    
    counts, edges, _ = plt.hist(data_high_resolution.loc[:,key], bins=bins)
    limits[key] = pd.DataFrame(edges.tolist(), columns=["edges"])
    limits[key]["counts"] = [np.nan]+ counts.tolist()
            
best_scoring_filter = best_scoring.copy(deep=True)
columns_to_check = best_scoring.loc[:,"dihedral_x1":"dist_cb1_ca2"].columns
for key in columns_to_check:
    best_scoring_filter[key] = best_scoring_filter[key].apply(lambda x: allowed_values(x, limits[key]["edges"], limits[key]["counts"]))
        
# Drop rows where any value is False
best_scoring_filter = best_scoring_filter[best_scoring_filter[columns_to_check].all(axis=1)]
best_scoring = best_scoring.loc[best_scoring_filter.index,: ]
 
########## Scoring analysis
columns = ["SS_rosetta","pack_res1","pack_res1_ctr","diff_pack_res1", "pack_res2","pack_res2_ctr","diff_pack_res2", "sasa_res1","sasa_res2","hydrophobic_res1","hydrophobic_res2","neighbor_CYS1","neighbor_CYS2","description"]
all_scoring_data = pd.DataFrame(index=best_scoring.index, columns = columns )

for idx in best_scoring.index:
    # General info
    all_scoring_data.loc[idx,"SS_rosetta" ] = best_scoring.loc[idx,"SS_rosetta"]
    all_scoring_data.loc[idx,"description" ] = best_scoring.loc[idx,"description"]   
    pdb = best_scoring.loc[idx,"description"]
    res1 = int(best_scoring.loc[idx,"SS_rosetta"].split('_')[0])
    res2 =  int(best_scoring.loc[idx,"SS_rosetta"].split('_')[1])  
    ## Controls    
    # Scoring
    pdb_control = pdb.split("_")[:4]
    pdb_control = "_".join(pdb_control)+".pdb"        
    pdb_control_path = os.path.join(ref_dir, pdb_control)
    control_pose =  pose_from_pdb(pdb_control_path)
    pose_control_monomer = control_pose.split_by_chain(1) #The control is a trimer while the designs are monomers
    total_res =  pose_control_monomer.total_residue()
    scorefxn(pose_control_monomer)
    # Neighbors
    neighbors1 = neighborhood(res1, pose_control_monomer, total_res, 7)
    neighbors2 = neighborhood(res2, pose_control_monomer, total_res, 7)
    # Packing 
    packing_ctr = per_resi_packing_score(pose_control_monomer, 5)
    all_scoring_data.loc[idx, "pack_res1_ctr"] = scores_in_neighborhood(packing_ctr, neighbors1)
    all_scoring_data.loc[idx, "pack_res2_ctr"] = scores_in_neighborhood(packing_ctr, neighbors2)
    # SASA
    sasas = get_sc_sasa(pose_control_monomer, [res1, res2], probe_size=2.2)
    all_scoring_data.loc[idx, "sasa_res1"] = sasas[0]
    all_scoring_data.loc[idx, "sasa_res2"] = sasas[1]
    # Hydrophobic?
    all_scoring_data.loc[idx, "hydrophobic_res1"] = hydrophobic(pose_control_monomer, res1)
    all_scoring_data.loc[idx, "hydrophobic_res2"] = hydrophobic(pose_control_monomer, res2)
    
    ## Design 
    # Scoring
    pdb_path = os.path.join(designs_dir, pdb)
    pdb_pose =  pose_from_pdb(pdb_path)
    scorefxn(pdb_pose) 
    # Packing 
    packing_design = per_resi_packing_score(pdb_pose, 5)
    all_scoring_data.loc[idx, "pack_res1"] = scores_in_neighborhood(packing_design, neighbors1)
    all_scoring_data.loc[idx, "pack_res2"] = scores_in_neighborhood(packing_design, neighbors2)
    # Packing difference compared to control
    all_scoring_data.loc[idx, "diff_pack_res1"] = all_scoring_data.loc[idx, "pack_res1"] - all_scoring_data.loc[idx, "pack_res1_ctr"]
    all_scoring_data.loc[idx, "diff_pack_res2"] = all_scoring_data.loc[idx, "pack_res2"] - all_scoring_data.loc[idx, "pack_res2_ctr"]
    # Neighbor CYS
    all_scoring_data.loc[idx, "neighbor_CYS1"] =  CYS_neighbors(res1, res2, pdb_pose, total_res, 7, format_res = "rosetta")
    all_scoring_data.loc[idx, "neighbor_CYS2"] =  CYS_neighbors(res2, res1, pdb_pose, total_res, 7, format_res = "rosetta")

##### Filtering
# Remove disulfides close to native disulfides to avoid unexpected interactions
all_scoring_data = all_scoring_data[(all_scoring_data["neighbor_CYS1"] == False) & (all_scoring_data["neighbor_CYS2"] == False)]
# Remove disulfides that disrupt the packing of buried hydrophobic residues
idx_to_remove = []
for idx in all_scoring_data.index:
    if all_scoring_data.loc[idx, "sasa_res1"] <= 20 and all_scoring_data.loc[idx, "hydrophobic_res1"] == True : # SASA threshold to identify buried residues
        if all_scoring_data.loc[idx, "diff_pack_res1"] <= -0.2: #if the packing got worse by more than 0.2 scoring units
            idx_to_remove.append(idx)
    if all_scoring_data.loc[idx, "sasa_res2"] <= 20 and all_scoring_data.loc[idx, "hydrophobic_res2"]  == True : # SASA threshold to identify buried residues
        if all_scoring_data.loc[idx, "diff_pack_res2"] <= -0.2: #if the packing got worse by more than 0.2 scoring units
            idx_to_remove.append(idx)
        
idx_to_remove = list(set(idx_to_remove))    
all_scoring_data.drop(index=idx_to_remove, inplace=True)  
best_scoring = best_scoring.loc[all_scoring_data.index,:]

###### Move best scoring candidates to a unique folder
out_path = os.path.join(home, "best_scoring_pdbs")
subprocess.call(["mkdir",out_path ])
for idx in best_scoring.index:
    pdb = best_scoring.loc[idx,"description"]
    old_path = os.path.join(designs_dir, pdb)
    new_path = os.path.join(home, out_path, pdb)
    subprocess.call(["cp",old_path, new_path ])

########## rotamer plots
##### angles_cb_s_s
bins = number_bins(data_high_resolution.loc[:,"angle_cb1_s1_s2"].tolist())  
customized_histogram(best_scoring, "angle_cb1_s1_s2", data_high_resolution, bins, "potential disulf", xticks= np.arange(96, 115, step=1), xlim_right=115, xlim_left=96)
##### dihedral1
bins = number_bins(data_high_resolution.loc[:,"dihedral_x1"].tolist())  
customized_histogram(best_scoring, "dihedral_x1", data_high_resolution, bins, "potential disulf", xticks= np.arange(-180, 185, step=30), xlim_right=185, xlim_left=-180)
##### dihedral2
bins = number_bins(data_high_resolution.loc[:,"dihedral_x2"].tolist())  
customized_histogram(best_scoring, "dihedral_x2", data_high_resolution, bins, "potential disulf", xticks= np.arange(-180, 190, step=30), xlim_right=190, xlim_left=-180)
##### dihedral3
bins = number_bins(data_high_resolution.loc[:,"dihedral_x3"].tolist())  
customized_histogram(best_scoring, "dihedral_x3", data_high_resolution, bins, "potential disulf", xticks= np.arange(-180, 185, step=30), xlim_right=185, xlim_left=-180)
##### dihedral3
bins = number_bins(data_high_resolution.loc[:,"dihedral_x3"].tolist())  
customized_histogram(best_scoring, "dihedral_x3", data_high_resolution, bins, "potential disulf", xticks= np.arange(-180, 185, step=30), xlim_right=185, xlim_left=-180)


################# RMSD analysis ################ 
# Obtain ca coordinates
pdb_post = arg_dict["pdb_post"]
monomer_chain_post = arg_dict["post_monomer_ch"]
coordinates_post = ca_coordinates(pdb_post, ch= monomer_chain_post)

pdb_pre = arg_dict["starting_unrelaxed_pdb_dir"] 
monomer_chain_pre = arg_dict["pre_monomer_ch"]
coordinates_pre = ca_coordinates(pdb_pre, ch= monomer_chain_pre)

# Use sequence alignment to identify residue positions that are shared between pre- and postfusion structures
alignment = read_alignment(arg_dict["alignment"])
alignment = alignment.rename(columns={arg_dict["pre_name_alignment"]: "prefusion", arg_dict["post_name_alignment"]: "postfusion"})
alignment = add_data_according_to_alignment("prefusion", alignment, coordinates_pre.loc[:,"rosetta_#"].tolist(), column_name="rosetta_#")
alignment = add_data_according_to_alignment("postfusion", alignment, coordinates_post.loc[:,"rosetta_#"].tolist(), column_name="rosetta_#")
                                                                                        
# Find the shared positions between pre- and postfusion                                                                                                                                                                                                                                                                         
idx_shared = []
for index, row in alignment.iterrows():
    if alignment.loc[index,"prefusion"] != "-" and alignment.loc[index,"postfusion"] != "-":
        idx_shared.append(index)
                                                    
# Trim coordinates to shared positions only
coordinates_pre_shared = coordinates_pre.loc[alignment.loc[idx_shared,"prefusion_rosetta_#"].tolist(), :]
coordinates_post_shared = coordinates_post.loc[alignment.loc[idx_shared,"postfusion_rosetta_#"].tolist(),:].set_index(pd.Index(coordinates_pre_shared.index.tolist()), inplace=False) # we need to have the two dataframes with the same indexes. The indexes of prefusion were chosen here  

# Calculate per residue root mean square deviation 
distance = []
for index, row in coordinates_pre_shared.iterrows():
    data1 = coordinates_pre_shared.loc[index, ["x","y","z"]].tolist()
    data2 = coordinates_post_shared.loc[index, ["x","y","z"]].tolist()
    rmsd = RMSD(data1,data2)
    distance.append(rmsd)
distance_df = pd.DataFrame(distance, index=coordinates_pre_shared.loc[:,"rosetta_#"], columns=["rmsd"])
                                                                      
# Add RMSD data to the alignment dataframe 
alignment.index = alignment.loc[:,"prefusion_rosetta_#"].tolist()                               
alignment.loc[distance_df.index.tolist(), "RMSD"] = distance

# Combine RMSD data with packing data
disulf = all_scoring_data.loc[:,"SS_rosetta"]
disulf_distance = []

for d in disulf:
    d = d.split("_")
    dist1 = alignment[alignment["prefusion_rosetta_#"] == int(d[0])].loc[:,"RMSD"].tolist()
    dist2 = alignment[alignment["prefusion_rosetta_#"] == int(d[1])].loc[:,"RMSD"].tolist()
    if np.isnan(dist1[0]):
        dist1[0] = 0
    if np.isnan(dist2[0]):
        dist2[0] = 0
                                          
    if dist1[0] > dist2[0]:
        disulf_distance.append(dist1[0] )           
    else:
        disulf_distance.append(dist2[0])   

all_scoring_data.insert(1, 'RMSD_pre_post', disulf_distance)

# Write output
summary = pd.concat([all_scoring_data.iloc[:,:-1], best_scoring.iloc[:,1:]], axis=1)    
summary.to_excel("all_scoring_data.xlsx", index=False) 

                                                                                       
                                                                                          

