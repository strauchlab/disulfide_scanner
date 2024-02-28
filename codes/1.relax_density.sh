#!/bin/bash

$Rosetta/main/source/bin/relax.linuxgccrelease \
-database $Rosetta/main/database \
-in:file:s INPUT_Monomer.pdb \ #the monomeric prefusion structure as symmetry is used to define the trimer
-relax:fast \
-relax:jump_move true \
-edensity:mapfile INPUT_map.ccp4 \ #density map of input PDB 
-edensity:mapreso x \ # x= map resolution of input PDB 
-edensity:fastdens_wt y \ # y= density weight according to map resolution. 
-symmetry_definition INPUT_symmetry_file.symm \ #symmetry definition file
-ignore_unrecognized_res \
-out::nstruct 200 \ #number of desired output (relaxed) structures
-ex1 -ex2 -use_input_sc \
-no_his_his_pairE \
-no_optH false \
-flip_HNQ \
-crystal_refine \
