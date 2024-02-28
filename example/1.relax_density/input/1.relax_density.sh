#!/bin/bash

$Rosetta/main/source/bin/relax.linuxgccrelease \
-database $Rosetta/main/database \
-in:file:s 7tn1_F_clean.pdb \
-relax:fast \
-relax:jump_move true \
-edensity:mapfile 7tn1_phases_2mFo-DFc.ccp4 \
-edensity:mapreso 3.19 \
-edensity:fastdens_wt 35.0 \
-symmetry_definition symmetry.symm \
-ignore_unrecognized_res \
-out::nstruct 200 \
-ex1 -ex2 -use_input_sc \
-no_his_his_pairE \
-no_optH false \
-flip_HNQ \
-crystal_refine \

