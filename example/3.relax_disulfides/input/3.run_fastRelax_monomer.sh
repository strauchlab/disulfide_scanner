#!/bin/bash

pdb=/example/2.build_disulfides/output/7tn1_F_clean1_0001_88_77.pdb
res1=88
res2=77


$Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease \
-database $Rosetta/main/database \
-in:file:s $pdb \
-parser:protocol 3.fastRelax_monomer.xml \
-ignore_unrecognized_res \
-nstruct 1 \
-ex1 -ex2 -use_input_sc \
-parser::script_vars res1=$res1 res2=$res2 \  







