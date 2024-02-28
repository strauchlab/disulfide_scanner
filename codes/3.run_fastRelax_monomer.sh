#!/bin/bash

$Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease \
-database $Rosetta/main/database \
-in:file:s $pdb \  #PDB model with newly introduced disulfide
-parser:protocol 3.fastRelax_monomer.xml \
-ignore_unrecognized_res \
-nstruct 1 \
-ex1 -ex2 -use_input_sc \
-parser::script_vars res1=$res1 res2=$res2 \  # res1 and res2 denote the residue positions involved in forming the disulfide bond (in Rosetta numbering). These positions can be found at the end of the PDB model's name. 







