#!/bin/bash

$Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease \
-database $Rosetta/main/database \
-in:file:s $pdb \ #input PDB
-parser:protocol fastRelax_highResolution.xml \
-ignore_unrecognized_res \
-nstruct 1 \
-use_input_sc \






