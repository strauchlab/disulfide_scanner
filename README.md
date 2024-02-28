# Disulfide design
## Requirements:

1. Rosetta version: 2020.10.post.dev+12.master.c7b9c3e c7b9c3e4aeb1febab211d63da2914b119622e69b  
   Instructions on how to install Rosetta can be found [here](https://new.rosettacommons.org/demos/latest/tutorials/install_build/install_build).
3. PyRosetta-4 2019 [Rosetta PyRosetta4.conda.linux.CentOS.python37.Release 2019.47+release.3d995b15922374726493453159734beedd7e28be 2019-11-20T17:52:20]  
   Instructions on how to install PyRosetta can be found [here](https://www.pyrosetta.org/downloads).
5. Python3.7
6. Python packages: 
	- os
 	- pandas (version 1.3.5)
 	- numpy
 	- subprocess
 	- argparse
 	- pyrosetta 
 	- rosetta
 	- math
  	- Bio.PDB
  	- dill
   	- json
   	- matplotlib

> [!NOTE]
> All scripts were tested in Ubuntu 18.04.5 LTS.

## Step-by-step scripts
--------------------- 
> [!TIP]
> - Additional information is available within each script in the codes folder.
> - Before initiating the design process, we advise using the script $Rosetta/tools/protein_tools/scripts/clean_pdb.py to preprocess the input PDB. Utilizing a clean PDB helps prevent potential confusion during the design, as all calculations in the subsequent processes rely on the Rosetta numbering system. 

### 1. Structural relaxation guided by density data  

As the R-1b protein is a homotrimer, structural relaxation can be streamlined using a symmetry definition file. Refer to [this tutorial](https://faculty.washington.edu/dimaio/files/_density_tutorial_aug18.pdf) for information on symmetry files.

To generate the symmetry file, use the following command line: 

  	perl $/main/source/src/apps/public/symmetry/make_symmdef_file.pl -m NCS -a A -i B -p INPUT.pdb > symmetry.symm

Once the symmetry file is obtained, the relaxation process is carried out with the script:

	1.relax_density.sh
 
To account for the dynamic nature of proteins, it is essential to generate a comprehensive set of relaxed models (~200). This approach ensures a thorough exploration of the protein's structure.


### 2. Disulfide scanning

Disulfide scanning is carried out for each model obtained through relaxation. For this, use the script:

	2.build_disulfides.py [-h] arg_file pdb_model

This script requires an external file (arg_file) containing various arguments and the PDB model for scanning. As mentioned, the script should be independently executed for each PDB generated during relaxation. An illustrative example of running the script is provided below:

	python3.7 2.build_disulfides.py arguments_file model_001.pdb

> [!NOTE]
> For a detailed illustration of the external file (arg_file), please consult /example/2.build_disulfides/input/2.arguments_file. This file contains all the necessary information needed as arguments; refer to it for guidance.

The script produces an output folder named 'disulfide_pdbs,' which includes PDB files showcasing potential disulfide bonds within each relaxed model. The generated PDBs are labeled based on their parent PDB, with the residue positions involved in forming the disulfide bond specified at the end of each PDB name. Notably, residue positions follow the Rosetta numbering system.

> [!IMPORTANT]
> This script is specifically designed for analyzing intraprotomer disulfides; consequently, output structures correspond to monomeric subunits.
	
### 3. Refinement of PDB structures incorporating newly formed disulfide bonds

After preliminary identification of potential disulfide bonds, performing energy minimization and refinement around the newly formed disulfide is essential to assess the geometric adequacy of the bond. This refinement process must be performed on all PDBs obtained in step 2.

The refinement is done with the script:  

	3.fastRelax_monomer.xml 
 
which can be executed using:

 	3.run_fastRelax_monomer.sh
   
> [!IMPORTANT]
> Ensure all residue positions provided as script input are in Rosetta numbering.


### 4. Selection of potential disulfide bonds

The selection of candidate disulfides can be done using the script:

	4.disulfides_selection.py [-h] arg_file
  
This script takes as input an external file (arg_file) containing various arguments. For a detailed illustration of the arg_file, please refer to /example/4.selection/input/4.arguments_file. 

To assess the geometric adequacy of potential disulfides, the script compares the modeled disulfides' geometric values with a reference set derived from crystallized high-resolution disulfides. You can locate this dataset at /example/4.selection/high_resolution_disulfides/relaxed_pdbs.<br><br>

> The protocol used to relax the high-resolution dataset is found at /codes/4.0.fastRelax_highResolution.xml with the execution script /codes/4.0.run_fastRelax_highResolution.sh.
>
> To streamline the comparison process, we have pre-calculated the high-resolution dataset's geometric parameters and provided this information in the file /example/4.selection/input/geometry_high_resolution.pkl. Alternatively, these values can be generated anew using the following script:
>
>		4.1.disulfides_analysis_highResolution.py [-h] pdbs_directory
>
> An example of how to run the above script is provided at /example/4.selection/high_resolution_disulfides/geometry_evaluation/input

<br>The disulfide selection script produces several outputs: a new folder containing PDB structures of the selected disulfide candidates after filtering, various plots illustrating the geometric parameters of the newly formed disulfides compared to the high-resolution dataset, and an Excel file summarizing all calculated parameters. Below, you will find a detailed explanation of each of these parameters.

- ***SS_rosetta***: Indicates the positions in Rosetta numbering of residues forming the disulfide bond.
  
- ***RMSD_pre_post***: Measures the root-mean-square deviation (RMSD) between the prefusion and postfusion structures, identifying regions with notable conformational dynamics. Disulfides in regions with high RMSD values may be more stabilizing of the prefusion state.
  
- ***Pack_res1***: Represents the packing score within a 7Å radius around the first residue (res 1) of the newly formed disulfide bond. Higher scores suggest better packing. Due to the stochastic nature of this calculation, slight variations in the computed values may occur with each run.

- ***Pack_res1_ctr***: Denotes the packing score within a 7Å radius around the original residue (res 1) before mutation to CYS in the new disulfide bond. Higher scores indicate better packing. Due to the stochastic nature of this calculation, slight variations in the computed values may occur with each run.
  
- ***Diff_pack_res1***: Measures the difference in packing score between areas surrounding the original residue and the CYS mutant (packing_mutant_residue1 - packing_control_residue1). Positive values indicate better packing in the mutant residue, while negative values suggest better packing for the original residue.
  
- ***Sasa_res1***: Reflects the solvent-accessible surface area of the first residue involved in the disulfide bond. Lower values indicate buried residues.
  
- ***Hydrophobic_res1***: Indicates whether the original residue (before mutation to CYS) from the new disulfide bond is hydrophobic. This parameter, combined with SASA and packing values, helps determine the potential disruption of protein packing.

- ***Neighbor_CYS1***: Counts the number of cysteines within a 7Å radius around the first residue forming the new disulfide bond. This parameter aids in identifying neighboring native disulfide bonds that could interact nonspecifically with the introduced disulfide.
  
- ***Fa13***: Represents the Rosetta-computed energy of the newly formed disulfide bond.
  
- ***Dihedral_x1, x2, x3***: Refers to the dihedral angles of the newly formed disulfide bond.
  
- ***Angle_cb1_s1_s2***: Measures the angle between the Cβ atom and sulfur atoms forming the disulfide bond.
  
- ***Dist_cb1_ca2***: Represents the distance between the Cβ atom from residue 1 and the Cα atom from residue 2 forming the disulfide bond.
   
- ***Description***: Provides the name of the PDB structure containing the analyzed disulfide bond.

Additionally, the script generates an optional "designed_disulfides_data.pkl" file, which contains the geometric analysis performed before the filtering step. This file can be beneficial for avoiding redundant calculations across all PDB structures with new disulfides, especially considering the computational demands involved.

To access this information, simply uncomment the following lines in the 4.disulfides_selection.py script:

	import dill
 
	with open('designed_disulfides_data.pkl', 'rb') as file:
	    data = dill.load(file)
	average_native_data = data['average_native_data']
	best_scoring = data['best_scoring']

