# BioHet
This programme takes a list of input biological PDBs (PDB1) and automates AlphaFold prediction for the orthologue of your complex. 
![all_multimers_tile](https://github.com/user-attachments/assets/947d6abd-df04-40d0-932b-c176a0691333)


- Prerequisites:
	- Programme is written in Python and run with Python 3.10.2
	- Programme was tested on a workstation running Ubuntu 20.04 with two Nvidia A5000 (2x 24 Gb VRAM) running AlphaFold 2.2 with the Kalinina lab install (https://github.com/kalininalab/alphafold_non_docker) 
	- A miniconda based local environment was made in Conda and BioPython installed 


- What you need:

	- Inputs should be placed in the 'input' directory 
	- Inputs:
		- List PDB IDs to be run in a text file (doesn't need a file extension). Each PDB should be on a new line 
		
		
- Running programme:
	- To run programme:	
		- bash BioHet_run.sh 
	- Programme will ask for:
		- Name of the file containing the list of input complexes 
		- Taxid for the organism of interest e.g. '562' for E. coli  
		- AlphaFold database date, is needed to run AlphaFold but not essential for this programme 
		- An email address in order for the BLAST submission to work, highly recommended as it seems that you get put to back of queue if the same email address submits many requests  


- Limitations:
	- The programme currently discards complexes with 6 or more chains, as these are assumed to be too large to be run with AlphaFold 
		- Will make this more intelligent 
	- BLAST is carried out with the BioBlast tool which is limited by requests so can be slow 
