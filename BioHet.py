import numpy as np
import sys
import os
import io
import requests
from Bio import SeqIO
from Bio import SearchIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
import argparse

def download_pdb(PDB_name, output_dir): # Setting download PDB function
	os.makedirs(output_dir, exist_ok=True)
	url = f"https://files.rcsb.org/download/{PDB_name}" #Download from RCSB
	output_filename = os.path.join(output_dir, f"{PDB_name}") # Setting output file path
	try:
		response = requests.get(url)
		with open(output_filename, 'w') as file:
			file.write(response.text)
		print(f"Successfully downloaded {PDB_name}")
	except requests.RequestException as error:
		print(f"Failed downloading {PDB_name}: {error}")

def AF_pred_file(input_file, output_path, command_output, cwd, org_taxid, AF_date, user_email):
	Entrez.email = user_email
	list = open(input_file, 'r')
	list = list.readlines()
	list_toUse = []
	global not_run_list # = []
	for blank in list:
		if blank != '\n':
			list_toUse.append(blank)	
#	print(list_toUse)
	for lone in list_toUse:
		PDB_name = lone[:4]+'.pdb1'
		print(f"Downloading first PDB: {PDB_name}")
		output_dir = "input/pdbs"

		download_pdb(PDB_name, output_dir)		
			
		PDB = open(f'input/pdbs/{PDB_name}')
		#PDB = PDB.readlines()
		cur_chain = -1
		chain_count = 0
		chain_id_list = []
		chain_good = False
		model_id = -1
		model_count = 0
		model_id_list = []
		model_good = False
		model_chain_id = -1
		model_chain_count = 0
		model_chain_list = [] 
		
		for ltwo in PDB:
			if ltwo[0:6] == "SEQRES":
				if ltwo[11] != cur_chain:
					cur_chain = ltwo[11]
					chain_count += 1
					chain_id_list.append(ltwo[11])
			if ltwo[0:5] == "MODEL":
				if ltwo[13] != model_id:
					model_id = ltwo[13]
					model_count += 1
					model_id_list.append(ltwo[13])	

#		print(chain_count)	
#		print(model_count)
		
		
		if chain_count > 1 and chain_count < 10:
			print(f"{chain_count} chain(s) found, does this match?")
			print(chain_id_list)
			chain_good = True 
		if model_count > 1 and model_count < 10:
			print(f"{model_count} model chain(s) found, does this match?")
			print(model_id_list)
			model_good = True
			

						
			
		if chain_count == 1 or chain_count > 10: # Modified this to 10 for paper structures
			if model_count == 1 or model_count > 10:
				print("Monomer or too many chains, not running")
				not_run_list.append(PDB_name)
			
		if chain_good == True or model_good == True:
			FASTA_temp = []
			#for lone in list_toUse: # Going through the input files one by one # Remove
			AF_FASTA = ""
			#chain_id_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'] # Assuming .pdb1s start at A, not more than 13 chains and are in order... 
			resDictionary = {"ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "GLN": "Q", "PHE": "F", "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N", "PRO": "P", "VAL": "V", "TRP": "W", "TYR": "Y", "ARG": "R", "SER": "S", "THR": "T", "   ": "", "MSE": "M"}
			histidine = "HIS"
			isoleucine = "ILE"
			lysine = "LYS"
			leucine = "LEU"
			methionine = "MET"
			asparagine = "ASN"
			proline = "PRO"
			valine = "VAL"
			phenylalanine = "PHE"
			tryptophan = "TRP"
			tyrosine = "TYR"
			arginine = "ARG" 
			serine = "SER"
			threonine = "THR"
			#chain_list = np.arange(chain_count)
	#		chain_list = chain_list + 1
			FASTA_checkForNext = ""
			
	#			num_model = len(model_id_list) # Since model is str, this gets number of models
	#			print(f'this number of models: {num_model}')
				
			for i in chain_id_list:
#				print(f'in loop {i}')
				chain_id = i
	#			print(chain_id)
				#FASTA = [f"> Chain id = {chain_id} \n"]
				FASTA = []
				
#				for lone in list_toUse:
#					#PDB_name = lone[:9]
				PDB_second = open(f'input/pdbs/{PDB_name}') 
				for lthree in PDB_second:
					if (lthree[:6] == "SEQRES" and lthree[11] == chain_id):
						FASTA.append(resDictionary[lthree[19:22]]) # getting first residue from each SEQRES line, looking up single letter code from dictionary and appending to FASTA 
						FASTA.append(resDictionary[lthree[23:26]])
						FASTA.append(resDictionary[lthree[27:30]])
						FASTA.append(resDictionary[lthree[31:34]])
						FASTA.append(resDictionary[lthree[35:38]])
						FASTA.append(resDictionary[lthree[39:42]])
						FASTA.append(resDictionary[lthree[43:46]])
						FASTA.append(resDictionary[lthree[47:50]])	
						FASTA.append(resDictionary[lthree[51:54]])
						FASTA.append(resDictionary[lthree[55:58]])
						FASTA.append(resDictionary[lthree[59:62]])
						FASTA.append(resDictionary[lthree[63:66]])
						FASTA.append(resDictionary[lthree[67:70]])
						#print(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10,res11,res12,res11)
						#FASTA.append(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10,res11,res12,res13)
						#' ,'.join(FASTA)
				print(f'Chain {chain_id} ')
#						print(FASTA)
				old_stdout = sys.stdout
				new_stdout = io.StringIO()
				sys.stdout = new_stdout
				print(*FASTA, sep='') # This is needed to clean up the FASTA sequence
				FASTA_clean = new_stdout.getvalue()
				sys.stdout = old_stdout
				if FASTA_checkForNext != FASTA_clean: # Attempting to prevent blasting of same subunit to speed up programme
					FASTA_checkForNext = FASTA_clean
					print(f'Sequence: {FASTA_clean}')
					result_handle = NCBIWWW.qblast("blastp", "nr", FASTA_clean, entrez_query=f'txid{org_taxid}[ORGN]', hitlist_size=1, alignments=1)
	#				print("Ran blast")
					blast_result = open(f'output/XMLs/{PDB_name}_chain{chain_id}.xml', "w")
	#				print("Saving file")
					blast_result.write(result_handle.read())
	#				print("extracting result")
	#				blast_read = SearchIO.read(f'{PDB_name}_chain{chain_id}.xml', 'blast-xml')
	#				print(blast_read)
	#				SeqIO.write(blast_read, f'{PDB_name}_chain{chain_id}.fasta', "fasta")
	#				for lthree in blast_result:
	#					if lthree[1:14] == "Hit_accession"
					
					blast_result = open(f'output/XMLs/{PDB_name}_chain{chain_id}.xml', 'r')
					blast_records = NCBIXML.parse(blast_result)
					
					
					for blast_record in blast_records:
						for alignment in blast_record.alignments:
							for hsp in alignment.hsps:
								search_ID = alignment.accession
								handle = Entrez.efetch(db="protein", id=search_ID, rettype="fasta", retmode="text")
								record = handle.read()
								record = (record.rstrip('\n'))
								AF_FASTA = (AF_FASTA+f'{record}\n')
								print(AF_FASTA)
				     #blast_result.close()
			   	     #result_handle.close()
					handle.close()
							    
				else:
					print(f'Chain {chain_id} same as previous chain, skipping blast')
					AF_FASTA = (AF_FASTA+f'{record}\n')
					print(AF_FASTA)

			
		if model_count > 1:
			AF_FASTA_model = ""
			for mod in model_id_list:
				AF_FASTA_model = AF_FASTA_model+AF_FASTA
			print(f"Model fasta: \n{AF_FASTA_model}")
		else:
			AF_FASTA_model = AF_FASTA
		
		output_dir_name = f'{output_path}/{PDB_name[0:4]}_homologue'
		os.makedirs(f'{output_dir_name}', exist_ok=True)
		FASTA_write = open(f'{output_dir_name}/{PDB_name[0:4]}_{org_taxid}homologue_FASTA.fasta', 'w')
		FASTA_write.write(AF_FASTA_model)
		FASTA_write.close()
		command_output.write(f'bash run_alphafold.sh -d {AF_db} -o {output_dir_name}/BioHet_AlphaFold/{PDB_name[0:4]}_{org_taxid}_homologue -f {output_dir_name}/{PDB_name[0:4]}_{org_taxid}homologue_FASTA.fasta -t {AF_date} -m multimer -e false -a 0,1 -l 1\n')
		
parser = argparse.ArgumentParser(description="")
parser.add_argument('input', metavar='name of input file', type=str, help='Program will look in input directory for this file')
parser.add_argument('AF_path', metavar='/path/to/alphafold', type=str, help='absolute path to main alphafold directory, for running alphafold')
parser.add_argument('AF_db', metavar='/path/to/alphafold/directories', type=str, help='absolute path to alphafold database')
parser.add_argument('AF_date', metavar='e.g., 2024-10-16', type=str, help='AlphaFold database time')
parser.add_argument('email', metavar='user@email.com', type=str, help='email needed for BLAST')
parser.add_argument('taxid', metavar='e.g. 562', type=str, help='NCBI Taxonomic identifier')

args = parser.parse_args()

input_name = args.input
AF_path = args.AF_path
AF_db = args.AF_db
AF_date = args.AF_date
email = args.email
org_taxid = args.taxid

cwd = os.getcwd()
#print(cwd)
#input_name = sys.argv[1]
input_file = os.path.join(cwd, 'input', f'{input_name}')
output_path = os.path.join(cwd, 'output')

#print("How much GPU memory is available (Gb)? (PDB files containing complexes which are too large for AlphaFold will be skipped, a list of these will be available afterwards in the log file) ")  # Incorporate system which counts number residues in chains from PDB file and then determines whether too large after user input of VRAM?
#VRAM = input()
#max_resid = VRAM*80 # Working on principle of ~ 80 residues per Gb


#print("Input AlphaFold database location ")
#AF_db = AF_db()
#print("Input AlphaFold database date (YYYY-MM-DD) ")
#AF_date = input()
#print("Enter email (blank uses default email, see README) ")
user_email = {email}
if user_email == "":
	user_email = ("wss205@exeter.ac.uk")
	

command_output = open(f'output/{input_name}_command', 'a')

not_run_list = []
AF_pred_file(input_file, output_path, command_output, cwd, org_taxid, AF_date, user_email)
command_output.close()

not_run_log = open(f'{output_path}/{input_name}.log', 'a')
not_run_log_join = "\n".join(not_run_list)
not_run_log.write(not_run_log_join)
not_run_log.close()
print(f'these were not run: {not_run_list} \nThese will be in the log file in the output directory')




