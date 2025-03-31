#echo Input P2RANK directory
#read P2RANK_dir
echo Input filename with list of PDB1s
read input

if [ -z "$input" ]; then  
	echo No input provided 
else
	echo Running BioHet on: $input
	echo Input path to AlphaFold working directory from '~/'
	read AF_path
	echo Input name '(base)' of AlphaFold conda environment
	read AF_env
	python BioHet.py $input
	varEnd="_command"
	input_dir=$(pwd)
	echo Outputs will be in: $input_dir/output 
	source ~/miniconda3/etc/profile.d/conda.sh
	conda activate $AF_env
	cd ~/
	cd $AF_path
	sh ${input_dir}/output/${input}${varEnd}
fi


echo BioHet Finished
