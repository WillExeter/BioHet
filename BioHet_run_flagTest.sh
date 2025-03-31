#echo Input P2RANK directory
#read P2RANK_dir

mkdir output
mkdir output/XMLs

input=""
AF_path=""
AF_env=""
AF_db=""
AF_date=""
email="wss205@exeter.ac.uk"
taxid=""

# Function to display help message
show_help() {
    echo "Usage: $0 [options]"
    echo
    echo "Options:"
    echo "  -i  <input>          Input file or directory with list of PDB1s"
    echo "  -p  <AF_path>        Absolute path to AlphaFold installation"
    echo "  -e  <AF_env>         Name of conda AlphaFold environment"
    echo "  -d  <AF_db>          Absolute AlphaFold database path"
    echo "  -t  <AF_date>        AlphaFold database time (e.g., 2024-10-16)"
    echo "  -e  <email>          Email address, needed to run BLAST"
    echo "  -o  <taxid>          Taxonomic identifier (e.g., 9606 for human, 562 for E. coli)"
    echo "  -h                   Show this help message"
    echo
    echo "Example:"
    echo "  $0 -i input.txt -p /path/to/af -e env -d /path/to/af/db -t 2024-10-16 -e user@example.com -o 9606"
}

# Parse the flags using getopts
while getopts "i:p:c:d:t:e:o:h" flag; do
    case "${flag}" in
        i) input=${OPTARG} ;;
        p) AF_path=${OPTARG} ;;
        c) AF_env=${OPTARG} ;;
        d) AF_db=${OPTARG} ;;
        t) AF_date=${OPTARG} ;;
        e) email=${OPTARG} ;;
        o) taxid=${OPTARG} ;;
        h) show_help; exit 0 ;;
        *) echo "Invalid option: -${OPTARG}" ;;  # Handling invalid flags
    esac
done


shift $((OPTIND-1))


#echo Input filename with list of PDB1s
#read input

if [ -z "$input" ]; then  
	echo No input provided 
else
	echo Running BioHet on: $input
	#read AF_path
	#echo Input name '(base)' of AlphaFold conda environment
	#read AF_env
	python BioHet.py $input $AF_path $AF_db $AF_date $email $taxid 
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



