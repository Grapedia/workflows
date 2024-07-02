#!/bin/bash

usage() {
	cat <<-__EOF__
		Usage:
			./edta.sh -g genome -n threads [-h]

		Description:
			Launch Extensive de-novo TE Annotator (EDTA) on genome

		Arguments:
			-g, --genome genome to annotate (fasta)
			-n, --threads number of CPUs
			-h, --help

		Exemple: ./edta.sh -g genome.fasta -n threads
		__EOF__
}

# Eval command line arguments given in input
ARGS=$(getopt -o "g:n:h" --long "genome:,threads:,help" -- $@ 2> /dev/null)

# Check if the return code of the previous command is not equal to 0 (command ...
# ... didn't work)
# >&2 send the message to standard error (stderr) instead of standard out (stdout)
if [ $? -ne 0 ]; then
	echo "Error in the argument list. Use -h or --help to display the help." >&2
	exit 1
fi

eval set -- ${ARGS}
while true
do
	case $1 in
		-g|--genome)
			GENOME=$2
			shift 2
			;;
		-n|--threads)
			THREADS=$2
			shift 2
			;;
		-h|--help)
			usage
			exit 0
			;;
		--) shift
			break
			;;
		*)  echo "Option $1 is not recognized. Use -h or --help to display the help." && \
			exit 1
			;;
	esac
done

eval "$(conda shell.bash hook)"
conda activate EDTA2

/EDTA/EDTA.pl --genome $GENOME --species others --step all --sensitive 1 --anno 1 --overwrite 1 --threads $THREADS
