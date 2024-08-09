#!/bin/bash

usage() {
	cat <<-__EOF__
		Usage:
			./Stringtie_merging.sh -o transcriptome_gtf -g gtf_files [-h]

		Description:
			Launch Stringtie on bam file to create a transcriptome

		Arguments:
			-o, --output output transcriptome filename
			-g, --gtf input gtf files
			-h, --help

		Exemple: ./Stringtie_merging.sh -o transcriptome_gtf -g gtf_files
		__EOF__
}

# Eval command line arguments given in input
ARGS=$(getopt -o "o:g:h" --long "output:,gtf:,help" -- $@ 2> /dev/null)

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
		-o|--output)
			OUTPUT=$2
			shift 2
			;;
		-g|--gtf)
			GTF=$2
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

echo $GTF | tr ',' '\n' > GTF_list.txt
stringtie --merge -o $OUTPUT GTF_list.txt
rm GTF_list.txt
