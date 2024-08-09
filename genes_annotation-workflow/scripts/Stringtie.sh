#!/bin/bash

usage() {
	cat <<-__EOF__
		Usage:
			./Stringtie.sh -t threads -o transcriptome_gtf -b bam_file -r read [-h]

		Description:
			Launch Stringtie on bam file to create a transcriptome

		Arguments:
			-t, --threads number of CPUs
			-o, --output output transcriptome filename
			-b, --bam input bam file
			-r, --read read type (short or long)
			-h, --help

		Exemple: ./Stringtie.sh -t threads -o transcriptome_gtf -b bam_file -r read
		__EOF__
}

# Eval command line arguments given in input
ARGS=$(getopt -o "t:o:b:r:h" --long "threads:,output:,bam:,read:,help" -- $@ 2> /dev/null)

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
		-t|--threads)
			THREADS=$2
			shift 2
			;;
		-o|--output)
			OUTPUT=$2
			shift 2
			;;
		-b|--bam)
			BAM=$2
			shift 2
			;;
		-r|--read)
			READ=$2
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

if [[ $READ == "short" ]]
then

	OUTPUT_FINAL=$( basename $OUTPUT | sed "s/_Aligned.sortedByCoord.out.bam//" )

	BAM_FINAL=$( echo $BAM | sed "s/.*\/work\//\/work\//" )

	stringtie -p $THREADS -o $OUTPUT_FINAL $BAM_FINAL

elif [[ $READ == "long" ||  ]]
then

	OUTPUT_FINAL=$( basename $OUTPUT | sed "s/_Aligned.sorted.bam//" )

	BAM_FINAL=$( echo $BAM | sed "s/.*\/work\//\/work\//" )

	stringtie -L -p $THREADS -o $OUTPUT_FINAL $BAM_FINAL

fi