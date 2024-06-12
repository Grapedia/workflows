#!/bin/bash

usage() {
  cat <<-__EOF__
    Usage:
      ./trainglimmerhmm.sh -a assembly -t transcriptome -d output_dir -i intermediate_files_outdir [-h]

    Description:
      GlimmerHMM training
      This script will test different number of transcripts extracted from the transcriptome
      given in input. The highest value will be kept

    Argument:
      -a, --assembly Assembly to annotate in FASTA format
      -t, --transcriptome Transcriptome in GFF3 format
      -d, --dir Output directory
      -i, --inter intermediate_files outdir
      -h, --help

    Exemple: ./trainglimmerhmm.sh -a assembly.fasta -t transcriptome.gff3 -d training_dir -i intermediate_files_outdir
		__EOF__
}

# Eval command line arguments given in input
ARGS=$(getopt -o "a:t:d:i:h" --long "assembly:,transcriptome:,dir:,inter:,help" -- $@ 2> /dev/null)

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
    -a|--assembly)
			ASSEMBLY=$2
			shift 2
			;;
		-t|--transcriptome)
			TRANSCRIPTOME=$2
			shift 2
			;;
    -d|--dir)
  		DIR=$2
  		shift 2
  		;;
    -i|--inter)
    	INTER=$2
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

transcript_number=$(grep "transcript" ${TRANSCRIPTOME} | wc -l)

for i in 1 1.2 1.5 1.7 2 3 4 5
do
  mkdir -p $INTER/GlimmerHMM/
  transcript_subset=$(bc <<< ${transcript_number}/${i}-1)

  # Extract random transcripts from transcriptome
  /scripts/extract_exons_random_transcripts.py -i ${TRANSCRIPTOME} -o $INTER/GlimmerHMM/exons.tsv -n ${transcript_subset}

  # Train GlimmerHMM
  trainGlimmerHMM ${ASSEMBLY} $INTER/GlimmerHMM/exons.tsv -d ${DIR}

  # If training was successful, break the loop. The training directory is done and it is ...
  # ... not necessary to test other values
  # If not, delete the directory to test another number of transcripts for the training
  if [ $? = 0 ]
  then
    echo "Success with ${i}"
    break
  else
    rm -r $INTER/GlimmerHMM/
  fi
done
