#!/bin/bash

usage() {
	cat <<-__EOF__
		Usage:
			./psiclass_gff3_formatting.sh -i gtf_input -o gff3_output -f fasta_output -s stranded_or_unstranded -g genome -d directory [-h]

		Description:
			Formatting of GTF file from PSIclass to GFF3 file for EVM.

		Arguments:
			-i, --gtf_input GTF file from PSIclass (*_vote.gtf)
			-o, --gff3_output GFF3 file formatted for EVM
			-f, --fasta_output fasta from the gff3 output (transcript sequence to give to PASA)
			-s, --stranded stranded or unstranded transcriptome
			-g, --genome genome assembly
			-d, --directory scriptdir
			-h, --help

		Exemple: ./psiclass_gff3_formatting.sh -i transcriptome.gff3 -o transcriptome_for_EVM.gff3 -s stranded_or_unstranded -g genome -f fasta_output -d scriptdir
		__EOF__
}

# Eval command line arguments given in input
ARGS=$(getopt -o "i:o:f:s:g:d:h" --long "gtf_input:,gff3_output:,fasta_output:,stranded:,genome:,directory:,help" -- $@ 2> /dev/null)

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
		-i|--gtf_input)
			INPUT_GFF3=$2
			shift 2
			;;
		-o|--gff3_output)
			OUTPUT_GFF3=$2
			shift 2
			;;
		-f|--fasta_output)
			OUTPUT_FASTA=$2
			shift 2
			;;
		-s|--stranded)
			STRANDED=$2
			shift 2
			;;
		-g|--genome)
			GENOME=$2
			shift 2
			;;
		-d|--directory)
			DIRECTORY=$2
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

gffread -E ${INPUT_GFF3} -o ${INPUT_GFF3}.gff3
${DIRECTORY}/scripts/add_IDs_to_exons.py -g ${INPUT_GFF3}.gff3 -o ${OUTPUT_GFF3}
sed -i "s/PsiCLASS/PsiCLASS_${STRANDED}/" ${OUTPUT_GFF3}
rm ${INPUT_GFF3}.gff3
gffread -E --keep-genes ${OUTPUT_GFF3} -o- > ${OUTPUT_GFF3}.tmp

gffread -w ${OUTPUT_FASTA} -g $GENOME ${OUTPUT_GFF3}.tmp
rm ${OUTPUT_GFF3}.tmp
