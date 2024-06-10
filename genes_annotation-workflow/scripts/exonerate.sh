#!/bin/bash

function usage {
	cat <<-__EOF__
		Usage:
			./exonerate.sh -g assembly -a alignments -q query -o outpref -d directory [-h]

		Description:
			Run exonerate alignment with exonerate-server

		Arguments:
			-a, --alignments pblat alignments result (psl)
			-g, --genome genome assembly to target (fasta)
			-q, --query Sequences you want to map
			-o, --outpref Output file prefix (organism for exemple)
			-d, --directory scriptdir
			-h, --help

		Exemple: ./exonerate.sh -g assembly.fasta -a alignments.psl -q query.fasta -o prefix -d scriptdir
		__EOF__
}

# Eval command line arguments given in input
ARGS=$(getopt -o "a:g:q:o:d:h" --long "alignments:,genome:,query:,outpref:,directory:,help" -- $@ 2> /dev/null)

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
		-a|--alignments)
			ALIGNMENTS=$2
			shift 2
			;;
		-q|--query)
			QUERY=$2
			shift 2
			;;
		-o|--outpref)
			OUT_PREFIX=$2
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

target_file=${TARGET##*/}
target_file_name=${target_file%%.*}

exonerate_subset_path="intermediate_files/evidence_data/protein/exonerate/${OUT_PREFIX}"
exonerate_tmp_path="intermediate_files/evidence_data/protein/exonerate/${OUT_PREFIX}/TMP"
exonerate_concat_path="intermediate_files/evidence_data/protein_final_alignments"
log_path="logs/exonerate/${OUT_PREFIX}"

mkdir -p ${exonerate_subset_path}
mkdir -p ${exonerate_tmp_path}

# Exonerate doesn't works with ambiguous or exceptional amino acid (like B, J, O, U, Z)
# so it generates an error for eudicotyledons_odb10 fasta file for example
# replace these letters by X
echo "sed '/^[^>]/s/[BJOUZ]/X/g' ${QUERY} > ${QUERY}.tmp" &> ${log_path}/exonerate_serveur.log
sed '/^[^>]/s/[BJOUZ]/X/g' ${QUERY} > ${QUERY}.tmp

for protein in `cut -f10 $ALIGNMENTS | sort -u`
do
	# retrieve all the psl results for the protein
	awk -v prot="$protein" '$10 == prot { print $0 }' $ALIGNMENTS > ${exonerate_tmp_path}/${protein}.psl
	# extract regions
	${DIRECTORY}/extract_mapped_regions_from_psl.py -i ${exonerate_tmp_path}/${protein}.psl -o ${exonerate_tmp_path}/${protein}.bed
	# sort, merge and expand each coordinates by 100
	bedtools sort -i ${exonerate_tmp_path}/${protein}.bed > ${exonerate_tmp_path}/${protein}.sorted.bed
	bedtools merge -i ${exonerate_tmp_path}/${protein}.sorted.bed > ${exonerate_tmp_path}/${protein}.merged.bed
	if [ ! -f ${GENOME}.fai ]; then
		samtools faidx ${GENOME}
	fi
	cut -f1,2 ${GENOME}.fai > ${exonerate_tmp_path}/genome.txt
	bedtools slop -i ${exonerate_tmp_path}/${protein}.merged.bed -g ${exonerate_tmp_path}/genome.txt -b 500 > ${exonerate_tmp_path}/${protein}.mergedfinal.bed
	# extract target region of the query for exonerate
	${DIRECTORY}/bed_to_coords_for_fasta_extraction.py -i ${exonerate_tmp_path}/${protein}.mergedfinal.bed -o ${exonerate_tmp_path}/${protein}.converted.txt
	samtools faidx --region-file ${exonerate_tmp_path}/${protein}.converted.txt ${GENOME} -o ${exonerate_tmp_path}/regions_mapped_on_new_assembly.${protein}.fasta
	# retrieve protin sequence (query)
	samtools faidx ${QUERY}.tmp "$protein" > ${exonerate_tmp_path}/${protein}.fasta
done

for protein in `cut -f10 $ALIGNMENTS | sort -u`
do
	if [[ ! -f "${exonerate_subset_path}/${protein}.exonerate.gff" && ! -s "${exonerate_subset_path}/${protein}.exonerate.gff" ]]; then
		exonerate --maxintron 20000 --model p2g --showtargetgff yes --showalignment no --showvulgar no \
		--ryo 'AveragePercentIdentity:%pi\nAveragePercentSimilarity:%ps\nScore:%s\nQueryLen:%ql\nAlignment:%qal\n' \
		-g --softmasktarget --geneseed 100 \
		--query ${exonerate_tmp_path}/${protein}.fasta \
		--target ${exonerate_tmp_path}/regions_mapped_on_new_assembly.${protein}.fasta > ${exonerate_subset_path}/${protein}.exonerate.gff || \
		exonerate --maxintron 20000 --model p2g --showtargetgff yes --showalignment no --showvulgar no \
		--ryo 'AveragePercentIdentity:%pi\nAveragePercentSimilarity:%ps\nScore:%s\nQueryLen:%ql\nAlignment:%qal\n' \
		-g --softmasktarget --geneseed 100 \
		--query ${exonerate_tmp_path}/${protein}.fasta \
		--target ${exonerate_tmp_path}/regions_mapped_on_new_assembly.${protein}.fasta > ${exonerate_subset_path}/${protein}.exonerate.gff
	fi
done

rm -fr $exonerate_tmp_path

mkdir -p ${exonerate_concat_path}

find ${exonerate_subset_path} -type f -name "*.exonerate.gff" -exec cat {} + > ${exonerate_concat_path}/${OUT_PREFIX}.gff
find ${exonerate_subset_path} -name "*.exonerate.gff" -type f -delete
rm ${QUERY}.tmp
