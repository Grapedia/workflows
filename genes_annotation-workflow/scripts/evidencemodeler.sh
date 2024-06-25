#!/bin/bash

usage() {
	cat <<-__EOF__
		Usage:
			./evidencemodeler.sh -e exec_dir -a genome -w weights -g glimmerhmm -i geneid -l liftoff -m maker -b braker2 -p protein_alignments -t unstranded_transcriptome -s stranded_transcriptome -r repeats -u annotations.txt -o annotations.gff3 -c 2_abinitio_output.gff3 -d 1_abinitio_output.gff3 -f evidencedata_only_output.gff3 -j 1_abinitio_output.proteins.fasta [-h]

		Description:
			Launch EvidenceModeler as a final step for genome annotation

		Arguments:
			-e, --exec_dir working directory
			-o, --output final GFF3 output file
			-u, --untreated_output final EVM output file
			-c, --two_abinitio EVM GFF3 file output with gene models supported by at least 2 abinitio predictors
			-d, --one_abinitio EVM GFF3 file output with gene models supported by at least 1 abinitio predictors
			-f, --evidencedata EVM GFF3 file output with gene models supported only by evidencedata
			-j, --proteins proteins fasta file with gene models supported by at least 1 abinitio predictors
			-a, --genome genome to annotate (fasta)
			-w, --weights weights file to pioritize some tools
      -g, --glimmerhmm glimmerhmm ab initio predictions (gff3)
      -i, --geneid geneid ab initio predictions (gff3)
      -l, --liftoff liftoff gene transfer (gff3)
      -m, --maker maker ab initio predictions (gff3)
      -b, --braker2 braker2 ab initio predictions (gff3)
			-p, --protein_alignments protein alignments (gff3)
			-t, --unstranded_transcriptome unstranded transcriptome (gff3)
      -s, --stranded_transcriptome stranded transcriptome (gff3)
			-r, --repeats repeat regions (gff)
			-h, --help

		Exemple: ./evidencemodeler.sh -e exec_dir -a genome -w weights -g glimmerhmm -i geneid -l liftoff -m maker -b braker2 -p protein_alignments -t unstranded_transcriptome -s stranded_transcriptome -r repeats -o GFF3_annotations -u EVM_annotations -c 2_abinitio_output -d 1_abinitio_output -f evidencedata_only_output -j 1_abinitio_output.proteins.fasta
		__EOF__
}

# Eval command line arguments given in input
ARGS=$(getopt -o "e:a:w:g:i:l:m:b:p:t:s:r:o:u:c:d:f:j:h" --long "exec_dir:,genome:,weights:,glimmerhmm:,geneid:,liftoff:,maker:,braker2:,gene_predictions:,protein_alignments:,unstranded_transcriptome:,stranded_transcriptome:,repeats:,output:,untreated_output:,two_abinitio:,one_abinitio:,evidencedata:,proteins:,help" -- $@ 2> /dev/null)

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
		-e|--exec_dir)
			EXEC_DIR=$2
			shift 2
			;;
		-o|--output)
			OUTPUT=$2
			shift 2
			;;
		-u|--untreated_output)
			EVM_OUTPUT=$2
			shift 2
			;;
		-c|--two_abinitio)
			TWO_ABINITIO=$2
			shift 2
			;;
		-d|--one_abinitio)
			ONE_ABINITIO=$2
			shift 2
			;;
		-f|--evidencedata)
			EVIDENCEDATA=$2
			shift 2
			;;
		-j|--proteins)
			PROTEINS=$2
			shift 2
			;;
		-a|--genome)
			GENOME=$2
			shift 2
			;;
		-w|--weights)
			WEIGHTS=$2
			shift 2
			;;
    -g|--glimmerhmm)
  		GLIMMERHMM=$2
  		shift 2
  		;;
    -i|--geneid)
    	GENEID=$2
    	shift 2
    	;;
    -l|--liftoff)
    	LIFTOFF=$2
    	shift 2
    	;;
    -m|--maker)
    	MAKER=$2
    	shift 2
    	;;
    -b|--braker2)
    	BRAKER2=$2
    	shift 2
    	;;
    -p|--protein_alignments)
    	PROTEIN=$2
    	shift 2
    	;;
    -t|--unstranded_transcriptome)
    	UNSTRANDED=$2
    	shift 2
    	;;
    -s|--stranded_transcriptome)
      STRANDED=$2
    	shift 2
    	;;
    -r|--repeats)
      REPEATS=$2
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

workdir=$PWD

mkdir -p ${EXEC_DIR}_tmp ${EXEC_DIR}
cat $UNSTRANDED $STRANDED > ${EXEC_DIR}_tmp/merged_transcriptomes.gff3

cp $WEIGHTS ${WEIGHTS}_tmp

if [[ -f $GENEID && -f $LIFTOFF ]]
then
	echo "cat $GLIMMERHMM $GENEID $LIFTOFF $MAKER $BRAKER2 > ${EXEC_DIR}_tmp/merged_ab_initio_predictions.gff"
	cat $GLIMMERHMM $GENEID $LIFTOFF $MAKER $BRAKER2 > ${EXEC_DIR}_tmp/merged_ab_initio_predictions.gff
elif [[ -f $GENEID && ! -f $LIFTOFF ]]
then
	sed -i "/Liftoff/d" ${WEIGHTS}_tmp
	echo "cat $GLIMMERHMM $GENEID $MAKER $BRAKER2 > ${EXEC_DIR}_tmp/merged_ab_initio_predictions.gff"
	cat $GLIMMERHMM $GENEID $MAKER $BRAKER2 > ${EXEC_DIR}_tmp/merged_ab_initio_predictions.gff
elif [[ ! -f $GENEID && -f $LIFTOFF ]]
then
	sed -i "/geneid/d" ${WEIGHTS}_tmp
	echo "cat $GLIMMERHMM $LIFTOFF $MAKER $BRAKER2 > ${EXEC_DIR}_tmp/merged_ab_initio_predictions.gff"
	cat $GLIMMERHMM $LIFTOFF $MAKER $BRAKER2 > ${EXEC_DIR}_tmp/merged_ab_initio_predictions.gff
elif [[ ! -f $GENEID && ! -f $LIFTOFF ]]
then
	sed -i "/geneid/d" ${WEIGHTS}_tmp
	sed -i "/Liftoff/d" ${WEIGHTS}_tmp
	echo "cat $GLIMMERHMM $MAKER $BRAKER2 > ${EXEC_DIR}_tmp/merged_ab_initio_predictions.gff"
	cat $GLIMMERHMM $MAKER $BRAKER2 > ${EXEC_DIR}_tmp/merged_ab_initio_predictions.gff
fi

echo "Weights file used : "
cat ${WEIGHTS}_tmp

sed "s/:.*-.*\texonerate\t/\t/" ${PROTEIN} > ${PROTEIN}_tmp

# we have to launch EVM per chromosome because then, in the EVM output, there is no information on the original chromosome, just the coordinates
# so we have to create a file per chromosome
# we retrieve the length of the bigger chromsome/contig to partition genome by chromosome/contig and not by segments inside each chromosome
bigger_chromosome=$( awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' $GENOME | sed "/>/d" | sort -n | tail -n1 )
cd ${EXEC_DIR}_tmp/

if [ -s $REPEATS ]
then
# partitionning by chromosome
echo "/usr/local/bin/partition_EVM_inputs.pl --genome $GENOME \
--gene_predictions ${EXEC_DIR}_tmp/merged_ab_initio_predictions.gff --protein_alignments ${PROTEIN}_tmp \
--transcript_alignments ${EXEC_DIR}_tmp/merged_transcriptomes.gff3 --repeats $REPEATS \
--segmentSize $bigger_chromosome --overlapSize 10000 --partition_listing ${EXEC_DIR}_tmp/partitions_list.out"

/usr/local/bin/partition_EVM_inputs.pl --genome $GENOME \
--gene_predictions ${EXEC_DIR}_tmp/merged_ab_initio_predictions.gff --protein_alignments ${PROTEIN}_tmp \
--transcript_alignments ${EXEC_DIR}_tmp/merged_transcriptomes.gff3 --repeats $REPEATS \
--segmentSize $bigger_chromosome --overlapSize 10000 --partition_listing ${EXEC_DIR}_tmp/partitions_list.out
else
	echo "/usr/local/bin/partition_EVM_inputs.pl --genome $GENOME \
	--gene_predictions ${EXEC_DIR}_tmp/merged_ab_initio_predictions.gff --protein_alignments ${PROTEIN}_tmp \
	--transcript_alignments ${EXEC_DIR}_tmp/merged_transcriptomes.gff3 \
	--segmentSize $bigger_chromosome --overlapSize 10000 --partition_listing ${EXEC_DIR}_tmp/partitions_list.out"

	/usr/local/bin/partition_EVM_inputs.pl --genome $GENOME \
	--gene_predictions ${EXEC_DIR}_tmp/merged_ab_initio_predictions.gff --protein_alignments ${PROTEIN}_tmp \
	--transcript_alignments ${EXEC_DIR}_tmp/merged_transcriptomes.gff3 \
	--segmentSize $bigger_chromosome --overlapSize 10000 --partition_listing ${EXEC_DIR}_tmp/partitions_list.out
fi

# create EVM command on each partition
echo "/usr/local/bin/write_EVM_commands.pl --genome $GENOME --weights ${WEIGHTS}_tmp \
--gene_predictions ${EXEC_DIR}_tmp/merged_ab_initio_predictions.gff --protein_alignments ${PROTEIN}_tmp \
--transcript_alignments ${EXEC_DIR}_tmp/merged_transcriptomes.gff3 \
--output_file_name evm.out  --partitions ${EXEC_DIR}_tmp/partitions_list.out >  ${EXEC_DIR}_tmp/commands.list"

/usr/local/bin/write_EVM_commands.pl --genome $GENOME --weights ${WEIGHTS}_tmp \
--gene_predictions ${EXEC_DIR}_tmp/merged_ab_initio_predictions.gff --protein_alignments ${PROTEIN}_tmp \
--transcript_alignments ${EXEC_DIR}_tmp/merged_transcriptomes.gff3 \
--output_file_name evm.out  --partitions ${EXEC_DIR}_tmp/partitions_list.out >  ${EXEC_DIR}_tmp/commands.list

echo "/usr/local/bin/execute_EVM_commands.pl ${EXEC_DIR}_tmp/commands.list | tee ${EXEC_DIR}_tmp/EVM_run.log"

/usr/local/bin/execute_EVM_commands.pl ${EXEC_DIR}_tmp/commands.list | tee ${EXEC_DIR}_tmp/EVM_run.log

echo "/usr/local/bin/convert_EVM_outputs_to_GFF3.pl --partitions ${EXEC_DIR}_tmp/partitions_list.out --output ${EXEC_DIR}_tmp/evm.out --genome $GENOME"

/usr/local/bin/convert_EVM_outputs_to_GFF3.pl --partitions ${EXEC_DIR}_tmp/partitions_list.out --output evm.out --genome $GENOME

cat ${EXEC_DIR}_tmp/*/evm.out > $workdir/$EVM_OUTPUT
cat ${EXEC_DIR}_tmp/*/evm.out.gff3 > $workdir/$OUTPUT

#---------------------------------------------------------
# filtering EVM output and generated the GFF3 files
#---------------------------------------------------------

# filter EVM output in 3 parts
# genes predicted by at least 2 ab initio predictors
# genes predicted by 1 ab initio predictor (and evidence data or not)
# genes predicted evidence data only
for chr in `grep ">" $GENOME | sed "s/>//"`
do
cd ${EXEC_DIR}_tmp/${chr}/
mkdir -p tmp
perl /scripts/evm-out.filter.pl ${EXEC_DIR}_tmp/${chr}/evm.out
echo "/usr/local/bin/EVM_to_GFF3.pl ${EXEC_DIR}_tmp/${chr}/tmp/evm.at_least_2_ABINITIO.out $chr > ${EXEC_DIR}/evm.2_abinitio.$chr.gff3"
/usr/local/bin/EVM_to_GFF3.pl ${EXEC_DIR}_tmp/${chr}/tmp/evm.at_least_2_ABINITIO.out $chr > ${EXEC_DIR}/evm.2_abinitio.$chr.gff3
echo "/usr/local/bin/EVM_to_GFF3.pl ${EXEC_DIR}_tmp/${chr}/tmp/evm.evidencedata_only.out $chr > ${EXEC_DIR}/evm.evidencedata.$chr.gff3"
/usr/local/bin/EVM_to_GFF3.pl ${EXEC_DIR}_tmp/${chr}/tmp/evm.evidencedata_only.out $chr > ${EXEC_DIR}/evm.evidencedata.$chr.gff3
echo "/usr/local/bin/EVM_to_GFF3.pl ${EXEC_DIR}_tmp/${chr}/tmp/evm.1_ABINITIO.out $chr > ${EXEC_DIR}/evm.1_abinitio.$chr.gff3"
/usr/local/bin/EVM_to_GFF3.pl ${EXEC_DIR}_tmp/${chr}/tmp/evm.1_ABINITIO.out $chr > ${EXEC_DIR}/evm.1_abinitio.$chr.gff3
done

cd ${EXEC_DIR}_tmp/

# 4 files generated
# tmp/evm.evidencedata_only.out
# tmp/evm.out.details
# tmp/evm.at_least_2_ABINITIO.out
# tmp/evm.1_ABINITIO.out

echo "Merging all EVM gff3 outputs ..."

echo "cat ${EXEC_DIR}/evm.2_abinitio.*.gff3 > ${workdir}/${TWO_ABINITIO}"
cat ${EXEC_DIR}/evm.2_abinitio.*.gff3 > ${workdir}/${TWO_ABINITIO}

echo "cat ${EXEC_DIR}/evm.1_abinitio.*.gff3 > ${workdir}/${ONE_ABINITIO}"
cat ${EXEC_DIR}/evm.1_abinitio.*.gff3 > ${workdir}/${ONE_ABINITIO}
gffread -y ${workdir}/${PROTEINS} -g $GENOME ${workdir}/${ONE_ABINITIO}
#gffread -x ${ONE_ABINITIO}.CDS.fasta -g $GENOME ${ONE_ABINITIO}

echo "cat ${EXEC_DIR}/evm.evidencedata.*.gff3 > ${workdir}/${EVIDENCEDATA}"
cat ${EXEC_DIR}/evm.evidencedata.*.gff3 > ${workdir}/${EVIDENCEDATA}

#rm -fr ${EXEC_DIR}_tmp ${WEIGHTS}_tmp ${PROTEIN}_tmp
