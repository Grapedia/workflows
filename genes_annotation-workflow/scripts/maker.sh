#!/bin/bash

usage() {
	cat <<-__EOF__
		Usage:
			./maker.sh -a assembly -t transcriptome -p prot1,prot2 -o organism -n nThreads -d outdir [-h]

		Description:
			Ab initio prediction with SNAP using MAKER

		Arguments:
			-a, --assembly Assembly of the genome to annotate in FASTA format
			-t, --transcriptome Transcriptome in GFF3 format
			-p, --protein Protein sequences files in FASTA format (each separated by a comma)
			-o, --organism Organism name (prefix name for intermediate and output files)
			-n, --nThreads Number of threads/cores to use
			-d, --outdir Output directory
			-h, --help

		Exemple: ./maker.sh -a assembly.fasta -t transcriptome.gff3 -p prot1.fasta,prot2.fasta -o PN40024_v4 -n 24 -d /outdir
		__EOF__
}

config_files () {
	round=$1
	hmm_file=$2

	# Generate configuration files and add path to genome assembly
	maker -CTL
	sed -i "s@^genome=@genome=${ASSEMBLY}@g" maker_opts.ctl

	# Change 'max_dna_len' and 'split_hit' parameters
	sed -i "s/^max_dna_len=100000/max_dna_len=300000/g" maker_opts.ctl
	sed -i "s/^split_hit=10000/split_hit=20000/g" maker_opts.ctl

	# Skip Repeat masking
	sed -i "s/^model_org=all/model_org=/g" maker_opts.ctl

	# Population variables with input files : RNA-Seq stranded transcriptome ...
	# ... and Viridiplantae and Eudicotyledones protein sequences
	sed -i "s@^est_gff=@est_gff=${TRANSCRIPTOME}@g" maker_opts.ctl
	sed -i "s@^protein=@protein=${PROTEIN}@g" maker_opts.ctl

	if [[ ${round} == 'init' ]]
	then
		# Pass parameters 'est2genome' and 'protein2genome' values to 1
		sed -i "s/^est2genome=0/est2genome=1/g" maker_opts.ctl
		sed -i "s/^protein2genome=0/protein2genome=1/g" maker_opts.ctl
	elif [[ ${round} == 'snap' ]]
	then
		# Populate variable snaphmm with input file
		sed -i "s@^snaphmm=@snaphmm=${hmm_file}@g" maker_opts.ctl
	fi
}

run_maker () {
	prefix=$1
	# Run maker
	maker -c ${THREADS} -base ${prefix} --ignore_nfs_tmp -TMP /tmp
}

snap_training () {
	prefix=$1

	# 1. Merge all GFF3 files into one - output: PN40024_v4_init_prediction.all.gff
	gff3_merge -d ${OUTDIR}/MAKER/${prefix}.maker.output/${prefix}_master_datastore_index.log

	# 2. Convert gff to ZFF format - output: genome.ann and genome.dna
	maker2zff -n ${prefix}.all.gff

	# 3. Filter input gene models, capture genomic sequences surrounding each ...
	# ... model locus and uses those segments to produce the HMM file
	fathom genome.ann genome.dna -categorize 1000
	fathom uni.ann uni.dna -export 1000 -plus

	# 4. Parameters estimation
	mkdir -p ${prefix}_params/
	cd ${prefix}_params/
	forge ../export.ann ../export.dna
	cd ..

	# 5. Build HMM file
	hmm-assembler.pl ${prefix}_snap_training ${prefix}_params > ${prefix}.hmm
}

# Eval command line arguments given in input
ARGS=$(getopt -o "a:t:p:o:n:d:h" --long "assembly:,transcriptome:,protein:,organism:,nThreads:,outdir:,help" -- $@ 2> /dev/null)

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
		-p|--protein)
			PROTEIN=$2
			shift 2
			;;
		-o|--organism)
			PREFIX=$2
			shift 2
			;;
		-n|--nThreads)
			THREADS=$2
			shift 2
			;;
		-d|--outdir)
			OUTDIR=$2
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

mkdir -p ${OUTDIR}/MAKER/
# cd ${OUTDIR}/MAKER/

# Initial prediction with MAKER using BLAST and Exonerate to align evidence data like proteins ...
# ... and RNA-Seq stranded transcriptomes
echo "Initial prediction"
config_files "init"
run_maker "${PREFIX}_init_prediction"

# SNAP
# First round
echo "First SNAP prediction"
mkdir -p SNAP_first_round/
cd SNAP_first_round/
snap_training "${PREFIX}_init_prediction"

cd ..
config_files "snap" SNAP_first_round/${PREFIX}_init_prediction.hmm
run_maker "${PREFIX}_snap_first_prediction"

# Second round
echo "Second SNAP prediction"
mkdir -p SNAP_second_round/
cd SNAP_second_round/
snap_training "${PREFIX}_snap_first_prediction"

cd ..
config_files "snap" SNAP_second_round/${PREFIX}_snap_first_prediction.hmm
run_maker "${PREFIX}_snap_second_prediction"

# Final merge
# The option '-n' is used to produce a GFF file without genome sequences
gff3_merge -n -d ${PREFIX}_snap_second_prediction.maker.output/${PREFIX}_snap_second_prediction_master_datastore_index.log
