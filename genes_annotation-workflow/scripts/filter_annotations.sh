#!/bin/bash

# Subsequently, the raw gene models were quality filtered. Gene models only supported by ab initio predictors were kept
# if at least two gene prediction programs predicted them, if the start and stop codon was present and if the gene length
# was equal or larger than 300 bp. However, the only ab initio supported gene models not matching these constraints were
# kept if they had a database hit with the UniProt/SwissProt or NCBI non-redundant database. For that, a blastp search of
# the protein sequences against the two databases was run allowing hits with an e-value <1e-6. Of the gene models only
# supported by evidence data or by the VCost.v3 lifted annotation, those gene models with missing start and stop and a
#  gene length smaller than 300 bp were discarded.

# Source for the AGAT scripts : https://github.com/NBISweden/AGAT/

function usage {
	cat <<-__EOF__
		Usage:
			./filter_annotations.sh -g genome -i gff3 -r EVM_annotations -a two_abinitio -b one_abinitio -c evidencedata -o gff3_filtered -d one_abinitio_proteins -n NR_proteins -u Uniprot_proteins -t threads [-h]

		Description:
      Ab initio supported gene models were kept if :
        * there are predicted by at least 2 ab initio predictors
        * if the start and stop condon was present
        * if gene length >= 300 bp
      Ab initio supported gene models not matching these constraints were kept if :
        * they had a database hit with the Uniprot/SwissProt or NCBI NR database (blastp hits with en evalue <1e-6)
      Gene models only supported by evidence data or by lifted annotation were kept if :
        * if the start and stop condon was present
        * if gene length >= 300 bp

		Arguments:
      -g, --genome assembly to annotate (fasta)
			-i, --input raw EVM gff3 file
			-r, --raw raw EVM output file
			-a, --two_abinitio gene models supported by at least 2 ab initio gff3 file
			-b, --one_abinitio gene models supported by at least 1 ab initio gff3 file
			-c, --evidencedata gene models supported by evidencedata only gff3 file
			-d, --EVM_proteins proteins fasta file with gene models supported by at least 1 abinitio predictors
			-o, --output filtered gff3 file
			-n, --NR_prot NR_proteins
			-u, --uniprot Uniprot_proteins
			-t, --threads Threads number
			-h, --help

		Exemple: ./filter_annotations.sh -g genome.fasta -i evm.annotations.gff3 -r EVM_annotations.out -a two_abinitio.gff3 -b one_abinitio.gff3 -c evidencedata.gff3 -o evm.annotations.filtered.gff3 -d one_abinitio_proteins.fasta -n NR_proteins.fasta  -u Uniprot_proteins.fasta -t threads
		__EOF__
}

# Eval command line arguments given in input
ARGS=$(getopt -o "g:i:r:o:a:d:b:c:n:u:t:h" --long "genome:,input:,raw:,output:,two_abinitio:,one_abinitio:,evidencedata:,EVM_proteins:,NR_prot:,uniprot:,threads:,help" -- $@ 2> /dev/null)

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
		-i|--input)
			INPUT=$2
			shift 2
			;;
		-r|--raw)
			EVM_OUTPUT=$2
			shift 2
			;;
		-o|--output)
			OUTPUT=$2
			shift 2
			;;
		-a|--two_abinitio)
			TWO_ABINITIO=$2
			shift 2
			;;
		-b|--one_abinitio)
			ONE_ABINITIO=$2
			shift 2
			;;
		-d|--proteins)
			EVM_PROTEINS=$2
			shift 2
			;;
		-c|--evidencedata)
			EVIDENCEDATA=$2
			shift 2
			;;
		-n|--NR_prot)
			NR_PROT=$2
			shift 2
			;;
		-u|--uniprot)
			UNIPROT=$2
			shift 2
			;;
		-t|--threads)
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

###########---------------
# 1. Ab initio supported gene models filtering
###########---------------

# discard genes less than 300 bp
echo "perl $( which agat_sp_filter_gene_by_length.pl ) --gff $TWO_ABINITIO --size 300 --test ">" -o ${TWO_ABINITIO}.length_filtered"
perl $( which agat_sp_filter_gene_by_length.pl ) --gff $TWO_ABINITIO --size 300 --test ">" -o ${TWO_ABINITIO}.length_filtered

# We then check that the gene models are complete (has a start and stop codon), and remove the incomplete ones.
echo "perl $( which agat_sp_filter_incomplete_gene_coding_models.pl ) -gff ${TWO_ABINITIO}.length_filtered --fasta $GENOME --output ${TWO_ABINITIO}.complete"
perl $( which agat_sp_filter_incomplete_gene_coding_models.pl ) -gff ${TWO_ABINITIO}.length_filtered --fasta $GENOME --output ${TWO_ABINITIO}.complete

# ###########---------------
# # All the genes Ab initio supported gene models not matching these constraints are aligned on databases
# ###########---------------

echo "diamond makedb --in $NR_PROT --threads $THREADS -d $NR_PROT"

diamond makedb --in $NR_PROT --threads $THREADS -d $NR_PROT

echo "diamond makedb --in $UNIPROT --threads $THREADS -d $UNIPROT"

diamond makedb --in $UNIPROT --threads $THREADS -d $UNIPROT

echo "diamond blastp -d $NR_PROT \
-q $EVM_PROTEINS --sensitive --top 1 -e 1e-6 \
-o FINAL_OUTPUT/EVM_prot_vs_NR_database.txt -f 6 --threads $THREADS"

# align $EVM_PROTEINS
diamond blastp -d $NR_PROT \
-q $EVM_PROTEINS --sensitive --top 1 -e 1e-6 \
-o FINAL_OUTPUT/EVM_prot_vs_NR_database.txt -f 6 --threads $THREADS

echo "diamond blastp -d $UNIPROT \
-q $EVM_PROTEINS --sensitive --top 1 -e 1e-6 \
-o FINAL_OUTPUT/EVM_prot_vs_Uniprot_database.txt -f 6 --threads $THREADS"

diamond blastp -d $UNIPROT \
-q $EVM_PROTEINS --sensitive --top 1 -e 1e-6 \
-o FINAL_OUTPUT/EVM_prot_vs_Uniprot_database.txt -f 6 --threads $THREADS

cut -f1 FINAL_OUTPUT/EVM_prot_vs_Uniprot_database.txt | sort -u > FINAL_OUTPUT/list_to_keep_from_Uniprot.txt
cut -f1 FINAL_OUTPUT/EVM_prot_vs_NR_database.txt | sort -u > FINAL_OUTPUT/list_to_keep_from_NR.txt

cat FINAL_OUTPUT/list_to_keep_from_Uniprot.txt FINAL_OUTPUT/list_to_keep_from_NR.txt | sort -u > FINAL_OUTPUT/list_to_keep.FINAL.txt
if [ -f ${ONE_ABINITIO}_report.txt ]
then
rm ${ONE_ABINITIO}_report.txt
fi

echo "perl $( which agat_sp_filter_feature_from_keep_list.pl ) --gff ${ONE_ABINITIO} --keep_list FINAL_OUTPUT/list_to_keep.txt  --output ${ONE_ABINITIO}.final"
perl $( which agat_sp_filter_feature_from_keep_list.pl ) --gff ${ONE_ABINITIO} --keep_list FINAL_OUTPUT/list_to_keep.FINAL.txt  --output ${ONE_ABINITIO}.final

if [ -f ${ONE_ABINITIO}.final ]
then
    if [ -s ${ONE_ABINITIO}.final ]
		then
			rm FINAL_OUTPUT/list_to_keep_from_Uniprot.txt FINAL_OUTPUT/list_to_keep_from_NR.txt FINAL_OUTPUT/EVM_prot_vs_Uniprot_database.txt FINAL_OUTPUT/EVM_prot_vs_NR_database.txt
		fi
fi

# ###########---------------
# # Gene models only supported by evidence data or by lifted annotation
# ###########---------------

if [ -f $EVIDENCEDATA ]; then
	# discard genes less than 300 bp
	echo "perl $( which agat_sp_filter_gene_by_length.pl ) --gff $EVIDENCEDATA --size 300 --test '>' -o ${EVIDENCEDATA}.length_filtered"
	perl $( which agat_sp_filter_gene_by_length.pl ) --gff $EVIDENCEDATA --size 300 --test ">" -o ${EVIDENCEDATA}.length_filtered

	# We then check that the gene models are complete (has a start and stop codon), and remove the incomplete ones.
	echo "perl $( which agat_sp_filter_incomplete_gene_coding_models.pl ) -gff ${EVIDENCEDATA}.length_filtered --fasta $GENOME --output ${EVIDENCEDATA}.complete"
	perl $( which agat_sp_filter_incomplete_gene_coding_models.pl ) -gff ${EVIDENCEDATA}.length_filtered --fasta $GENOME --output ${EVIDENCEDATA}.complete
fi

# ###########---------------
# # Merge des 3 GFF3
# ###########---------------

cat ${TWO_ABINITIO}.complete ${EVIDENCEDATA}.complete ${ONE_ABINITIO}.final > $OUTPUT

if [ -f $OUTPUT ]
then
    if [ -s $OUTPUT ]
		then
			rm ${TWO_ABINITIO}.complete ${EVIDENCEDATA}.complete ${ONE_ABINITIO}.final ${EVIDENCEDATA}.length_filtered ${TWO_ABINITIO}.length_filtered
		fi
fi
