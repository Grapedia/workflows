# Markerlift
Markerlift is heavily based on SNPlift: Lift over SNP positions to match a new reference genome

Source: https://github.com/enormandeau/snplift

Coords input file is a tab delimited file with 4 columns
<chr> <start> <end> <name>

if <start> is equal to <end> it just run as SNPlift, otherwise SNPlift runs on <start> and <end> separetely. Results are mergedat the end.
```
----- INPUT
old genome  : <FASTA_FILE>
new genome  : <FASTA_FILE>
coords file : <TAB_DELIMITED_FILE>
working directory: <PATH>
window length: <NUMBER>
---- OUTPUT
output file : <TAB_DELIMITED_FILE>
```

## Example:

```
nextflow run main.nf --genome_old /path/to/Vitis12x.fa --genome_new /path/to/PN12Xv2_chloro_mito.fa --coords_old /path/to/SNPs_GrapeReseq.vcf --coords_new /path/to/coords_new.vcf --working_dir /path/to/
```

```nextflow.config``` should contains the ```--working_dir``` path as follow

``` 
docker {
    enabled = true
    runOptions = '-w /snplift -v /path/to:/path/to'
}
```
