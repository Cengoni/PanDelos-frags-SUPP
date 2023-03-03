#!/bin/bash

# Check arguments
if [ $# -ne 1 ]
  then
    echo 'Usage: ./panprova_blast_scripts/run_extract_survival_genes.sh synth_paeru|synth_ecoli|synth_myco'
    echo 'Choose one of the 3 species e.g: ./panprova_blast_scripts/run_extract_survival_genes.sh synth_myco'
    exit 1
fi


basedir=${PWD}
species=$1

# For each fragmentation level
for fragmentation in 0.5 0.6 0.7 0.8 0.9 1; do
	
	# Create directory to store fragmented genes, will be used to build blastdb of each fragmentation level
	mkdir -p ${basedir}/PANPROVA_${species}/${fragmentation}/blastdb
	
	# Extract genes at each fragmentation level
	for i in `ls  ${basedir}/PANPROVA_${species}/ori/*fna`;do
		genome=$(basename $i .fna)
		number=${genome#*_}
		echo $genome
		echo $i
		# create fragmented genes such that blast can be run on them
		python3 ${basedir}/panprova_blast_scripts/extract_survival_genes.py ${basedir}/PANPROVA_${species}/${fragmentation}/${genome}_fr.log ${basedir}/PANPROVA_${species}/ori/${genome}.fna ${basedir}/PANPROVA_${species}/ori/${genome}.gff ${number} ${basedir}/PANPROVA_${species}/${fragmentation}/blastdb/${genome}_fr_genes.fna

	done
done
