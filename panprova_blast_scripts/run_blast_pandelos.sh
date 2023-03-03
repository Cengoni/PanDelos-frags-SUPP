#!/bin/bash

# Prepare files
basedir=${PWD}

# Check arguments
if [ $# -ne 1 ]
  then
    echo 'Usage: ./panprova_blast_scripts/run_blast_pandelos.sh synth_paeru|synth_ecoli|synth_myco'
    echo 'Choose one of the 3 species e.g: ./panprova_blast_scripts/run_blast_pandelos.sh synth_paeru'
    exit 1
fi

species=$1
filedir=${basedir}/${species}
dbbase=${basedir}/PANPROVA_${species}
coverage=50 

# create blastdb of complete PANPROVA genomes
for d in `ls  ${dbbase}/1/blastdb/*fna`;do
	makeblastdb -in ${d} -out ${d}.blastdb -dbtype nucl
	        
done

# For each fragmentation level
for fragmentation in 0.5 0.6 0.7 0.8 0.9 1 ; do

	pf=${filedir}/${fragmentation}/PanDelos/output/fragmented #this is a directory of fastas dirs
	dbdir=${dbbase}/1 # db to map against: original complete genomes 

	# if it doesn't exist create directory
	mkdir -p ${dbbase}/mappings/qcov_${coverage}

	# make sure file does not exist already 
	rm -f ${dbbase}/mappings/qcov_${coverage}/PanDelos_${fragmentation}_output.delos2ori
	rm -f ${dbbase}/mappings/qcov_${coverage}/PanDelos_${fragmentation}_output.ori2delos

	# for each genome
	for d in `ls ${pf}`; do
		wholegenome=$d  # whole genome file name
		ogenome=`echo $d | sed s/\.fasta//g` # without extension
		d="${pf}/${wholegenome}/results/gene_sequences.fna"

		# prepare for blast (change headers) and create blastdb
		python3 panprova_blast_scripts/make_blastdbin_pan.py $d ${d}_fixed.fna 
		makeblastdb -in ${d}_fixed.fna -out $pf/$wholegenome/results/gene_sequences.fna.blastdb -dbtype nucl

		# map in both directions
		blastn -query $dbdir/blastdb/${ogenome}_genes.fna -db $pf/$wholegenome/results/gene_sequences.fna.blastdb -perc_identity 0 -qcov_hsp_perc ${coverage} -outfmt 6 > $pf/$wholegenome/ori2delos.out 
		blastn -query $pf/$wholegenome/results/gene_sequences.fna_fixed.fna  -db $dbdir/blastdb/${ogenome}_genes.fna.blastdb -perc_identity 0 -qcov_hsp_perc ${coverage} -outfmt 6 > $pf/$wholegenome/delos2ori.out 

		# put all genome information into one file
		cat $pf/$wholegenome/ori2delos.out >> ${dbbase}/mappings/qcov_${coverage}/PanDelos_${fragmentation}_output.ori2delos
		cat $pf/$wholegenome/delos2ori.out >> ${dbbase}/mappings/qcov_${coverage}/PanDelos_${fragmentation}_output.delos2ori
	done

	# create file with bidirectional best hit
	python3 panprova_blast_scripts/bidirectional_best_blast.py  ${dbbase}/mappings/qcov_${coverage}/PanDelos_${fragmentation}_output.ori2delos ${dbbase}/mappings/qcov_${coverage}/PanDelos_${fragmentation}_output.delos2ori > ${dbbase}/mappings/qcov_${coverage}/PanDelos_${fragmentation}_bbmapping.csv

done
