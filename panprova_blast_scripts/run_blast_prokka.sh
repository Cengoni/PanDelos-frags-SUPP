#!/bin/bash

# Prepare files
basedir=${PWD}

# Check arguments
if [ $# -ne 1 ]
  then
    echo 'Usage: ./panprova_blast_scripts/run_blast_prokka.sh synth_paeru|synth_ecoli|synth_myco'
    echo 'Choose one of the 3 species e.g: ./panprova_blast_scripts/run_blast_prokka.sh synth_paeru'
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
for fragmentation in 0.5 0.6 0.7 0.8 0.9 1; do
	pf=${filedir}/${fragmentation}/Prokka/fasta 
		
        if [[ -d ${pf} ]];then
		rm -r ${pf}
	fi
		
	mkdir -p ${pf}
	for x in $(ls -d ${filedir}/${fragmentation}/Prokka/*.fasta); do 
		bs=$(basename $x)
		bs=`echo $bs | sed s/prokka_//g` 
	 	bs=`echo $bs | sed s/.fasta//g` 
		cp ${x}/*ffn  ${pf}/$bs.fa
	done
		
	dbdir=${dbbase}/1 # db to map against: original complete genomes 	 

	# if it doesn't exist create directory 
	mkdir -p ${dbbase}/mappings/qcov_${coverage}

	# make sure file does not exist already 
	rm -f ${dbbase}/mappings/qcov_${coverage}/Prokka_${fragmentation}_output.delos2ori
	rm -f ${dbbase}/mappings/qcov_${coverage}/Prokka_${fragmentation}_output.ori2delos

	# for each genome
	for d in `ls ${pf}`;do
		#ogenome=`echo $d | sed s/\.fa//g` 
		ogenome=`echo "$d" | cut -f 1 -d '.'`
		echo $ogenome
		d="${pf}/$d"
		echo $d

		#prepare for blast (change headers) and create blastdb
		mkdir -p $pf/$ogenome/blastdb
		python3 panprova_blast_scripts/make_blastdbin_prokka.py $pf/$ogenome.fa $pf/$ogenome/blastdb/gene_sequences.fna
		makeblastdb -in $pf/$ogenome/blastdb/gene_sequences.fna -out $pf/$ogenome/blastdb/gene_sequences.fna.blastdb -dbtype nucl

		# map in both directions
		blastn -query $dbdir/blastdb/${ogenome}_genes.fna -db $pf/$ogenome/blastdb/gene_sequences.fna.blastdb  -perc_identity 0 -qcov_hsp_perc ${coverage} -outfmt 6 > $pf/$ogenome/ori2delos.out
		blastn -query $pf/$ogenome/blastdb/gene_sequences.fna  -db $dbdir/blastdb/${ogenome}_genes.fna.blastdb -perc_identity 0 -qcov_hsp_perc ${coverage} -outfmt 6 > $pf/$ogenome/delos2ori.out  

		# put all genomes into one file
		cat $pf/$ogenome/ori2delos.out >> ${dbbase}/mappings/qcov_${coverage}/Prokka_${fragmentation}_output.ori2delos
		cat $pf/$ogenome/delos2ori.out >> ${dbbase}/mappings/qcov_${coverage}/Prokka_${fragmentation}_output.delos2ori
	done

	# create file with bidirectional best hit
	python3 panprova_blast_scripts/bidirectional_best_blast.py ${dbbase}/mappings/qcov_${coverage}/Prokka_${fragmentation}_output.ori2delos ${dbbase}/mappings/qcov_${coverage}/Prokka_${fragmentation}_output.delos2ori > ${dbbase}/mappings/qcov_${coverage}/Prokka_${fragmentation}_bbmapping.csv

done
