#!/bin/bash

basedir=${PWD}


# Check arguments
if [ $# -ne 1 ]
  then
    echo 'Usage: ./metagenome.sh synth_paeru|synth_ecoli|synth_myco'
    echo 'Choose one of the 3 species e.g: ./synthetic.sh synth_myco'
    exit 1
fi

species=$1 

# Compute for the available fragmentation levels
for frag in 0.5 0.6 0.7 0.8 0.9 1
do

	# Prepare output directory
	mkdir -p ${basedir}/${species}/${frag}/input/fragmented/
	rm -f ${basedir}/${species}/${frag}/input/fragmented/*

	# Get data from PANPROVA and put in th expected directory
	panprova_dir= ${basedir}/PANPROVA_${species}/${frag}/ # Location of PANPROVA data
	cp ${panprova_dir}/*fasta ${basedir}/${species}/${frag}/input/fragmented/

	# Run Prokka
	echo 'Running Prokka...'
	${basedir}/analysis_scripts/run_prokka.sh ${species}/${frag}

	# Run Roary
	echo 'Running Roary...'
	${basedir}/analysis_scripts/run_roary.sh ${species}/${frag}


	# Run Genapi
	echo 'Running GenAPI...'
	${basedir}/analysis_scripts/run_genapi.sh ${species}/${frag}

	# Run Panaroo
	echo 'Running Panaroo...'
	${basedir}/analysis_scripts/run_panaroo.sh ${species}/${frag}


	# Run PanDelos
	echo 'Running PanDelos-frags...'
	${basedir}/analysis_scripts/run_pandelos.sh ${species}
done

