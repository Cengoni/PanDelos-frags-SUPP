#!/bin/bash

basedir=${PWD}

# Check arguments
if [ $# -ne 1 ]
  then
    echo 'Usage: ./metagenome.sh abiotrophia_defectiva|bacteroides_nordii|pseudomonas_aeruginosa'
    echo 'Choose one of the 3 species e.g: ./metagenome.sh abiotrophia_defectiva'
    exit 1
fi

species=$1 

# Download metagenome files
./analysis_scripts/download_frag.sh ${basedir}/${species}/sequence_url_opendata.txt ${basedir}/${species}/input/fragmented/

# Download reference file
./analysis_scripts/download_ref.sh ${basedir}/${species}/reference.txt ${basedir}/${species}/input/complete/

# Run Prokka (required by Roary, Panaroo, and GenAPI)
echo 'Running Prokka...'
${basedir}/analysis_scripts/run_prokka.sh ${species}

# Run Roary
echo 'Running Roary...'
${basedir}/analysis_scripts/run_roary.sh ${species}

# Run Genapi
${basedir}/analysis_scripts/run_genapi.sh ${species}

# Run Panaroo
${basedir}/analysis_scripts/run_panaroo.sh ${species}

# Run PanDelos-frags
echo 'Running PanDelos-frags...'
${basedir}/analysis_scripts/run_pandelos.sh ${species}
