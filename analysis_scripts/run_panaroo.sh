#!/bin/bash

# Read in arguments
basedir=${PWD}
species=$1 

# Prepare output folder
mkdir -p ${basedir}/${species}/Panaroo

# Run tool
panaroo -i ${basedir}/${species}/Prokka/prokka_gffs/*.gff -o ${basedir}/${species}/Panaroo --clean-mode strict # expects panaroo installed and in path
