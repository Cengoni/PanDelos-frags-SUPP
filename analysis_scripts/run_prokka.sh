#!/bin/bash

# Read in arguments
basedir=${PWD}
species=$1 

# Prepare output directory
mkdir -p ${basedir}/${species}/Prokka

# Run prokka on fragmented genomes
for x in $(ls ${basedir}/${species}/input/fragmented); do  
	prokka --kingdom Bacteria --outdir ${basedir}/${species}/Prokka/prokka_$x ${basedir}/${species}/input/fragmented/$x
done

# Run prokka on complete genomes (if they exist)
if [[ -d ${basedir}/${species}/input/complete ]]
then
    for x in $(ls ${basedir}/${species}/input/complete); do
    	prokka --kingdom Bacteria --outdir ${basedir}/${species}/Prokka/prokka_$x ${basedir}/${species}/input/complete/$x
    done
fi

# Put all gff files together
mkdir ${basedir}/${species}/Prokka/prokka_gffs

for x in $(ls ${basedir}/${species}/Prokka); do
	bs=$(basename $x)
	mv ${basedir}/${species}/Prokka/$x/*gff  ${basedir}/${species}/Prokka/prokka_gffs/$bs.gff
done
