#!/bin/bash

# Read in arguments
ifile=$1 # input file with links
odir=$2 # output dir where genomes will be stored

# Prepare output directory
mkdir -p $odir
rm -f $odir/*

# Download reference genome
for e in `cat $ifile`; do
    wget -P "$odir" $e
    e=`basename $e`
    gunzip $odir/$e
done
