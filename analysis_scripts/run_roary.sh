#!/bin/bash

# Read in arguments
basedir=${PWD}
species=$1

# Prepare output
mkdir -p ${basedir}/${species}/Roary

# Run Roary
roary -f ${basedir}/${species}/Roary/output -e -n -v ${basedir}/${species}/Prokka/prokka_gffs/*.gff
