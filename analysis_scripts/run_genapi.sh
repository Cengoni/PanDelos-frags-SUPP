#!/bin/bash

# Read in arguments
basedir=${PWD}
species=$1

# Prepare input/output folder
mkdir -p ${basedir}/${species}/GenAPI

# Put all gff files produced by Prokka together
cp ${basedir}/${species}/Prokka/prokka_gffs/* ${basedir}/${species}/GenAPI

# Run tool
cd ${basedir}/${species}/GenAPI
genapi -n genapi --threads 10 # expects genapi installed and in path
cd ${basedir}
