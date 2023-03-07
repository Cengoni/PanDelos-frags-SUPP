#!/bin/bash

basedir=${PWD}
species=$1
faa=$2
out=$3

echo ${basedir}
echo $species
echo ${faa}

mafft ${faa} > ${out}

