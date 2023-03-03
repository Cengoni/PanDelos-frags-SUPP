#!/bin/bash

# Read in arguments
basedir=${PWD}
ilist="ilist.csv"   
species=$1 

#tool=${basedir}/../pandelos-frags

# Prepare output folders
mkdir -p ${basedir}/${species}/PanDelos/output
rm -rf ${basedir}/${species}/PanDelos/output/*

# Prepare input file
if [ -f ${basedir}/${species}/PanDelos/${ilist} ]; then
	rm ${basedir}/${species}/PanDelos/${ilist}
fi
 
# Write fragmented genome paths to ilist file
if [ -d ${basedir}/${species}/input/fragmented ]; then
	for e in $(ls ${basedir}/${species}/input/fragmented); do 
		echo "${basedir}/${species}/input/fragmented/${e},fragmented" >> ${basedir}/${species}/PanDelos/${ilist}
	done
fi

# Write complete genome paths to ilist file
if [ -d ${basedir}/${species}/input/complete ]; then
	rm -rf ${basedir}/${species}/input/complete/*fai #make sure .fai files have not been created already
	for e in $(ls ${basedir}/${species}/input/complete); do 
		echo "${basedir}/${species}/input/complete/${e},complete" >> ${basedir}/${species}/PanDelos/${ilist}
	done
fi
 
# Run tool
pandelos-frags ${basedir}/${species}/PanDelos/${ilist} ${basedir}/${species}/PanDelos/output
