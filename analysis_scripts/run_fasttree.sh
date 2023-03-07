#!/bin/bash

basedir=${PWD}
species=$1
faa=$2
out=$3 

./../FastTree  -gtr -nt ${faa} > ${out}

