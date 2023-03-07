#!/bin/bash
basedir=${PWD}
for sp in synth_myco synth_ecoli synth_paeru; do
	for frag in 0.5 0.8 1; do 
		species=${sp}/${frag}
		echo ${species}
		#conda activate ete3
		ete3 compare -t ${basedir}/${species}/Comparison/GFs/CG/msa_delos/concat.tree \
		${basedir}/${species}/Roary/output/core_gene_alignment_renamed.tree \
		${basedir}/${species}/Comparison/GFs/CG/msa_genapi/concat.tree \
		${basedir}/${species}/Comparison/GFs/CG/msa_panaroo/concat.tree  \
		-r ${basedir}/${sp}/1/Comparison/GFs/CG/msa_prova/${sp}.tree --unrooted > ${basedir}/${species}/Comparison/ete3.log
	done
done
