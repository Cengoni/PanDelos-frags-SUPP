# PanDelos-frags-SUPP
This repository contains scripts to run the analyses performed in the methodological paper PanDelos-frags (tool available at [InfOmics/PanDelos-frags](https://github.com/InfOmics/PanDelos-frags))

These scripts assume that the tools used are installed and accessible in the current working directory.
The tools being compared are [PanDelos-frags](https://github.com/InfOmics/PanDelos-frags), [Roary](), [GenAPI](), and [Panaroo]().

Other tools being used in the analysis are [Prokka](), [Diamond](), [ete3](), [FastTree](), [mafft]().

## Analysis of synthetic data
### Data
In order to run the pangenomic tools, it is required that you have the synthetic data generated by [PANPROVA](https://doi.org/10.1093/bioinformatics/btac158). 
These are available in the directories `PANPROVA_synth_ecoli/`, `PANPROVA_synth_myco/`, and `/PANPROVA_synth_paeru`. 
Within these directories there are the directories `ori` which contains info on complete genomes generated by PANPROVA. 
The other directories are `0.5.`, `0.6`, `0.7`, `0.8`, `0.9`, and `1`, which represent the respective percentage of genome being 
virtually sequenced when artificially fragmenting the genomes, and contain the respective fragmented genomes (e.g. `genome_9_fr.fasta`)
and `.log` files generated by PANPROVA that are required to extract the genes.
Additionally the `survival_families` and `survival_gene.list` files containing information on  gene families and genes that are 
retained after fragmentation are present.These files correspond to the ground truth gene family composition at different levels of fragmentation.

### Analysis
The script `synthetic.sh` runs the four tools on these genomes.
Downstream analysis on the resulting pangenomes are performed in the first section of `tools_comparison.ipynb` and by running the script synthetic_trees.py
The notebook includes analyses on homology relationship, diffusivity, and core genes while the script contains the phylogenetic inference steps.

## Analysis of metagenomic data
### Data
In our study we consider metagenomic assembled genomes from three species belonging to different phyla that were reconstructed, annotated, and made available online by [Pasolli et al. 2019](10.1016/j.cell.2019.01.001).
Data is downloaded within the analysis scripts, however within the directories `abiotrophia_defectiva`, `bacteroides_nordii`, and `pseudomonas_aeruginosa` there are files containing links 
to download the genomes, namely `reference.txt` for the NCBI RefSeq genome relative to the species and `sequence_url_opendata.txt` for metagenomic assembled genomes.

### Analysis
The script `metagenome.sh` runs the four tools on these genomes.
Downstream analysis on the resulting pangenomes are performed in the second section of `tools_comparison.ipynb` and in `process_metagenomes.ipynb`
In the first notebook, diffusivity of the four tools is compared, in the second one we show the comparison between PanDelos-frags and Roary.

The directory `analysis_scripts/` contains all the internal scripts called by the main scripts.
