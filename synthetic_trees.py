#!/usr/bin/python3
import subprocess
import pandas as pd
import os
import re

# Function used to read in PANPROVA GFs of complete genomes
# It returns a set of core gene families (represented as tuples)
def read_panprova_CG(species):
    basedir=os.getcwd()
    ipanprova = basedir+'/PANPROVA_'+species+'/1/survival_families'  

    gf_cg_prova = set()
    for line in open(ipanprova,'r'):
        cc = line.strip().split(' ')[1:]
            
        genomes = set()
        cc_filt=set()
        for c in cc:
            genome = c.split(',')[0].lstrip('(')
            if genome not in genomes:
                genomes.add(genome)
                cc_filt.add((c))

        size = len(genomes)

        # core GF
        if size == 10:
            gf_cg_prova.add(tuple(cc_filt))

    return gf_cg_prova


# Function used to read in PanDelos GFs of complete genomes
# It returns a set of core gene families (represented as tuples)
def read_pandelos_CG(species, frag):
    basedir=os.getcwd()

    species= species+'/'+frag
    delosdir = basedir+'/'+species+'/PanDelos'

    with open(delosdir+"/output/output.clus") as f:
        lines = f.readlines()

    cg_delos=set()
    for l in lines: 
        l = l.rstrip()
        l = l.split(' ')
        genome_set = set()
        selected_genes = []
        for gene in l:
            genome = gene.split(':')[0]     
            if genome in genome_set:
                continue
            else:
                genome_set.add(genome)
                selected_genes.append(gene)

        if len(selected_genes)==10:
            cg_delos.add(tuple(selected_genes)) 

    return cg_delos


# Function used to read in GenAPI core gene families
# It returns a dictionary gene family 2 genes
def read_genapi_CG(species,frag):
    basedir=os.getcwd()
    species= species+'/'+frag
    itool = basedir+'/'+species+'/GenAPI/clustered_genes_genapi.ffn.clstr'
    gf2genes=dict()
    with open(itool,'r') as f:
        lines = f.readlines()
        GF_name=lines[0].lstrip('>').split(' ')[1].rstrip()
        GF_list=[]
        genomes=set()
        CG=False

        for l in lines[1:]:
            if l.startswith('>'):
                if len(GF_list) == 10 and CG == True:
                    gf2genes[GF_name]=GF_list
                GF_name=l.lstrip('>').split(' ')[1].rstrip()
                CG = True
                genomes=set()
                GF_list=[]
            else:
                gene= l.split('>')[1].split(':')[0]
                genome=gene.split('_')[0]
                if genome in genomes:
                    CG=False
                genomes.add(genome)
                GF_list.append(gene)

    return gf2genes


# Function used to read in Panaroo core gene families
# It returns a dictionary gene family 2 genes
def read_panaroo_CG(species,frag):
    basedir=os.getcwd()
    species= species+'/'+frag
    itool = basedir+'/'+species+'Panaroo/gene_presence_absence.csv'
    gf2genes=dict()
    with open(itool,'r') as f:
        lines = f.readlines()[1:]  
        for line in lines:
            CG=True
            name=line.split(',')[0]
            cc = line.strip().split(',')[3:]
            cc=[g for g in cc if g != '']
            if len(cc)==10:
                for c in cc:
                    if len(c.split(';'))>1 or len(c.split('_'))>2:
                        CG=False
                if CG:
                    gf2genes[name]=cc
           
    return gf2genes


# Function used to read in PANPROVA complete gene sequences
# It returns a dictionary gene name 2 sequence
def read_panprova_genes(species):
    basedir=os.getcwd()

    # read all genes
    gene2seq=dict()
    ipanprova_blastdb = basedir+'/PANPROVA_'+species+'/1/blastdb/'
    for file in os.listdir(ipanprova_blastdb):
        if file.endswith('fna'):
            with open(ipanprova_blastdb+'/'+file) as f:
                lines = f.readlines()
                for l in lines:
                    if l.startswith('>'):
                        gene_name = l[1:].rstrip()
                    else: 
                        gene2seq[gene_name]=l.rstrip()
    
    return gene2seq


# Function used to read in PanDelos-frag gene sequences
# It returns a dictionary gene name 2 sequence
def read_pandelos_genes(species,frag):
    basedir=os.getcwd()
    species= species+'/'+frag

    delosdir = basedir+'/'+species+'/PanDelos'

    frag_dir = delosdir+'/output/fragmented_coordinates/'
    records_pan = {}
    for d in os.listdir(frag_dir):
        with open(frag_dir+d+'/'+"coordinates_frag.sam") as f:
            lines = f.readlines()
            for l in lines:
                if l.startswith('@'):
                    continue
                l = l.split('\t')
                name = l[0]
                seq = l[9]
                records_pan[name] = seq
    return records_pan

# Function used to read in GenAPI and Panaroo gene sequences (obtained from Prokka)
# It returns a dictionary gene name 2 sequence
def read_prokka_genes(species, frag):
    basedir = os.getcwd()
    species = species+'/'+frag
    prokkadir = basedir+'/'+species+'/Prokka' 
    gene2seq = dict()

    for file in os.listdir(prokkadir):
        all_files = os.listdir(prokkadir+'/'+file)
        ffn_file = [f for f in all_files if f.endswith('ffn')][0]
        with open(prokkadir+'/'+file+'/'+ffn_file) as f:
            lines=f.readlines()
            gene=lines[0].split(' ')[0]
            seq=''
            for l in lines[1:]:
                if l.startswith('>'):
                    gene2seq[gene.lstrip('>')]=seq
                    gene=l.split(' ')[0]
                    seq=''
                else:
                    seq+=l.strip()
            gene2seq[gene.lstrip('>')]=seq
    
    return gene2seq

# This function writes to fasta files the sequences of genes belonging to core gene families
# Each file represents a gene family.
def write_fasta_to_file(tool,gf_cg,gene2seq,species,frag):
    basedir=os.getcwd()
    CG_dir = basedir+'/'+species+'/'+frag+'/Comparison/CG'

    os.makedirs(CG_dir,exist_ok=True)
    os.makedirs(CG_dir+'/fasta_'+tool,exist_ok=True)

    # write to file
    i = 0 
    for family in gf_cg:
        with open(CG_dir+'/fasta_'+tool+'/fam_'+str(i)+'.fa','w') as f:
            i += 1
            for gene in family: 
                f.write('>'+str(gene)+'\n')
                f.write(str(gene2seq[gene])+'\n')


# Run Mafft on fasta files
def run_mafft(tool,species, frag):
    basedir=os.getcwd()
    CG_dir = basedir+'/'+species+'/'+frag+'/Comparison/CG'

    cg_fa = CG_dir+'/fasta_'+tool+'/'
    msa = CG_dir+'/msa_'+tool+'/'
    os.makedirs(msa,exist_ok=True)

    for f in os.listdir(CG_dir+'/fasta_'+tool+'/'):
        if f not in os.listdir(CG_dir+'/msa_'+tool+'/'):
            ret=subprocess.call(['bash',basedir+"/analysis_scripts/run_mafft.sh",species,cg_fa+f,msa+f],stdout=subprocess.DEVNULL)
    
# Concatenate alignments of multiple gene families 
def concatenate_alignments(tool,species, frag):
    basedir=os.getcwd()
    CG_dir = basedir+'/'+species+'/'+frag+'/Comparison/CG'
    msa = CG_dir+'/msa_'+tool+'/'

    genome2seqs=dict()
    for file in os.listdir(msa):
        with open(msa+file) as f:
            lines=f.readlines()

            if tool=='prova':
                genome=lines[0].split(',')[0].lstrip('>(')
            elif tool=='delos':
                genome=lines[0].split(':')[0]
            elif tool=='genapi' or tool=='panaroo':
                genome=lines[0].split('_')[0]

            seq=''
            for l in lines[1:]:
                if l.startswith('>'):

                    if genome in genome2seqs:
                        genome2seqs[genome]=genome2seqs[genome]+seq 
                    else:
                        genome2seqs[genome]=seq
                    seq=''

                    if tool=='prova':
                        genome=l.split(',')[0].lstrip('>(')
                    elif tool=='delos':
                        genome=l.split(':')[0]
                    elif tool=='genapi' or tool=='panaroo':
                        genome=l.split('_')[0]

                else:
                    seq+=l.strip()
                

            if genome in genome2seqs:
                genome2seqs[genome]=genome2seqs[genome]+seq 
            else:
                genome2seqs[genome]=seq
    
    #check all went well
    for g in genome2seqs:
        print(len(genome2seqs[g]))

    return genome2seqs

# Function used to map the genome names of the tools to those of PANPROVA
def translate_genome(tool,species, frag):
    basedir = os.getcwd()
    species = species+'/'+frag+'/'
    # read pandelos file to associate pandelos aliases to genome name
    if tool == 'delos':
        delosdir = basedir+'/'+species+'/PanDelos'

        genome_names = os.listdir(delosdir+'/output/fragmented')
        matchfile='predictedCDSs_filtered_only_genes.bed'
        dict_list_names = []

        for g in genome_names:
            files = os.listdir(delosdir+'/output/fragmented/' + g + '/artifacts')
            afile = [f for f in files if f.endswith(matchfile)][0]
            with open(delosdir+'/output/fragmented/' + g + '/artifacts/' +afile) as f:
                delos_name = f.readline().split('\t')[0]
                
            dict_row = {'genome':g, 'pandelos_name':delos_name}
            dict_list_names.append(dict_row)
        
        names_df=pd.DataFrame(dict_list_names)
        names_df

        delos2prova = dict()
        for r in names_df.index:
            delos2prova[names_df.loc[r,'pandelos_name']]=names_df.loc[r,'genome']

        return delos2prova
    else:
        # Change names of genomes to correposnd to PANPROVA genomes
        prokkadir= basedir+'/'+species+'/Prokka' 
        delosdir = basedir+'/'+species+'/PanDelos'

        roary2genome = dict()
        genome_names = os.listdir(delosdir+'/output/fragmented')
        matchfile='predictedCDSs_filtered_only_genes.bed'
        dict_list_names = [] 
        prokka_dirs=os.listdir(prokkadir)
        for genome_full in prokka_dirs:
            genome=genome_full.split('_')[2]
            files=os.listdir(prokkadir+'/'+genome_full)
            index_tsv = [i for i,name in enumerate(files) if re.search("tsv", name)][0]

            with open(prokkadir+'/'+genome_full+'/'+files[index_tsv]) as f:
                code=f.readlines()[1].split('\t')[0].split('_')[0]
         
            roary2genome[code]=genome

        return roary2genome


# Write to file concatenate alignments of multiple gene families 
def concatenated_to_file(tool,species,frag,genome2seqs,translationdict=None):
    basedir=os.getcwd()
    CG_dir = basedir+'/'+species+'/'+frag+'/Comparison/CG'
    msa = CG_dir+'/msa_'+tool+'/'
    #write to file
    with open(msa+'concat.msa','w') as f:
        for g in genome2seqs:
            if tool=='prova':
                f.write('>'+g+'\n')
            elif tool=='delos':
                f.write('>'+translationdict[g.lstrip('>')].split('_')[1]+'\n')
            elif tool=='genapi':
                f.write('>'+translationdict[g.lstrip('>')]+'\n')

            f.write(str(genome2seqs[g])+'\n')

# Build tree with FastTree
def build_tree(tool,species, frag):
    basedir = os.getcwd()
    CG_dir = basedir+'/'+species+'/'+frag+'/Comparison/CG'
    msa = CG_dir+'/msa_'+tool+'/'
    ret = subprocess.call(['bash',basedir+"/analysis_scripts/run_fasttree.sh",species,msa+'concat.msa',msa+'concat.tree'],stdout=subprocess.DEVNULL)


# Function to build tree and rename it on Roary core genes
def run_roary(species,frag):
    basedir=os.getcwd()
    species=species+'/'+frag
    alndir = basedir+'/'+species+'/Roary/output/' # core gene alignment is already provided by Roary
    ret=subprocess.call(['bash',basedir+"/analysis_scripts/run_fasttree.sh",species,alndir+'core_gene_alignment.aln',alndir+'core_gene_alignment.tree'],stdout=subprocess.DEVNULL)

    with open(alndir+'core_gene_alignment.tree') as infile, open(alndir+'core_gene_alignment_renamed.tree', 'w') as outfile:
        for line in infile:
            line = line.replace('prokka_genome_','')
            line = line.replace('_fr.fasta','')

            outfile.write(line)


def main():
    basedir=os.getcwd()
    for species in ['synth_myco','synth_ecoli','synth_paeru']:
        # PANPROVA 
        gf_cg = read_panprova_CG(species) 
        gene2seq = read_panprova_genes(species)
        write_fasta_to_file('prova',gf_cg, gene2seq,species,1)
        run_mafft('prova',species, 1)
        genome2seqs=concatenate_alignments('prova',species, 1)
        concatenated_to_file('prova',species,frag,genome2seqs)
        build_tree('prova',species, 1)

        for frag in ['0.5','0.8','1']:
            # PanDelos-frags
            gf_cg = read_pandelos_CG(species, frag)
            gene2seq = read_pandelos_genes(species,frag)
            write_fasta_to_file('delos',gf_cg, gene2seq,species,frag)
            run_mafft('delos',species, frag)
            concatenate_alignments('delos',species, frag)
            genome2prova = translate_genome('delos',species,frag)
            concatenated_to_file('delos',species,frag,genome2seqs,genome2prova)
            build_tree('delos',species, frag)

            # Roary
            run_roary(species, frag) # Roary already provides aligned core gene families, so only tree building is required

            # GenAPI
            gf_cg = read_genapi_CG(species, frag)
            gene2seq = read_prokka_genes(species,frag)
            write_fasta_to_file('genapi',gf_cg, gene2seq,species,frag)
            run_mafft('genapi',species, frag)
            concatenate_alignments('genapi',species, frag)
            genome2prova = translate_genome('genapi',species,frag)
            concatenated_to_file('genapi',species,frag,genome2seqs,genome2prova)
            build_tree('genapi',species, frag)

            # Panaroo
            gf_cg = read_panaroo_CG(species, frag)
            gene2seq = read_prokka_genes(species,frag)
            write_fasta_to_file('panaroo',gf_cg, gene2seq,species,frag)
            run_mafft('panaroo',species, frag)
            concatenate_alignments('panaroo',species, frag)
            genome2prova = translate_genome('panaroo',species,frag)
            concatenated_to_file('panaroo',species,frag,genome2seqs,genome2prova)
            build_tree('panaroo',species, frag)

    subprocess.call(['bash',basedir+"/analysis_scripts/run_ete3.sh"],stdout=subprocess.DEVNULL)
    
    # Show results
    for species in ['synth_myco','synth_ecoli','synth_paeru']:
        for frag in ['0.5','0.8','1']:
            print(species+'/'+frag)
            subprocess.call(['cat',basedir+'/'+species+'/'+frag+'/Comparison/ete3.log'])

if __name__ == "__main__":
    main()