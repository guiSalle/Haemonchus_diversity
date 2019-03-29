# Haemonchus_diversity
Code and data for global diversity study. 

These scripts are distributed in the hope that they will be useful,
but WITHOUT ANY WARRANTY.

# Main script
GlobalDiversityScript.R : R script used to analyse the data and produce figures and supplementary figures.

# Helper scripts
These are a few tools used to create intermediate files needed in GlobalDiversityScript.R.

* ``` convert_hits_to_gene.sh ```:
A shell script needed to convert v3 assembly coordinates into v0 assembly gene identifier. This script uses a python script (```split_ref_by_contigs.py```) and another shell for interecting bed file with gtf (```retrieve_gene_gtf_from_blast.sh```).

* ``` create_GradientForest_input_ANGSD.py ``` :
A python script to combine minor allele frequency estimates from ANGSD with environmental data into an input file needed for gradientForest analysis. 

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
