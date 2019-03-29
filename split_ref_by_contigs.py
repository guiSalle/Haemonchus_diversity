#!/bin/python3
#Usage: python ./scripts/split_ref_by_contig.py ./reference.fa

import os
from os import path
import sys

infile = open(sys.argv[1])

ref_file = sys.argv[1].split('/')[1]

print('\n' + "reference file: " + ref_file + '\n')
num_contigs = 'grep "^>" %s  | wc -l' % (ref_file)
print("The number of contigs in reference fasta:")
os.system(num_contigs)

opened = False # Assume outfile is not open

for line_ref in infile:
    line_ref = line_ref.split(' ')[0]
    if line_ref[0] == ">": # If line begins with ">"
        if(opened): 
            outfile.close() # Will close the outfile if it is open (see below and follow loop)
        opened = True # Set opened to True to represent an opened outfile
        contig_name = line_ref[1:].rstrip() # Extract contig name
        print("contig: " + contig_name)
        print("Name of outfile: " + str(ref_file.split(".")[0]) + "_" + str(contig_name) + ".fasta") 
        outfile=open(str(ref_file.split(".")[0]) + "_" + str(contig_name) + ".fasta", 'w')
        outfile.write(line_ref) ## Output
    else:
        outfile.write(line_ref) # Output in case line does not begin with '>'
outfile.close()
print "Done"
