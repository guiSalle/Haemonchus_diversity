#!/usr/bin/env python
###--- 24.05.2018
###--- Takes x mafs.gz files from x POP + env variables
###--- Output SNPs shared across POP
###--- Output a matrix: each row = 1 pop / columns = envVar + SNP


import sys
import os
import gzip
#from collections import defaultdict

def store_env(envfile,poplist):
	ENV = {}
	SIZE = {}
	with open(envfile,'r') as qf:
		for line in qf:
			for line in qf:
				line = line.strip().split("\t")
				iso = line[5]
				if iso in SIZE.keys():
					SIZE[iso] = SIZE[iso] + 1
				else:
					SIZE[iso] = 1
				if iso in poplist and not(iso in ENV.keys()):
					ENV[iso] = line[9:len(line)-1]

	return ENV,SIZE

###-- Store SNP of each pop
def read_pop_freq(poplist,SIZE,chrom):
#	SNP = {} ##- Dict to store pop with given SNP
	SNPOP_dicts = {}
	for pop in poplist:
		SNPOP_dicts[pop] = {}
		nm = 0
		print("%s\t%s" % ("Processing Pop:",pop))
		frqfile = pop + '.' + str(chrom) + '.mafs.gz'
		with gzip.open(frqfile,'r') as qf:
			for line in qf:
				for line in qf:
					line = line.strip().split("\t")
					if int(line[6]) >= int(0.9*int(SIZE[pop])): ##-- 90% call rate filter
						maf = line[4]
						mrk = line[0] + '.' + line[1]
						SNPOP_dicts[pop][mrk] = maf
						nm +=1
		print ("%d\t%s" % (nm,"SNPs stored"))
						
	return SNPOP_dicts

###-- trim file
def trim(SNPOP,poplist):
	trimmedSNP = []
	keyset = []
	n = 0
	for k in SNPOP.keys():
		d = SNPOP[k]
		keyset.append(set(d.keys()))
	u = set.intersection(*keyset) # retained 
	for k in u:
		trimmedSNP.append(k)
		n+=1
	print("%d\t%s" % (n,"SNPs retained across populations"))

	return trimmedSNP

###-- Output file
def outGDF(ENV,trimmedSNP,dicoSNPbypop,poplist,chrom):
	out = os.getcwd() + "/GdF."+str(chrom)+".in"
	with open(out, 'w') as o:
		for p in poplist:
			line = ''
			for i in range(len(ENV[p])):
				line = line + '\t' + ENV[p][i]
			for i in range(len(trimmedSNP)):
				snp = trimmedSNP[i]
				line = line + '\t' + str(dicoSNPbypop[p][snp])
			o.write("%s\t%s\n" % (p,line))
	print("%s" % "---=== Output file ready for use ===---")

#############---- RUN
def main():
	import argparse
	## Get arguments
	parser = argparse.ArgumentParser(description='Given a list of pops, output GdF input file')
	parser.add_argument('envfile', type=str, help='File with env variables')
	parser.add_argument('nchrom', type=str, help='Number of chromosomes')
	parser.add_argument('querypop', metavar='queryFile', type=str, help='Pop list separated by comma')
	args = parser.parse_args() #gets the arguments

	## Poplist
	poplist = args.querypop.split(',')
	
	## Store env variables
	print("%s" % ("---=== Store env variables ===---"))
	env = store_env(args.envfile,poplist)[0]
	size = store_env(args.envfile,poplist)[1]

	## Loop over every chromosome
	for chromosome in range(int(args.nchrom)+1)[1:]:
	## Mark common SNPs + store MAF by pop
		print("%s%d%s" % ("---=== Preprocessing SNPs CHROM",chromosome ,"===---"))
		dicoSNPbypop = read_pop_freq(poplist,size,chromosome)

	## Retain common SNPs 
		print("%s" % ("---=== Trimming SNPs ===---"))
		trimmedSNP = trim(dicoSNPbypop,poplist)
	
	## Output GdF file
		print("%s" % ("---=== Output GradientForest input file ===---"))
		outGDF(env,trimmedSNP,dicoSNPbypop,poplist,chromosome)

if __name__ == '__main__':
	main()			      






