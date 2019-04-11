#!/usr/binenv python 
#=================================================
#= Compute parameter uncertainty
#=================================================
import numpy 
import dadi
import Models_2D
import Models_ms
import os
import sys

wk = os.getcwd() + '/'

##-- arguments
popc = str(sys.argv[2])
bestmod = str(sys.argv[1])

print(popc,bestmod)

##--- Pop Couple 
#popc = 'FRA.1-NAM'

p1 = popc.split("-")[0]
p2 = popc.split("-")[1]

print(p1,p2)

## SFS file
fsfile = wk + 'dadi_sfs/'+ popc +'/2d_wg_'+ p1 + '.' + p2 + '.fs'
data = dadi.Spectrum.from_file(fsfile)
ns = data.sample_sizes

print(ns)

## Dadi output
resfile = wk + 'dadi_sfs/' + popc + '/Results_Summary_Short.txt'
print("%s\t%s" % ("resfile found: ",resfile))

##--- Retrieve nSites
gensize = 237372181
with open(wk + 'Summary_nSites.txt','r') as s:
    for line in s:
        line = line.strip().split("\t")
        if line[0] == p1 and line[1] == p2:
            fracSites = int(line[3])
            totSites = int(line[2])
            chrL = float(gensize)*float(fracSites)/float(totSites) ## L = fraction of Sites considered * gensize

print(chrL)

##--- Retrieve best model and parameters
popt = []
#bestmod = ''

with open(resfile,'r') as res:
    for l in res:
        for l in res:
##-- for automatic best mod retrieval
#            if 'inf' not in l:
#                while bestmod =='':
#                    l = l.strip().split('\t')
#                    bestmod = l[0]
#                    poptt = l[len(l)-1].split(',')
#                    for i in range(len(poptt)):
#                        popt.append(float(poptt[i]))
##-- from file with best model
            l = l.strip().split('\t')
            if l[0] == bestmod:
                poptt = l[len(l)-1].split(',')
                scaledtheta = float(l[5])
                for i in range(len(poptt)):
                    popt.append(float(poptt[i]))
print(bestmod,scaledtheta,popt)

# Data generation using ms
demofun = getattr(Models_ms, bestmod+'_ms')
func_ex = dadi.Numerics.make_extrap_log_func(getattr(Models_2D, bestmod))

##-- ms command
mscore = demofun(popt) 
print(mscore)

# call ms
datasize = ns
reps = 1
rec = 1.68

# ms command to generate replicates
mscommand = dadi.Misc.ms_command(scaledtheta, datasize, mscore, reps, recomb = rec, rsites = chrL)

# create output directory
bsdir = wk + 'boot_dadi/' + popc + '/boostraps/'
if not os.path.exists(bsdir):
    os.makedirs(bsdir)

##-- Create 100 fs bs files
print("Start bootstrapping")
for ii in range(100):
    msfs = dadi.Spectrum.from_ms_file(os.popen(mscommand))
    ## fold simulated data because original data is also folded
    msfold = msfs.fold()
    msfold.to_file(bsdir + '{0:02d}.fs'.format(ii)) 
print "Bootstraps done"

# Estimate parameter uncertainties using the Godambe Information Matrix, to
# account for linkage in the data. From bootstrapped data.

all_boot = [dadi.Spectrum.from_file(bsdir + '{0:02d}.fs'.format(ii)) 
            for ii in range(100)]

# uncert contains the estimated standard deviations of each parameter, with
# theta as the final entry in the list.
pts_l = [30,40,50]
uncerts = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot, popt, data,
                                  multinom = True)

print('Estimated parameter standard deviations from GIM: {0}'.format(uncerts))
