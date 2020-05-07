"""
Created on Fri May 1 2020

@author: landengozashti
"""

"""
Associate high parsimony mutations with submitting and originating labs



Inputs: parsimony BED file and corresponding VCF



Outputs: 

flagged_snps_by_lab.out --- Flagged snps by submitting and originating lab

flagged_snps_summary.tsv --- flagged mutations and respective reasoning

""" 

import scipy.stats as stats
import os
import numpy as np
from statistics import mean 
import argparse
import time

dir = 'tests'

"""

Arguments

metafilter = False
if True, program will check metadata for possible errors and return file containing metadata flags

output = ''
User specified output directory

Nextstrain Metadata = ''
Path to metadata file

Nextstrain Parsimony Bed File = ''
Path to nextrain parsimony bed file

Nextstrain VCF = ''
Path to nextrain VCF file



"""


parser = argparse.ArgumentParser(description='Flag suspicious high parsimony mutations (Talk to me nice)')
parser.add_argument('-m', nargs='?', required=True,
                    help='Path to metadata file')
parser.add_argument('-p', nargs='?', required=True,
                    help='Path to nextrain parsimony file')
parser.add_argument('-v', nargs='?', required=True,
                    help='Path to nextrain VCF file')
parser.add_argument('-o', nargs='?', required=True,
                    help='Path to output directory')
parser.add_argument('-b', nargs='?',const='T', default=None,
                    help='Program will also flag borderline suspicious variants that exhibit low minor allele frequency and are significantly assosiated with 1 or more particular lab')
parser.add_argument('-min', nargs='?',
                    help='minimum parsimony (must be an integer); default = 6')
parser.add_argument('-include_others', nargs='?',const='T', default=None,
                    help='Include all mutations with parsimony > min regardless of flags') 

args = vars(parser.parse_args())
if args['min'] == None:
    args['min'] = 6



parsimonyDic = {}
sourceList = []
sourceDic = {}
totalParsimoniousmutationCount = 0
parsimonySourceDic = {}

print('Reading parsimony file')
with open('{0}'.format(args['p']), 'r') as f:
    for line in f:
        sp = line.split('\t')
        if int(sp[-1]) >= int(args['min']):
            parsimonyDic[sp[2]] = sp[-1]

#add VCF parser, perhaps use dataframe rather than dic of lists...
start_time = time.time()
print('Reading VCF')
with open('{0}'.format(args['v']), 'r') as vcf:
    for line in vcf:
        if '##' in line:
            pass
        elif '#' in line:
            sourceList = line.split('\t')
        elif line.split('\t')[1] in parsimonyDic:
            parsimonySourceDic['{0}'.format(line.split('\t')[2])] = {}
            for i in range(9,len(line.split('\t'))):
                parsimonySourceDic['{0}'.format(line.split('\t')[2])][sourceList[i].split('|')[0]] = line.split('\t')[i].split(':')[0]
print('VCF took {0} seconds, thats annoying'.format(int(round((time.time() - start_time)/float(60),0))))
oriParsCountDic = {}
subParsCountDic = {}
subAccessionDic = {}
oriAccessionDic = {}
start_time2 = time.time()
print('Reading metadata')
with open('{0}'.format(args['m']), 'r') as meta:
    for line in meta:
        sp = line.split('\t')
        if 'strain' == sp[0]:
            indexPosOri = sp.index('originating_lab')
            indexPosSub = sp.index('submitting_lab')

        else:
            try:
                if sp[indexPosOri] not in oriAccessionDic:
                    oriAccessionDic[sp[indexPosOri]] = [sp[2]]
                else:
                    oriAccessionDic[sp[indexPosOri]].append(sp[2])

                if sp[indexPosSub] not in subAccessionDic:
                    subAccessionDic[sp[indexPosSub]] = [sp[2]]

                else:
                    subAccessionDic[sp[indexPosSub]].append(sp[2])
            except IndexError:
                pass
print('metadata took {0} seconds'.format((time.time() - start_time2)))


print('Writing output files')
with open('{0}/flagged_snps_by_lab.tsv'.format(args['o']), 'w') as f:
    f.write('bin\tsource\tlab\tsnp\tref\talt\tglobal ref\tglobal alt\tfisher exact\tproportion of calls\tMAF\n')

altDic = {}
refDic = {}
globalRefCountDic = {}
globalAltCountDic = {}
for ori in subAccessionDic:
    refDic[ori] = {}
    altDic[ori] = {}
    for snp in parsimonySourceDic:
        refCount = 0
        altCount = 0
        for accession in subAccessionDic[ori]:
            try:
                if int(parsimonySourceDic[snp][accession]) == 0:
                    refCount += 1
                    if snp not in globalRefCountDic:
                        globalRefCountDic[snp] = 1
                    else:
                        globalRefCountDic[snp] += 1

                elif int(parsimonySourceDic[snp][accession]) > 0:
                    altCount += 1
                    if snp not in globalAltCountDic:
                        globalAltCountDic[snp] = 1
                    else:
                        globalAltCountDic[snp] += 1
           
            except KeyError:
                pass

        if refCount > 0 or altCount > 0:
            refDic[ori][snp] = refCount
            altDic[ori][snp] = altCount

trashDic = {}
susDic = {}
compDic = {}
otherDic = {}


with open('{0}/flagged_snps_by_lab.tsv'.format(args['o']), 'a') as f:
    for ori in refDic:
        for snp in refDic[ori]:
            oddsratio, pvalue = stats.fisher_exact(np.array([[globalRefCountDic[snp]-refDic[ori][snp],refDic[ori][snp]],[globalAltCountDic[snp]-altDic[ori][snp],altDic[ori][snp]]]))
            MAF = float(globalAltCountDic[snp])/(float(globalAltCountDic[snp])+ float(globalRefCountDic[snp]))
            if (float(altDic[ori][snp])/float(globalAltCountDic[snp]) > 0.5):
                if snp not in compDic or float(altDic[ori][snp])/float(globalAltCountDic[snp]) > compDic[snp]:
                    compDic[snp] = float(altDic[ori][snp])/float(globalAltCountDic[snp])

                    trashDic[snp] = '{0}% of alternate allele calls stem from {1}'.format(round(float(altDic[ori][snp])/float(globalAltCountDic[snp])*100,2), ori)
                f.write('trash\tsubmitting lab\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(ori, snp, refDic[ori][snp], altDic[ori][snp], globalRefCountDic[snp], globalAltCountDic[snp],pvalue,float(altDic[ori][snp])/float(globalAltCountDic[snp]),float(globalAltCountDic[snp])/(float(globalAltCountDic[snp])+ float(globalRefCountDic[snp])) ))

            elif (pvalue < 0.05 and MAF < 0.01 and float(altDic[ori][snp])/float(globalAltCountDic[snp]) > 0.1 and args['b']=='T' and snp not in trashDic):
                f.write('suspicious\tsubmitting lab\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(ori, snp, refDic[ori][snp], altDic[ori][snp], globalRefCountDic[snp], globalAltCountDic[snp],pvalue,float(altDic[ori][snp])/float(globalAltCountDic[snp]),float(globalAltCountDic[snp])/(float(globalAltCountDic[snp])+ float(globalRefCountDic[snp])) ))
                susDic[snp] = 'Minor allele frequency is {0} and {1} contributes a suspiciously large proportion of minor allele calls (p = {2})'.format(MAF, ori, pvalue)
            elif args['include_others']=='T' and snp not in trashDic and snp not in susDic:
                otherDic[snp] = 'parsimony > {0}'.format(int(args['min'])-1)
                f.write('no flag\tsubmitting lab\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(ori, snp, refDic[ori][snp], altDic[ori][snp], globalRefCountDic[snp], globalAltCountDic[snp],pvalue,float(altDic[ori][snp])/float(globalAltCountDic[snp]),float(globalAltCountDic[snp])/(float(globalAltCountDic[snp])+ float(globalRefCountDic[snp])) ))
altDic = {}
refDic = {}
globalRefCountDic = {}
globalAltCountDic = {}

for ori in oriAccessionDic:
    refDic[ori] = {}
    altDic[ori] = {}
    for snp in parsimonySourceDic:
        refCount = 0
        altCount = 0
        for accession in oriAccessionDic[ori]:
            try:
                if int(parsimonySourceDic[snp][accession]) == 0:
                    refCount += 1
                    if snp not in globalRefCountDic:
                        globalRefCountDic[snp] = 1
                    else:
                        globalRefCountDic[snp] += 1

                elif int(parsimonySourceDic[snp][accession]) > 0:
                    altCount += 1
                    if snp not in globalAltCountDic:
                        globalAltCountDic[snp] = 1
                    else:
                        globalAltCountDic[snp] += 1
            except KeyError:
                pass

        if refCount > 0 or altCount > 0:
            refDic[ori][snp] = refCount
            altDic[ori][snp] = altCount


with open('{0}/flagged_snps_by_lab.tsv'.format(args['o']), 'a') as f:
    for ori in refDic:
        for snp in refDic[ori]:
            MAF = float(globalAltCountDic[snp])/(float(globalAltCountDic[snp])+ float(globalRefCountDic[snp]))
            oddsratio, pvalue = stats.fisher_exact(np.array([[globalRefCountDic[snp]-refDic[ori][snp],refDic[ori][snp]],[globalAltCountDic[snp]-altDic[ori][snp],altDic[ori][snp]]]))

            if (float(altDic[ori][snp])/float(globalAltCountDic[snp]) > 0.5):
                if snp not in compDic or float(altDic[ori][snp])/float(globalAltCountDic[snp]) > compDic[snp]:
                    compDic[snp] = float(altDic[ori][snp])/float(globalAltCountDic[snp])

                    trashDic[snp] = '{0}% of alternate allele calls stem from {1}'.format(round(float(altDic[ori][snp])/float(globalAltCountDic[snp])*100,2), ori)
                f.write('trash\toriginating lab\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(ori, snp, refDic[ori][snp], altDic[ori][snp], globalRefCountDic[snp], globalAltCountDic[snp],pvalue,float(altDic[ori][snp])/float(globalAltCountDic[snp]),float(globalAltCountDic[snp])/(float(globalAltCountDic[snp])+ float(globalRefCountDic[snp])) ))

            elif (pvalue < 0.05 and MAF < 0.01 and float(altDic[ori][snp])/float(globalAltCountDic[snp]) > 0.1 and args['b']=='T' and snp not in trashDic):
                f.write('suspicious\toriginating lab\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(ori, snp, refDic[ori][snp], altDic[ori][snp], globalRefCountDic[snp], globalAltCountDic[snp],pvalue,float(altDic[ori][snp])/float(globalAltCountDic[snp]),float(globalAltCountDic[snp])/(float(globalAltCountDic[snp])+ float(globalRefCountDic[snp])) ))
                susDic[snp] = 'Minor allele frequency is {0} and {1} contributes a suspiciously large proportion of minor allele calls (p = {2})'.format(MAF, ori, pvalue)
            elif args['include_others']=='T' and snp not in trashDic and snp not in susDic:
                otherDic[snp] = 'parsimony > {0}'.format(int(args['min'])-1)
                f.write('no flag\toriginating lab\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(ori, snp, refDic[ori][snp], altDic[ori][snp], globalRefCountDic[snp], globalAltCountDic[snp],pvalue,float(altDic[ori][snp])/float(globalAltCountDic[snp]),float(globalAltCountDic[snp])/(float(globalAltCountDic[snp])+ float(globalRefCountDic[snp])) ))

with open('{0}/flagged_snps_summary.tsv'.format(args['o']), 'w') as f:
    f.write('snp\tbin\treasoning\n')

    for snp in trashDic:
        f.write('{0}\ttrash\t{1}\n'.format(snp, trashDic[snp]))
    if args['b']=='T':
        for snp in susDic:
            f.write('{0}\tsuspicious\t{1}\n'.format(snp, susDic[snp]))
    if args['include_others']=='T':
        for snp in otherDic:
            if snp not in susDic and snp not in trashDic:
                f.write('{0}\tno flag\t{1}\n'.format(snp, otherDic[snp]))


print('That all took like {0} minutes... I should probably write better code'.format(int(round((time.time() - start_time)/float(60),0))))

