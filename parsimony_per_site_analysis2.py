"""
Created on Fri May 1 2020

@author: landengozashti
"""

"""
Associate high parsimony mutations with submitting and originating labs



Inputs: parsimony BED file, corresponding VCF, GISAID Metadata file



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
import pandas as pd

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


parser = argparse.ArgumentParser(description='Flag suspicious variants in COVID-19 genomic data (Talk to Russ and I nice)')
parser.add_argument('-m', nargs='?', required=True,
                    help='Path to metadata file')
parser.add_argument('-p', nargs='?', required=True,
                    help='Path to nextrain parsimony file')
parser.add_argument('-v', nargs='?', required=True,
                    help='Path to nextrain VCF file')
parser.add_argument('-o', nargs='?', required=True,
                    help='Path to output directory')
parser.add_argument('-b', nargs='?',const='T', default=None,
                    help='Program will also flag suspicious variants that exhibit low minor allele frequency and are significantly assosiated with 1 or more particular lab')
parser.add_argument('-min_parsimony', nargs='?',
                    help='Minimum parsimony; must be an integer; default = 4 (Beware that reducing this parameter from the default can obscure accuracy)')
parser.add_argument('-include', nargs='?',const='T', default=None,
                    help='Include all mutations with parsimony > min_parsimony regardless of flags')
parser.add_argument('-min_contribution', nargs='?',
                    help='Minimum percent of alternate alleles contributed by a specific lab for variants to be considered highly suspect; default = 80 (beware that reducing this parameter from the default can obscure accuracy)') 
parser.add_argument('-track', nargs='?',const='T', default=None,
                    help='Program will output a track annotating highly suspect lab assosiated mutations and a track annotating Artic primers that overlap or are within 10bp of lab-associated mutations in BED Detail format. Both can be directly uploaded to the UCSC Genome Browser.')


args = vars(parser.parse_args())
if args['min_parsimony'] == None:
    args['min_parsimony'] = 4
if args['min_contribution'] == None:
    args['min_contribution'] = .8
else:
    args['min_contribution'] = float(args['min_contribution'])/100.0

refChrom = 'NC_045512v2'
parsimonyDic = {}
sourceList = []
sourceDic = {}
totalParsimoniousmutationCount = 0
parsimonySourceDic = {}

print('Reading parsimony file')
with open('{0}'.format(args['p']), 'r') as f:
    for line in f:
        sp = line.split('\t')
        if int(sp[-1]) >= int(args['min_parsimony']):
            parsimonyDic[sp[2]] = sp[-1]

#add VCF parser, perhaps use dataframe rather than dic of lists...
start_time = time.time()
print('Reading VCF')

"""
vcf = pd.DataFrame.read_csv(args['v'], sep='\t', header = 2)
index = 0
for snp in vcf['POS']:
    index +=1 
    if snp in parsimonyDic:
        snpList.append(vcf.loc[[index]])
for i in range(9,len(line.split('\t'))):
                parsimonySourceDic['{0}'.format(line.split('\t')[2])][sourceList[i].split('|')[0]] = line.split('\t')[i].split(':')[0]
"""

with open(args['v'], 'r') as vcf:
    for line in vcf:
        if '##' in line:
            pass
        elif '#' in line:
            sourceList = line.split('\t')
        else:
            cols = line.split('\t')
            (pos, mut) = (cols[1], cols[2])
            if pos in parsimonyDic:
                parsimonySourceDic[mut] = {}
                for i in range(9,len(cols)):
                    epiId = sourceList[i].split('|')[0]
                    parsimonySourceDic[mut][epiId] = cols[i].split(':')[0]


print('VCF took {0} seconds'.format(time.time() - start_time))
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
mafDic = {}

with open('{0}/flagged_snps_by_lab.tsv'.format(args['o']), 'a') as f:
    for ori in refDic:
        for snp in refDic[ori]:
            oddsratio, pvalue = stats.fisher_exact(np.array([[globalRefCountDic[snp]-refDic[ori][snp],refDic[ori][snp]],[globalAltCountDic[snp]-altDic[ori][snp],altDic[ori][snp]]]))
            MAF = float(globalAltCountDic[snp])/(float(globalAltCountDic[snp])+ float(globalRefCountDic[snp]))
            mafDic[snp] = MAF

            if (float(altDic[ori][snp])/float(globalAltCountDic[snp]) > args['min_contribution']):
                if snp not in compDic or float(altDic[ori][snp])/float(globalAltCountDic[snp]) > compDic[snp]:
                    compDic[snp] = float(altDic[ori][snp])/float(globalAltCountDic[snp])

                    trashDic[snp] = '{0}% of alternate allele calls stem from {1}'.format(round(float(altDic[ori][snp])/float(globalAltCountDic[snp])*100,2), ori)
                f.write('highly suspect\tsubmitting lab\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(ori, snp, refDic[ori][snp], altDic[ori][snp], globalRefCountDic[snp], globalAltCountDic[snp],pvalue,float(altDic[ori][snp])/float(globalAltCountDic[snp]),float(globalAltCountDic[snp])/(float(globalAltCountDic[snp])+ float(globalRefCountDic[snp])) ))

            elif (pvalue < 0.05 and MAF < 0.01 and float(altDic[ori][snp])/float(globalAltCountDic[snp]) > 0.1 and args['b']=='T' and snp not in trashDic):
                f.write('suspect\tsubmitting lab\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(ori, snp, refDic[ori][snp], altDic[ori][snp], globalRefCountDic[snp], globalAltCountDic[snp],pvalue,float(altDic[ori][snp])/float(globalAltCountDic[snp]),float(globalAltCountDic[snp])/(float(globalAltCountDic[snp])+ float(globalRefCountDic[snp])) ))
                susDic[snp] = 'Minor allele frequency is {0} and {1} contributes {3}% of minor allele calls (p = {2})'.format(MAF, ori, pvalue,round(float(altDic[ori][snp])/float(globalAltCountDic[snp])*100,2))
            elif args['include']=='T' and snp not in trashDic and snp not in susDic:
                otherDic[snp] = 'parsimony > {0}'.format(int(args['min_parsimony'])-1)
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

            if (float(altDic[ori][snp])/float(globalAltCountDic[snp]) > args['min_contribution']):
                if snp not in compDic or float(altDic[ori][snp])/float(globalAltCountDic[snp]) > compDic[snp]:
                    compDic[snp] = float(altDic[ori][snp])/float(globalAltCountDic[snp])

                    trashDic[snp] = '{0}% of alternate allele calls stem from {1}'.format(round(float(altDic[ori][snp])/float(globalAltCountDic[snp])*100,2), ori)
                f.write('highly suspect\toriginating lab\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(ori, snp, refDic[ori][snp], altDic[ori][snp], globalRefCountDic[snp], globalAltCountDic[snp],pvalue,float(altDic[ori][snp])/float(globalAltCountDic[snp]),float(globalAltCountDic[snp])/(float(globalAltCountDic[snp])+ float(globalRefCountDic[snp])) ))

            elif (pvalue < 0.05 and MAF < 0.01 and float(altDic[ori][snp])/float(globalAltCountDic[snp]) > 0.1 and args['b']=='T' and snp not in trashDic):
                f.write('suspect\toriginating lab\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(ori, snp, refDic[ori][snp], altDic[ori][snp], globalRefCountDic[snp], globalAltCountDic[snp],pvalue,float(altDic[ori][snp])/float(globalAltCountDic[snp]),float(globalAltCountDic[snp])/(float(globalAltCountDic[snp])+ float(globalRefCountDic[snp])) ))
                susDic[snp] = 'Minor allele frequency is {0} and {1} contributes {3}% of minor allele calls (p = {2})'.format(MAF, ori, pvalue,round(float(altDic[ori][snp])/float(globalAltCountDic[snp])*100,2))
            elif args['include']=='T' and snp not in trashDic and snp not in susDic:
                otherDic[snp] = 'parsimony > {0}'.format(int(args['min_parsimony'])-1)
                f.write('no flag\toriginating lab\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(ori, snp, refDic[ori][snp], altDic[ori][snp], globalRefCountDic[snp], globalAltCountDic[snp],pvalue,float(altDic[ori][snp])/float(globalAltCountDic[snp]),float(globalAltCountDic[snp])/(float(globalAltCountDic[snp])+ float(globalRefCountDic[snp])) ))

with open('{0}/flagged_snps_summary.tsv'.format(args['o']), 'w') as f:
    f.write('snp\tbin\treasoning\n')

    for snp in trashDic:
        f.write('{0}\thighly suspect\t{1}\n'.format(snp, trashDic[snp]))
    if args['b']=='T':
        for snp in susDic:
            f.write('{0}\tsuspect\t{1}\n'.format(snp, susDic[snp]))
    if args['include']=='T':
        for snp in otherDic:
            if snp not in susDic and snp not in trashDic:
                f.write('{0}\tno flag\t{1}\n'.format(snp, otherDic[snp]))

if args['track'] != None:
    primerTrackList = []
    primerDic = {}
    with open('primers.txt', 'r') as primeFile:
        for line in primeFile:

            for snp in trashDic:

                sp = line.split('\t')
                if int(snp[1:-1]) > (int(sp[1]) - 10) and int(snp[1:-1]) < (int(sp[2]) + 10):
                    primerTrackList.append(line)
                    if snp not in primerDic:
                        primerDic[snp] = sp[3]
                    else:
                        primerDic[snp] += (',' + sp[3])

    with open('{0}/lab_associated_error.bed'.format(args['o']), 'w') as f:
        for snp in trashDic:
            if snp in primerDic:
                f.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{8}\t{7}\t{6}\n'.format(refChrom, str(int(snp[1:-1])-1),snp[1:-1],snp.replace('T','U'), parsimonyDic[snp[1:-1]].strip('\n'),globalAltCountDic[snp],trashDic[snp],primerDic[snp],mafDic[snp]))
            else:
                f.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{8}\t{7}\t{6}\n'.format(refChrom, str(int(snp[1:-1])-1),snp[1:-1],snp.replace('T','U'), parsimonyDic[snp[1:-1]].strip('\n'),globalAltCountDic[snp],trashDic[snp],'NA',mafDic[snp]))
    with open('{0}/Artic_primers.bed'.format(args['o']), 'w') as f:
        for line in primerTrackList:
            f.write('{0}\t{1}\t{2}\t{3}\n'.format(refChrom, line.split('\t')[1], line.split('\t')[2], line.split('\t')[3]))
    os.system('sort -k1,1 -k2,2n {0}/Artic_primers.bed > {0}/Artic_primers_sorted.bed'.format(args['o']))
    os.system('sort -k1,1 -k2,2n {0}/lab_associated_error.bed > {0}/lab_associated_error_sorted.bed'.format(args['o']))

    os.system('rm -r {0}/lab_associated_error.bed'.format(args['o']))
    os.system('rm -r {0}/Artic_primers.bed'.format(args['o']))


    with open('{0}/lab_associated_error_sorted.bed'.format(args['o']),'r') as sorted:
        lineList = []
        for line in sorted:
            lineList.append(line)

    with open('{0}/lab_associated_error_final.bedDetail'.format(args['o']),'w') as final:
        final.write('browser position NC_045512v2:0-29903\ntrack name=lab type=bedDetail description="Lab-associated Mutations" visibility=pack\n#chrom\tchromStart\tchromStop\tname\tparsimony score\tnumber of alt alleles\tminor allele frequency\tPrimer within 10bp\tcomment\n')
        for line in lineList:
            cols = line.split("\t")
            bedDetailCols = cols[0:4] + [cols[3]] + [cols[8]]
            final.write('\t'.join(bedDetailCols))

    with open('{0}/lab_associated_error_final.bedForBigBed'.format(args['o']),'w') as final:
        for line in lineList:
            final.write(line)

    os.system('rm -r {0}/lab_associated_error_sorted.bed'.format(args['o']))

    with open('{0}/Artic_primers_sorted.bed'.format(args['o']),'r') as sorted:
        lineList = []
        for line in sorted:
            lineList.append(line)
    with open('{0}/Artic_primers_final.bed'.format(args['o']),'w') as final:
        final.write('browser position NC_045512v2:0-29903\ntrack name=primers type=bed description="Artic Primers" visibility=pack\n#chrom\tchromStart\tchromStop\tname\n')
        for line in lineList:
            final.write(line)
    os.system('rm -r {0}/Artic_primers_sorted.bed'.format(args['o']))




print('That all took like {0} minutes... sorry for the wait'.format(int(round((time.time() - start_time)/float(60),0))))

