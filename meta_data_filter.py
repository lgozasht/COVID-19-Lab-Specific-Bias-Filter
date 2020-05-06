"""

Created on Fri May 1 2020

@author: landengozashti

"""


oriDic = {}
subDic = {}
title  = ""
from difflib import SequenceMatcher
import argparse

"""

Input:  GISAID metadata filter

Output: Merged metadata file, log file

"""


parser = argparse.ArgumentParser(description='Filter GISAID metadata for errors in originating and submitting lab names. Merge suspected errors.')
parser.add_argument('-m', nargs='?', required=True,
                    help='Path to metadata file')
parser.add_argument('-o', nargs='?', required=True,
                    help='Path to output directory')
parser.add_argument('-c', nargs='?',
                    help='If TRUE only merge high confidence errors')


args = vars(parser.parse_args())

with open(args['m'], 'r') as f:
    for line in f:
        
        sp = line.split('\t')
        if 'strain' == sp[0]:
            indexPosOri = sp.index('originating_lab')
            indexPosSub = sp.index('submitting_lab')
            title  = line
        else:
            try:

                if sp[indexPosOri] not in oriDic.keys():
                    oriDic[sp[indexPosOri]] = [line]
                else: 
                    oriDic[sp[indexPosOri]].append(line)
            except IndexError:
                pass

oriflagScore = {}
subflagScore = {}

def compare(a, b):
    count = 0
    for x, y in zip(a, b):
        if x == y:
            count += 1
    return count

def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

#with open('{0}/name_flags.tsv'.format(args['o']), 'w') as  f:
#    pass

oriflagDic = {}
with open('{0}/metadata_filter.log'.format(args['o']), 'w') as logFile:
    pass

def filter(oriDic):
    with open('{0}/metadata_filter.log'.format(args['o']), 'a') as logFile:
        for ori1 in oriDic:
            for ori2 in oriDic:
                if ori1 != ori2:
                    if float(len(ori2))/float(len(ori1)) > 0.7:
                        if float(compare(ori1, ori2))/float(len(ori1)) > 0.8:
                            if ori1 not in oriflagDic.values():  
                                oriflagDic[ori1] = ori2
                                oriflagScore[ori1] = float(compare(ori1, ori2 ))/float(len(ori1))
                                if float(compare(ori1, ori2 ))/float(len(ori1)) > 0.95:
                                    logFile.write('Merging {0} and {1}\n\n'.format(ori1, ori2))
                                    newList = []
                                    for line in oriDic[ori2]:
                                        line = line.replace(ori2, ori1)
                                        newList.append(line)

                                    oriDic[ori2] = newList

                                
                                elif args['c'] != 'TRUE':
                                    logFile.write('Merging {0} and {1}, but with less confidence... you should check this one\n\n'.format(ori1, ori2))
                                    newList = []
                                    for line in oriDic[ori2]:
                                        line = line.replace(ori2, ori1)
                                        newList.append(line)

                                    oriDic[ori2] = newList


                        elif similar(ori1, ori2) > 0.9:
                            if ori1 not in oriflagDic.values():
                                oriflagScore[ori1] = similar(ori1, ori2)
                                oriflagDic[ori1] = ori2

                                if similar(ori1, ori2) > 0.95:
                                    newList = []
                                    for line in oriDic[ori2]:
                                        line = line.replace(ori2, ori1)
                                        newList.append(line)

                                    oriDic[ori2] = newList

                                    logFile.write('Merging {0} and {1}\n\n'.format(ori1, ori2))
                                elif args['c'] != 'TRUE':
                                    logFile.write('Merging {0} and {1}, but with less confidence... you should check this one\n\n'.format(ori1, ori2))
                                
                                    newList = []
                                    for line in oriDic[ori2]:
                                        line = line.replace(ori2, ori1)
                                        newList.append(line)

                                    oriDic[ori2] = newList





    return oriflagDic, oriflagScore, oriDic

oriflagDic, oriflagScore, oriDic = filter(oriDic)



#with open('{0}/name_flags.tsv'.format(args['o']), 'a') as  f:
#    for ori in oriflagDic:
#        f.write('Originating Lab\t{0}\t{1}\t{2}\n'.format(ori, oriflagDic[ori], oriflagScore[ori]))

for ori in oriDic:
    for line in oriDic[ori]:
        sp = line.split('\t')
        try:
            if sp[indexPosSub] not in subDic.keys(): 
                subDic[sp[indexPosSub]] = [line]
            else:
                subDic[sp[indexPosSub]].append(line)
        except IndexError:
            pass

subflagDic, subflagScore, subDic = filter(subDic)



#with open('{0}/name_flags.tsv'.format(args['o']), 'a') as  f:
#    for sub in subflagDic:
#        f.write('Submitting Lab\t{0}\t{1}\t{2}\n'.format(sub, subflagDic[sub], subflagScore[sub]))

with open('{0}/{1}_merged.tsv'.format(args['o'], args['m'].split('/')[-1].split('.')[0]), 'w') as f:
    f.write(title)
    for sub in subDic:
        for line in subDic[sub]:
            f.write(line)


