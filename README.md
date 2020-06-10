# COVID-19-Lab-Specific-Bias-Filter
We developed a program in python 3 to systematically flag COVID-19 genomes for variants resulting from possible lab specific biases. The program requires a concurrent VCF file, Nextstrain parsimony file and GISAID metadata file as input. It outputs two tsv files: one with each flagged snp and the respective reasoning underlying the flag, and another providing more detailed information on each flagged lab.

We also added a GISAID metadata filter, which filters errors in "submitting lab" and "originating lab" names and generates a merged metadata file.

##########################################

Lab Specific Bias Filter Usage:

Make sure that parsimony_per_site_analysis.py and sequenceAnalyzer.py are in the same directory. We provide the files we used when running the program.

python3 parsimony_per_site_analysis2.py [options] -m [Path to GISAID metadata file] -p [Path to Nextrain parsimony file] -v [Path to VCF file] -o [Path to output directory]

Input:

GISAID metadata file, Nextrain parsimony file and VCF file

Options:

-b                   Program will also flag borderline suspicious variants that exhibit low minor allele frequency and are                          significantly assosiated with 1 or more particular lab.

-include             Snps with parsimony > min will be included in the output regardless of flags

-track               Program will output a track annotating highly suspect lab assosiated mutations and a track annotating                          Artic primers that overlap or are within 10bp of lab-associated mutations in BED Detail format. Both                          tracks can be directly uploaded to the UCSC Genome Browser. The user must also download "primers.txt"                          (available in this repository) and store it in the current working directory.

-min_parsimony       Minimum parsimony (must be an integer) default = 4 (Beware that reducing this parameter from the default
                     can obscure accuracy)

-min_contribution    Minimum percent of alternate alleles contributed by a
                     specific lab for variants to be considered highly
                     suspect; default = 80 (beware that reducing this
                     parameter from the default can obscure accuracy)

Output:

flagged_snps_by_lab.tsv    A tab delimitted file displaying bin, source, lab, snp, ref, alt, global ref, global alt, fisher                              exact, proportion of calls, MAF for each flagged snp where:

bin = "highly suspect", "suspect," or "no flag" 

source = "originating lab" or "submitting lab"

lab = lab name

snp = snp

ref = reference allele count

alt = alternate allele count

global ref = global reference allele count

global alt = global alternate allele count

fisher exact = p value for fishers exact test assosiating lab specific with global allele frequency

proportion of calls = proportion of minor allele calls attributed to respective lab

MAF = minor allele frequency

flagged_snps_summary.tsv   A tab delimitted file displaying each flagged snp, bin and the reasoning underlying the flag

lab_associated_error_final.bed     A BED Detail file annotating lab-associated mutations compatable with the UCSC Genome                                        Browser

Artic_primers_final.bed     A BED Detail file annotating Artic Primers compatable with the UCSC Genome Browser

##########################################

GISAID metadata filter Usage:

We provide our metadata_1_merged.tsv output along with the our log file when we used metadata_1.tsv as input. We hand checked the log file by hand to ensure accuracy and commented out merger errors.

python3 meta_data_filter.py [options] -m [Path to GISAID metadata file] -o [Path to output directory]

Input:

GISAID metadata file

Options:

-c     If "TRUE" only merge high confidence errors.

Output:

metadata_filter.log     Log of changes made to the metadata set and warnings if any changes are made with lower confidence

metadata_merged.tsv     Merged metadata file


