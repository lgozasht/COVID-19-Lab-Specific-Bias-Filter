# COVID-19-Lab-Specific-Bias-Filter
See also https://github.com/lgozasht/SARS-COV-2_COLLECTIVE_ANALYSIS


We developed a program in python 3 to systematically flag COVID-19 genomes for variants resulting from possible lab specific biases. The program requires a concurrent VCF file, Nextstrain parsimony file and GISAID metadata file as input. It outputs two tsv files: one with each flagged snp and the respective reasoning underlying the flag, and another providing more detailed information on each flagged lab.

We also added a GISAID metadata filter, which filters errors in "submitting lab" and "originating lab" names and generates a merged metadata file.



## Lab Specific Bias Filter Usage:

Make sure that parsimony_per_site_analysis.py and sequenceAnalyzer.py are in the same directory. We provide the files we used when running the program.

```
python3 parsimony_per_site_analysis3.py [options] -m [Path to GISAID metadata file] -p [Path to Nextrain parsimony file] -v [Path to VCF file] -o [Path to output directory]
```

### Input:

GISAID metadata file, Nextrain parsimony file and VCF file; for -track option, also primers.txt.

### Options:

**-b T**:            Program will also flag borderline suspicious variants that exhibit low minor allele frequency and are significantly associated with 1 or more particular lab.

**-include T**:      Snps with parsimony > min will be included in the output regardless of flags

**-track**:          Program will output a track annotating highly suspect lab associated mutations and a track annotating                          Artic primers that overlap or are within 10bp of lab-associated mutations.
                     See [Viewing tracks in the UCSC Genome Browser](#viewing-tracks-in-the-UCSC-Genome-Browser) below for instructions.
                     The file "primers.txt" (available in this repository) must be present in the current working directory.

**-min_parsimony *I***: Minimum parsimony (must be an integer) default = 4 (Beware that reducing this parameter from the default
                     can obscure accuracy)

**-min_contribution *N***: Minimum percent of alternate alleles contributed by a
                     specific lab for variants to be considered highly
                     suspect; default = 80 (beware that reducing this
                     parameter from the default can obscure accuracy)
                    
**-fast**:           (highly recommended) Program will only consider highly suspicious lab-associated mutations and will not                         output flagged_snps_by_lab.tsv. However, it will run much faster!

### Output:

**flagged_snps_by_lab.tsv**:    A tab delimitted file displaying bin, source, lab, snp, ref, alt, global ref, global alt, fisher                              exact, proportion of calls, MAF for each flagged snp where:

| Column | Description |
| ------ | ----------- |
| bin | "highly suspect", "suspect," or "no flag"  |
| source | "originating lab" or "submitting lab" |
| lab | lab name |
| snp | snp |
| ref | reference allele count |
| alt | alternate allele count |
| global ref | global reference allele count |
| global alt | global alternate allele count |
| fisher exact | p value for fishers exact test associating lab specific with global allele frequency |
| proportion of calls | proportion of minor allele calls attributed to respective lab |
| MAF | minor allele frequency |

**flagged_snps_summary.tsv**:   A tab delimitted file displaying each flagged snp, bin and the reasoning underlying the flag

If the -track option is added:

**lab_associated_error_final.bedDetail**: A [bedDetail](https://genome.ucsc.edu/FAQ/FAQformat.html#format1.7) file annotating lab-associated mutations.  Contents can be pasted directly into the UCSC Genome Browser as a custom track.

**Artic_primers_final.bed**: A [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) custom track file annotating Artic Primers.  Contents can be pasted directly into the UCSC Genome Browser as a custom track.

**lab_associated_error_final.bedForBigBed**:     A [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file with extra data columns fully annotating lab-associated mutations.  Can be converted to [bigBed](https://genome.ucsc.edu/goldenPath/help/bigBed.html) format for viewing in the UCSC Genome Browser.

##########################################

## Viewing tracks in the UCSC Genome Browser

If you ran parsimony_per_site_analysis2.py with the `-track` option, then the script generated three extra output files suitable for viewing in the UCSC Genome Browser: lab_associated_error_final.bedDetail, Artic_primers_final.bed, and lab_associated_error_final.bedForBigBed.

The contents of lab_associated_error_final.bedDetail and Artic_primers_final.bed can be copy-pasted directly into the UCSC Genome Browser's [Custom Track](https://genome.ucsc.edu/cgi-bin/hgCustom) input form.

lab_associated_error_final.bedForBigBed has additional columns for a more complete set of annotation details to be displayed when you click on an error in the UCSC Genome Browser, but first it must be converted to the [bigBed](https://genome.ucsc.edu/goldenPath/help/bigBed.html) format as follows.

* Download the [linux](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed) or [Mac OS X](http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/bedToBigBed) version of `bedToBigBed` program, depending on your computer's platform.  After downloading, you may need to make the file executable using the command `chmod a+x bedToBigBed`.

* Run bedToBigBed like this (the labAssocMuts.as file is included in this repository):

```
./bedToBigBed -type=bed4+5 -as=labAssocMuts.as -tab lab_associated_error_final.bedForBigBed http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/bigZips/wuhCor1.chrom.sizes lab_associated_error_final.bb
```

* Place the resulting file `lab_associated_error_final.bb` on a website or FTP server (see UCSC's suggestions for [hosting](https://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html#Hosting)).

* Copy and paste the following custom track definition line into the UCSC Genome Browser's [Custom Track](https://genome.ucsc.edu/cgi-bin/hgCustom) input form, but replace the example link `http://mylab.org/lab_associated_error_final.bb` with the actual link to lab_associated_error_final.bb on your web/FTP server:

```
track name=lab_associated_errors description="Lab-associated Mutations" type=bigBed mouseOverField=comment visibility=pack bigDataUrl=http://mylab.org/lab_associated_error_final.bb"
```

##########################################

## GISAID metadata filter Usage:

We provide our metadata_1_merged.tsv output along with the our log file when we used metadata_1.tsv as input. We hand checked the log file by hand to ensure accuracy and commented out merger errors.

```
python3 meta_data_filter.py [options] -m [Path to GISAID metadata file] -o [Path to output directory]
```

### Input:

GISAID metadata file

### Options:

-c     If "TRUE" only merge high confidence errors.

### Output:

**metadata_filter.log**:     Log of changes made to the metadata set and warnings if any changes are made with lower confidence

**metadata_merged.tsv**:     Merged metadata file


##########################################

### Reference
* Yatish Turakhia, Bryan Thornlow, Landen Gozashti, Angie S. Hinrichs, Jason D. Fernandes, David Haussler, and Russell Corbett-Detig, "Stability of SARS-CoV-2 Phylogenies", bioRxiv [pre-print](https://www.biorxiv.org/content/10.1101/2020.06.08.141127v1) 2020.
