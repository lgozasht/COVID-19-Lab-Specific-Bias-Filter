# COVID-19-Lab-Specific-Bias-Filter
We developed a tool to systematically flag COVID-19 genomes for variants resulting from possible lab specific biases. The program requires a concurrent nextstrain VCF file, nextstrain parsimony file and nextstrain metadata file as input. It outputs two tsv files: one with each flagged snp and the respective reasoning underlying the flag, and another providing more detailed information on each flagged lab.
