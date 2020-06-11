table labAssocMuts
"Lab-associated mutations identified using https://github.com/lgozasht/COVID-19-Lab-Specific-Bias-Filter"
    (
    string chrom;       "Reference sequence name"
    uint   chromStart;  "Start position in reference sequence"
    uint   chromEnd;    "End position in reference sequence"
    string name;        "Mutation name (reference base value, position, alternate base value)"
    int parsimony;      "Number of base value changes in most parsimonious labeling of tree"
    int altAlleleCount; "Alternate allele count"
    float maf;          "Minor Allele Frequency"
    string articPrimer; "ARTIC primer(s) within 10bp of mutation or 'NA' if not applicable"
    string aaChange;    "Amino acid change caused by mutation or 'NA' if not applicable"
    string comment;     "Explanation for inclusion in lab-associated category"
    )
