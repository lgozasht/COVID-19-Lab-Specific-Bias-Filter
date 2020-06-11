table labAssocMuts
"Lab-associated mutations identified using https://github.com/lgozasht/COVID-19-Lab-Specific-Bias-Filter"
    (
    string chrom;       "Reference sequence chromosome or scaffold"
    uint   chromStart;  "Start position in chromosome"
    uint   chromEnd;    "End position in chromosome"
    string name;        "Mutation name"
    int parsimony;      "Number of base value changes in most parsimonious labeling of tree"
    int altAlleleCount; "Alternate allele count"
    float maf;          "Minor Allele Frequency"
    string articPrimer; "ARTIC primer(s) within 10bp of mutation, if applicable"
    string comment;     "Explanation for inclusion in lab-associated category"
    )
