# Breast cAncer Subtype AnaLysis

### Overview

Runs PAM50 and SCMGENE subtype methods on a log2 expression data

### Requirements

A file with samplenames as column IDs and gene symbols as row IDs

```
Hugo_symbol	EFM192A	EVSAT	MDAMB453	MDAMB231	MDAMB468	BT549	CAL851	EFM19
A1BG	8.315511	5.454009	7.822566	6.140037	5.111359	7.599171	5.173046	6.038426
A2M	4.107914	3.714508	3.907026	4.059009	7.406219	3.973131	4.741874	3.801728
NAT1	8.914029	8.343696	9.054734	8.80809	8.708712	7.315499	8.805189	6.833612
NAT2	6.564016	4.953733	4.763946	4.644254	4.364252	4.228607	5.419792	4.522074
SERPINA3	8.339674	5.361677	4.422622	4.641157	9.53551	5.258785	11.69059	8.127704
AADAC	4.326528	4.037927	3.680285	3.766857	3.94583	3.694712	4.694278	3.882659
AAMP	9.266102	10.43423	9.745142	8.790493	9.474694	8.564627	8.914895	9.234852
AANAT	4.118585	4.168364	3.84862	3.839901	4.013059	4.515908	3.721725	3.903788
AARS	10.73831	10.82067	10.40384	9.457173	10.26427	10.38802	10.48933	9.772005
``` 

### Installation

```
install.packages("devtools")
devtools::install_github("elswob/BASAL")
```

### To run

Need three arguements:

1. The name of the directory containing the expression data
2. The name of the expression data file
3. A short name for the data

```
run_basal("ccle_data","CCLE.tsv","ccle")
```

### Output

1. A combined tsv file for all subtype methods (subtype_summary.tsv)
2. Detailed output from PAM50 and SCMGENE