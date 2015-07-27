# Breast cAncer Subtype AnaLysis Tool

### Overview

Runs PAM50 and SCMGENE subtype methods on log2 expression data

### Requirements

A tab delimited file with samplenames as column IDs, gene symbols as row IDs and log2 expression data. Any gene symbol can be present in the file as long as ESR1 and ERBB2 are there, otherwise PAM50 won't work. The more gene symbols matching the PAM50 set the better.

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

Install some libraries

```
install.packages("ggplot2")
install.packages("gplots")
install.packages("org.Hs.eg.db")
install.packages("heatmap.plus")
install.packages("reshape")
install.packages("RColorBrewer")
install.packages("devtools")
source("http://bioconductor.org/biocLite.R")
biocLite("genefu")
biocLite("ctc")
```

Load them up

```
library(ggplot2)
library(gplots)
library(org.Hs.eg.db)
library(heatmap.plus)
library(reshape)
library(RColorBrewer)
library(genefu)
library(ctc)
library(devtools)
```

Install and load BASAL
```
install_github("elswob/BASAL")
libary(BASAL)
```

### To run

Need three arguements:

1. The name of the directory for the output
2. The full path and name of the expression data file
3. A short name for the data

Example breast cancer data from the Cancer Cell Line Encyclopedia (CCLE) is available for testing.

```
s=system.file("extdata","CCLE.tsv",package="BASAL")
run_basal("~/subtype_test",s,"ccle")
```

### Output

1. A combined tsv file for all subtype methods (subtype_summary.tsv)
2. Detailed output from PAM50 and SCMGENE


### Session Info

```
sessionInfo()
R version 3.2.0 (2015-04-16)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.9.5 (Mavericks)

locale:
[1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] impute_1.42.0        ctc_1.42.0           amap_0.8-14          RColorBrewer_1.1-2   reshape_0.8.5        heatmap.plus_1.3     genefu_1.18.0        biomaRt_2.24.0      
 [9] mclust_5.0.1         survcomp_1.18.0      prodlim_1.5.1        survival_2.38-3      org.Hs.eg.db_3.1.2   RSQLite_1.0.0        DBI_0.3.1            AnnotationDbi_1.30.1
[17] GenomeInfoDb_1.4.1   IRanges_2.2.5        S4Vectors_0.6.1      Biobase_2.28.0       BiocGenerics_0.14.0  gplots_2.17.0        ggplot2_1.0.1        BASAL_0.1           

loaded via a namespace (and not attached):
 [1] Rcpp_0.11.6        plyr_1.8.3         bitops_1.0-6       tools_3.2.0        digest_0.6.8       gtable_0.1.2       proto_0.3-10       bootstrap_2015.2   stringr_1.0.0     
[10] SuppDists_1.1-9.1  gtools_3.5.0       caTools_1.17.1     XML_3.98-1.3       gdata_2.17.0       lava_1.4.1         rmeta_2.16         reshape2_1.4.1     magrittr_1.5      
[19] scales_0.2.5       survivalROC_1.0.3  MASS_7.3-42        splines_3.2.0      colorspace_1.2-6   labeling_0.3       KernSmooth_2.23-15 stringi_0.5-5      RCurl_1.95-4.7    
[28] munsell_0.4.2   
```
