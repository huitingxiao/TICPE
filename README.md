# TICPE
TICPE is a cancer-specific qualitative method based on the relative expression orderings (REOs) of gene pairs within a sample to estimate the proportion of immune cells.

# Install
To install the TICPE, install from github using devtools
```
library(devtools)
install_github("huitingxiao/TICPE")
```
Or you can download the .ZIP file and unzip it.
```
install.packages("TICPE",repos = NULL,type="source")
#The "TICPE" should be combined with the absolute path.
```
# Usage
```
select_siggene=SelectSiggene(cancer_exp,immune_exp,marker_geneset,0.05)
stablepairs=StablePairs(cancer_exp,geneid,0.99)
parameter=Stimulatedmodel(combat_edata,select_siggene,stablepairs)
estimated_proportion=TICPEScores(expr,select_siggene,stablepairs,parameter,0.5)
```
Arguments|Description
:--|:---
cancer_exp|the gene expression data set of cancer cell lines, samples in columns, genes in rows, rownames corresponding to Entrez IDs.
immune_exp|the gene expression data set of different types of immune cells, samples in columns, genes in rows, rownames corresponding to Entrez IDs.
marker_geneset|a list of marker genes for different immune cell subsets.
combat_edata| the gene expression data including immune cells,cancer cells,and normal tissues using ComBat to remove batch effects,samples in columns, genes in rows, rownames corresponding to Entrez IDs.
expr| the gene expression data set of tumor samples,samples in columns, genes in rows, rownames corresponding to Entrez IDs.

This function performs all three steps in TICPE, which can be performed seperately as well:

1.SelectSiggene

2.Stimulatedmodel

3.TICPEScores

# Example
```
library(TICPE)
data(example)
estimated_proportion=TICPEScores(GSE39582_IHC_exp,CRC_select_siggene,CRC_stable_pairs,CRC_parameter,0.5)
```
The gene expression data are from the database Gene Expression Omnibus (GEO), **CRC_select_siggene**,**CRC_stable_pairs** and **CRC_parameter** are calculated for colorectal cancer.  
The **GSE39582_IHC_exp** is the colorectal cancer samples who had immunohistochemistry data of CD3, CD8, and CD68.
In addition, the relevant parameters of melanoma and lung cancer were added in the example.Rdata.

# Contact email
Please don't hesitate to address comments/questions/suggestions regarding this R package to: Huiting Xiao xht_fuj@hotmail.com
