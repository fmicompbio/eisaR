# `eisaR`: Exon-Intron Split Analysis (EISA) in R

<br>

## Overview

Exon-intron split analysis (`EISA`) uses ordinary RNA-seq data to measure
changes in mature RNA and pre-mRNA reads across different experimental
conditions to quantify transcriptional and post-transcriptional regulation
of gene expression.

For details see Gaidatzis et al., Nat Biotechnol 2015. doi: 10.1038/nbt.3269.  
eisaR implements the major steps of EISA in R.

Current contributors include:

- [Michael Stadler](https://github.com/mbstadler)
- Dimos Gaidatzis
- [Lukas Burger](https://github.com/LukasBurger)
- [Charlotte Soneson](https://github.com/csoneson)

## Functionality

All you need is RNA-seq data from at least two conditions (e.g. wildtype and
mutant). The `eisaR` package contains convenience functions to facilitate the
steps in an exon-intron split analysis, which consists of:  

1. preparing the annotation (exonic and gene body coordinate ranges)  
2. quantifying RNA-seq alignments in exons and introns  
3. calculating and comparing exonic and intronic changes across conditions  
4. visualizing the results  

For the steps 1. and 2. above, this `eisaR` vignette makes use of Bioconductor
annotation and the [QuasR](https://bioconductor.org/packages/QuasR/) package.
It is also possible to obtain count tables for exons and introns using some
other pipeline or approach, and directly start with step 3.


## Reference
`EISA` has been described in:  

"Analysis of intronic and exonic reads in RNA-seq data characterizes
transcriptional and post-transcriptional regulation."  
Gaidatzis D., Burger L., Florescu M. and Stadler, M.B.  
*Nat Biotechnol.* **2015**; 33(7):722-9.
[PubMed: 26098447](https://www.ncbi.nlm.nih.gov/pubmed/26098447), [doi: 10.1038/nbt.3269](https://doi.org/10.1038/nbt.3269)

<!--
## Download from Bioconductor
[QuasR download page](https://bioconductor.org/packages/QuasR/)

## Software status

| Platforms |  OS  | R CMD check | Coverage | 
|:----------------:|:----------------:|:----------------:|:----------------:|
| Travis CI | Linux | [![Travis CI build status](https://travis-ci.com/fmicompbio/QuasR.svg?branch=master)](https://travis-ci.com/fmicompbio/QuasR) | [![Codecov.io coverage status](https://codecov.io/github/fmicompbio/QuasR/coverage.svg?branch=master)](https://codecov.io/github/fmicompbio/QuasR) |
| Bioc ([_devel_](http://bioconductor.org/packages/devel/bioc/html/QuasR.html)) | Multiple | [![Bioconductor-devel Build Status](http://bioconductor.org/shields/build/devel/bioc/QuasR.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/QuasR) | `NA` |
| Bioc ([_release_](http://bioconductor.org/packages/release/bioc/html/QuasR.html)) | Multiple | [![Bioconductor-release Build Status](http://bioconductor.org/shields/build/release/bioc/QuasR.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/QuasR) | `NA` |
-->
