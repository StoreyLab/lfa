# lfa

Logistic factor analysis

## Installation

To install latest version on Bioconductor, open R and type:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("lfa")
```

You can also install development version from GitHub this way:
```R
install.packages("devtools")
library("devtools")
install_github("Storeylab/lfa")
```
Apple OS X users, see FAQ below.

## Data input

We recommend using the `genio` or `BEDMatrix` packages to read genotype data into an R matrix.

Be warned that genotype matrices from `genio` and some `lfa` functions require a lot of memory.
As a rule of thumb, the in memory sizes of a few relevant genotype matrices:

- 431345 SNPs by 940 individuals: 1.5 GB needed for genotype matrix, about 9 GB to run `lfa`.
- 339100 SNPs by 1500 individuals: 1.9 GB needed for genotype matrix, about 11.5 GB to run `lfa`.

`BEDMatrix` inputs consume much less memory but can be slower otherwise.

## FAQ

Apple OS X users may experience a problem due to Fortran code that is included in this package.  This gfortran issue is discussed here: 

http://www.thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error. 

A solution that has worked for us is to follow the advice given above. Specifically, open a Terminal and type:

```
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```

If you know what causes this, let us know!

## Citations

Hao, Wei, Minsun Song, and John D. Storey. "Probabilistic Models of Genetic Variation in Structured Populations Applied to Global Human Studies." Bioinformatics (Oxford, England) 32, no. 5 (March 1, 2016): 713–21. [doi:10.1093/bioinformatics/btv641](https://doi.org/10.1093/bioinformatics/btv641). [arXiv](http://arxiv.org/abs/1312.2041).

