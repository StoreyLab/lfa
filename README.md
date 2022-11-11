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
Apple OS X users, see Troubleshooting below.

## Data input

We recommend using the `genio` or `BEDMatrix` packages to read genotype data into an R matrix.

Be warned that genotype matrices from `genio` and some `lfa` functions require a lot of memory.
As a rule of thumb, the in memory sizes of a few relevant genotype matrices:

- 431345 SNPs by 940 individuals: 1.5 GB needed for genotype matrix, about 9 GB to run `lfa`.
- 339100 SNPs by 1500 individuals: 1.9 GB needed for genotype matrix, about 11.5 GB to run `lfa`.

`BEDMatrix` inputs consume much less memory but can be slower otherwise.

## Troubleshooting

Apple OS X users may experience a problem due to Fortran code that is included in this package. You must install the X code command line tools (XCode CLI) and `gfortran`.  Try the following commands on terminal:

```
xcode-select --install 
brew install gcc
```

If XCode installation fails, you may have to sign up on Apple Developer: https://www.ics.uci.edu/~pattis/common/handouts/macmingweclipse/allexperimental/macxcodecommandlinetools.html

Alternatively, this Installer Package for macOS R toolchain may work https://github.com/rmacoslib/r-macos-rtools/

## Citations

Hao, Wei, Minsun Song, and John D. Storey. "Probabilistic Models of Genetic Variation in Structured Populations Applied to Global Human Studies." Bioinformatics 32, no. 5 (March 1, 2016): 713–21. [doi:10.1093/bioinformatics/btv641](https://doi.org/10.1093/bioinformatics/btv641). [arXiv](https://arxiv.org/abs/1312.2041).

