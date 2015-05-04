lfa
===

logistic factor analysis

Pre-print available: http://arxiv.org/abs/1312.2041

Installation
===

Apple OS X users, see FAQ below.  To install, open R and type:

```R
install.packages("devtools")
library("devtools")
install_github("Storeylab/lfa")
```

Data input
===

Recommend using `read.bed` to read in a PLINK binary genotype file, or manually reading a text file of 0, 1, and 2 via `read.table` with the `colClasses` option set to `rep("integer", n)`, where `n` is the number of individuals.

Be warned that both the `lfa` function and the genotype matrices require a lot of memory, since we're just storing the genotype matrices as integer arrays (which take up a lot of space in R!) and doing normal math operations on them. As a rule of thumb, the in memory sizes of a few relevant genotype matrices:

- 431345 SNPs by 940 individuals: 1.5 Gb needed for genotype matrix, about 9 Gb to run `lfa`.
- 339100 SNPs by 1500 individuals: 1.9 Gb needed for genotype matrix, about 11.5 Gb to run `lfa`.

If you want to run something but don't have the computing power for it, shoot me an email at whao@princeton.edu. There are number of workarounds that can trade off between memory and time that I've used in the past.


FAQ
===

Apple OS X users may experience a problem due to Fortran code that is included in this package.  This gfortran issue is discussed here: 

http://www.thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error. 

A solution that has worked for us is to follow the advice given above. Specifically, open a Terminal and type:

```
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```

If you know what causes this, let us know!

