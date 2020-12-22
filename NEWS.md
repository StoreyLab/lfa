# 2020-09-18 - lfa 2.0.0.9000

Major overhaul from last version (1.9.0, last updated 2018-02-11).
Overall, added unit testing to all functions, which resulted in the identification and fixing of several edge-case bugs, and also minor improvements.

- User-facing changes
  - Removed redundant/obsolete exported functions:
    - `model.gof`: entirely redundant with `sHWE` in this same package
    - `read.tped.recode`: obsolete file format; instead, use `plink` for file conversions!
    - `read.bed`: instead, use `read_bed` from the `genio` package!
	- `center`: only worked for matrices without missingness (useless for real data), no dependencies in code, plus centering is trivial in R
  - Renamed remaining functions, replacing periods with underscores.
    - Only specific change: `trunc.svd` -> `trunc_svd`
	- NOTE `af_snp` and `pca_af` already had underscores instead of periods (unchanged).
	  All other functions unchanged.
  - Function `trunc_svd`
    - Debugged `d=1` case (output matrices were dropped to vectors, which was a fatal error in `lfa` after it called `trunc_svd`).
    - Added option `force` that, when `TRUE`, forces the Lanczos algorithm to be used in all cases (most useful for testing purposes).
  - Function `lfa`
    - Improved documentation.
  - Functions `af_snp` and `af`
	- Fixed a bug in which large logistic factor coefficients resulted in `NaN` allele frequencies instead of 1 as they should be in the limit.
    - Improved code to "impute" allele frequencies for `NA` genotypes.
      Original version preserved `NA` values (a genotype that was `NA` for a particular individual and locus combination resulted in an `NA` in the corresponding allele frequency only, and conversely non-`NA` genotypes never resulted in `NA` allele frequencies).
	  The new version takes advantage of knowing the LFs of all individuals (regardless of genotype missingness), and LFs and their coefficients are never `NA`, permitting allele frequencies to always be imputed into non-`NA` values!
  - Function `pca_af`
    - Similarly imputes allele frequencies for `NA` genotypes (internal change was trivial)
	- Debugged `d=1` case, which incorrectly returned an intercept column matrix instead of an allele frequency matrix.
  - Function `check_geno`
    - Debugged testing for matrix class (original code when run in R 4.0 generated warning "the condition has length > 1 and only the first element will be used")
  - Function `sHWE` (through internal `inverse_2x2`)
    - When a test was "singular" at a single SNP, function used to die; now that SNP gets an `NA` p-value.
	- Other previous `NA` cases here are avoided now that `af` never returns `NA` values.

- Internal changes
  - Separated R functions into one source file each.
  - Added `.gitignore` files from another project.
    - Removed `src/lfa.so` from version control tracking.
  - Added unit tests for all functions using `testthat`.
  - Updates to C code solely to pass latest `R CMD check` requirements.

# 2020-11-11 - lfa 2.0.1.9000

- Function `lfa` added support for BEDMatrix objects for the genotype matrix `X`.
  - This consumes lower memory when the number of loci `m` is very large, so it enables analysis of larger datasets.
  - Algorithm for BEDMatrix case is different: instead of Lanczos truncated SVD, covariance matrices are computed explicitly and truncated eigendecomposition performed.  This means runtime and memory complexity are very different here as the number of individuals `n` gets larger.
  - Added `RSpectra` package dependency (for fast truncated eigendecomposition).

# 2020-11-12 - lfa 2.0.2.9000

- More functions updated to support BEDMatrix inputs for the genotype matrix `X`.  Although BEDMatrix is supported, in these cases there are minimal memory reduction advantages as outputs or intermediate matrices are necessarily as large as the input genotype data.
  - Function `af`.  Although there is memory saving by not loading `X` entirely into memory, the output individual-specific allele frequency matrix `P` has the same dimensions so memory usage may still be excessive for in large datasets, negating the BEDMatrix advantage.
  - Function `pca_af`.  Note same memory usage issue as `af`.
  - Function `sHWE`.  A worse memory problem is present, since internal code calculates the entire individual-specific allele frequency matrix `P`, then simulates `B` random genotype matrices of the same dimensions as input (each in memory) from which `LF` and ultimately HWE statistics are obtained.

# 2020-12-16 - lfa 2.0.3.9000

- Fixed an integer overflow error that occurred in `sHWE` (in internal function `compute_nulls`), which only occurred if the number of individuals times the number of loci exceeded the maximum integer size in R (the value of `.Machine$integer.max`, which is 2,147,483,647 in my 64-bit machine).
- Function `lfa` added `rspectra` option (`FALSE` by default), as an alternative way of calculating SVD internally (for experimental use only).
- Function `trunc_svd` is now exported.
- Minor, user-imperceptible changes in other functions.

# 2020-12-22 - lfa 2.0.4.9000

- Function `sHWE` fixed bug: an error could occur when internal statistics vector included `NA` values.
  - Original error gave this obscure message, which occurred because an index went out of bounds due to a discrepancy in vector lengths due to the presence of `NA` values:
```
Error in while ((i2 <= B0) & (obs_stat[i1] >= stat0[i2])) { : 
  missing value where TRUE/FALSE needed
```
  - Now empirical p-value code is separated into new internal function `pvals_empir`, and its tested against a new naive version `pvals_empir_brute` (slower brute-force algorithm, used to validate outputs only) in unit tests including simulated data with `NA` values.
  - Also refactored other internal `sHWE` code into a new internal function `gof_stat`, which by itself better handles BEDMatrix files (though overall memory savings are not yet there on the broader `sHWE`).
- Spell-checked this news file (edited earlier entries).
