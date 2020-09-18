# 2020-09-18 - lfa 2.0.0.9000

Major overhaul from last version (1.9.0, last updated 2018-02-11).
Overall, added unit testing to all functions, which resulted in the identification and fixing of several edge-case bugs, and also minor improvements.

- User-facing changes
  - Removed redundant/obsolete exported functions:
    - `model.gof`: entirely redundant with `sHWE` in this same package
    - `read.tped.recode`: obsolete file format; instead, use `plink` for file conversions!
    - `read.bed`: instead, use `read_bed` from the `genio` package!
	- `center`: only worked for matrices without missingness (useless for real data), no depenencies in code, plus centering is trivial in R
  - Renamed remaining functions, replacing periods with underscores.
    - Only specific change: `trunc.svd` -> `trunc_svd`
	- NOTE `af_snp` and `pca_af` already had undescores instead of periods (unchanged).
	  All other functions unchanged.
  - Function `trunc_svd`
    - Debugged `d=1` case (output matrices were dropped to vectors, which was a fatal error in `lfa` after it called `trunc_svd`).
    - Added option `force` that, when `TRUE`, forces the Lanczos algorithm to be used in all cases (most useful for testing purposes).
  - Function `lfa`
    - Improved documentation.
  - Functions `af_snp` and `af`
	- Fixed a bug in which large logistic factor coefficients resulted in `NaN` allele frequencies instead of 1 as they should be in the limit.
    - Improved code to "impute" allele frequencies for `NA` genotypes.
      Original version preserved `NA` values (a genotype that was `NA` for a particular individual and locus combination resulted in an `NA` in the corresponding allele frequency only, and conversely non-`NA` genotypes never resulted in `NA` alele frequencies).
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
