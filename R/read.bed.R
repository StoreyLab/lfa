#' @title File input: .bed
#' @description Reads in genotypes in .bed format with corresponding bim
#' and fam files
#' @details Use plink with --make-bed
#' @return Genotype matrix
#' @param bed.prefix Path leading to the bed, bim, and fam files.
#' @name read.bed-deprecated
#' @usage read.bed(bed.prefix)
#' @seealso [lfa-deprecated()]
#' @keywords internal
NULL

#' @rdname lfa-deprecated
#' @section `read.bed`:
#' For `read.bed`, use [genio::read_plink()].
#' @export
read.bed <- function(bed.prefix) {
    .Deprecated("genio::read_plink")
    bed.filename <- paste(bed.prefix, ".bed", sep = "")
    bim.filename <- paste(bed.prefix, ".bim", sep = "")
    fam.filename <- paste(bed.prefix, ".fam", sep = "")
    if (!file.exists(bed.filename))
        stop("need .bed file")
    if (!file.exists(bim.filename))
        stop("need .bim file")
    if (!file.exists(fam.filename))
        stop("need .fam file")
    buffer <- utils::read.table(fam.filename, colClasses = "character")
    n <- nrow(buffer)
    buffer <- utils::read.table(bim.filename, colClasses = "character")
    m <- nrow(buffer)
    rm(buffer)
    X <- matrix(0, m, n)
    snp.map <- binary.genotype.map()
    bed <- file(bed.filename, "rb")
    if (readBin(bed, what = "integer", n = 1, size = 1) != 108)
        stop("not valid bed file (magic number fail)")
    if (readBin(bed, what = "integer", n = 1, size = 1) != 27)
        stop("not valid bed file (magic number fail)")
    buffer <- readBin(bed, what = "integer", n = 1, size = 1)
    if (buffer == 0) {
        stop("individual major mode not yet supported")
    } else if (buffer == 1) {
        print("snp major mode")
    } else {
        stop("bed mode problem")
    }
    numbytes <- ceiling(n/4)
    for (i in seq_len(m)) {
        indices <- readBin(bed, what = "int", n = numbytes, size = 1,
            signed = FALSE) + 1
        snp.in <- snp.map[, indices]
        X[i, ] <- as.vector(snp.in[seq_len(n)])
    }
    close(bed)
    X
}

binary.genotype.map <- function() {
    combinations <- as.matrix(expand.grid(0:3, 0:3, 0:3, 0:3))
    snp.map <- matrix(0, 4, 256)
    colnames(combinations) <- NULL
    bitstring <- list()
    bitstring[[1]] <- "00"
    bitstring[[2]] <- "01"
    bitstring[[3]] <- "10"
    bitstring[[4]] <- "11"
    indices <- apply(combinations, 1, function(x) {
        strtoi(paste(bitstring[[x[1] + 1]], bitstring[[x[2] + 1]],
            bitstring[[x[3] + 1]], bitstring[[x[4] + 1]], sep = ""), base = 2)
    })
    indices <- indices + 1
    combinations[combinations == 1] <- NA  #PLINK IS BACKWARDS
    combinations[combinations == 2] <- 1  #PLINK IS BACKWARDS
    combinations[combinations == 0] <- 2  #PLINK IS BACKWARDS
    combinations[combinations == 3] <- 0  #PLINK IS BACKWARDS
    snp.map <- apply(combinations, 1, rev)
    snp.map[, indices] <- snp.map
    snp.map
}
