#' @title Read .tped
#' @description Reads a .tped format genotype matrix and returns the R
#' object needed by \code{\link{lfa}}.
#' @details Use --transpose and --recode12 on your plink formatted genotypes
#' to generate the proper tped file. This is a pretty terrible function
#' that uses a growing matrix for the genotypes so it is to your 
#' benefit to have as large a \code{buffer.size} as possible. 
#' @param tped.filename Path to your .tped file after tranposing and recoding.
#' @param buffer.size Number of characters to keep in the buffer 
#' @examples
#' #assuming you have a .tped file in the right directory
#' x = NULL
#' \dontrun{x = read.tped.recode('file.tped')}
#' @return genotype matrix with elements 0, 1, 2, and NA.
#' @name read.tped.recode-deprecated
#' @usage read.tped.recode(tped.filename, buffer.size=5e8)
#' @seealso [lfa-deprecated()]
#' @keywords internal
NULL

#' @rdname lfa-deprecated
#' @section `read.tped.recode`:
#' For `read.tped.recode`, use `plink` (external binary) to convert to
#' BED/BIM/FAM, then parse with
#' [genio::read_plink()].
#' @export
read.tped.recode <- function(tped.filename, buffer.size = 5e+08) {
    .Deprecated(msg = "Use `plink` (external binary) for file conversions!")
    tped.line <- readLines(tped.filename, n = 1)
    if (nchar(tped.line) > buffer.size/10)
        warning("recommend increasing buffer")
    tped.line <- strsplit(tped.line, " ")[[1]]
    if (length(tped.line) <= 4)
        stop("expecting SNPs in tped (line length <= 4)")
    if (!(as.integer(tped.line[5]) %in% 0:2))
        stop("expecting -recode12")
    n <- (length(tped.line) - 4)/2
    message("reading in", n, "individuals")
    X <- NULL
    buffer <- NULL
    con <- file(tped.filename, "r")
    m <- 0
    while (TRUE) {
        buffer <- paste(buffer, readChar(con, buffer.size), sep = "")
        if (identical(buffer, character(0)))
            break
        in.lines <- strsplit(buffer, "\n")[[1]]
        new.m <- length(in.lines) - 1
        if (new.m < 2)
            stop("probably should increase buffer")
        if (substr(buffer, nchar(buffer), nchar(buffer)) == "\n") {
            new.m <- new.m + 1
            snps <- in.lines
            buffer <- NULL
        } else {
            snps <- in.lines[seq_len(new.m)]
            buffer <- in.lines[new.m + 1]
        }
        geno.tmp <- matrix(0, new.m, n)
        for (i in seq_len(new.m)) geno.tmp[i, ] <- .tped_line(in.lines[i])
        X <- rbind(X, geno.tmp)
        m <- m + new.m
        message("finished snp ", m)
    }

    close(con)
    X
}

.tped_line <- function(tped.line) {
    snps = strsplit(tped.line, " ")[[1]]
    if (length(snps) <= 4)
        stop("invalid tped (line length <= 4)")
    snps <- as.integer(snps[5:length(snps)])
    if (length(snps)%%2 == 1)
        stop("snp length error")
    even <- seq(2, length(snps), 2)
    odds <- seq(1, length(snps), 2)
    ret <- snps[even] + snps[odds] - 2
    ret[ret < 0] <- NA
    ret
}
