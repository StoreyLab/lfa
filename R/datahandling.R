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
#' \dontrun{x = read.tped.recode("file.tped")}
#' @return genotype matrix with elements 0, 1, 2, and NA.
#' @export
read.tped.recode <- function(tped.filename, buffer.size=5e8){
    
    # check first line of tped to set up the parameters
    tped.line <- readLines(tped.filename, n=1)
    if(nchar(tped.line) > buffer.size/10){
        warning("recommend increasing buffer")
    }
    tped.line <- strsplit(tped.line, " ")[[1]]
    if(length(tped.line) <=4){
        stop("expecting SNPs in tped (line length <= 4)")
    }
    if(!(as.integer(tped.line[5]) %in% 0:2)) {
        stop("expecting -recode12")
    }
    
    n <- (length(tped.line)-4) / 2
    
    print(paste("reading in", n, "individuals"))
    
    # prepare genotype matrix and buffer
    X <- NULL
    buffer <- NULL
    con <- file(tped.filename, "r")
    m <- 0
    
    # use readChar() to read in the genotype matrix to fill buffer
    while(TRUE) {
        buffer <- paste(buffer, readChar(con, buffer.size), sep="")
        if(identical(buffer, character(0)))
            break
        
        # count the number of newlines in tped.line
        in.lines <- strsplit(buffer, "\n")[[1]]
        new.m <- length(in.lines) - 1
        if(new.m < 2){
            stop("probably should increase buffer")
        }
        
        # check if last element of buffer is a newline, then set new buffer
        if(substr(buffer, nchar(buffer), nchar(buffer)) == "\n") {
            new.m <- new.m + 1
            snps <- in.lines
            buffer <- NULL
        } else{
            snps <- in.lines[1:new.m]
            buffer <- in.lines[new.m+1]
        }
        
        geno.tmp <- matrix(0, new.m, n)
        for(i in 1:new.m){
            geno.tmp[i,] <- process.tped.recode.line(in.lines[i])
        }
        
        X <- rbind(X, geno.tmp)
        m <- m + new.m
        print(paste("finished snp", m))
    }
    
    close(con)
    X
}

process.tped.recode.line <- function(tped.line){
    snps = strsplit(tped.line, " ")[[1]]
    if (length(snps) <= 4) {
        stop("error in process.tped.recode.line: invalid tped (line length <= 4)")
    }
    
    snps <- as.integer(snps[5:length(snps)])
    
    if(length(snps) %% 2 == 1){
        stop("snp length error")
    }
    
    even <- seq(2, length(snps), 2)
    odds <- seq(1, length(snps), 2)

    ret <- snps[even] + snps[odds] - 2
    ret[ret<0] <- NA
    
    ret
}

#' @title File input: .bed
#' @description Reads in genotypes in .bed format with corresponding bim
#' and fam files
#' @details Use plink with --make-bed
#' @return Genotype matrix
#' @param bed.prefix Path leading to the bed, bim, and fam files.
#' @examples
#' # assuming you have PLINK format HapMap data from: http://pngu.mgh.harvard.edu/~purcell/plink/res.shtml
#' # run this in the unpacked folder
#' x = NULL
#' \dontrun{x = read.bed("hapmap_r23a")}
#' @export
read.bed <- function(bed.prefix){
    bed.filename <- paste(bed.prefix, ".bed", sep="")
    bim.filename <- paste(bed.prefix, ".bim", sep="")
    fam.filename <- paste(bed.prefix, ".fam", sep="")
    
    if(!file.exists(bed.filename))
        stop("need .bed file")
    if(!file.exists(bim.filename))
        stop("need .bim file")
    if(!file.exists(fam.filename))
        stop("need .fam file")
    
    #figure out number of individuals by counting newlines in fam
    buffer <- read.table(fam.filename, stringsAsFactors=FALSE, colClasses="character")
    n <- nrow(buffer)
    
    print(paste("reading in", n, "individuals"))
        
    #figure out number of SNPs by counting newlines in bim
    buffer <- read.table(bim.filename, stringsAsFactors=FALSE, colClasses="character")
    m <- nrow(buffer)
    rm(buffer)

    print(paste("reading in", m, "snps"))
    
    #initialize genotype matrix
    X <- matrix(0, m, n)
    snp.map <- binary.genotype.map()
    
    #open stream
    bed <- file(bed.filename, "rb")
    
    #check that beginning of bim satisfies magic numbers
    if(readBin(bed, what="integer", n=1, size=1) != 108)
        stop("not valid bed file (magic number fail)")
    if(readBin(bed, what="integer", n=1, size=1) != 27)
        stop("not valid bed file (magic number fail)")
    buffer <- readBin(bed, what="integer", n=1, size=1)
    if(buffer == 0) {
        stop("individual major mode not yet supported")
    } else if(buffer==1) {
        print("snp major mode")
    } else{
        stop("bed mode problem")
    }
    
    #calculate the blocksize
    numbytes <- ceiling(n/4)

    
    #read in SNPs!
    for(i in 1:m){
        indices <- readBin(bed, what="int", n=numbytes, size=1, 
	  signed=FALSE)+1
        snp.in <- snp.map[,indices]
        X[i,] <- as.vector(snp.in[1:n])
        
        if(i %% 20000 == 0)
            print(paste("reading snp", i))
    }
    
    close(bed)
    X
}

#generates signed integer to genotypes conversion matrix
#columns of matrix coorespond to unsigned integers PLUS ONE
binary.genotype.map <- function(){    
    combinations <- as.matrix(expand.grid(0:3,0:3,0:3,0:3))
    snp.map <- matrix(0, 4, 256)
    colnames(combinations) <- NULL
    
    #again offset by 1
    bitstring <- list()
    bitstring[[1]] <- "00"
    bitstring[[2]] <- "01"
    bitstring[[3]] <- "10"
    bitstring[[4]] <- "11"
    
    #generate the indices
    indices <- apply(combinations, 1, function(x){
            strtoi(paste(bitstring[[x[1]+1]],
                         bitstring[[x[2]+1]],
                         bitstring[[x[3]+1]],
                         bitstring[[x[4]+1]], sep=""), base=2)})
    
    indices <- indices+1
    
    #this order matters...
    combinations[combinations==1] <- NA #PLINK IS BACKWARDS
    combinations[combinations==2] <- 1  #PLINK IS BACKWARDS
    combinations[combinations==0] <- 2  #PLINK IS BACKWARDS
    combinations[combinations==3] <- 0  #PLINK IS BACKWARDS
    
    snp.map <- apply(combinations, 1, rev)
    snp.map[,indices] <- snp.map
    
    snp.map
}

