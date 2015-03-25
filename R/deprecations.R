# Deprecated functions

#was deprecated 2015-03-24 
impBamGRL <- function()
{
	    .Deprecated("impBamGAL")
    ## use new function, or remainder of myOldFunc
}

#in 2014
getAlleleCount <- function() {
    .Deprecated("getAlleleCounts")
    ## use new function, or remainder of myOldFunc
}

#' Import Bam-2
#' 
#' Imports bla bal bal a specified genomic region from a bam file using a GenomicRanges
#' object as search area.
#' 
#' These functions are right  on tahea wrappers to import bam files into R and store them into
#' either GRanges, GAlignments or GappedAlignmentpairs objects.
#' 
#' It is recommended to use the impBamGAL() which takes information of gaps
#' into account. It is also possible to use the other variants as well, but
#' then pre-filtering becomes important keps to understand because gapped, intron-spanning reads
#' will cause problems. This is because the GRanges objects can not handle if
#' gaps are present and will then give a wrong result when calculating the
#' allele (SNP) count table.
#' 
#' @name import-bam-2
#' @rdname import-bam-2
#' @aliases import-bam-2 impBamGRL
#' @param UserDir The relative or full path of folder containing bam files.
#' @param searchArea A \code{GenomicRanges object} that contains the regions of
#' interest
#' @param verbose Setting \code{verbose=TRUE} gives details of procedure during
#' function run.
#' @return \code{impBamGRL} returns a GRangesList object containing the RNA-seq
#' reads in the region defined by the \code{searchArea} argument.
#' \code{impBamGAL} returns a list with GAlignments objects containing the
#' RNA-seq reads in the region defined by the \code{searchArea} argument.
#' \code{funImpBamGAPL} returns a list with GappedAlignmentPairs object
#' containing the RNA-seq reads in the region defined by the \code{searchArea}
#' argument.
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords bam import
#' @examples
#' 
#' #Declare searchArea
#' searchArea <- GRanges(seqnames=c('17'), ranges=IRanges(79478301,79478361))
#' 
#' #Relative or full path  
#' pathToFiles <- system.file('extdata/ERP000101_subset', package='AllelicImbalance')
#' 
#' 
#' @export impBamGRL
NULL


#' @rdname import-bam-2
impBamGRL <- function(UserDir, searchArea, verbose = TRUE) {
    # Set parameters
    which <- searchArea  #A GRanges, RangesList, RangedData, or missing object, from which a IRangesList instance will be constructed.
    what <- scanBamWhat()  #A character vector naming the fields to return. scanBamWhat() returns a vector of available fields. Fields are described on the scanBam help page.
    flag <- scanBamFlag(isUnmappedQuery = FALSE)
    param <- ScanBamParam(flag = flag, which = which, what = what)  #store ScanBamParam in param.
    
    # Point to correct directory and create a BamFileList object
    bamDir <- normalizePath(UserDir)
    allFiles <- list.files(bamDir, full.names = TRUE)
    bamFiles <- allFiles[grep(".bam$", allFiles)]
    if (length(bamFiles) == 0) {
        stop(paste("No bam files found in", bamDir))
    }
    if (!all(file.exists(paste(bamFiles, ".bai", sep = "")))) {
        if (verbose) {
            cat(paste("The bam files in UserDir are required to also have", ".bam.bai index files.", 
                " Trying to run indexBam function on each", "\n"), )
        }
        indexBam(bamFiles)
        if (!all(file.exists(paste(bamFiles, ".bai", sep = "")))) {
            stop("The bam files in UserDir are required to also have", ".bam.bai index files.")
        } else {
            if (verbose) {
                cat(paste("Succesfully indexed all bamFiles in UserDir", UserDir, 
                  "\n"))
            }
        }
    }
    # store all the .bam paths in a BamFile.
    bamFilesList <- BamFileList(bamFiles)
    
    # check that sequences in searchArea are actually found in the bam files
    header <- scanBamHeader(bamFiles)
    checkSeqNameExists <- function(bamHeader, requestedSeqNames) {
        as.character(requestedSeqNames) %in% names(bamHeader[["targets"]])
    }
    if (!all(unlist(lapply(header, checkSeqNameExists, seqnames(searchArea))))) {
        # not all searchArea requested seq-names found in bam files. Create nice error
        # report and stop
        seqNotFoundErrors <- lapply(header, checkSeqNameExists, seqnames(searchArea))
        seqNotFounds <- vector()
        for (sampleName in names(seqNotFoundErrors)) {
            seqNotFounds <- c(seqNotFounds, as.character(seqnames(searchArea)[!seqNotFoundErrors[[sampleName]]]))
        }
        stop(paste("The following seq name(s) not found in the bam files:", paste(sort(unique(seqNotFounds)), 
            collapse = ", ")))
    }
    
    # Loop through, open scanBam, store in GRList and then close each object in the
    # BamFileList object.
    i <- 1
    BamGRL <- GRangesList()
    for (bamName in names(bamFilesList)) {
        # Description
        bf <- bamFilesList[[bamName]]
        open(bf)
        if (verbose) {
            cat(paste("Reading bam file", i, "with filename", basename(bamName)), 
                "\n")
        }
        bam <- scanBam(bf, param = param)
        # Description
        for (rangeName in names(bam)) {
            
            # if NA values your in trouble. That means the read didnt map
            ranges <- IRanges(start = bam[[rangeName]][["pos"]], width = cigarWidthAlongReferenceSpace(bam[[rangeName]][["cigar"]]))
            GRangeBam <- GRanges(seqnames = as.character(bam[[rangeName]][["rname"]]), 
                ranges = ranges, strand = bam[[rangeName]][["strand"]], names = bam[[rangeName]][["qname"]], 
                flag = bam[[rangeName]][["flag"]], cigar = bam[[rangeName]][["cigar"]], 
                mapq = bam[[rangeName]][["mapq"]], mpos = bam[[rangeName]][["mpos"]], 
                isize = bam[[rangeName]][["isize"]], seq = bam[[rangeName]][["seq"]], 
                qual = bam[[rangeName]][["qual"]])
            # This way of merging the different chromosomes to the same GRangeObject is maybe
            # not the best way. Later try to store them in a separate list, and then unlist
            # before importing to GrangeBam Store GRangeBam in BamGRL (which is the GRange
            # List object)
            if (basename(bamName) %in% names(BamGRL)) {
                BamGRL[[basename(bamName)]] <- c(BamGRL[[basename(bamName)]], GRangeBam)
            } else {
                BamGRL[[basename(bamName)]] <- GRangeBam
            }
        }
        if (verbose) {
            cat(paste("stored", basename(bamName), "in BamGRL"), "\n")
        }
        i <- 1 + i
        gc()
        close(bf)
    }
    return(BamGRL)
}
