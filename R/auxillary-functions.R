#'@include ASEset-class.R
NULL


#' ASEset from bam file
#' 
#' count alleles and create an ASEset direct from bam file instead of reading into R first.
#' 
#' counts the alleles in a bam file based on GRanges positions. 
#' 
#' 
#' @param gr GenomicRanges of SNPs to create ASEset for
#' @param PE if paired end or not (default: TRUE)
#' @param pathToDir Directory of bam files with index in same directory
#' @param strandUnknown default: FALSE
#' @param ... passed on to countAllelesFromBam function
#' @param flagsMinusStrand flags that mark reads coming from minus strand
#' @param flagsPlusStrand flags that mark reads coming from plus strand
#' @author Jesper R. Gadin
#' @keywords ASEset
#' @examples
#'
#' data(GRvariants)
#' gr <- GRvariants
#'
#' ##no execution at the moment
#' #pathToDir <- system.file('inst/extdata/ERP000101_subset', package='AllelicImbalance')
#' #a <- ASEsetFromBam(gr, pathToDir)
#'  
#' @export ASEsetFromBam

ASEsetFromBam <- function(gr, pathToDir,PE=TRUE, flagsMinusStrand=c(83,163), flagsPlusStrand=c(99,147), strandUnknown=FALSE, ...) {

	if(!PE){
		stop("no support for SE atm")
	}

	if(PE==TRUE){
		#minus strand
		arm1 <- countAllelesFromBam(gr, pathToDir, flag=83)
		arm2 <- countAllelesFromBam(gr, pathToDir, flag=163)
		arm <- arm1 + arm2

		#plus strand
		arp1 <- countAllelesFromBam(gr, pathToDir, flag=99)
		arp2 <- countAllelesFromBam(gr, pathToDir, flag=147)
		arp <- arp1 + arp2
	}
	
	#ASEsetFromArray
	if(!strandUnknown){
		a <- ASEsetFromArrays(gr, countsPlus = arp, 
			countsMinus = arm)

	}else{	
		a <- ASEsetFromArrays(gr, countsUnknown = arp+arm) 
	}
	a
}


#' makes masked fasta reference
#' 
#' Replaces all selected positions in a fasta file with the character N
#' 
#' @param fastaIn character string of the path for the fasta file to be used
#' @param fastaOut character string of the path for the masked fasta file (no extension)
#' @param posToReplace GRanges object with the genomic ranges to replace 
#' @param splitOnSeqlevels write on file for each seqlevel to save memory
#' @param verbose makes function more talkative
#' @author Jesper R. Gadin
#' @keywords masked fasta reference
#' @examples
#'
#' data(ASEset.sim)
#' gr <- rowRanges(ASEset.sim) 
#' fastaIn <- system.file('extdata/hg19.chr17.subset.fa', package='AllelicImbalance')
#' makeMaskedFasta(fastaIn=fastaIn, fastaOut="fastaOut",posToReplace=gr) 
#' 
#'  
#' @export makeMaskedFasta
makeMaskedFasta <- function(fastaIn, fastaOut, posToReplace, splitOnSeqlevels=TRUE, verbose=TRUE){

	#does the inFasta file exist?
	if(!file.exists(fastaIn)){
		stop("fasta infile doesnt exist")
	}

	#make FaFile
	fl <- FaFile(fastaIn)

	#check if the index file is present otherwise tell the user to use the
	#indexFa(FaFile("pathToReference")) command
	if(!file.exists(paste(fastaIn,".fai",sep=""))){
		cat("could not find index file\n")
		cat("creates a new index file")
		indexFa(FaFile(fastaIn))
		cat("finished creating new index file")
	}
	#indexFa(fl) #creates a new file as index in the same directory but
	#with extension *.fai 

	#open,scan,close file
	open(fl)
	#check if seqleveles are present
	fa.info <- scanFaIndex(fl)
	if(!(all(seqlevels(posToReplace) %in% seqlevels(fa.info)))){
		close(fl)
		stop("seqlevels in object x are not in fasta index file")
	}

	###############LOOOP
	chrs <- seqlevels(posToReplace)
	for (chr in chrs){
	#	chr <- "chr1"

		searchArea <- fa.info[seqnames(fa.info)==chr]
		seq <- scanFa(fl,param=searchArea)

		#replace the SNPs with N
		toReplace <- unique(as.integer(ranges(posToReplace)))
		seq <- seq[[1]]
		seq[toReplace] <- "N"

		if(verbose){cat("replaced", length(toReplace),"instances with N\n")}

		#write new file
		#library("seqinr")
		outfile <- paste(fastaOut,chr,".fa",sep="")

		write.fasta(as.character(seq),names=chr, file.out=outfile, nbchar=80)	
		if(verbose){cat("wrote chr",chr,"to file\n")}

		gc()
	}
	close(fl)
	if(verbose){cat("all chromosomes written to file\n")}
}



#' global analysis wrapper
#' 
#' A wrapper to make a global analysis based on paths for BAM, VCF and GFF files
#' 
#' @param pathBam path to bam file
#' @param pathVcf path to vcf file
#' @param pathGFF path to gff file
#' @param verbose makes function more talkative
#' @author Jesper R. Gadin
#' @keywords global wrapper
#' @examples
#'
#' #empty as function doesn't exist
#' 
#' @export 
gba <- function(pathBam,pathVcf,pathGFF=NULL, verbose){

	#summarize counts
	
	#detectAI

	
}



