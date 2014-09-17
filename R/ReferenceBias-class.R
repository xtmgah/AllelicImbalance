#'@include initialize-methods.R
NULL

#' ReferenceBias class
#' 
#' Object that holds Reference bias information
#'
#' Reference bias-class contain annotated fractions of the reference allele
#'
#' Shortly the ReferenceBias-class most prominent purpose is quikly to be able to dispatch on
#' methods like, 'plot', 'summary' and similar.
#' 
#' @name ReferenceBias-class
#' @rdname ReferenceBias-class
#' @aliases ReferenceBias-class ReferenceBias ReferenceBias-method
#' @docType class
#' @param x ReferenceBias object
#' @param verbose makes function more talkative
#' @return An object of class ReferenceBias storing reference fractions.

#' @section Constructor: refBias(x = ASEset)
#' 
#' \describe{
#' 
#' Arguments: \item{x}{an ASEset object  }
#'
#' }
#'
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords class ASEset
#' @examples
#'
#' data(ASEset)
#' a <- ASEset
#' genotype(a) <- inferGenotypes(a)
#' a <- refAllele(a,
#'  	fasta=system.file('extdata/hg19.chr17.fa', 
#'  	package='AllelicImbalance'))	
#' refbiasObject <- refBias(a)
#' 
#' @exportClass ReferenceBias
#' @export

setClass("ReferenceBias", contains = "SummarizedExperiment")
	

#setClass("ReferenceBias", representation(
#			refFraction="array"))
#
#.valid_RefBias_object <- function(object) {
#
#	#check dimensions
#	#if(!length(dim(object))%in% c(1,2,3)){
#	#	stop("maximum size of dimension is 3")
#	#}
#
#	return(TRUE)
#}
#setValidity("ReferenceBias", .valid_RefBias_object)
#

# @rdname ReferenceBias-class
#setGeneric("table")
#setGeneric("table", function(x, strand = "*", sortBy="none", ...) {
#    standardGeneric("table")
#})

#setMethod("table", signature(... = "ReferenceBias"), function(...) {
#
#	args <- list(...)
#	if (length(args) > 1)
#	  stop("Only one argument in '...' supported")
#	x <- args[[1L]]
#
#	#because the generis of table is rubbish we have to return a list for each strand
#	retList <- list()
#
#	for(strand in c("+","-","*")){
#
#		lst <- list()
#		y <- x@refFraction[,,strand]
#		mean.na.rm <- function(x){mean(x,na.rm=TRUE)}
#		lst[["samples"]] <- apply(y,2,mean.na.rm)
#		lst[["SNPs"]] <- apply(y,1,mean.na.rm)
#		lst[["all"]] <- mean(y,na.rm=TRUE)
#		
#		retList[[strand]] <- lst
#
#	}
#	return(SimpleList(retList))
#
#})	
#
