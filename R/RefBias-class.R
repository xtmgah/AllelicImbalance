#'@include initialize-methods.R
NULL

#' RefBias class
#' 
#' Object that holds RefBias information
#'
#' RefBias-class is in principle an array with the fractions of the reference allele
#' in the dimensions SNP, sample, and strand
#'
#' Shortly the RefBias-class most prominent purpose is to be able to dispatch on
#' methods like, 'plot', 'table', 'summary' and similar.
#' 
#' @name RefBias-class
#' @rdname RefBias-class
#' @aliases RefBias-class RefBias 
#' @docType class
#' @param x RefBias object
#' @param strand which strand of '+', '-' or '*'
#' @param verbose makes function more talkative
#' @return An object of class RefBias storing reference fractions.

#' @section Constructor: RefBias(x = ASEset)
#' 
#' \describe{
#' 
#' Arguments: \item{x}{an ASEset object  }
#'
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords class ASEset
#' @examples
#'
#' data(ASEset)
#' refbiasObject <- refBiad(ASEset)
#'
#' @exportClass RefBias
#' 
#' @export frequency

setClass("ReferenceBias", representation(
			refFraction="array"))

.valid_RefBias_object <- function(object) {

	#check dimensions
	#if(!length(dim(object))%in% c(1,2,3)){
	#	stop("maximum size of dimension is 3")
	#}

	return(TRUE)
}
setValidity("ReferenceBias", .valid_RefBias_object)


#' @rdname ASEset-class
setGeneric("table")
#setGeneric("table", function(x, strand = "*", sortBy="none", ...) {
#    standardGeneric("table")
#})

setMethod("table", signature(... = "ReferenceBias"), function(...) {

	args <- list(...)
	if (length(args) > 1)
	  stop("Only one argument in '...' supported")
	x <- args[[1L]]

	#because the generis of table is rubbish we have to return a list for each strand
	retList <- list()

	for(strand in c("+","-","*")){

		lst <- list()
		y <- x@refFraction[,,strand]
		mean.na.rm <- function(x){mean(x,na.rm=TRUE)}
		lst[["samples"]] <- apply(y,2,mean.na.rm)
		lst[["SNPs"]] <- apply(y,1,mean.na.rm)
		lst[["all"]] <- mean(y,na.rm=TRUE)
		
		retList[[strand]] <- lst

	}
	return(SimpleList(retList))

})	

