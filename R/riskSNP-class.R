#'@include initialize-methods.R
NULL

#' riskVariant class
#' 
#' Object that holds results from AI detection.
#'
#' The riskVariant-class contains 
#'
#' @name riskVariant-class
#' @rdname riskVariant-class
#' @aliases riskVariant-class riskVariant riskVariant-method
#' @docType class
#' @param x riskVariant object or list of ASEsets
#' @param return.class type of class returned eg. "list or ""array".
#' @param ... pass arguments to internal functions
#'
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords class ASEset
#' @examples
#'
#' #some code
#'
#' @exportClass riskVariant
NULL

#' @rdname riskVariant-class
#' @exportClass riskVariant
setClass("riskVariant", contains = "SummarizedExperiment",
	representation(
		meta = "list"
	)
)


