#'@include initialize-methods.R
NULL

#' DetectedAI class
#' 
#' Object that holds results from a global AI analysis including reference bias
#' estimations and AI detection.
#'
#' The DetectedAI-class contains summaries needed to create an AI-report
#'
#' @name DetectedAI-class
#' @rdname DetectedAI-class
#' @aliases DetectedAI-class DetectedAI DetectedAI-method
#' @docType class
#' @param x ASEset object or list of ASEsets
#' @param TxDb A \code{transcriptDb} object
#' @param ... pass arguments to internal functions
#' @return An object of class DetectedAI containing all data to make report.
#'
#'
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords class ASEset
#' @examples
#'
#' data(ASEset)
#' a <- ASEset
#' gba <- gba(a)
#' 
#' #summary(gba)
#' #write.tables(gba)
#'
#' @exportClass DetectedAI
NULL

#' @rdname DetectedAI-class
#' @exportClass DetectedAI
setClass("DetectedAI", contains = "SummarizedExperiment",
	representation(
		strand = "character"))


