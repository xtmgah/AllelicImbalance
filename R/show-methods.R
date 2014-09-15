#'@include ASEset-class.R
NULL

#' genotype filter methods 
#' 
#' useful genotype filters
#' 
#' These filters are called upon ASEset objects 
#' 
#' @name RefBias-show
#' @rdname RefBias-show
#' @aliases RefBias-show RefBias-show,RefBias-method 
#' @docType methods
#' @param object \code{RefBias} object
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords show
#' @examples
#' 
#' #load example data
#' cat("there is no example data yet for this method")
#' 
#' importFrom(methods, show)
#' exportMethods(show)

setMethod("show","ReferenceBias",
	function(object)
		{

		#Header
		cat("Class:","ReferenceBias\n")
		cat("Contains: RefFraction array\n")
		}
)



