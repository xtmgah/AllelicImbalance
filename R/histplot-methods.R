#' histogram plots
#' 
#' uses base graphics hist plot
#' 
#' The histogram will show the density over frequencies for each sample
#' 
#' @name histplot
#' @rdname histplot
#' @aliases hist hist,ReferenceBias-method 
#' @docType methods
#' @param x \code{ReferenceBias} object
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords plot hist
#' @examples
#' 
#' #load example data
#' 
#' data("ReferenceBias")
#' hist(ReferenceBias)
#' 
NULL

#' @rdname histplot-class
setMethod("hist", signature(x = "ReferenceBias"), function(x, strand="*", ...){
	hi <- hist(frequency(x,strand=strand),breaks = 40, freq=TRUE, ...)

	#add red line for 0.5
	abline(v=0.5, col="red")

	invisible(hi)
})

#' @rdname histplot-class
setMethod("hist", signature(x = "ASEset"), function(x, strand="*", type="mean", log=1, ...){

		if(log==1){
			counts <- countsPerSnp(x, strand=strand, return.type=type, return.class="vector")
		}else{
			counts <- log(countsPerSnp(x, strand=strand, return.type=type, return.class="vector"),log)
		}

		hi <- hist(counts,breaks = 40, freq=TRUE, ...)
		invisible(hi)
})


