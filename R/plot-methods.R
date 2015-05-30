#' DetectedAI plot
#' 
#' plot functions for the DetectedAI-class
#' 
#' plot helper functions. The documentation will
#' be improved before next release.
#' 
#' @name DetectedAI-plot
#' @rdname DetectedAI-plot
#' @aliases frequency_vs_threshold_variable_plot frequency_vs_threshold_variable_plot,DetectedAI-class
#' @param x detectedAI object
#' @param var string, see details for available options
#' @param smoothscatter boolean, smoothscatter over the means
#' @param ... pass on variables internally
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords list
#' @examples
#' 
#' #some example code here
#' generate example
#' data(ASEset)
#' a <- ASEset
#' object <- detectAI(a, 
#' 			threshold.count.sample=1:50,
#' 			threshold.frequency=seq(0,0.5,by=0.01),
#' 			threshold.delta.frequency=seq(0,0.5,by=0.01),
#' 			threshold.pvalue=rev(seq(0.001,0.05, by=0.005))
#' )
#' 
#' frequency_vs_threshold_variable_plot(object)
#' detectedAI_vs_threshold_variable_plot(object)
#' detectedAI_vs_threshold_variable_multigraph_plot(object)
#' frequency_vs_threshold_variable_multigraph_plot(object)
#'
NULL

#' @rdname DetectedAI-plot
#' @export
setGeneric("frequency_vs_threshold_variable_plot", function(x, ... 
	){
    standardGeneric("frequency_vs_threshold_variable_plot")
})

#' @rdname DetectedAI-plot
#' @export
setMethod("frequency_vs_threshold_variable_plot", signature(x = "DetectedAI"),
	function(x,var="threshold.count.sample",smoothscatter=FALSE){

	mat <- frequency_vs_threshold_variable_summary(x, var)

	#mean for every variable
	vec <- apply(mat,2,mean,na.rm=TRUE)
	labels <- slot(x, paste(var,".names", sep=""))
	labels <- as.numeric(labels)
	if(var=="threshold.pvalue"){
		labels <- -log10(labels)
	}

	data <- data.frame(refFreq=vec, variable=labels)

	if(smoothscatter){
		#include points for all samples
		data2 <- data.frame(refFreq=as.vector(t(mat)))
		data2$variable <- rep(labels,nrow(mat))			

		colnames(data2) <- c("refFreq","variable")

		xyplot(refFreq ~ variable, data=data2, grid=T, panel=function(x,y, ...){
			panel.smoothScatter(x, y, ...)
			panel.linejoin(x, y, fun=function(y, ...){
								   mean(y,  na.rm=TRUE)
						   },
						   horizontal=FALSE, lwd=2, col="red", ...)
		})
		
	}else{

		xyplot(refFreq ~ variable,data,grid=TRUE, lwd=2, pch=16,type=c("p","l"),
				xlab=var,ylab="Reference Frequency",
		)
	}
})

#' @rdname DetectedAI-plot
#' @export
setGeneric("detectedAI_vs_threshold_variable_plot", function(x, ... 
	){
    standardGeneric("detectedAI_vs_threshold_variable_plot")
})

#' @rdname DetectedAI-plot
#' @export
setMethod("detectedAI_vs_threshold_variable_plot", signature(x = "DetectedAI"),
	function(x, var="threshold.count.sample", summaryOverSamples="sum", smoothscatter=FALSE){

	mat <- detectedAI_vs_threshold_variable_summary(x, var)

	summaryOverSamplesFun <- function(mat,opt=summaryOverSamples){
		if(opt=="sum"){
			return(apply(mat,2,sum,na.rm=TRUE))
		}else if(opt=="mean"){
			return(apply(mat,2,mean,na.rm=TRUE))
		}
	}
	vec <- summaryOverSamplesFun(mat, summaryOverSamples)

	labels <- slot(x, paste(var,".names", sep=""))
	labels <- as.numeric(labels)

	if(var=="threshold.pvalue"){
		labels <- -log10(labels)
	}

	data <- data.frame(detectedAI=vec, variable=labels)

	if(smoothscatter){
		#include points for all samples
		data2 <- data.frame(detectedAI=as.vector(t(mat)))
		data2$variable <- rep(labels, nrow(mat))			

		colnames(data2) <- c("detectedAI","variable")

		xyplot(detectedAI ~ variable, data=data2, grid=T, panel=function(x,y, ...){
			panel.smoothScatter(x, y, ...)
			panel.linejoin(x, y, fun=function(y, ...){
								   mean(y,  na.rm=TRUE)
						   },
						   horizontal=FALSE, lwd=2, col="red", ...)
		})

	}else{

		xyplot(detectedAI ~ variable, data,grid=TRUE, lwd=2, pch=16,type=c("p","l"),
				xlab=var,ylab="Detected AI"
		)
	}
})

#' @rdname DetectedAI-plot
#' @export
setGeneric("detectedAI_vs_threshold_variable_multigraph_plot", function(x, ... 
	){
    standardGeneric("detectedAI_vs_threshold_variable_multigraph_plot")
})

#' @rdname DetectedAI-plot
#' @export
setMethod("detectedAI_vs_threshold_variable_multigraph_plot", signature(x = "DetectedAI"),
	function(x, ncol=2, ...){
	
	require(gridExtra)

	graphs <- c("threshold.frequency", "threshold.count.sample",
	"threshold.delta.frequency", "threshold.pvalue" ) 

	pl <- lapply(graphs, function(x,object){
			detectedAI_vs_threshold_variable_plot(object, var=x, ...)
		}, object
	)
	do.call(grid.arrange, c(pl, ncol=2))

})

#' @rdname DetectedAI-plot
#' @export
setGeneric("frequency_vs_threshold_variable_multigraph_plot", function(x, ... 
	){
    standardGeneric("frequency_vs_threshold_variable_multigraph_plot")
})

#' @rdname DetectedAI-plot
#' @export
setMethod("frequency_vs_threshold_variable_multigraph_plot", signature(x = "DetectedAI"),
	function(x, ncol=2){

	require(gridExtra)
	
	graphs <- c("threshold.frequency", "threshold.count.sample",
	"threshold.delta.frequency", "threshold.pvalue" ) 

	pl <- lapply(graphs, function(x,object){
			frequency_vs_threshold_variable_plot(object, var=x)
		}, object
	)
	do.call(grid.arrange, c(pl, ncol=2))

})


