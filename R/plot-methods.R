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
#' @param hetOverlay logical, if TRUE show nr of het SNPs used to calculate the reference allele frequency mean
#' @param summaryOverSamples 'mean' or 'sum'
#' @param ncol nr of columns for multiplots
#' @param ... pass on variables internally
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords list
#' @examples
#' 
#' #some example code here
#' #generate example
#' data(ASEset)
#' a <- ASEset
#'dai <- detectAI(a, 
#' 			threshold.count.sample=1:50,
#' 			threshold.frequency=seq(0,0.5,by=0.01),
#' 			threshold.delta.frequency=seq(0,0.5,by=0.01),
#' 			threshold.pvalue=rev(seq(0.001,0.05, by=0.005))
#' )
#' 
#' frequency_vs_threshold_variable_plot(dai)
#' detectedAI_vs_threshold_variable_plot(dai)
#' detectedAI_vs_threshold_variable_multigraph_plot(dai)
#' frequency_vs_threshold_variable_multigraph_plot(dai)
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
	function(x,var="threshold.count.sample", hetOverlay=TRUE, smoothscatter=FALSE){

	mat <- frequency_vs_threshold_variable_summary(x, var)
	#mean.het <- apply(mat,2,function(x){sum(!is.na(x))})
	mean.het <- apply(usedSNPs_vs_threshold_variable_summary(x, var), 2, mean)

	#mean for every variable
	vec <- apply(mat,2,mean,na.rm=TRUE)
	labels <- slot(x, paste(var,".names", sep=""))
	labels <- as.numeric(labels)
	if(var=="threshold.pvalue"){
		labels <- -log10(labels)
	}

	data <- data.frame(refFreq=vec, variable=labels)
	data2 <- data.frame(meanHetSNPs=mean.het, variable=labels)

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

		a <- xyplot(refFreq ~ variable,data,grid=TRUE, lwd=2, pch=16,type=c("p","l"),
				xlab=var,ylab="Reference Frequency Mean", col="blue"
		)

		if(hetOverlay==TRUE){

			library(latticeExtra)
			b <- xyplot(meanHetSNPs ~ variable, data2,grid=TRUE, lwd=2, pch=16,type=c("p","l"),
				xlab="",ylab="Nr. of Het SNPs Left", add=TRUE, col="green")

			doubleYScale(a, b, style1 = 0, style2 = 3, add.ylab2 = TRUE)
			    
		}
		
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
	function(x, var="threshold.count.sample", summaryOverSamples="sum", hetOverlay=TRUE, smoothscatter=FALSE){

	mat <- detectedAI_vs_threshold_variable_summary(x, var)
	mean.het <- apply(usedSNPs_vs_threshold_variable_summary(x, var), 2, mean)

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
	data2 <- data.frame(meanHetSNPs=mean.het, variable=labels)

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

		a <- xyplot(detectedAI ~ variable, data,grid=TRUE, lwd=2, pch=16,type=c("p","l"),
				xlab=var,ylab="Detected AI"
		)

		if(hetOverlay){

			library(latticeExtra)
			b <- xyplot(meanHetSNPs ~ variable, data2,grid=TRUE, lwd=2, pch=16,type=c("h"),
				xlab=var,ylab="", add=TRUE, col="green")

			doubleYScale(a, b, style1 = 0, style2 = 3, add.ylab2 = TRUE)
			    
		}
	}
})

#' @rdname DetectedAI-plot
#' @export
setGeneric("reference_frequency_density_vs_threshold_variable_plot", function(x, ... 
	){
    standardGeneric("reference_frequency_density_vs_threshold_variable_plot")
})

#' @rdname DetectedAI-plot
#' @export
setMethod("reference_frequency_density_vs_threshold_variable_plot", signature(x = "DetectedAI"),
	function(x, var="threshold.count.sample"){

	fr <- frequency_vs_threshold_variable_summary(x, var, return.class="array")

	labels <- slot(x, paste(var,".names", sep=""))
	labels <- as.numeric(labels)
	if(var=="threshold.pvalue"){
		labels <- -log10(labels)
	}

	labels <- aperm(array(labels, dim=dim(fr)[c(3,2,1)]), c(3,2,1))

	#calculate density
	nbins <- dim(fr)[3] 
	x.bin <- seq(min(labels), max(labels), length=nbins)
	y.bin <- seq(0, 1, length=nbins)

	mat.i <- matrix(NA, ncol=2, nrow=length(fr))
	mat.i[,1] <- findInterval(labels, x.bin)
	mat.i[,2] <- findInterval(fr, y.bin)
	mat.i <- mat.i[!is.na(mat.i[,2]),]

	df.long <- as.data.frame(table(mat.i[,2],mat.i[,1]))
	colnames(df.long) <- c("refFreq","threshold", "z")
	
	#data to long format
	#df.long <- data.frame(refFreq=as.vector(fr))
	#df.long[["threshold"]] <- as.vector()


	levelplot(z ~ refFreq*threshold, data = df.long,
			    xlab = var, ylab = "frequency",
				  main = "density of ref allele frequency",
				    col.regions = rev(heat.colors(100))
				)
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

	pl <- lapply(graphs, function(y,x ){
			detectedAI_vs_threshold_variable_plot(x, var=y, ...)
		}, x
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
	function(x, ncol=2, ...){

	require(gridExtra)
	
	graphs <- c("threshold.frequency", "threshold.count.sample",
	"threshold.delta.frequency", "threshold.pvalue" ) 

	pl <- lapply(graphs, function(y,x){
			frequency_vs_threshold_variable_plot(x, var=y)
		}, x
	)
	do.call(grid.arrange, c(pl, ncol=2))

})

#' @rdname DetectedAI-plot
#' @export
setGeneric("reference_frequency_density_vs_threshold_variable_multigraph_plot", function(x, ... 
	){
    standardGeneric("reference_frequency_density_vs_threshold_variable_multigraph_plot")
})

#' @rdname DetectedAI-plot
#' @export
setMethod("reference_frequency_density_vs_threshold_variable_multigraph_plot", signature(x = "DetectedAI"),
	function(x, ncol=2, ...){

	require(gridExtra)
	
	graphs <- c("threshold.frequency", "threshold.count.sample",
	"threshold.delta.frequency", "threshold.pvalue" ) 

	pl <- lapply(graphs, function(y,x){
			reference_frequency_density_vs_threshold_variable_plot(x, var=y)
		}, x
	)
	do.call(grid.arrange, c(pl, ncol=2))

})


