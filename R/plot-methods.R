#generate example
#data(ASEset)
#object <- detectAI(ASEset, 
#			threshold.count.sample=1:50,
#			threshold.frequency=seq(0,0.5,by=0.01),
#			threshold.delta.frequency=seq(0,0.5,by=0.01),
#			threshold.pvalue=rev(seq(0.001,0.05, by=0.005))
#)

frequency_vs_threshold_variable_plot <- 
	function(object,var="threshold.count.sample",smoothscatter=FALSE){

	mat <- frequency_vs_threshold_variable_summary(object, var)

	#mean for every variable
	vec <- apply(mat,2,mean,na.rm=TRUE)
	labels <- slot(object, paste(var,".names", sep=""))
	labels <- as.numeric(labels)
	if(var=="threshold.pvalue"){
		labels <- -log10(labels)
	}

	data <- data.frame(refFreq=vec, variable=labels)

	if(smoothscatter){
		#include points for all samples
		data2 <- data.frame(matrix(c(mat,rep(labels,nrow(mat))), 
									nrow=ncol(mat)*nrow(mat), ncol=2, byrow=FALSE))

		colnames(data2) <- c("refFreq","variable")
		groups <- var
		
		data$groups <- var

		#test
		p2 <- xyplot(refFreq ~ variable, data=data2, grid=T, panel=function(x,y,...){
			panel.smoothScatter(x,y,...)
			panel.linejoin(x,y, fun=function(y){mean(y,na.rm=TRUE)}, horizontal=FALSE, col="green", ...)
		})
		
	}else{

		xyplot(refFreq ~ variable,data,grid=TRUE, lwd=2, pch=16,type=c("p","l"),
				xlab=var,ylab="Reference Frequency",
		)
	}
}

