\name{ASEset-gviztrack}
\docType{methods}
\alias{ASEset-gviztrack}
\alias{locationplot}
\alias{locationplot,ASEset-method}


\title{locationplot ASEset objects}
\description{
	plotting ASE effects over a specific genomic region
}
\usage{
\S4method{locationplot}{ASEset}(
	x,
	GR=rowData(x),
	type="fraction",
	strand="+",
	mainVec=vector(),
	verbose=TRUE,
	...)

}
\arguments{
	\item{x}{ an ASEset object.}
	\item{type}{"fraction" or "count"}
	\item{strand}{"+","-". This argument determines which strand is plotted. }
	\item{verbose}{Setting \code{verbose=TRUE} gives details of procedure during function run}
	\item{...}{arguments passed on to barplot function}
}	


\details{
	For information of how to use these tracks in more ways, visit the Gviz package manual.			
} 
 
\author{Jesper R. Gadin}

\seealso{
	\itemize{
		\item The \code{\link{ASEset}} class which the functions can be called up on.
	}
}

\examples{

	data(ASEset)


	GR <- GRanges(seqnames=seqlevels(x),ranges=IRanges(start=min(start(x)),end=max(end(x))),strand="+", genome=genome(x))

	deTrack <- ASEDAnnotationTrack(x, GR=GR, type="fraction",strand="+")
	covTracks <- CoverageDataTrack(x,BamList=r,strand="+") 

	lst <- c(deTrack,covTracks)

	sizes <- c(0.5,rep(0.5/length(covTracks),length(covTracks)))
	plotTracks(lst, from=min(start(x)), to=max(end(x)),sizes=sizes, col.line = NULL, showId = FALSE, main="mainText", cex.main=1, title.width=1, type="histogram")



}
\keyword{locationplot}


