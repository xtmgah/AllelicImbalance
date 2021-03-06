% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/GvizTrack-methods.R
\docType{methods}
\name{ASEset-gviztrack}
\alias{ASEDAnnotationTrack}
\alias{ASEDAnnotationTrack,ASEset-method}
\alias{ASEset-gviztrack}
\alias{CoverageDataTrack}
\alias{CoverageDataTrack,ASEset-method}
\title{ASEset-gviztrack ASEset objects}
\usage{
ASEDAnnotationTrack(x, GR = rowRanges(x), type = "fraction", strand = "*",
  trackName = paste("deTrack", type), verbose = TRUE, ...)

\S4method{ASEDAnnotationTrack}{ASEset}(x, GR = rowRanges(x),
  type = "fraction", strand = "*", trackName = paste("deTrack", type),
  verbose = TRUE, ...)

CoverageDataTrack(x, GR = rowRanges(x), BamList = NULL, strand = NULL,
  start = NULL, end = NULL, trackNameVec = NULL, meanCoverage = FALSE,
  verbose = TRUE, ...)
}
\arguments{
\item{x}{an ASEset object.}

\item{GR}{genomic range of plotting}

\item{type}{'fraction' or 'count'}

\item{strand}{'+','-'. This argument determines which strand is plotted.}

\item{trackName}{name of track (ASEDAnnotationTrack)}

\item{verbose}{Setting \code{verbose=TRUE} gives details of procedure during
function run}

\item{...}{arguments passed on to barplot function}

\item{BamList}{GAlignmnentsList object of reads from the same genomic region
as the ASEset}

\item{start}{start position of reads to be plotted}

\item{end}{end position of reads to be plotted}

\item{trackNameVec}{names of tracks (CoverageDataTrack)}

\item{meanCoverage}{mean of coverage over samples (CoverageGataTrack)}
}
\description{
plotting ASE effects over a specific genomic region
}
\details{
For information of how to use these tracks in more ways, visit the Gviz
package manual.
}
\examples{
data(ASEset)
x <- ASEset[,1:2]
r <- reads[1:2]
genome(x) <- 'hg19'
seqlevels(r) <- seqlevels(x)

GR <- GRanges(seqnames=seqlevels(x),
			ranges=IRanges(start=min(start(x)),end=max(end(x))),
			strand='+', genome=genome(x))

deTrack <- ASEDAnnotationTrack(x, GR=GR, type='fraction',strand='+')
covTracks <- CoverageDataTrack(x,BamList=r,strand='+')

lst <- c(deTrack,covTracks)

sizes <- c(0.5,rep(0.5/length(covTracks),length(covTracks)))
#temporarily do not run this function
#plotTracks(lst, from=min(start(x)), to=max(end(x)),
#sizes=sizes, col.line = NULL, showId = FALSE, main='mainText',
#cex.main=1, title.width=1, type='histogram')
}
\author{
Jesper R. Gadin
}
\seealso{
\itemize{ \item The \code{\link{ASEset}} class which the functions
can be called up on.}
}
\keyword{ASEDAnnotationTrack}
\keyword{CoverageDataTrack}

