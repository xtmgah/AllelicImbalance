% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/auxillary-methods.R
\docType{methods}
\name{coverageMatrixListFromGAL}
\alias{coverageMatrixListFromGAL}
\alias{coverageMatrixListFromGAL,GAlignmentsList-method}
\title{coverage matrix of GAlignmentsList}
\usage{
coverageMatrixListFromGAL(BamList, ...)

\S4method{coverageMatrixListFromGAL}{GAlignmentsList}(BamList, strand = "*",
  ignore.empty.bam.row = TRUE)
}
\arguments{
\item{BamList}{GAlignmentsList containing reads over the region to calculate
coverage}

\item{...}{arguments to pass on}

\item{strand}{strand has to be '+' or '-'}

\item{ignore.empty.bam.row}{argument not in use atm}
}
\description{
Get coverage per nucleotide for reads covering a region
}
\details{
a convenience function to get the coverage from a list of reads stored in
GAlignmnetsList, and returns by default a list with one matrix, and
information about the genomic start and stop positions.
}
\examples{
r <- reads
seqlevels(r) <- '17'
covMatList <- coverageMatrixListFromGAL(BamList=r, strand='+')
}
\author{
Jesper R. Gadin
}
\keyword{coverage}

