% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/filter-methods.R
\docType{methods}
\name{genofilters}
\alias{genofilters}
\alias{hetFilt}
\alias{hetFilt,ASEset-method}
\alias{multiAllelicFilt}
\alias{multiAllelicFilt,ASEset-method}
\title{genotype filter methods}
\arguments{
\item{x}{ASEset object}

\item{source}{'genotype' or 'alleleCounts'}

\item{x}{\code{ASEset} object}

\item{strand}{strand to infer from}

\item{inferOver}{'eachSample' or 'allSamples'}

\item{threshold.count.sample}{least amount of counts to try to infer allele}

\item{threshold.frequency}{least fraction to classify (see details)}
}
\description{
useful genotype filters

filter on multiallelic snps
}
\details{
hetFilt returns TRUE if the samples is heterozygote, based on stored genotype information
present in the phase data.

based on the allele counts for all four variants A, T, G and C and returns true
if there is counts enough suggesting a third or more alleles. The sensitivity can
be specified using 'threshold.count.sample' and 'threshold.frequency'.
}
\examples{
#load example data
data(ASEset)
a <- ASEset

genotype(a) <- inferGenotypes(a)
hets <- hetFilt(a)
#load example data
data(ASEset)
a <- ASEset

multiAllelicFilt(a)
}
\author{
Jesper R. Gadin, Lasse Folkersen

Jesper R. Gadin, Lasse Folkersen
}
\keyword{filter}

