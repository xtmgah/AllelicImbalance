\name{inferVariants-utilities}
\alias{inferVariants-utilities}
\alias{inferGenotypes}
\alias{inferSnp}


\title{infer snps and genotypes}
\description{
    inference method that uses count data from the two alleles to guess genotype  
}
\usage{
    inferGenotypes(x, min.allele.fraction=0.05, min.infer.count=10,return.type="list", biallelic.only=TRUE)
    inferSnp(x,strand="nonStranded",biallelic.only=TRUE)
}

\arguments{
  \item{x}{ASEset}
  \item{min.allele.fraction}{the least fraction of reads that must come from each allele}
  \item{min.infer.count}{the least amount of reads that must come from both alleles in total}
  \item{return.type}{return a 'list' or 'array'}
  \item{biallelic.only}{If TRUE we only infer biallelic alleles}
  \item{strand}{the strand to use counts from in the inference}
}

\details{
    InferGenotypes' is a very simple variant caller. But might be sufficient if the read depth is high.	We recommend strongly to use for example samtools/vcftools in the general case.
}

\value{
    \code{inferGenotypes} returns inferred genotypes
    \code{inferSnp} returns a matrix of biallelic SNPs
}

\author{Jesper R. Gadin}
\examples{
    data(ASEset)
    mat <- inferSnp(ASEset)
    lst <- inferGenotypes(ASEset)
}

\keyword{ infer }


