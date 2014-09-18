#' detectAI
#' 
#' detection of AllelicImbalance 
#' 
#' threshold.frequency is the least fraction needed to classify as bi tri or
#' quad allelic SNPs. If 'all' then all of bi tri and quad allelic SNPs will use the same
#' threshold. Everything under the treshold will be regarded as noise. 'all' will return 
#' a matrix with snps as rows and uni bi tri and quad will be columns. For this function
#' Anything that will return TRUE for tri-allelicwill also return TRUE for uni and bi-allelic
#' for the same SNP an Sample.
#'
#' @name detectAI
#' @rdname detectAI
#' @aliases detectAI
#' detectAI,ASEset-method 
#' @docType methods
#' @param x ASEset
#' @param strand strand to infer from
#' @param threshold.count.sample least amount of counts to try to infer allele
#' @param threshold.frequency least fraction to classify (see details)
#' @author Jesper R. Gadin
#' @keywords infer
#' @examples
#' 
#' data(ASEset)
#' i <- detectAI(ASEset)
#' 
#' @exportMethod inferAlleles

#' @rdname detectAI
setGeneric("detectAI", function(x, ...){
    standardGeneric("detectAI")
})

setMethod("detectAI", signature(x = "ASEset"), function(x, 
	return.class = "vector", return.type="all", strand = "*",
	threshold.frequency=0.2, threshold.count.sample=5,
	min.delta.frequency=0.1, max.pvalue=0.05,
	function.test="binom.test") {

	fr <- refFraction(x, strand=strand,
			  threshold.count.sample=threshold.count.sample)
	fr[fr < threshold.frequency] <- NaN
	fr[fr > (1-threshold.frequency)] <- NaN

	#survive min delta freq
	fr2 <- fr
	fr2[is.na(fr2)] <- 0.5

	if(any(fr2<0.5)){
		fr[(0.5-fr[fr2<0.5]) < min.delta.frequency] <- NaN
	}
	if(any(fr2>0.5)){
		fr[(fr[fr>0.5]-0.5) < min.delta.frequency] <- NaN
	}

	#select return type
	if(return.type=="higher"){fr[fr>0.5]<- NaN}
	if(return.type=="lower"){fr[fr<0.5]<- NaN}

	#survive stat test max p-value
	if(function.test=="binom.test" & !max.pvalue==1){
		idx <- which(!is.na(fr))
		arn <- arank(x,return.type="names",return.class="matrix") 
		ac <- alleleCounts(x, strand=strand, return.class="array")

	mat1 <- matrix(ac[aperm(array(matrix(x@variants, ncol=length(x@variants),
			 nrow=nrow(x), byrow=TRUE)==arn[,1]
		 ,dim=c(nrow(x), length(x@variants), ncol=ncol(x))),c(1,3,2))
		],ncol=ncol(x),nrow=nrow(x),dimnames=list(rownames(x),colnames(x)))


	mat2 <- matrix(ac[aperm(array(matrix(x@variants, ncol=length(x@variants),
			 nrow=nrow(x), byrow=TRUE)==arn[,2]
		 ,dim=c(nrow(x), length(x@variants), ncol=ncol(x))),c(1,3,2))
		],ncol=ncol(x),nrow=nrow(x),dimnames=list(rownames(x),colnames(x)))




		mat1 <- alleleCounts(x, strand=strand, return.type="array")
		allele1 <- matrix 

		allele2 <- 

		lapply(idx,function(x){
			#test the two most expressed alleles
			ar <- alleleCounts(x, strand=strand, return.type="array")
			binom.test(fr[x],1-fr[x] )
		}
	}else stop("function.test must be binom.test")}
	if(return.class=="logical"){
		
	}else if(return.class=="vector"){
	
	

	}

})

