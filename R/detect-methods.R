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
	threshold.frequency=0.05, threshold.count.sample=5,
	min.delta.frequency=0.05, max.pvalue=0.05,
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

		mat1 <- matrix(aperm(ac,c(3,2,1))[aperm(array(matrix(x@variants, ncol=length(x@variants),
			 nrow=nrow(x), byrow=TRUE)==arn[,1]
		 ,dim=c(nrow(x), length(x@variants), ncol=ncol(x))),c(2,3,1))
		],ncol=ncol(x),nrow=nrow(x), byrow=TRUE, dimnames=list(rownames(x),colnames(x)))

		mat2 <- matrix(aperm(ac,c(3,2,1))[aperm(array(matrix(x@variants, ncol=length(x@variants),
			 nrow=nrow(x), byrow=TRUE)==arn[,2]
		 ,dim=c(nrow(x), length(x@variants), ncol=ncol(x))),c(2,3,1))
		],ncol=ncol(x),nrow=nrow(x), byrow=TRUE, dimnames=list(rownames(x),colnames(x)))

		mat3 <- matrix(aperm(ac,c(3,2,1))[aperm(array(matrix(x@variants, ncol=length(x@variants),
			 nrow=nrow(x), byrow=TRUE)==arn[,2]
		 ,dim=c(nrow(x), length(x@variants), ncol=ncol(x))),c(2,3,1))
		],ncol=ncol(x),nrow=nrow(x), byrow=TRUE, dimnames=list(rownames(x),colnames(x)))

		mat4 <- matrix(aperm(ac,c(3,2,1))[aperm(array(matrix(x@variants, ncol=length(x@variants),
			 nrow=nrow(x), byrow=TRUE)==arn[,2]
		 ,dim=c(nrow(x), length(x@variants), ncol=ncol(x))),c(2,3,1))
		],ncol=ncol(x),nrow=nrow(x), byrow=TRUE, dimnames=list(rownames(x),colnames(x)))

		allele1 <- mat1[idx] 
		allele2 <- mat2[idx]
		allele3 <- mat3[idx]
		allele4 <- mat4[idx]

		#test the two most expressed alleles
		biasmat1 <- matrix(aperm(mapBias(x, return.class="array"),c(3,2,1))[aperm(array(matrix(
			x@variants, ncol=length(x@variants),
			 nrow=nrow(x), byrow=TRUE)==arn[,1]
		 ,dim=c(nrow(x), length(x@variants), ncol=ncol(x))),c(2,3,1))
		],ncol=ncol(x),nrow=nrow(x), byrow=TRUE, dimnames=list(rownames(x),colnames(x)))

		biasmat2 <- matrix(aperm(mapBias(x, return.class="array"),c(3,2,1))[aperm(array(matrix(
			x@variants, ncol=length(x@variants),
			 nrow=nrow(x), byrow=TRUE)==arn[,1]
		 ,dim=c(nrow(x), length(x@variants), ncol=ncol(x))),c(2,3,1))
		],ncol=ncol(x),nrow=nrow(x), byrow=TRUE, dimnames=list(rownames(x),colnames(x)))

		biasAllele1 <- biasmat1[idx] 
		biasAllele2 <- biasmat2[idx]
		
		#maybe put a check here later that bias1 and bias2 sum to 1

		ml <- mapply(allele1,allele2,biasAllele1,FUN=function(x,y,z){
			binom.test(c(x,y), p = z)[[3]]
		})

		#cerate matrix in same size as fr	
		pv <- fr	
		pv[idx] <- ml

	}else{stop("function.test must be binom.test")}

	if(return.class=="logical"){
		
	}else if(return.class=="vector"){
		stop("return.class vector as option is not available atm")
	}

})

