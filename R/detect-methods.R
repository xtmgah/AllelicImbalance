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
#' return.type 'ref' return only AI when reference allele is more expressed. 'alt' return only
#' AI when alternative allele is more expressed or 'all' for both 'ref' and 'alt' alleles.
#' Reference allele is the one present in the reference genome on the forward strand.
#'
#' min.delta.frequency and function.test will use the value in mapBias(x) as expected value. 
#' 
#' function.test will use the two most expressed alleles for testing. Make therefore sure there
#' are no tri-allelic SNPs or somatic mutations among the SNPs in the ASEset. 
#' 
#' inferGenotype(), set TRUE it should be used with as much samples as possible. If you split up the
#' samples and run detectAI() on each sample separately, please make sure you have inferred 
#' the genotypes in before hand, alternatively used the genotypes detected by another variantCaller
#' or chip-genotypes.
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
#' @param return.type 'ref' ,'alt' or default: 'all'
#' @param return.class class to return (atm only class 'logical')
#' @param min.delta.frequency minimum of frequency difference from 0.5 (or mapbias adjusted value)
#' @param max.pvalue pvalue over this number will be filtered out
#' @param function.test At the moment the only available option is 'binomial.test'
#' @param ... internal arguments 
#' @author Jesper R. Gadin
#' @keywords infer
#' @examples
#' 
#' #load example data
#' data(ASEset)
#' a <- ASEset
#'
#' dai <- detectAI(a)
#' 
#'
#' @rdname detectAI
#' @export
setGeneric("detectAI", function(x, ...){
    standardGeneric("detectAI")
})

#' @rdname detectAI
#' @export
setMethod("detectAI", signature(x = "ASEset"), function(x, 
	return.class = "DetectedAI", return.type="all", strand = "*",
	threshold.frequency=0, threshold.count.sample=1,
	min.delta.frequency=0, max.pvalue=0.05,
	inferGenotype=FALSE,
	random.ref=FALSE,
	function.test="binom.test",
	gc=TRUE) {

	if(inferGenotype){
		genotype(x) <- inferGenotypes(x,threshold.frequency = 0.05, return.allele.allowed ="bi")
	}
	#check for presence of genotype data
    if (!("genotype" %in% names(assays(x)))) {
		stop(paste("genotype matrix is not present as assay in",
				   " ASEset object, see '?inferGenotypes' ",
				   "or set 'inferGenotypes=TRUE'"))
    }

	if(!random.ref){
		if(!("ref" %in% colnames(mcols(x)))){
			stop("column name 'ref' in mcols(x) is required")
		}
	}else{
		mcols(x)[,"ref"] <- randomRef(x)
	}

	#main return class should be DetectedAI, but should also be able to return arrays or lists
	#return.type may be logical or .. 

	#make refFreq array (third dim has length 1 and softest condition)
	fr <- refFraction(x, strand=strand,
			  threshold.count.sample= 1)[,,1]

	#make t.c.s array
	acounts <- alleleCounts(x, strand=strand, return.class="array")
	acounts[array(!hetFilt(x), dim=c(nrow(x), ncol(x), 4))] <- 0
	allele.count.tot <- apply(acounts, c(1,2), sum)
	t.c.s  <- array(allele.count.tot,dim=c(nrow(x),ncol(x),length(threshold.count.sample))) < 
					aperm(array(threshold.count.sample,
						dim=c(length(threshold.count.sample),nrow(x),ncol(x) )),
					c(2,3,1))
	dimnames(t.c.s) <- list(rownames(fr),colnames(fr),paste(threshold.count.sample))

	#make thr.req array (if TRUE it fullfills the condition)
	newDims <- length(threshold.frequency)
	thr.fr.freq <- array(fr,dim=c(nrow(x),ncol(x),newDims))
	thr.freq <- aperm(array(threshold.frequency, dim=c(newDims,nrow(x),ncol(x))),c(2,3,1))

	thr.freq.ret <- !(thr.fr.freq < thr.freq | thr.fr.freq > (1-thr.freq))
	dimnames(thr.freq.ret) <- list(rownames(fr),colnames(fr),paste(threshold.frequency))

	#delta freq array (move to mapbias methods)
	biasmatRef <- matrix(aperm(mapBias(x, return.class="array"),c(3,2,1))[aperm(array(matrix(
		x@variants, ncol=length(x@variants),
		 nrow=nrow(x), byrow=TRUE) == mcols(x)[,"ref"]
	 ,dim=c(nrow(x), length(x@variants), ncol=ncol(x))),c(2,3,1))
	],ncol=ncol(x),nrow=nrow(x), byrow=TRUE, dimnames=list(rownames(x),colnames(x)))

	fr2 <- fr
	fr2[is.na(fr2)] <- biasmatRef[is.na(fr2)] 

	#make fr2 and biasmatRef to array dim 3
	newDims <- length(min.delta.frequency)
	#fr5 <- array(fr2,dim=c(nrow(x),ncol(x),newDims))
	#biasmatRef2 <- array(biasmatRef,dim=c(nrow(x),ncol(x),newDims))

	#to keep1 (fullfills cond min.delta.freq)
	if(any(fr2<biasmatRef)){
		fr3 <- fr2
		fr3[fr2<biasmatRef] <- biasmatRef[fr2<biasmatRef] - fr[fr2<biasmatRef]
		#reset the NA again
		fr3[is.na(fr)] <- NA
		#make 3d array return object 
		fr3 <- array(fr3,dim=c(nrow(x),ncol(x),newDims))
		tf.keep1 <- !(fr3 < min.delta.frequency)
	}else{
		#set all FALSE (NAs will be kept )
		tf.keep1 <- matrix(TRUE,ncol=ncol(fr),nrow=nrow(fr))
		tf.keep1[is.na(fr)] <- NA
		#make 3d array return object 
		tf.keep1 <- array(tf.keep1,dim=c(nrow(x),ncol(x),newDims))
	}

	#to keep2 (fullfills cond min.delta.freq)
	if(any(fr2>biasmatRef)){
		fr3 <- fr2
		fr3[fr2>biasmatRef] <- fr[fr2>biasmatRef] - biasmatRef[fr2>biasmatRef] 
		#reset the NA again
		fr3[is.na(fr)] <- NA
		#make 3d array return object 
		fr3 <- array(fr3,dim=c(nrow(x),ncol(x),newDims))
		tf.keep2 <- !(fr3 < min.delta.frequency)
	}else{
		#set all FALSE (NAs will be kept )
		tf.keep2 <- matrix(TRUE,ncol=ncol(fr),nrow=nrow(fr))
		tf.keep2[is.na(fr)] <- NA
		#make 3d array return object 
		tf.keep2 <- array(tf.keep2,dim=c(nrow(x),ncol(x),newDims))
	}
	
	delta.freq <- tf.keep1 | tf.keep2
	dimnames(delta.freq) <- list(rownames(fr),colnames(fr),paste(min.delta.frequency))

	#select return type (better put this as method on the detectedAI class)
	#if(return.type=="ref"){fr[fr<biasmatRef]<- NaN}
	#if(return.type=="alt"){fr[fr>biasmatRef]<- NaN}

	# p-value array
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


		allele1 <- mat1[idx] 
		allele2 <- mat2[idx]

		#test the two most expressed alleles
		biasmat1 <- matrix(aperm(mapBias(x, return.class="array"),c(3,2,1))[aperm(array(matrix(
			x@variants, ncol=length(x@variants),
			 nrow=nrow(x), byrow=TRUE)==arn[,1]
		 ,dim=c(nrow(x), length(x@variants), ncol=ncol(x))),c(2,3,1))
		],ncol=ncol(x),nrow=nrow(x), byrow=TRUE, dimnames=list(rownames(x),colnames(x)))


		biasAllele1 <- biasmat1[idx] 
		
		#maybe put a check here later that bias1 and bias2 sum to 1

		ml <- mapply(allele1,allele2,biasAllele1,FUN=function(x,y,z){
			binom.test(c(x,y), p = z)[[3]]
		})

		#cerate matrix in same size as fr	
		pv <- fr	
		pv[idx] <- ml

		#make pv array
		pv <- array(pv,dim=c(nrow(fr), ncol(fr),length(max.pvalue)),
					dimnames=list(rownames(fr),colnames(fr),paste(max.pvalue)))

		#filter
		newDims <- length(max.pvalue)
		pv.thr <- aperm(array(max.pvalue, dim=c(newDims,nrow(x),ncol(x))), c(2,3,1))
		pv.thr.ret <- pv < pv.thr
	
	}else{stop("function.test must be binom.test")}

	if(return.class=="DetectedAI"){
		#make DetectedAI object
		DetectedAIFromArray(
			x, 
			strand=strand,
			reference.frequency=fr,
			threshold.frequency=thr.freq.ret,
			threshold.count.sample=t.c.s,
			threshold.delta.frequency=delta.freq,
			threshold.pvalue=pv.thr.ret
		)

	}else if(return.class=="list"){
		stop("return.class list as option is not available atm")
	}else{
		stop(paste("return.class",return.class,"as option is not available"))
	}

})

