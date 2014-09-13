#' mapbias effect
#' 
#' measures the effect of mapbias over heterozygous SNPs
#' 
#' The distribution of fractions for the ref allele and the 
#' alternative are compared
#' 
#' @name mapbias
#' @rdname mapbias
#' @aliases mapbias,ASEset-method
#' @docType methods
#' @param x \code{ASEset} object
#' @param strand strand option
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords mapbias
#' @examples
#' 
#' #load example data
#' data(ASEset)
#' a <- ASEset
#'
#' #this example data contains to few SNPs to actually 
#' #measure the effect reliably, but as example it serves
#' #the purpose.
#'	
#' #prepare ASEset
#' genotype(a) <- inferGenotypes(a)
#' a <- refAllele(a,
#'		fasta=system.file('extdata/hg19.chr17.fa', 
#'		package='AllelicImbalance'))	
#'
#' m <- mapbias(a)
#' 
#' @exportMethod mapbiasEffect
NULL

#' @rdname ASEset-class
setGeneric("mapbiasEffect", function(x, strand="*"){
    standardGeneric("mapbiasEffect")
})

setMethod("mapbiasEffect", signature(x = "ASEset"), function(x, strand="*"){
	
	#check for presence of genotype data
	if(!("genotype" %in% assays(x))){
		stop(paste("genotype assay does not exist, add genotypes or infer them",
			 "eg. '?inferGenotypes'"),sep="")
	}
	#check for presence of reference allele
	if(!("ref" %in% colnames(mcols(x)))){
		stop("column name 'ref' in mcols(x) is required")
	}
	



})


#' Reference allele
#' 
#' Extract the allele based on SNP location from the reference fasta file
#' 
#' The alleles will be placed in the rowData() meta column 'ref'
#' 
#' 
#' @name refAllele 
#' @rdname refAllele
#' @aliases refAllele,ASEset-method
#' @docType methods
#' @param x \code{ASEset} object
#' @param strand strand option
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords reference mapbias
#' @examples
#' 
#' #load example data
#' data(ASEset)
#'
#' fasta <- system.file('extdata/hg19.chr17.fa', package='AllelicImbalance')
#' refAllele(ASEset,fasta=fasta)
#' a <- refAllele(ASEset,fasta=fasta) 
#'
#' @exportMethod refAllele
NULL

#' @rdname refAllele
setGeneric("refAllele", function(x, fasta ){
    standardGeneric("refAllele")
})

setMethod("refAllele", signature(x = "ASEset"), function(x, fasta){

	#does the fasta file exist?
	if(!file.exists(fasta)){
		stop("fasta file doesnt exist")
	}

	#make FaFile
	fl <- FaFile(fasta)

	#check if the index file is present otherwise tell the user to use the
	#indexFa(FaFile("pathToReference")) command
	if(!file.exists(paste(fasta,".fai",sep=""))){
		cat("could not find index file\n")
		cat("creates a new index file")
		indexFa(FaFile(fasta))
		cat("finished creating new index file")
	}
	#IMPORTANT! The  index command only needs to be executed once
	#indexFa(fl) #creates a new file as index in the same directory but
	#with extension *.fai 

	#open,scan,close file
	open(fl)
	#check if seqleveles are present
	fa.info <- scanFaIndex(fl)
	if(!all(seqlevels(x) %in% seqlevels(fa.info))){
		stop("seqlevels in object x are not in fasta index file")
	}
	
	ref <- scanFa(fl, param=rowData(x))
	close(fl)
	
	mcols(x)[["ref"]] <- as.vector(ref)
	x
})



