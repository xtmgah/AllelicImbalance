
#test_mapBiasReads <- function() {
#
#	#load example data
#	data(ASEset)
#
#	#add plus and minus counts
#	x2 <- ASEset
#	assays(x2)$countsPlus <- assays(x2)$countsNonStranded 
#	assays(x2)$countsMinus <- assays(x2)$countsNonStranded 
#	#path to refgenome
#	"hg19.fa"
#	
#	#add genotypes to the object
#	#assays(ASEset)$genotypes <- array(c("C","G",""),dim=c(3,2),dimnames=c("SNPs","Alleles"))
#
#	#use one sample and one snp
#	x <- ASEset[1,1]
#	x <- renameSeqlevels(x,"chr17")
#
#	#use all samples
#	x <- ASEset
#	x <- renameSeqlevels(x,"chr17")
#
#	#here check something
#	checkEquals(funPOS,readPOS)
#}
#
#
#
