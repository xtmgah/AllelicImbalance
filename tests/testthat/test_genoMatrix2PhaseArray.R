
context("internal functions for genoMatrix2PhaseArray")

test_that(paste("checking .splitGenotypeMatrix"), {

	#####################
	#prepare testdata
	#####################
	#A matrix with both alleles present in all samples (simplest case)
	#Dim: 3 SNPs 4 Samples
	a1 <- c("A", "A", "T", "T",
			"G", "G", "C", "C",
			"C", "G", "C", "G")
	a2 <- c("T", "A", "A", "T",
			"C", "C", "C", "G",
			"G", "G", "C", "C")
	mat1 <- matrix(paste(a1,"/",a2,sep=""), nrow=3, ncol=4, byrow=TRUE)
	exp1 <- aperm(array(c(a1, a2), c(4,3,2)),c(2,1,3))
			
	#run tests
	res1 <- .splitGenotypeMatrix(mat1)
	
	#test equality
    expect_that(exp1, equals(res1))

})

test_that(paste("checking .splitGenotypeRank"), {

	#####################
	#prepare testdata
	#####################
	#A matrix with both alleles present in all samples (simplest case)
	#Dim: 3 SNPs 4 Samples
	a1 <- c("A", "A", "T", "T",
			"G", "G", "C", "C",
			"C", "G", "C", "G")
	a2 <- c("T", "A", "A", "T",
			"C", "C", "C", "G",
			"G", "G", "C", "C")
	mat1 <- matrix(paste(a1,"/",a2,sep=""), nrow=3, ncol=4, byrow=TRUE)
	exp1 <- c("A", "T", "C", "G",
			"C", "G", "C", "G",
			"G", "G", "C", "C")
	#run tests
	res1 <- .splitGenotypeMatrix(mat1, return.ranknames=FALSE)
	
	#test equality
    expect_that(exp1, equals(res1))


    expect_that(colnames(plus), equals(colnames(arp)))
    expect_that(rownames(minus), equals(rownames(arm)))
    expect_that(colnames(minus), equals(colnames(arm)))
	expect_identical(plus, arp)
	expect_identical(minus, arm)
})



