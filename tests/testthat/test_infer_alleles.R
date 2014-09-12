context("inferAlleles ASEset")

test_that("correct inference of alleles from ASEset return.type='all'", {

	#prepare testdata
	data(ASEset)	
	x <- ASEset[1:3,1:3]
		
	res <- inferAlleles(x, return.type = "all",return.type.criteria = 0.05 )

    expect_that(rownames(res), equals(rownames(x)))
    expect_that(colnames(res), equals(c("bi","tri","quad")))
    expect_that(matrix(res,ncol=3), 
		equals(matrix(c(TRUE,TRUE,rep(FALSE,7)),ncol=3,byrow=FALSE)))



})



