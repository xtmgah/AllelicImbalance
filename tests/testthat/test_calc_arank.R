context("calculate rank from ASEset")

test_that(paste("calculate arank from ASEset return.class='matrix' and",
				"return.type='names'"), {

	#prepare testdata
	data(ASEset)	
	x <- ASEset
		
	res <- arank(x, return.type = "names",
		return.class="matrix")

    expect_that(rownames(res), equals(rownames(x)))
    expect_that(colnames(res), equals(c("1","2","3","4")))
    expect_that(matrix(res,ncol=4), 
		equals(matrix(c(rep(TRUE,5),rep(FALSE,7)),ncol=4,byrow=FALSE)))



})



