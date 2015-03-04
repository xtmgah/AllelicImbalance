
context("test fraction functionality")

test_that(paste("test fraction default"), {

	#prepare testdata
	data(ASEset)	
	x <- ASEset

	res <- fraction(x)

	#check 
    expect_that(colnames(res), equals(rownames(x)))
    expect_that(rownames(res), equals(colnames(x)))
    expect_that(as.vector(res)[1:3], equals(c(0.8, 1.0, 1.0)))
    expect_that(as.vector(res)[57:60], equals(c(NaN,1,1,1)))

})


test_that(paste("test fraction usePhase=TRUE"), {

	#prepare testdata
	data(ASEset)	
	x <- ASEset

	set.seed(1)
	pha <- matrix(sample(c(1,0), ncol(x)*nrow(x), replace=TRUE), ncol=ncol(x), nrow=nrow(x))

	#store
	phase(x) <- pha
	#access
	res <- phase(x)
	#check 

	res <- fraction(x, usePhase=TRUE)

	#check 
    expect_that(colnames(res), equals(rownames(x)))
    expect_that(rownames(res), equals(colnames(x)))
    expect_that(as.vector(res)[1:5], equals(c(0.8, 0, 0, 1, 0)))
    expect_that(as.vector(res)[57:60], equals(c(NaN,1,1,1)))

})

