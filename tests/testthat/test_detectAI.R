context("detect AI from ASEset")

test_that(paste("detectAI from ASEset return.class='DetectedAI',", 
				"multiple samples, multipliple SNPs"), {

	#prepare testdata
	data(ASEset)	
	x <- ASEset


	#make this when the new reference fasta file is in place
	res <- detectAI(x,strand="*", threshold.count.sample=1) 

	#check row and colnames
    expect_that(rownames(res), equals(rownames(x)))
    expect_that(colnames(res), equals(colnames(x)))

	#check strand
    expect_that(rownames(res), equals(rownames(x)))

	#check reference.frequency
    expect_that(as.vector(res)[1:3], equals(c(0.8000000, 0.8076923, NaN)))
    expect_that(round(as.vector(res)[57:60],7), equals(c(NaN, 0.7656250, 0.8163265, NaN)))

	#check threshold.frequency
    expect_that(as.vector(res)[1:3], equals(c(0.8000000, 0.8076923, NaN)))
    expect_that(round(as.vector(res)[57:60],7), equals(c(NaN, 0.7656250, 0.8163265, NaN)))
	
	#check threshold.count.sample
    expect_that(as.vector(res)[1:3], equals(c(0.8000000, 0.8076923, NaN)))
    expect_that(round(as.vector(res)[57:60],7), equals(c(NaN, 0.7656250, 0.8163265, NaN)))

	#check threshold.delta.frequency
    expect_that(as.vector(res)[1:3], equals(c(0.8000000, 0.8076923, NaN)))
    expect_that(round(as.vector(res)[57:60],7), equals(c(NaN, 0.7656250, 0.8163265, NaN)))

	#check threshold.pvalue
    expect_that(as.vector(res)[1:3], equals(c(0.8000000, 0.8076923, NaN)))
    expect_that(round(as.vector(res)[57:60],7), equals(c(NaN, 0.7656250, 0.8163265, NaN)))

})



