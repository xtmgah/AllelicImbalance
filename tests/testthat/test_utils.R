context("test utils functions")

test_that(paste("test trailing slash remover"), {


	dir1 <- system.file('extdata/ERP000101_subset', package='AllelicImbalance')
	file1 <- "ERR009113.bam" 

	dir2 <- paste(dir,"/",sep="")	
	file2 <- paste(file,"/",sep="")	
	file3 <- paste("/", file,sep="")	


	#prepare result string
	res <- paste("/mnt/kelewan/pappewaio/Documents/PHD/repos/AllelicImbalance",
				 "/bioCgit/AllelicImbalance/inst/extdata/ERP000101_subset/",
				 "ERR009113.bam",
				 sep="")

	#check 
    expect_that(res, equals(.mergeDirAndFilename(dir1, file1)))
    expect_that(res, equals(.mergeDirAndFilename(dir1, file2)))
    expect_that(res, equals(.mergeDirAndFilename(dir1, file3)))
    expect_that(res, equals(.mergeDirAndFilename(dir2, file1)))
    expect_that(res, equals(.mergeDirAndFilename(dir2, file2)))
    expect_that(res, equals(.mergeDirAndFilename(dir2, file3)))

})


