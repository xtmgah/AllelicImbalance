### =========================================================================
### Helper functions not exported
### =========================================================================

#supposed to merge paths and file irrespective OS and presence of trailing slash
.mergeDirAndFilename <- function(dir,file){
	#check for presence of / in filename in that case remove
	file <- sub("/","",file)
	paste(normalizePath(dir),"/",file, sep="")
}

