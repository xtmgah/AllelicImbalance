###################
#
#  internal List methods
#  
###################

#' general list methods
#' 
#' general list methods that supports flatteing and indexing of lists
#' 
#' This is a collection of supportive list methods. The documentation will
#' be improved before next release.
#' 
#' @name general-list
#' @rdname general-list
#' @aliases general-list list.depth multiUnlist multiUnlist.index multiUlist.index.names
#' @param this a list
#' @param thisdepth, start value for recursion
#' @param lst, the single snp name to plot
#' @param expand.lowest.level logical
#' @param ar array to use for repopulation
#' @param idx.mat array to use for repopulation
#' @param idx.mat.names array to use for repopulation
#' @param ... pass on variables internally
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords list
#' @examples
#' 
#' #some example code here
#' 
#' 
NULL

#' @rdname general-list
#' @export
list.depth <- function(this, thisdepth = 0) {
	  if(!is.list(this)) {
		return(thisdepth)
	  }else {
		return(list.depth(this[[1]], thisdepth = thisdepth+1))
	  }
}
			

#' @rdname general-list
#' @export
multiUnlist <- function(lst, ...){
	 if(!is.list(lst)){
		 return(lst)
	 }else{
		  multiUnlist(do.call(c, unname(lst)))
	 }
}

#' @rdname general-list
#' @export
multiUnlist.index <- function(lst, expand.lowest.level=FALSE){
	
	list.idx.vec <- function(this, i=vector(),vec = vector()) {
		  if(class(this) == "GRanges") {
			  return(c(vec,i,length(this)))
	  } else {
			return(
				unlist(lapply(seq_along(this), 
					   function(y, z, i) { list.idx.vec(y[[i]], i, z)}, y=this, z=c(vec,i))
				))
		  }
	}

	idx.mat <- matrix(list.idx.vec(lst),nrow=list.depth(lst)+2)
					  
	if(expand.lowest.level){
		matrix(inverse.rle(structure(list(lengths = rep(idx.mat[nrow(idx.mat),],nrow(idx.mat)-1), 
								   values  = as.vector(t(idx.mat[-nrow(idx.mat),]))), class = "rle"))
						   ,nrow=nrow(idx.mat)-1, byrow=TRUE)
	}else{
		idx.mat
	}
}

#' @rdname general-list
#' @export
multiUnlist.index.names <- function(lst, expand.lowest.level=FALSE){
	
	list.idx.vec <- function(this, i=vector(),vec = vector(), nms=names(this)) {
		
		  if(class(this) == "GRanges") {
			  return(c(vec, nms[i], length(this)))
		  }else {
			return(
				unlist(
					lapply(seq_along(this), 
					   function(y, i, z, n) { 
						   list.idx.vec(y[[i]], i=i, vec=z, nms=n)
					   },
					   y=this, z=c(vec, nms[i]), n=names(this)
					)
				)
			)
		  }
	}

	list.idx.vec(lst)

	idx.mat <- matrix(list.idx.vec(lst),nrow=list.depth(lst)+2)
					  
	if(expand.lowest.level){
		matrix(inverse.rle(structure(list(lengths = rep(idx.mat[nrow(idx.mat),],nrow(idx.mat)-1), 
								   values  = as.vector(t(idx.mat[-nrow(idx.mat),]))), class = "rle"))
						   ,nrow=nrow(idx.mat)-1, byrow=TRUE)
	}else{
		idx.mat
	}
}

#' @rdname general-list
#' @export
region.list.populate <- function(ar, idx.mat, idx.mat.names ){

	if(!class(idx.mat) == "matrix") {
		l <- lapply(unique(idx.mat), 
			   function(i, a, m){ 
				   a[,,m==i] 
			   },
			   a=ar, m=idx.mat)
		names(l) <- unique(idx.mat.names) 
		l
	}else{
		l <- lapply(unique(idx.mat[1,]), 
		   function(i, a, m, m2){ 

			   region.list.populate(a[,,m[1,]==i], m[-1, m[1, ]==i], m2[-1, m[1, ]==i]) 

		   },
		   a=ar, m=idx.mat, m2=idx.mat.names)
		names(l) <- unique(idx.mat.names[1,]) 
		l
	}
}	

