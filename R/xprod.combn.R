#' All Combinations as Strings
#' 
#' Find all n-wise combinations of elements from n input sets.
#' 
#' @param ... Two or more character vectors (or any vector type coercible to character class). Vectors can be any non-zero length. If input vectors are very long this function may crash. Similarly (and moreso), a large number of input vectors will likely crash this function.
#' @param unique.mode One of types "unordered", "ordered", or "index" (default is "index").
#'  If unique.mode="unordered", output sets (character strings with elements separate by " ") are considered equivalent if they contain the same elements regardless of order, i.e. set (A,B) = set (B,A), and only the former will be returned.
#'  If unique.mode="ordered", output sets are considered equivalent if they contain the same elements in the same order.
#'  If unique.mode="index". Output sets are retained even if they contain the same elements in the same order, which occurs when at least one input vector contains a repeated value.
#' @return A character vector. Each vector element is a character string of the input vector elements separated by a space " ". Definition of unique combination is controlled by the unique.mode parameter.
#' @export xprod.combn
xprod.combn <- function(..., unique.mode="index") {
	list.of.vects <- list(...)
	nloops        <- (length(list.of.vects)-1)
	for(i in 1:nloops){
		if(i==1){
			res.temp <- c(list.of.vects[[1]])
		}
		res.temp <- c(outer(res.temp, c(list.of.vects[[i+1]]), FUN=paste))
	}
	res.temp <- strsplit(res.temp,split=" ")
	res.mat  <- matrix(unlist(res.temp), ncol=length(list.of.vects), byrow=T)
	
	if(unique.mode == "unordered"){
	   res.df  <- data.frame(matrix(apply(X=res.mat,MARGIN=1,FUN=sort),ncol=ncol(res.mat),byrow=T))
	   res.mat <- as.matrix(dplyr::distinct(res.df))
	}
	if(unique.mode=="ordered"){
	 res.mat <- as.matrix(dplyr::distinct(data.frame(res.mat)))
	}
	result    <- apply(X=res.mat,MARGIN=1,FUN=paste,collapse=" ")
	result
}
