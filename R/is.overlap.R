#' @title Test if Two Ranges Overlap
#'
#' This function takes as input two integer vectors. The min and max values of each vector represent the lower and upper limits of each range, respectively. Returns TRUE if input ranges overlap, otherwise FALSE.
#'
#' @param r1 Integer vector with min and max values representing limits of first range
#' @param r2 Integer vector with min and max values representing limits of second range
#' @return Logical indicating whether or not r1 and r2 overlap
#' @export is.overlap
is.overlap <- function(r1,r2){
	 min(r1) <= max(r2) && min(r2) <= max(r1)
}

#' @title Return range of merged intervals if they overlap
#'
#' This function takes as input two integer vectors. The min and max values of each vector represent the lower and upper limits of each range, respectively. Returns the range of the merged interval if r1 and r2 overrlap, otherwise NA.
#'
#' @param r1 Integer vector with min and max values representing limits of first range
#' @param r2 Integer vector with min and max values representing limits of second range
#' @param ifFALSE What to return when r1 and r2 are nonoverlapping. Default NA; ifFALSE=1 returns r1; ifFALSE=2 returns r2; ifFALSE=12 returns a length 2 list with r1 and r2
#' @return range of merged r1 and r2 if is.overlap(r1,r2) is true, otherwise NA.
#' @export mergedRange
mergedRange <- function(r1,r2,ifFALSE=NA){
	if(is.overlap(r1,r2)) {
		range(c(r1,r2))
	} else {
		if(all(is.na(ifFALSE))) {
			return(NA)
		}
		if(all(ifFALSE==1)) {
			return(r1)
		}
		if(all(ifFALSE==2)) {
			return(r2)
		}
		if(all(ifFALSE==12)) {
			return(list(r1,r2))
		}

	}
}

#' @title Clusters intervals 
#'
#' Merges overlapping intervals and returns a data frame with names of merged intervals and start and end of merged intervals.
#'
#' @param data two column data frame with limits of each interval held on a separate row
#' @return Three column data frame with the set of nonoverlapping intervals formed after merging overlapping intervals of data; the first column holds the comma-separated rownames of data that were merged; columns two and three hold the start and end positions of the post-merge interval
#' @export clusterIntervals
clusterIntervals <- function(data){
	# sorts values in rows such that column 1 values are less those in column 2
	df1 <- rbind(data[c(data[,1] < data[,2]),c(1,2)],data[c(data[,1] > data[,2]),c(2,1)])
	# sorts rows by first column and then by successive columns to break ties
	df2 <- df1[do.call(order, as.list(df1)),]
	res=list(); length(res) <- nrow(df2)
	names(res)[1] <- rownames(df2)[1]
	res[[1]] <- unname(unlist(df2[1,]))
	for(i in 2:nrow(df2)){
		res[[i]] <- mergedRange(res[[(i-1)]],df2[i,],ifFALSE=2)
		if(is.overlap(res[[(i-1)]],df2[i,])) {
			names(res)[[i]] <- paste(rownames(df2)[i],names(res)[[i-1]],sep=",")
			res[[(i-1)]] <- NA
		} else {
			names(res)[[i]] <- rownames(df2)[i]
		}
	}
	intervalSets        <- do.call(rbind,res[!sapply(res,function(x){all(is.na(x))})])
	result.df           <- as.data.frame(unname(as.matrix(cbind(rownames(intervalSets),intervalSets))))
	colnames(result.df) <- c("merged","start","end")
	result.df
}

#' @title Test if subset
#' 
#' Test if a range of numbers is a subset of another range of numbers.
#' 
#' @param r1 Integer vector with min and max values representing limits of first range
#' @param r2 Integer vector with min and max values representing limits of second range
#' @return Logical indicating whether or not r1 is a non-equivalent subset of r2.
#' @export is.subset
is.subset <- function(r1,r2){
	test1  <- (min(r1) >= min(r2) && max(r1) <= max(r2))
	test2  <- any(c(min(r1),max(r1)) != c(min(r2),max(r2)))
	result <- all(test1,test2)
	result
}
