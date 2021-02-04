#' Find First N Matches
#' 
#' This function returns the indices of the first n matches of each value in a vector x present in a vector y.
#'
#' @param x Vector of values to look for. x should not contain NA values.
#' @param y Vector to search within.
#' @param n Number of matches to return.
#' @return when simplify = FALSE, this function returns a list of numeric vectors. The ith vector contains the locations of (at most) the first n matches of x[i] in y.
#'  If simplify = TRUE, a single numeric vector is returned that corresponds to concatenating the list of vectors that would hav ebeen returned if simplify=F.
#' @export nmatch
nmatch <- function(x,y,n){
	queries  <- x
	subjects <- y
	#subjects.mat <- cbind(1:length(subjects),subjects)
	result.temp <- NULL
	for(i in 1:n){
		#if(i==1){
			## Find ith match of each x in y
			matches.temp         <- match(queries,subjects)
			## Remove any NA values generated for x values with fewer than i matches in y.
			matches.temp         <- matches.temp[!is.na(matches.temp)]
			## Break out of the loop if matches.temp is empty because all elements in x have fewer than i matches in y
			if(length(matches.temp)==0){
				break
			}
			#stats.table3.temp   <- stats.table2[to.keep2.temp,]
			# Updates subjects to block out the ith matches, so that the match function can find the i+1 matches in the next loop
			subjects[matches.temp] <- NA
			# Add matches.temp to the set of earlier matches
			result.temp            <- c(result.temp,matches.temp)
			#remainder.temp       <- subjects[-matches.temp,]
		#} else {
		#	matches.temp         <- match(queries,remainder.temp)
		#	#stats.table3.temp   <- cbind(stats.table3.temp,remainder.table.temp[to.keep2.temp,])
		#	remainder.table.temp <- remainder.table.temp[-to.keep2.temp,]
		#}
	}
	result <- result.temp
	result
}

#' Find First N Matches
#' 
#' This function returns the indices of the first n matches of each value in a vector x present in a vector y.
#' This does the nearly the same thing as the function nmatch, but uses base::lapply on base::grep instead of looping base::match
#' Unlike nmatch, nmatch.v2 has the option of returning a list of vectors holding match indices, rather than a concatenated list of vectors holding match indices.
#' This method is probably slower than nmatch, but I havent verified this.
#'
#' @param x Vector of values to look for.
#' @param y Vector to search within.
#' @param n Number of matches to return.
#' @param simplify Should the result be coerced to single numeric vector. Default is FALSE.
#' @return when simplify = FALSE, this function returns a list of numeric vectors. The ith vector contains the locations of (at most) the first n matches of x[i] in y.
#'  If simplify = TRUE, a single numeric vector is returned that corresponds to concatenating the list of vectors that would hav ebeen returned if simplify=F.
#' @export nmatch.v2
nmatch.v2 <- function(x,y,n,simplify=F){
	#queries     <- c("A","B","C")
	queries     <- x
	queries2    <- paste0("^",queries,"$")
	#subjects    <- c("AB","A","A","A","B","C","BC","CC")
	subjects <- y
	#all.matches <- lapply(queries2,FUN=grep,x=subjects)
	all.matches <- lapply(queries2,FUN=function(pat,sub,num){res.temp=grep(pat,sub);res=res.temp[1:min(length(res.temp),n)]},sub=subjects,num=n)
	if(simplify){
		result <- unlist(all.matches)
	} else {
		result <- all.matches
	}
	result
}




