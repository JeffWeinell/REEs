#' Filename from Filepath
#'
#' Returns the name of a file given the full or relative path to the file, optionally dropping the file extension.
#'
#' @param path Path to file
#' @param include.extension If the filename extension should be included in the result (default is TRUE).
#' @return Name of file
#' @export
nameFromPath <- function(path,include.extension=TRUE){
	pathLength <- nchar(path)
	all.locs   <- stringr::str_locate_all(path,"/")[[1]][,1]
	if(all(max(all.locs)==pathLength)){
		start.loc <- (all.locs[(length(all.locs)-1)]+1)
		end.loc   <- (pathLength-1)
		name <- substring(path,first=start.loc,last=end.loc)
	} else {
		start.loc <- (max(all.locs)+1)
		name      <- substring(path,first=start.loc)
	}
	if(!include.extension){
		name <- substr(x=name,start=1,stop=(str_locate_last(name,"\\.")-1))
	}
	name
}
