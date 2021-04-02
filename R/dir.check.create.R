#### Neither of these is needed anymore, because you can use the base function dir.create(dir.temp,recursive=T)


#' Create directory unless it already exists
#' 
#' Checks if a directory exists, and if not, creates it. Parent directories are also created if they do not already exist.
#' Avoids calling 'mkdir -p dir'
#'
#' @param directory name or path to the directory to check/create
#' @return A new directory unless it already existed
#' @export dir.check.create.v2
dir.check.create.v2 <- function(directory){
	dir.parts <- unlist(strsplit(directory,"/"))
	for(i in 1:(length(dir.parts)-1)){
		dir.temp <- paste(dir.parts[1:(i+1)],collapse="/")
		if (!file.exists(dir.temp)){
			dir.create(dir.temp)
		}
	}
}

#' Create directory unless it already exists
#' 
#' Checks if a directory exists, and if not, creates it. Parent directories are also created if they do not already exist.
#' 
#' 
#' @param directory name or path to the directory to check and create
#' @return A new directory unless it already existed
#' @export dir.check.create
dir.check.create <- function(directory){
	dirs.exist <- sapply(X=BuildPath.dir(directory),FUN=dir.exists)
	if(any(!dirs.exist)){
		result <- sapply(names(dirs.exist[which(!dirs.exist)]),FUN=dir.create)
	}
}
