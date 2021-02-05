#' Build Path of Directory
#' 
#' This function returns the full path to each directory along the input path. The input path does not need to exist, but if it does, it must lead to a directory, not a file.
#' 
#' @param path Character string with path to a directory
#' @return Character vector containing paths to each directory on the input path.
#' @export BuildPath.dir
BuildPath.dir <- function(path){
    if(file.exists(path)==T & dir.exists(path)==F){
        stop("input path is to a file, but should be a directory")
    } else {
        subpath <- path
        result  <- subpath
        while(!(dirname(subpath) %in% c("/",""))){
            subpath <- dirname(subpath)
            result  <- c(subpath,result)
        }
    }
    result
}
