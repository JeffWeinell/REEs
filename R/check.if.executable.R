#' @title Check Executable
#' 
#' Checks if the input path is executable.
#' 
#' @param exe.path Character string with input path to check.
#' @return Integer 0 if exe.path is executable, otherwise 1.
#' @export check.if.executable
check.if.executable <- function(exe.path){
	system(paste("command -v",paste0("'",exe.path,"'"),">/dev/null 2>&1 || { echo >&2 ''; exit 1; }"))
}