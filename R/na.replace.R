#' Replace missing values
#' 
#' Replace missing values. This function is identical to the na.replace function in the gtools package.
#' Convenience function that is the same as x[is.na(x)] <- replace
#' 
#' @param x Vector possibly containing missing (NA) values.
#' @param replace Either a scalar replacement value, or a function returning a scalar value.
#' @param ... Optional arguments to be passed to replace.
#' @return Vector with missing values (NA) replaced by the value of replace.
#' @export na.replace
na.replace <- function (x, replace, ...) {
	if (is.function(replace)){ 
		replace <- replace(x, ...)
	}
	x[is.na(x)] <- replace
	x
}