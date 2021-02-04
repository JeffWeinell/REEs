#' Whole number test
#' 
#' Test if an input number is a whole number
#' 
#' @param x input number
#' @param tol tolerance?
#' @return TRUE if x is a whole number, otherwise FALSE
#' @export
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
	abs(x - round(x)) < tol
}
