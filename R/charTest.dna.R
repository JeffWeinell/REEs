#' Check if Nucleotides Present
#'
#' Check if a character vector contains either at least one, or all, of the four DNA nucleotides (A C G T).
#'
#' @param ls Input DNA sequence held in character vector
#' @param method Either "any" or "all". If "any", then the function checks if any of the four nucleotides are present. If "all", then the function checks if all of the nucleotides are present.
#' @return TRUE or FALSE
#' @export
charTest.dna <- function(ls,method=c("any","all")){
		ls <- tolower(ls)
		if(method=="any"){
			result <- any(c("a","c","g","t") %in% ls)
		}
		if(method=="all"){
			result <- all(c("a","c","g","t") %in% ls)
		}
		result
}
