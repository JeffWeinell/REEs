#' Any AAs Present
#'
#' Check that at least one amino acid is present in a putative amino acid sequence
#'
#' @param ls Input amino acid sequence held in character vector
#' @return TRUE or FALSE
#' @export
charTest.aa <- function(ls){
		#ls <- tolower(ls)
		#all(c("a","c","g","t") %in% ls)
		any(GENETIC_CODE %in% ls)
}
