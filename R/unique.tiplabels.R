#' All Unique Tree Tips
#' 
#' Return the unique set of tip labels for an input set of trees
#'
#' @param trees.multiPhylo Input set of trees held in an object of class multiPhylo (see ape Package for more details on multiPhylo class objects)
#' @return Character vector of tree tip labels occurring in at least one tree in the 
#' @export unique.tiplabels
unique.tiplabels <- function(trees.multiPhylo){
	label.search   <- stringr::str_locate(names(unlist(trees.multiPhylo)),pattern="tip.label")
	tip.labels.all <- unlist(trees.multiPhylo)[which(!is.na(label.search[,1]))]
	result         <- unique(tip.labels.all)
	result
}
