#' Read Multiple Newick Trees
#' 
#' Reads in newick trees from a vector of filepaths.
#' 
#' @param filepaths Vector of filepaths to the newick tree files to read.
#' @return An object of class multiPhylo that holds the trees.
#' @export 
read.tree.multipleFiles <- function(filepaths){
	treenames.local <- paste0("tree",seq(1:length(filepaths)))
	for(i in 1:length(filepaths)){
		assign(treenames.local[i],ape::read.tree(filepaths[i]))
		if(i==1){
			all.trees <- get(treenames.local[i])
		} else {
			all.trees <- c(all.trees,get(treenames.local[i]))
		}
	}
	all.trees
}
