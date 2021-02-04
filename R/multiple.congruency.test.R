#' Multiple Congruency Test
#' 
#' Pairwise comparison of input trees and if topologies of each pair is congruent at highly supported nodes.
#' 
#' @param ... Two or more input trees (object names separated by commas) of class 
#' @param min.support Nodes with support less than this value are collapse into a polytomy. Default is 80, but care should be taken for this value. 
#' @return Returns a three-column character matrix. Each row is a comparison of two input trees; first two columns indicate trees compared; third column indicates if the two trees are congruent.
#' @export 
multiple.congruency.test <- function(...,min.support=80){
	list.of.trees <- list(...)
	if(length(list.of.trees)==1){
		list.of.trees <- list.of.trees[[1]]
	} else {
		for(i in seq(length(list.of.trees))){
			names(list.of.trees)[i]  <-  as.character(sys.call()[[i+1]])
		}
	}
	if(any(is.null(names(list.of.trees)))){
		names(list.of.trees) <- paste0("tree.name",seq(1:length(list.of.trees)))
	}
	
	list.of.treeNames       <- names(list.of.trees)
	for(i in 1:length(list.of.treeNames)){
		assign(list.of.treeNames[i],list.of.trees[i])
		#class(get(list.of.treeNames[i]))
	}
	pairwise.treesNames.mat <- xprod.combn.mat(list.of.treeNames,list.of.treeNames)
	pairwise.treesNames.mat <- pairwise.treesNames.mat[!duplicated(t(apply(pairwise.treesNames.mat,1,sort))),]
	congruence.result       <- apply(X=pairwise.treesNames.mat,MARGIN=1,FUN=function(input,support.thresh=min.support){tree1=get(input[1]);tree2=get(input[2]);congruency.test(tree1,tree2,min.support=support.thresh)})
	result                  <- cbind(pairwise.treesNames.mat,congruence.result)
	#result                  <- result[!duplicated(t(apply(result,1,sort))),]
	colnames(result) <- c("tree1","tree2","congruent")
	result
}
