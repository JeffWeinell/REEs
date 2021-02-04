#' Collapse Low Supported Nodes
#' 
#' Given an input phylogeny, returns a phylogeny with low support nodes (below a threshold value of support) collapsed into polytomies.
#' 
#' @param rooted.tree Rooted input phylogeny.
#' @param support.threshold Nodes with support less than this value are collapsed into a polytomy. Default is 80, but care should be taken for this value. 
#' @return Returns a phylogeny with the low support nodes collapsed into polytomies.
#' @export collapse.low.support.nodes
collapse.low.support.nodes <- function(rooted.tree,support.threshold=80){
	tree1.rooted <- rooted.tree
	### This value is half the length of the shortest edge, and is used as the cuttoff for removing branches that will get set to zero-length
	tolerance                 <- min(tree1.rooted$edge.length)/2
	ntips                     <- length(tree1.rooted$tip.label)
	### Also may want to drop nodes without a support value other than the root?
	tree1.nodes.to.drop       <- (ntips + which(as.numeric(tree1.rooted$node.label) <= support.threshold))
	### The edges having the descendent node being a weakly supported node that we want to drop
	tree1.crown.edges.to.drop <- which(tree1.rooted$edge[,2] %in% tree1.nodes.to.drop)
	tree1.rooted2             <- tree1.rooted
	### Sets edge length to zero for the edges that you want to drop
	tree1.rooted2$edge.length[tree1.crown.edges.to.drop] <- 0
	### Deletes the zero-length edges deleted (this is the tree with low-support nodes collapsed!!!).
	### Switched from using phytools::di2multi to ape::di2multi, because it appears that the former is the same as the latter.
	tree1.rooted3             <- ape::di2multi(tree1.rooted2,tol=tolerance)
	tree1.rooted3
}
