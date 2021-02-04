#' Get All Quartets
#' 
#' Get the set of all 4-taxon sub-trees (quartets) of an input tree
#' 
#' @param ref.tree Newick format phylogenetic tree
#' @param alt.tips Alternative set of tips (characterized by tip labels) to create quartets for. Default is NULL, meaning that all tips in the input tree are all allowed in quartets and all possible quartets are included in output. If alt.tips is set to a value other than NULL, any quartets containing tips not in the alt.tips are filtered from the output. It is often useful to set alt.tips = unique.tiplabels(multiPhylo), where multiPhylo is an object of class multiPhylo that holds a set of gene trees.
#' @param support.scaler Value to multiply the support values by (default is 1). It might be useful to scale the support values of some trees if comparing trees with different methods of support (e.g., if comparing Bayesian and ML trees). I haven't experimented with this yet.
#' @return multiPhylo object containing the set of all quartets for an input tree and a set of tip labels.
#' @export 
get.all.quartets <- function(ref.tree,alt.tips=NULL,support.scaler=1){
	best.tree   <- ref.tree
	if(is.null(alt.tips)){
		taxa = best.tree$tip.label
	} else {
		taxa = alt.tips
	}
	group.index   <- sample(x=c(1:4),size=length(taxa),replace=T)                                                                                 ### assigns an index number to each taxon (i.e., tip label) in taxa vector
	all.quartets     <- xprod.combn(list(taxa[group.index==1],taxa[group.index==2],taxa[group.index==3],taxa[group.index==4]))                       ### returns a list of character vectors (each length four), and each is a unique combination of the tip labels
	all.quartets.mat <- mat.strsplit(all.quartets)                                                                                                      ### each row is a unique 4-taxon set (same as all.quartets list, but held in a four column matrix)
	best.tree$edge.length <- rep(1,length(best.tree$edge.length))                                                                                 ### sets all edge lengths equal to 1, because this is a test of topology
	best.tree$node.label[which(best.tree$node.label !="")] <- as.numeric(best.tree$node.label[which(best.tree$node.label !="")])*(support.scaler) ### useful to scale the support values if planning to compare the quartet.trees to a set of gene trees with different range of support values compared to ref.tree
	all.quartets.counter = 0                                                                                                                         ### empty vector that will contain the set of quartet trees within the best tree
	for(i in 1:nrow(all.quartets.mat)){
		if(all(all.quartets.mat[i,] %in% best.tree$tip.label)){
			drop.temp       <- best.tree$tip.label[!(best.tree$tip.label %in% all.quartets.mat[i,])] #### all tips except for those in all.quartets.mat[i,]
		} else {
			next
		}
		quartet.temp       <- drop.tip(best.tree,drop.temp)
		all.quartets.counter = (all.quartets.counter +1)        ### increase the counter by 1 if didnt skip to next cycle of loop
		if(all.quartets.counter==1){
			all.quartets.trees <- c(quartet.temp)
		} else {
			all.quartets.trees <- c(all.quartets.trees,quartet.temp)
		}
	}
	all.quartets.trees
}
