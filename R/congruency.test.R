#' Congruency Test
#' 
#' Pairwise comparison of two input trees testing if their topologies are congruent at highly supported nodes.
#' 
#' @param tree1 First input tree. Should be class Phylo or multiPhylo.
#' @param tree2 Second input tree. Should be class Phylo or multiPhylo.
#' @param min.support Threshold for a node to be considered strongly supported. Nodes with support less than this value are collapsed into a polytomy. Default is 0. Ignored for trees without support values.
#' @return Logical, TRUE if the two trees are congruent, and otherwise FALSE
#' @export congruency.test
congruency.test <-function(tree1,tree2,min.support=0){
	if(class(tree1)=="multiPhylo"){
		tree1 <- tree1[[1]]
	}
	if(class(tree2)=="multiPhylo"){
		tree2 <- tree2[[1]]
	}
	
	shared.tips    <- intersect(tree1$tip.label,tree2$tip.label)
	
	if(length(shared.tips)<4){
		isCongruent <- NA
	} else {
		### just picks the first shared terminal branch as the root
		arbitrary.root <- intersect(tree1$tip.label,tree2$tip.label)[1]
		tree1.rooted <- ape::root(tree1,arbitrary.root)
		tree2.rooted <- ape::root(tree2,arbitrary.root)

		# if node labels exist, collapse low supported nodes
		if(is.null(tree1.rooted$node.label)){
			tree1.collapsed <- tree1.rooted
		} else {
			#is.numeric(tree1$node.label)
			tree1.rooted$node.label <- as.numeric(gsub(".*/","",tree1.rooted$node.label))
			tree1.collapsed <- REEs::collapse.low.support.nodes(tree1.rooted,support.threshold=min.support)
		}
		if(is.null(tree2.rooted$node.label)){
			tree2.collapsed <- tree2.rooted
		} else {
			tree2.rooted$node.label <- as.numeric(gsub(".*/","",tree2.rooted$node.label))
			tree2.collapsed <- REEs::collapse.low.support.nodes(tree2.rooted,support.threshold=min.support)
		}
		#if (FALSE){
		#### collapsing low-support nodes
		#	tree1.collapsed <- collapse.low.support.nodes(tree1.rooted,support.threshold=min.support)
		#	tree2.collapsed <- collapse.low.support.nodes(tree2.rooted,support.threshold=min.support)
		#} else {
		#	tree1.collapsed <- tree1.rooted
		#	tree2.collapsed <- tree2.rooted
		#}
		
		### dropping tips that are not shared between trees
		shared.tips      <- intersect(tree1.collapsed$tip.label,tree2.collapsed$tip.label)
		tree1.drop.names <- setdiff(tree1.collapsed$tip.label,shared.tips)
		tree2.drop.names <- setdiff(tree2.collapsed$tip.label,shared.tips)
		
		tree1.trimmed <- drop.tip(tree1.collapsed,tree1.drop.names)
		tree2.trimmed <- drop.tip(tree2.collapsed,tree2.drop.names)
		
		### internal node numbers for the collapsed and trimmed trees
		tree1.NodesInternal <- (length(tree1.trimmed$tip.label)+1):max(tree1.trimmed$edge)
		tree2.NodesInternal <- (length(tree2.trimmed$tip.label)+1):max(tree2.trimmed$edge)
	
		### finding the set of clade-sets for each collapsed and trimmed trees
		tree1.cladeSets <- phangorn::Descendants(tree1.trimmed,node=tree1.NodesInternal,type="tips")
		tree2.cladeSets <- phangorn::Descendants(tree2.trimmed,node=tree2.NodesInternal,type="tips")
	
		tree1.cladeSets.names <- lapply(X=tree1.cladeSets,FUN=function(input){tree1.trimmed$tip.label[input]})
		tree2.cladeSets.names <- lapply(X=tree2.cladeSets,FUN=function(input){tree2.trimmed$tip.label[input]})
	
		### pairwise comparison of cladeSets to check if the trees are congruent
		pairwise.cladeSets <- list(); length(pairwise.cladeSets) <- length(tree1.cladeSets.names)
		for(i in 1:length(tree1.cladeSets)){
			 ### checks if every cladeSet of tree1 exists as either (1) a matching set, subset, or does not intersect with each cladeSet of tree2.
			pairwise.cladeSets[i]  <- all(unlist(lapply(X=tree2.cladeSets.names,FUN=function(input){all(input %in% tree1.cladeSets.names[[i]]) | all(tree1.cladeSets.names[[i]] %in% input) | all(! input %in% tree1.cladeSets.names[[i]])})))
		}
		pairwise.cladeSets <- unlist(pairwise.cladeSets)
		isCongruent        <- all(pairwise.cladeSets)
	}
	isCongruent
}


