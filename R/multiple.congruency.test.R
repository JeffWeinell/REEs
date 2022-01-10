#' Multiple Congruency Test
#' 
#' Pairwise comparison of input trees and if topologies of each pair is congruent at highly supported nodes.
#' 
#' @param ... Two or more input trees (object names separated by commas) of class multiPhylo
#' @param min.support Nodes with support less than this value are collapse into a polytomy. Default is 80, but care should be taken for this value. 
#' @param save.as Path to output file. Defult NULL. When NULL, runtimes should be faster, but crashes during long runs will have to be restarted. If an output path is provided then analysis can be resumed.
#' @return Returns a three-column character matrix. Each row is a comparison of two input trees; first two columns indicate trees compared; third column indicates if the two trees are congruent.
#' @export multiple.congruency.test
multiple.congruency.test <- function(...,min.support=80,save.as=NULL){
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
	#pairwise.treesNames.mat <- xprod.combn.mat(list.of.treeNames,list.of.treeNames)
	pairwise.treesNames.mat <- as.matrix(expand.grid(list.of.treeNames,list.of.treeNames))
	pairwise.treesNames.mat <- pairwise.treesNames.mat[!duplicated(t(apply(pairwise.treesNames.mat,1,sort))),]
	#tmpmatpath <- tempfile()
	#write.table(pairwise.treesNames.mat,tmpmatpath,col.names=F,row.names=F,sep='\t',quote=FALSE)
	#system(sprintf("cat '%s' | sort -V > '%s' ",tmpmatpath,tmpmatpath))
	if(is.null(save.as)){
		congruence.result       <- apply(X=pairwise.treesNames.mat,MARGIN=1,FUN=function(input,support.thresh=min.support){tree1=get(input[1]);tree2=get(input[2]);REEs::congruency.test(tree1,tree2,min.support=support.thresh)})
		result <- data.frame(tree1=pairwise.treesNames.mat[,1],tree2=pairwise.treesNames.mat[,2],congruent=congruence.result)
	} else {
		if(!dir.exists(dirname(save.as))){
			return(paste0("output directory does not exist:",dirname(save.as)))
		} else {
			if(file.exists(save.as)){
				result.last <- read.table(save.as,header=T,sep='\t')
				#tree1.last  <- result.last[nrow(result.last),1]
				#STARTi=match(tree1.last,list.of.treeNames) + 1
				pairwise.treesNames.df  <- data.frame(tree1=pairwise.treesNames.mat[,1],tree2=pairwise.treesNames.mat[,2])
				#print(dim(pairwise.treesNames.df))
				#print(colnames(pairwise.treesNames.df))
				#print(colnames(result.last))
				#stop()
				pairwise.treesNames.df2 <- data.table::setDT(pairwise.treesNames.df)[!result.last[,c('tree1','tree2')], on = names(pairwise.treesNames.df)]
				if(nrow(pairwise.treesNames.df2)==0){
					return(print("No new comparisons."))
				}
				pairwise.treesNames.mat <- as.matrix(pairwise.treesNames.df2)
				list.of.treeNames <- unique(pairwise.treesNames.mat[,1])
				for(i in 1:length(list.of.treeNames)){
					print(sprintf("%s/%s",i,length(list.of.treeNames)))
					pairwise.treesNames.mat.i <- pairwise.treesNames.mat[pairwise.treesNames.mat[,1] == list.of.treeNames[i],,drop=F]
					#pairwise.treesNames.mat.i <- pairwise.treesNames.mat[grep(list.of.treeNames[i],pairwise.treesNames.mat[,1],fixed=T),,drop=F]
					congruence.result.i       <- apply(X=pairwise.treesNames.mat.i,MARGIN=1,FUN=function(input,support.thresh=min.support){tree1=get(input[1]);tree2=get(input[2]);REEs::congruency.test(tree1,tree2,min.support=support.thresh)})
					result.i <- data.frame(tree1=pairwise.treesNames.mat.i[,1],tree2=pairwise.treesNames.mat.i[,2],congruent=congruence.result.i)
					write.table(result.i,save.as,append=T,col.names=F,row.names=F,quote=F,sep='\t')
				}
			} else {
				for(i in 1:length(list.of.treeNames)){
					print(sprintf("%s/%s/%s",STARTi,i,length(list.of.treeNames)))
					pairwise.treesNames.mat.i <- pairwise.treesNames.mat[grep(list.of.treeNames[i],pairwise.treesNames.mat[,1],fixed=T),,drop=F]
					congruence.result.i       <- apply(X=pairwise.treesNames.mat.i,MARGIN=1,FUN=function(input,support.thresh=min.support){tree1=get(input[1]);tree2=get(input[2]);REEs::congruency.test(tree1,tree2,min.support=support.thresh)})
					result.i <- data.frame(tree1=pairwise.treesNames.mat.i[,1],tree2=pairwise.treesNames.mat.i[,2],congruent=congruence.result.i)
					write.table(result.i,save.as,append=F,col.names=T,row.names=F,quote=F,sep='\t')
				}
			}
			congruence.result <- read.table(save.as,header=T,sep='\t')[,'congruent']
			result <- read.table(save.as,header=T,sep='\t')
		}
	result
	}
}


