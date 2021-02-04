#' Bait Coverage
#' 
#' For an input set of sequences and baits targetting those sequences, this function returns either the fraction of sites or number of sites with at least 1X bait coverage for each input target sequence.
#' 
#' @param targets Set of target sequences.
#' @param baits Bait sequences. Multiple baits can be associated with one target sequence, but each bait is specific to one target.
#' @param distance Linear distance (number of nucleotides) allowed between target sites and bait sequence sites (i.e., interval on contig between the target and bait) for a target site to be considered "covered" by the bait. Default is 0.
#' @param min.depth This argument is not currently implemented. Future versions of this code might allow for calculating bait coverage at a user-specified minimum depth.
#' @param bait.length Nucleotide length of bait sequences.
#' @param fraction Report coverage as a fraction of the target sequence covered for each target sequence. Default is TRUE. If set to FALSE, the absolute number of covered target sites is returned for each input target sequence.
#' @return A vector containing coverage (fraction or number of sites) for each input target sequence.
#' @export
bait.coverage <- function(targets,baits,distance=0,min.depth=1,bait.length=120,fraction=T) {
	bait.matrix           <- do.call(rbind,strsplit(names(baits),split="_"))
	bait.matrix[,2]       <- as.numeric(bait.matrix[,2])+1
	bait.matrix           <- cbind(bait.matrix,as.numeric(bait.matrix[,2])+(bait.length-1))
	rownames(bait.matrix) <- names(baits)
	colnames(bait.matrix) <- c("target","probe.start.in.target","probe.end.in.target")
	result                <- vector(mode="numeric",length=length(targets))
	buffer                <- distance
	for(i in 1:length(targets)){
		bait.matrix.temp       <- bait.matrix[which(bait.matrix[,"target"]==names(targets[i])),,drop=F]
		ranges.matrix.temp     <- cbind(as.numeric(bait.matrix.temp[,"probe.start.in.target"]),as.numeric(bait.matrix.temp[,"probe.end.in.target"]))
		buffered.ranges.matrix <- cbind((ranges.matrix.temp[,1]-buffer),(ranges.matrix.temp[,2]+buffer))
		covered.temp           <- unique(as.numeric(apply(X=buffered.ranges.matrix,MARGIN=1,FUN=function(X){seq(from=X[1],to=X[2])})))
		target.temp.sites      <- 1:width(targets[i])
		covered                <- intersect(target.temp.sites,covered.temp)
		if(fraction){
			bait.coverage      <- length(covered)/width(targets[i])
		} else {
			bait.coverage      <- length(covered)
		}
		result[i]              <- bait.coverage
	}
	result
}
