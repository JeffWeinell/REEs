#' @title Make Partitioned Alignments
#' 
#' Given an input DNA alignment and a reference CDS sequence for the locus in the input alignment, this function generates up to nine output alignments, including three partitioned DNA alignments, five DNA alignments containing a subset of the input DNA sites, and one amino acid alignment for the CDS region(s) of the input alignment.
#' 
#' @param input.path Folder containing fasta-formatted input DNA alignments.
#' @param output.dir Where to save alignments (alignments for each type of data will be in different subdirectory)
#' @param TargetCDS.path Full path to the fasta file containing only the CDS sequences of target loci (from which probes were designed). Sequence names must have the following format: "GeneName_TargetCDS_of_TargetLocusName_AnyAdditionalIformation", where GeneName and TargetLocusName are replaced with the actual names, and AnyAdditionalIformation can be a string of any characters
##' @param ref.type Indicates that the input alignments are DNA sequences (default is "DNA").
#' @param steps Character string with twelve characters, each 0 or 1, indicating which steps of the function to do (0) or skip (1). Used for debugging.
#' @param mafft.params Character string with parameters passed to MAFFT. Default is " --auto --adjustdirection --nwildcard --op 3 --ep 0.123 --quiet "
#' @param old.names Character vector of old sample names that should be changed to those defined in new.names parameter (default is NA, meaning don't change names).
#' @param new.names Character vector of new sample names to use instead of corresponding name in old.names vector (default is NA, meaning don't change names).
#' @param drop.reference Whether or not to include the reference sequence (from which baits were designed for a particular locus) in the alignment written to file (default is FALSE).
#' @param ith.locus.start First locus in input.path to process (default is 1).
#' @param ith.locus.end Last locus in input.path to process (default is "all", when means process until no more loci to process).
#' @param locus.names.omit Names of loci to skip (default is NULL).
#' @param remove.outliers Logical indivating if odseq::odseq should be used to detect and remove outlier sequences. Default FALSE.
#' @param AA.pdist.drop.thresh Maximum mean pairwise p-distance for AA sequence of an individual (translated from CDS regions of input DNA sequences) to keep the individual in AA or CDS-containing DNA alignments.
#' @param trimto Optional character string with names of sequences whose combined range should be used to trim the alignments. Default is NULL.
#' @return NULL; writes up to nine different fasta-formatted sequence alignments for each input DNA alignment in input.path, and partition files partitioned alignments. The nine output alignment types are:
#' (1) Partitioned version of the input DNA alignment. The alignment is partitioned into (if present): upstream noncoding DNA, each codon position (within CDS regions), and downstream noncoding DNA.
#' (2) Alignment containing only the CDS regions, which are partitioned by codon position.
#' (3) Alignment containing only the upstream noncoding DNA.
#' (4) Alignment containing only the downstream noncoding DNA.
#' (5) Alignment containing only noncoding DNA, partitioned into upstream and downstream regions.
#' (6) Alignment containing only the first codon positions of CDS regions.
#' (7) Alignment containing only the second codon positions of CDS regions.
#' (8) Alignment containing only the third codon positions of CDS regions.
#' (9) Alignment containing the amino acid sequence for the translated CDS region.
#' @export make.partitioned.alignment
make.partitioned.alignment  <- function(input.path,output.dir,TargetCDS.path,steps="000000000000",mafft.params=" --auto --adjustdirection --nwildcard --op 3 --ep 0.123 --quiet ",old.names=NA,new.names=NA,drop.reference=F,ith.locus.start=1,ith.locus.end="all",locus.names.omit=NULL,remove.outliers=FALSE,AA.pdist.drop.thresh=0.5,trimto=NULL){
	
	### makes necessary output subdirectories if they dont exist ###
	dir1  <- file.path(output.dir,"All_parts/alignmentFiles/")
	dir2  <- file.path(output.dir,"All_parts/partitionFiles/")
	dir3  <- file.path(output.dir,"CDS_only/alignmentFiles/")
	dir4  <- file.path(output.dir,"CDS_only/partitionFiles/")
	dir5  <- file.path(output.dir,"Upstream_noncoding/alignmentFiles/")
	dir6  <- file.path(output.dir,"Downstream_noncoding/alignmentFiles/")
	dir7  <- file.path(output.dir,"All_noncoding/alignmentFiles/")
	dir8  <- file.path(output.dir,"All_noncoding/partitionFiles/")
	dir9  <- file.path(output.dir,"CDS_FirstCodonPosition/alignmentFiles/")
	dir10 <- file.path(output.dir,"CDS_SecondCodonPosition/alignmentFiles/")
	dir11 <- file.path(output.dir,"CDS_ThirdCodonPosition/alignmentFiles/")
	dir12 <- file.path(output.dir,"AminoAcids/alignmentFiles/")

	subdirectories  <- c(dir1,dir2,dir3,dir4,dir5,dir6,dir7,dir8,dir9,dir10,dir11,dir12)
	makeDirectories <- lapply(X=subdirectories,FUN=dir.check.create)
	
	steps <- as.numeric(REEs::s2v(steps))
	mafft.params2 <- gsub("--adjustdirection"," ",mafft.params)
	
	###################
	### Input files ###
	###################
	### Input fasta file(s) with sequences, each with DNA sequences for a locus
	if(dir.exists(input.path)){
		input.seqs.filenames <- gtools::mixedsort(list.files(input.path,full.names=T))
	} else {
		input.seqs.filenames <- input.path
		ith.locus.start <- 1
		ith.locus.end   <- 1
	}
	### Read fasta file with DNA for the CDS regions of target references
	TargetDNA_CDS.regions       <- Biostrings::readDNAStringSet(TargetCDS.path,format="fasta")
	### Read the table specifying which species the probes were designed from for each locus 
#	bait.species.table          <- data.table::fread(bait.species.filename)
	
	### Extracting the names of target reference sequences from the names used in the file TargetDNA_CDS.regions
#	name.string.start           <- (str_locate_X(strings=names(TargetDNA_CDS.regions),pattern="_",X=3)+1)
#	name.string.end             <- (str_locate_X(strings=names(TargetDNA_CDS.regions),pattern="_",X=4)-1)
	# target loci names of each sequence in TargetDNA_CDS.regions
#	targetCDS.names   <- substring(names(TargetDNA_CDS.regions),first=name.string.start,last=name.string.end)
	targetCDS.names   <- gsub("_.*","",gsub(".*_TargetCDS_of_","",names(TargetDNA_CDS.regions)))

	### Next line assumes that files are named according to WeinelEntry names
	input.seqs.shortnames  <- gsub(".phy|.fa|.fasta","",basename(input.seqs.filenames))

	### Names of loci that have a CDS region included in TargetDNA_CDS.regions
	shared.names                <- Reduce(intersect, list(input.seqs.shortnames,targetCDS.names))
	### names of loci to ignore
	if(!is.null(locus.names.omit)){
		shared.names  <- setdiff(shared.names,locus.names.omit)
	}
	if(file.exists(file.path(output.dir,"partitioned_alignments_made.txt"))){
		alignments.made <- utils::read.table(file=file.path(output.dir,"partitioned_alignments_made.txt"), header=T, colClasses="character")
		if(!!length(setdiff(shared.names, rownames(alignments.made)))){
			newrows <- matrix(data="no",nrow=length(setdiff(shared.names, rownames(alignments.made))),ncol=9,dimnames=list(setdiff(shared.names, rownames(alignments.made)),colnames(alignments.made)))
			alignments.made <- rbind(alignments.made,newrows)
		}
	} else {
		alignments.made <- matrix(data="no", nrow=length(shared.names), ncol=9)
		rownames(alignments.made)   <- shared.names
		colnames(alignments.made)   <- c("all.data","CDS","5'NC","3'NC","Noncoding","AA","FirstCodon","SecondCodon","ThirdCodon")
	}
	if(ith.locus.end=="all"){
		last.locus.process          <- length(shared.names)
	} else {
		last.locus.process          <- ith.locus.end
	}
	for(i in ith.locus.start:last.locus.process){
		ALIGNMENTSi <- list(); length(ALIGNMENTSi) <- 9
		PARTSi <- ALIGNMENTSi
		#if(is.wholenumber(i/50)){
		#	print(i)
		#}
		print(i)
		locus.name.temp   <- shared.names[i]
		print(locus.name.temp)
		FILEPATHi         <- input.seqs.filenames[input.seqs.shortnames==locus.name.temp]
		### reading and unaligning input seqs
		novel             <- REEs::trimXN(Biostrings::DNAStringSet(x=gsub("-|\\?","",Biostrings::readDNAStringSet(FILEPATHi))))
		reference.cds     <- TargetDNA_CDS.regions[grep(paste0("_",locus.name.temp,"_"),names(TargetDNA_CDS.regions))]
		### Skip locus if 
		if(length(reference.cds)>1){
			print(sprintf("Multiple CDS regions for locus %s, skipping this locus", locus.name.temp))
			next
		}
		SEQSi             <- c(reference.cds, novel)
		# Skips locus if <4 sequences other than the reference CDS
		# if(FALSE) {
		#	if(length(SEQSi)<5){
		#		next
		#	}
		#}
		alignment.names  <- names(SEQSi)
		if(!all(is.na(old.names))){
			alignment.names  <- mgsub(old.names,new.names,alignment.names)
		}
		### Calling MAFFT
		#alignment   <- REEs::mafft(SEQSi, param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6")
		if(!steps[1]){
			print("Beginning MAFFT run 1/7")
			alignment   <- REEs::mafft(REEs::trimXN(SEQSi), param=mafft.params)
			alignment_A <- alignment
			### Trim alignment to the maximum width region that includes at least one individual in 'trimto' at the first and last base
			if(!is.null(trimto)){
				alignment   <- REEs::trimTo(aln=alignment,nam=trimto)
				alignment_B <- alignment
			} else {
				alignment_B <- NULL
			}
		}
		# plot1 <- REEs::plotAlignment(alignment_B,title=paste(locus.name.temp,"alignment_B"))
		### Use odseq package to remove outlier individuals from alignment, and then rerun mafft
		# print("Removing outlier individuals")
		if(remove.outliers){
			outliers         <- odseq::odseq(DNAMultipleAlignment(alignment), threshold = 0.01, distance_metric = "affine", B = 1000)
			if(!steps[2]){
				if(!!length(which(outliers))) {
					print(sprintf("%s outlier sequences removed",length(which(outliers))))
					print("Beginning MAFFT run 2/7")
					alignment <- REEs::mafft(REEs::trimXN(Biostrings::DNAStringSet(x=gsub("-","",alignment[!outliers]))), param=mafft.params2)
					#plot2     <- REEs::plotAlignment(alignment, title=paste(locus.name.temp,"MAFFT run 2"))
				} else {
					print("No outlier sequences, skippng MAFFT run 2/7")
				}
			}
		} else {
			print("skipping odseq outlier detection and MAFFT run 2/7")
		}
		### Define "upstream.noncoding", "CDS", and "downstream.noncoding" regions by comparing to the CDS of the reference target
		begin.exon       <- stringr::str_locate(string=alignment[1],pattern=as.character(subseq(reference.cds,start=1,end=1)))[1]
		end.exon         <- REEs::str_locate_last(string=alignment[1],pattern=as.character(subseq(reference.cds,start=(width(reference.cds)-1),end=width(reference.cds))))
		starts           <- c(1,begin.exon,begin.exon+1,begin.exon+2,end.exon+1)
		ends             <- c(begin.exon-1,end.exon,end.exon,end.exon,width(alignment[1]))
		widths           <- ends-(starts-1)
		extra.text       <- c("","\\3","\\3","\\3","")
		ranges           <- cbind(starts,ends,widths,extra.text)
		rownames(ranges) <- c("upstream.noncoding","CDS.1","CDS.2","CDS.3","downstream.noncoding")
		
		### logical vector indicating if upstream noncoding, CDS, and downstream noncoding data are present
		partition.groups        <- c(widths[1]>0, all(widths[2:4]>0), widths[5]>0)
		names(partition.groups) <- c("upstream.noncoding","CDS","downstream.noncoding")
		print(partition.groups)
		if(!steps[3]){
			#### CDS-only alignment (preliminary)
			if(partition.groups[2]){
				if(drop.reference==T){
						cds.alignment     <- Biostrings::subseq(x=alignment,start=begin.exon,end=end.exon)[-1]
					} else {
						cds.alignment     <- Biostrings::subseq(x=alignment,start=begin.exon,end=end.exon)
				}
				### specifies individuals with no data in the region
				drop.nodata               <- stringr::str_count(cds.alignment,"-") == width(cds.alignment)
				if(length(which(!drop.nodata))>3){
					cds.alignment         <- cds.alignment[!drop.nodata]
					print("Beginning MAFFT run 3/7")
					cds.alignment2        <- REEs::mafft(REEs::trimXN(Biostrings::DNAStringSet(x=gsub("-","", cds.alignment))), param=mafft.params2)
					if(!is.null(trimto)){
						cds.alignment2    <- REEs::trimTo(aln=cds.alignment2, nam=trimto)
					}
					alignment_C           <- cds.alignment2
				} else {
					### updates decision from TRUE to FALSE on whether to write the CDS alignment, because too few individuals with CDS data
					partition.groups[2] <- F
					alignment_C <- NULL
				}
			}
			if(!partition.groups[2]){
				print("No CDS region, skipping MAFFT run 3/7")
			}
		}
		if(!steps[4]){
			#### Upstream noncoding alignment
			if(partition.groups[1]){
				if(begin.exon!=1){
					if(drop.reference==T){
						upstream.alignment    <- Biostrings::subseq(x=alignment,start=1,end=(begin.exon-1))[-1]
					} else {
						upstream.alignment    <- Biostrings::subseq(x=alignment, start=1, end=(begin.exon-1))
					}
					drop.nodata               <- stringr::str_count(upstream.alignment,"-") == width(upstream.alignment)
					if(length(which(!drop.nodata))>3){
						upstream.alignment         <- upstream.alignment[!drop.nodata]
						print("Beginning MAFFT run 4/7")
						upstream.alignment2        <- REEs::mafft(REEs::trimXN(Biostrings::DNAStringSet(x=gsub("-","",upstream.alignment))), param=mafft.params2)
						if(!is.null(trimto)){
							upstream.alignment2    <- REEs::trimTo(aln=upstream.alignment2, nam=trimto)
						}
						alignment_D <- upstream.alignment2
					} else {
						partition.groups[1] <- F   ### updates decision from TRUE to FALSE on whether to write the upstream alignment, because too few individuals with upstream noncoding data
						alignment_D <- NULL
					}
				}
			}
		}
		if(!steps[5]){
			#### Downstream noncoding alignment
			if(partition.groups[3]){
				if(end.exon!=width(alignment[1])){
					if(drop.reference==T){
						downstream.alignment    <- Biostrings::subseq(x=alignment, start=(end.exon+1), end=width(alignment[1]))[-1]
					} else {
						downstream.alignment    <- Biostrings::subseq(x=alignment, start=(end.exon+1), end=width(alignment[1]))
					}				
					drop.nodata                 <- stringr::str_count(downstream.alignment,"-") == width(downstream.alignment)
					### Skips downstream noncoding alignment if 3 or less individuals
					if(length(which(!drop.nodata))>3){
						downstream.alignment    <- downstream.alignment[!drop.nodata]
						length(downstream.alignment)
						print("Beginning MAFFT run 5/7")
						downstream.alignment2   <- REEs::mafft(REEs::trimXN(Biostrings::DNAStringSet(x=gsub("-","",downstream.alignment))), param=mafft.params2)
						if(!is.null(trimto)){
							downstream.alignment2    <- REEs::trimTo(aln=downstream.alignment2,nam=trimto)
						}
						alignment_E <- downstream.alignment2
					} else {
						### Updates decision from TRUE to FALSE on whether to write the downstream alignment, because too few individuals with downstream noncoding data.
						partition.groups[3] <- F
						alignment_E <- NULL
					}
				}
			}
		}
		##########
		##### Remainder of Code below was added to include parts of the code in the R script "make_codon_alignments.R"
		##########
		if(partition.groups[2]){
			if(!steps[6]){
				cds.temp2   <- cds.alignment2
				p1 <- which(is.wholenumber(c(0:(width(cds.temp2[1])-1))/3)) ### numerical vector of first codon position sites
				p2 <- intersect((p1+1),c(2:width(cds.temp2[1])))            ### numerical vector of second codon position sites
				p3 <- intersect((p2+1),c(3:width(cds.temp2[1])))            ### numerical vector of third codon position sites
				### CDS sequences aligned (extracted from whole alignment) but not yet re-aligned
				cds.temp3 <- Biostrings::subseq(cds.temp2,start=min(p1), end=max(p3))
				### CDS sequences not aligned
				cds.temp4 <- Biostrings::DNAStringSet(x=gsub("-","",cds.temp3))   ### unaligns sequences (after removing gaps) so that the sequences can be translated 
				### Translates CDS sequences, then does MAFFT multiple sequence alignment to make an alignment of amino acids
				print("Translating CDS region")
				aa.temp             <- REEs::trimXN(suppressWarnings(Biostrings::translate(cds.temp4, if.fuzzy.codon="solve",no.init.codon=T)))
				print("Beginning MAFFT run 6/7")
				aa.alignment.temp   <- REEs::mafft(aa.temp, param=mafft.params2)
				# Nothing filtered by default
				aa.alignment.temp2  <- REEs::filter.alignment(aa.alignment.temp)
				# This alignment has no missing or ambiguous data. This alignment is used for calculating median pairwise p-distances, but is not the alignment written to file.
		##		aa.alignment.mdt0   <- REEs::filter.alignment(aa.alignment.temp, mdt=0)
				print("Calculating median pairwise distances for each individual in the protein alignment")
				### Median pairwise distances for the AA alignment that can contain any amount of missing or ambiguous data. This info not currently used.
				if(TRUE) {
				##	distAA.mdt1          <- Biostrings::stringDist(aa.alignment.temp2, upper=T,diag=T)
					distAA.mdt1          <- Biostrings::stringDist(aa.alignment.temp2, method="substitutionMatrix",substitutionMatrix=REEs::submat("AA"), upper=T,diag=T)
					distAA.mdt1.mat      <- matrix(nrow=length(aa.alignment.temp2), ncol=length(aa.alignment.temp2))
					distAA.mdt1.mat[lower.tri(distAA.mdt1.mat)] <- as.numeric(distAA.mdt1)
					distAA.mdt1.mat[upper.tri(distAA.mdt1.mat)] <- t(distAA.mdt1.mat)[upper.tri(distAA.mdt1.mat)]
					rownames(distAA.mdt1.mat) <- labels(aa.alignment.temp2)
					### Number of AAs comparable between sequences
					distAA.mdt1_compared <- Biostrings::stringDist(aa.alignment.temp2, method="substitutionMatrix",substitutionMatrix=REEs::submat("AA_comparable"), upper=T,diag=T)
					distAA.mdt1_compared.mat  <- matrix(nrow=length(aa.alignment.temp2), ncol=length(aa.alignment.temp2))
					distAA.mdt1_compared.mat[lower.tri(distAA.mdt1_compared.mat)] <- as.numeric(distAA.mdt1_compared)
					distAA.mdt1_compared.mat[upper.tri(distAA.mdt1_compared.mat)] <- t(distAA.mdt1_compared.mat)[upper.tri(distAA.mdt1_compared.mat)]
					rownames(distAA.mdt1_compared.mat) <- labels(aa.alignment.temp2)
					### fraction pairwise differences = number of differences/number of comparisons
					median.distAA.mdt1   <- apply(X=(distAA.mdt1.mat/distAA.mdt1_compared.mat),MARGIN=1,FUN=median,na.rm=T)
					#median.distAA.mdt1  <- apply(distAA.mdt1,1,median,na.rm=T)/width(aa.alignment.temp2[1])
					#median.distAA.mdt1  <- apply(X=distAA.mdt1,MARGIN=1,FUN=median)/width(aa.alignment.temp2[1])
				}
				### Median pairwise p-distances for the AA alignment containing no missing or ambiguous data
				if(FALSE){
					if(!all(width(aa.alignment.mdt0)==0)){
						distAA.mdt0        <- Biostrings::stringDist(aa.alignment.mdt0, upper=T,diag=T)
						# Move distance data to a matrix
						distAA.mdt0.mat    <- matrix(nrow=length(aa.alignment.mdt0), ncol=length(aa.alignment.mdt0))
						distAA.mdt0.mat[lower.tri(distAA.mdt0.mat)] <- as.numeric(distAA.mdt0)
						# reflect lower triangle to make symmetric, then add row labels
						distAA.mdt0.mat[upper.tri(distAA.mdt0.mat)] <- t(distAA.mdt0.mat)[upper.tri(distAA.mdt0.mat)]
						rownames(distAA.mdt0.mat) <- labels(aa.alignment.mdt0)
						# median percent distance for each individual
						median.distAA.mdt0 <- apply(distAA.mdt0.mat,1,median,na.rm=T)/width(aa.alignment.mdt0[1])
						# median percent pairwise distance for each sample
						# colnames(distAA.mdt0.mat) <- rownames(distAA.mdt0.mat) <- labels(aa.alignment.mdt0)
						# median.distAA.mdt0 <- apply(X=distAA.mdt0,MARGIN=1,FUN=median)/width(aa.alignment.mdt0[1])
						# median.distAA.mdt0 <- median(distAA.mdt0)/width(aa.alignment.mdt0[1])
					} else {
						distAA.mdt0        <- NA
						median.distAA.mdt0 <- NA
					}
				}
				print(sprintf("Removing sequences with median percent pairwise distance > %s in protein alignment", AA.pdist.drop.thresh))
			#	toDrop <- which(median.distAA.mdt0 > AA.pdist.drop.thresh) ###| individuals with median p-distance greater than AA.pdist.drop.thresh should be dropped
				toDrop <- which(median.distAA.mdt1 > AA.pdist.drop.thresh)
				print(sprintf("%s individuals removed from alignments", length(toDrop)))
				if(!!length(toDrop)){                                      ###| creates a character vector containing the names of individuals
			#		toDrop.names    <- names(median.distAA.mdt0)[toDrop]   ###| that should be dropped from the alignment because they are highly
					toDrop.names    <- names(median.distAA.mdt1)[toDrop]   ###| that should be dropped from the alignment because they are highly
				} else {                                                   ###| diverged from most other individuals
					toDrop.names <- NULL                                   ###|
				}                                                          ###|
				### Drops individuals with unusually high amount of AA divergence unless all individuals dropped
				if(!!length(toDrop)){
					if(length(toDrop)==length(aa.alignment.temp2)){
						next
					}
					aa.alignment.temp3 <- aa.alignment.temp2[-toDrop]                    ###| sites-unfiltered AA alignment with highly diverged individuals removed
					aa.alignment.temp3 <- Biostrings::AAStringSet(x=gsub("-","",aa.alignment.temp3)) ###|
					names.aa.temp3     <- names(aa.alignment.temp3)                      ###|
					### rerun alignment algorithm
					print("Beginning MAFFT run 7/7")
					aa.alignment.temp3 <- REEs::mafft(REEs::trimXN(aa.alignment.temp3), param = mafft.params2)
				} else {                                               ###| if no individuals needed to be removed, then aa.alignment.temp3 is the same as aa.alignment.temp2
					aa.alignment.temp3 <- aa.alignment.temp2           ###|
				}                                                      ###|
				if(!length(aa.alignment.temp3)){  ##|Skips to next locus if no sites are retained
					next                          ##|
				}                                 ##|
			}
			if(!steps[7]){
				print("Saving protein alignment")
				#### Writes target AA-only alignment
				alignments.made[locus.name.temp,6] <- "yes"
				Biostrings::writeXStringSet(x = aa.alignment.temp3, filepath=paste0(dir12,locus.name.temp,".fa"), append=FALSE,compress=FALSE, compression_level=NA, format="fasta")
				ALIGNMENTSi[[6]]=aa.alignment.temp3
				
				#### Updates CDS alignment to remove individuals with highly diverged AA sequences
				if(!!length(toDrop)){                          #|Drops individuals with unusually high amount of AA divergence
					cds.temp5 <- cds.temp2[-toDrop]            #|The difference between cds.temp5 and cds.temp6 is that cds.temp6
					cds.temp6 <- cds.temp3[-toDrop]            #|may have dropped one or two bases from the CDS region to ensure that
					cds.temp7 <- REEs::filter.alignment(cds.temp5)   #|the region sequence length is a multiple of three (i.e., a sequence
				} else {                                       #|of codons)
					cds.temp5 <- cds.temp2                     #|cds.temp7 removes columns with only missing data if they exist
					cds.temp6 <- cds.temp3                     #|
					cds.temp7 <- REEs::filter.alignment(cds.temp5)   #|
				}                                              #|
				if(length(cds.temp6)<2){
					next
				}
				print("Extracting bases separately at first, second, and third codon positions")
				cds.p1 <- REEs::filter.alignment(cds.temp5,keep.specific=p1) ### alignment containing only first codon position sites (nodata sites dropped)
				cds.p2 <- REEs::filter.alignment(cds.temp5,keep.specific=p2) ### alignment containing only second codon position sites (nodata sites dropped)
				cds.p3 <- REEs::filter.alignment(cds.temp5,keep.specific=p3) ### alignment containing only third codon position sites (nodata sites dropped)
				
				print("Writing the first, second, and third codon position alignments to file")
				alignments.made[locus.name.temp,7:9] <- "yes"
				Biostrings::writeXStringSet(x = cds.p1, file=paste0(dir9,locus.name.temp,".fa"))
				Biostrings::writeXStringSet(x = cds.p2, file=paste0(dir10,locus.name.temp,".fa"))
				Biostrings::writeXStringSet(x = cds.p3, file=paste0(dir11,locus.name.temp,".fa"))
				
				ALIGNMENTSi[[7]]=cds.p1
				ALIGNMENTSi[[8]]=cds.p2
				ALIGNMENTSi[[9]]=cds.p3
				
				###########
				#### Checking if CDS frame preserved after dropping individuals with high AA divergence
				###########
			
				cds.p1.v2 <- REEs::filter.alignment(cds.temp7,keep.specific=p1)  ###|Like cds.p1,cds.p2,cds.p3, except that missing data columns
				cds.p2.v2 <- REEs::filter.alignment(cds.temp7,keep.specific=p2)  ###|were removed prior to extracting first, second, and third codon
				cds.p3.v2 <- REEs::filter.alignment(cds.temp7,keep.specific=p3)  ###|positions
				
				### The next steps write DNA alignments after dropping individuals with unusually high AA divergence,
				### but only if first, second, and third codon positions are exactly the same when no data columns are 
				### removed before or after extracting each codon position. This is a sanity check before writing
				### alignments.
				
				print("Saving CDS-only alignment")
				if(all(cds.p1==cds.p1.v2) & all(cds.p2==cds.p2.v2) & all(cds.p3==cds.p3.v2)){
					if(length(cds.temp7)>3){
						#ape::write.dna(x= cds.temp7,file=paste0(dir3,locus.name.temp,".phy"),format="sequential",append=F,nbcol=1,colw=100000000)
						Biostrings::writeXStringSet(x= cds.temp7,file=paste0(dir3,locus.name.temp,".fa"))
						alignments.made[locus.name.temp,2]  <- "yes"
						ALIGNMENTSi[[2]]=cds.temp7
						cds.starts            <- c(1:3)
						cds.ends              <- rep(width(cds.temp7[1]),3)
						cds.ranges            <- cbind(cds.starts,cds.ends)
						cds.partition.line    <- NULL
						for(j in 1:nrow(cds.ranges)){
							cds.partition.line[j] <- paste0("DNA, ", rownames(cds.ranges)[j]," = ",cds.ranges[j,1],"-",cds.ranges[j,2],"\\3")
						}
						print("Saving partition file for CDS-only alignment")
						write(cds.partition.line,file=paste0(dir4,locus.name.temp,"_parts.txt"))
						PARTSi[[2]] <- cds.partition.line
					} else {
						partition.groups[2] <- F
					}
				} else {
					partition.groups[2] <- F
				}
			}
		}
		if(!steps[8]){
			### Writing upstream noncoding alignment
			if(partition.groups[1]){
				if(!!length(toDrop.names)){
					if(any(toDrop.names %in% names(upstream.alignment2))){
						toDrop.names.upstream <- toDrop.names[toDrop.names %in% names(upstream.alignment2)]
						upstream.alignment3   <- upstream.alignment2[-match(toDrop.names.upstream,names(upstream.alignment2))]
						upstream.alignment3   <- REEs::filter.alignment(upstream.alignment3)
					} else {
						### columns are required to have some non-gap, non-ambiguous data
						upstream.alignment3 <- REEs::filter.alignment(upstream.alignment2)
					}
				}	else {
					### columns are required to have some non-gap, non-ambiguous data
					upstream.alignment3 <- REEs::filter.alignment(upstream.alignment2)
				}
				if(length(upstream.alignment3)>3){
					#ape::write.dna(x = upstream.alignment3,file=paste0(dir5,locus.name.temp,".phy"),format="sequential",append=F,nbcol=1,colw=100000000)
					print("Saving upstream non-coding DNA alignment")
					Biostrings::writeXStringSet(x = upstream.alignment3,file=paste0(dir5,locus.name.temp,".fa"))
					alignments.made[locus.name.temp,3]  <- "yes"
					ALIGNMENTSi[[3]]=upstream.alignment3
				} else {
					partition.groups[1] <- F
				}
			}
		}
		if(!steps[9]){
			### Writing downstream noncoding alignment
			if(partition.groups[3]){
				if(!!length(toDrop.names)){
					if(any(toDrop.names %in% names(downstream.alignment2))){
						toDrop.names.downstream <- toDrop.names[toDrop.names %in% names(downstream.alignment2)]
						downstream.alignment3   <- downstream.alignment2[-match(toDrop.names.downstream,names(downstream.alignment2))]
						downstream.alignment3   <- REEs::filter.alignment(downstream.alignment3)
					} else {
						downstream.alignment3   <- REEs::filter.alignment(downstream.alignment2)
					}
				}	else {
					downstream.alignment3 <- REEs::filter.alignment(downstream.alignment2)
				}
				if(length(downstream.alignment3)>3){
					#ape::write.dna(x = downstream.alignment3,file=paste0(dir6,locus.name.temp,".phy"),format="sequential",append=F,nbcol=1,colw=100000000)
					print("Saving downstream non-coding DNA alignment")
					Biostrings::writeXStringSet(x = downstream.alignment3,file=paste0(dir6,locus.name.temp,".fa"))
					alignments.made[locus.name.temp,4]  <- "yes"
					ALIGNMENTSi[[4]]=downstream.alignment3
				} else {
					partition.groups[3] <- F
				}
			}
		}
		if(!steps[10]){
			### Writing all-noncoding alignment (upstream and downstream non-coding regions concatenated)
			if(partition.groups[1] & partition.groups[3]){
				upstream.dataTable          <- data.table::data.table(as.matrix(upstream.alignment3),keep.rownames = TRUE)
				downstream.dataTable        <- data.table::data.table(as.matrix(downstream.alignment3),keep.rownames = TRUE)
				noncoding.alignment3        <- merge(upstream.dataTable, downstream.dataTable, by="rn", all=TRUE)
				noncoding.alignment3        <- REEs::na.replace(noncoding.alignment3,"-")
				rn.temp                     <- noncoding.alignment3$rn
				noncoding.alignment3        <- apply(noncoding.alignment3[ ,!"rn"], 1, paste, collapse="")
				noncoding.alignment3        <- Biostrings::DNAStringSet(noncoding.alignment3)
				names(noncoding.alignment3) <- rn.temp
				#ape::write.dna(x = noncoding.alignment3,file=paste0(dir7,locus.name.temp,".phy"),format="sequential",append=F,nbcol=1,colw=100000000)
				print("Saving alignment with upstream and downstream non-coding DNA concatenated")
				Biostrings::writeXStringSet(x = noncoding.alignment3,file=paste0(dir7,locus.name.temp,".fa"))
				### updates the alignments.made matrix
				alignments.made[locus.name.temp,5]      <- "yes"
				ALIGNMENTSi[[5]]=noncoding.alignment3
				### preparing to write the partition file for the locus
				noncoding.ranges          <- matrix(nrow=2,ncol=2)
				noncoding.ranges[1,1]     <- 1
				noncoding.ranges[1,2]     <- width(upstream.alignment3)[1]
				noncoding.ranges[2,1]     <- (noncoding.ranges[1,2]+1)
				noncoding.ranges[2,2]     <- width(noncoding.alignment3)[1]
				noncoding.partition.line  <- NULL
				for(j in 1:2){
					noncoding.partition.line[j] <- paste0("DNA, ", rownames(noncoding.ranges)[j]," = ",noncoding.ranges[j,1],"-",noncoding.ranges[j,2])
				}
				print("Saving partition file for noncoding DNA alignment (upstream and downstream non-coding regions concatenated)")
				write(noncoding.partition.line,file=paste0(dir8,locus.name.temp,"_parts.txt"))
				PARTSi[[5]] <- noncoding.partition.line
			}
		}
		if(!steps[11]){
			### Writing the all-data DNA alignment (which contains upstream non-coding, CDS, and downstream non-coding regions)
			if(all(partition.groups)){
				if(partition.groups[1]){
					upstream.dataTable.all    <- data.table::data.table(as.matrix(upstream.alignment3),keep.rownames = TRUE)
				} else {
					upstream.dataTable.all    <- NULL
				}
				if(partition.groups[2]){
					cds.datatable.all         <- data.table::data.table(as.matrix(cds.temp7),keep.rownames = TRUE)
				} else {
					cds.datatable.all         <- NULL
				}
				if(partition.groups[3]){
					downstream.dataTable.all  <- data.table::data.table(as.matrix(downstream.alignment3),keep.rownames = TRUE)
				} else {
					downstream.dataTable.all  <- NULL
				}
				dat.list         <- list(upstream.dataTable.all,cds.datatable.all,downstream.dataTable.all)
				dat              <- dat.list[[which(partition.groups)[1]]]
				
	#			ranges2      <- matrix(nrow=5,ncol=3)
	#			ranges2[,3]  <- c(1,2,2,2,3)
				ranges2      <- matrix(nrow=3,ncol=3)
				ranges2[,3]  <- c(1,2,3)
				ranges2[which(partition.groups)[1],1] <- 1
				ranges2[which(partition.groups)[1],2] <- (ncol(dat)-1)
				
				row.index             <- 1
				
				if(length(which(partition.groups)>1)){
					for(k in which(partition.groups)[-1]){
						row.index <- row.index+1
						dat       <- merge(dat, dat.list[[k]],by="rn", all=TRUE)
						ranges2[k,1] <- (ranges2[which(partition.groups)[row.index-1],2]+1)
						ranges2[k,2] <- (ncol(dat)-1)
					}
				}
				all.alignment2        <- REEs::na.replace(dat,"-")
				rn.temp               <- all.alignment2$rn
				all.alignment2        <- apply(all.alignment2[ ,!"rn"], 1, paste, collapse="")
				all.alignment2        <- Biostrings::DNAStringSet(all.alignment2)
				names(all.alignment2) <- rn.temp
				if(!all(is.na(ranges2[2,c(1,2)]))){
					ranges2           <- rbind(ranges2[1:2,],c((ranges2[2,1]+1),(ranges2[2,2]),2),c((ranges2[2,1]+2),(ranges2[2,2]),2),ranges2[3,])
					# ranges2           <- rbind(ranges2[1:2,],c(NA,NA,2),c(NA,NA,2),ranges2[3,])
					# ranges2[3,1]      <- ranges2[2,1]+1
					# ranges2[4,1]      <- ranges2[2,1]+2
					# ranges2[c(3:4),2] <- ranges2[2,2]
				}
				extra.text2       <- c("","\\3","\\3","\\3","")
				ranges3           <- cbind(ranges2,extra.text2)
				rownames(ranges3) <- c("upstream.noncoding","CDS.1","CDS.2","CDS.3","downstream.noncoding")
				ranges3           <- ranges3[c(1,2,2,2,3) %in% which(partition.groups),]
				#ape::write.dna(x= all.alignment2,file=paste0(dir1,locus.name.temp,".phy"),format="sequential",append=F,nbcol=1,colw=100000000)
				print("Saving the all-DNA-data alignment")
				Biostrings::writeXStringSet(x= all.alignment2,file=paste0(dir1,locus.name.temp,".fa"))
				### updates the alignments.made matrix
				alignments.made[locus.name.temp,1]  <- "yes"
				ALIGNMENTSi[[1]]=all.alignment2
				### preparing to write the partition file
				#partition.line        <- NULL
				#for(j in 1:nrow(ranges3)){
				#	partition.line[j] <- paste("DNA, ", rownames(ranges3)[j]," = ",ranges3[j,1],"-",ranges3[j,2],ranges3[j,4],sep="")
				#}
				print("Saving partition file of all-DNA-data alignment")
				partition.line <- paste0("DNA, ", rownames(ranges3)," = ",ranges3[,1],"-",ranges3[,2],ranges3[,4])
				write(partition.line,file=paste0(dir2,locus.name.temp,"_parts.txt"))
				PARTSi[[1]] <- partition.line
			}
		}
		if(!steps[12]) {
			print("Saving table summarizing which alignments were written")
			write.table(x=alignments.made,file=file.path(output.dir,"partitioned_alignments_made.txt"),sep="\t",quote=F)
		}
		
	} # End i for loop
	print("Done")
	# Returns a list with alignments and partitions of the last locus
	RESULTSi <- list(ALIGNMENTSi,PARTSi)
	RESULTSi
}
#' @examples
#' input.alignments.directory <- "~/Immune/unpartitioned/"
#' output.directory           <- "~/Immune/partitioned/"
#' target.loci                <- "~/Weinell_TargetLoci_Snakes_Final_targetCDS_v3.fa"
#'
#' make.partitioned.alignment(input.path=input.alignments.directory, output.dir=output.directory, TargetCDS.path=target.loci)

#' @title trimXN
#' 
#' Remove or replace 5' and 3' polyNs and polyXs
#' 
#' @param xstrings Object of class XStringSet (Biostrings package; DNAStringSet, AAStringSet, or RNAStringSet)
#' @param repl Either an empty string (default; polyNs and polyXs), or a single character that will replace each N or X of 5' and 3' polyNs and polyXs
#' @return XStringSet object with same class as object supplied to 'xstrings', with leading and trailing polyNs/polyXs removed or replaced depending on the value supplied to 'repl'
#' @export trimXN
trimXN <- function(xstrings,repl=''){
	pos   <- stringr::str_locate_all(xstrings,"^[N,-]+|[N,-]+$|^[X,-]+|[X,-]+$")
	if(all(!lengths(pos))){
		return(xstrings)
	}
	strings <- xstrings
	if(repl==""){
		return(as(gsub("^[N,-]+|[N,-]+$|^[X,-]+|[X,-]+$","",strings),class(xstrings)))
	}
	pos2    <- lapply(1:length(pos),function(i){cbind(pos[[i]],repl=apply(pos[[i]],1,function(x){paste(rep(repl,length(x[1]:x[2])),collapse="")}))})
	for(i in which(lengths(pos)>0)){
		for(j in 1:nrow(pos[[i]])){
			stringr::str_sub(strings[i],pos[[i]][j,1],pos[[i]][j,2]) <- pos2[[i]][j,3]
		}
	}
	as(strings,class(xstrings))
}

#' @title Trim alignment to covered by spanned by indidividual
#' 
#' Trims an alignment to the largest window containing all non-gap characters of a subset of individuals. See examples.
#' 
#' @param aln Object of class XStringSet (Biostrings package; DNAStringSet, AAStringSet, or RNAStringSet)
#' @param nam Character string vector with one or more names of sequences in aln
#' @return DNAStringSet
#' @export trimTo
trimTo <- function(aln, nam){
	if(length(intersect(nam, names(aln)))>0){
		firstbase      <- as.character(subseq(Biostrings::DNAStringSet(x=gsub("-","",aln)),start=1, end=1))
		lastbase       <- as.character(subseq(Biostrings::reverse(Biostrings::DNAStringSet(x=gsub("-","",aln))), start=1, end=1))
		firstbases_pos <- sapply(1:length(aln),function(x){stringr::str_locate(string=aln[x],pattern=firstbase[x])[1]})
		lastbases_pos  <- sapply(1:length(aln),function(x){REEs::str_locate_last(string=aln[x],pattern=lastbase[x])})
		aln_out        <- subseq(aln,start=min(firstbases_pos[names(aln) %in% nam]),end=max(lastbases_pos[names(aln) %in% nam]))
		
	} else {
		aln_out <- NA
	}
	aln_out
}
#' @examples
#' 
#' Condsider the alignment:
#' > alignment
#'     width  seq                       names
#' [1]    25  -----ACGTACGTAC-ACTG-----  Seq1
#' [2]    25  ACGTAACGTACGTAC-AC-------  Seq2
#' [3]    25  CGTACA--TACGTAC-ACTGACTGN  Seq3
#' 
#' # Trimming alignment to Seq1:
#' > trimTo(aln=alignment,c("Seq1"))
#'     width  seq             names
#' [1]    15  ACGTACGTAC-ACTG  Seq1
#' [2]    15  ACGTACGTAC-AC--  Seq2
#' [3]    15  A--TACGTAC-ACTG  Seq3
#' 
#' # Trimming alignment to Seq2:
#' > trimTo(aln=alignment,c("Seq2"))
#'     width  seq                names
#' [1]    18  -----ACGTACGTAC-AC  Seq1
#' [2]    18  ACGTAACGTACGTAC-AC  Seq2
#' [3]    18  CGTACA--TACGTAC-AC  Seq3
#' 
#' # Trimming alignment to Seq1 and Seq2:
#' > trimTo(aln=alignment,c("Seq1","Seq2"))
#'     width  seq                  names
#' [1]    20  -----ACGTACGTAC-ACTG  Seq1
#' [2]    20  ACGTAACGTACGTAC-AC--  Seq2
#' [3]    20  CGTACA--TACGTAC-ACTG  Seq3


#' @title submat
#' 
#' Return a substitution matrix for AA or DNA sequences to be used with the Biostrings::stringDist function.
#' Ambiguous nucleotides or amino acids are treated as missing data and are ignored.
#' 
#' @param type One of: "DNA","AA", "DNA_comparable", "AA_comparable"
#' @return XStringSet object with same class as object supplied to 'xstrings', with leading and trailing polyNs/polyXs removed or replaced depending on the value supplied to 'repl'
#' @export submat
submat <- function(type="DNA") {
	if(type=="DNA"){
		CODE_MAP <- c(names(IUPAC_CODE_MAP),"-")
		submat <- matrix(data=(-1),nrow=16,ncol=16,dimnames=list(CODE_MAP,CODE_MAP))
		submat[!(lower.tri(submat) | upper.tri(submat))] <- 0
		submat[5:16,] <- submat[,5:16] <- 0
		return(submat)
	}
	if(type=="AA"){
		CODE_MAP <- c("*",names(AMINO_ACID_CODE),"-")
		submat <- matrix(data=(-1),nrow=28,ncol=28,dimnames=list(CODE_MAP,CODE_MAP))
		submat[!(lower.tri(submat) | upper.tri(submat))] <- 0
		submat[27:28,] <- submat[,27:28] <- 0
		return(submat)
	}
	if(type=="DNA_comparable"){
		CODE_MAP <- c(names(IUPAC_CODE_MAP),"-")
		submat <- matrix(data=(-1),nrow=16,ncol=16,dimnames=list(CODE_MAP,CODE_MAP))
		#submat[!(lower.tri(submat) | upper.tri(submat))] <- 0
		submat[5:16,] <- submat[,5:16] <- 0
		return(submat)
	}
	if(type=="AA_comparable"){
		CODE_MAP <- c("*",names(AMINO_ACID_CODE),"-")
		submat <- matrix(data=(-1),nrow=28,ncol=28,dimnames=list(CODE_MAP,CODE_MAP))
		#submat[!(lower.tri(submat) | upper.tri(submat))] <- 0
		submat[27:28,] <- submat[,27:28] <- 0
		return(submat)
	}
}
