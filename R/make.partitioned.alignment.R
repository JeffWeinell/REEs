#' Make Partitioned Alignments
#' 
#' Given an input DNA alignment and a reference CDS sequence for the locus in the input alignment, this function generates up to nine output alignments, including three partitioned DNA alignments, five DNA alignments containing a subset of the input DNA sites, and one amino acid alignment for the CDS region(s) of the input alignment.
#' 
#' @param InputAlignmentFolder Folder containing fasta-formatted input DNA alignments.
#' @param output.dir Where to save alignments (alignments for each type of data will be in different subdirectory)
#' @param TargetCDS.path Full path to the fasta file containing only the CDS sequences of target loci (from which probes were designed). Sequence names must have the following format: "GeneName_TargetCDS_of_TargetLocusName_AnyAdditionalIformation", where GeneName and TargetLocusName are replaced with the actual names, and AnyAdditionalIformation can be a string of any characters
#' @param bait.species.filename Full path to the file containing a two column table in which the first column contains name of locus and the second column contains name of species that the probes for that locus were designed from. This is needed so that the the function knows which genome to use as a reference.
#' I will update this function to remove the need for this argument. Write now this table is just used for renaming the reference sequence.
#' @param ref.type Indicates that the input alignments are DNA sequences (default is "DNA").
#' @param old.names Character vector of old sample names that should be changed to those defined in new.names parameter (default is NA, meaning don't change names).
#' @param new.names Character vector of new sample names to use instead of corresponding name in old.names vector (default is NA, meaning don't change names).
#' @param drop.reference Whether or not to include the reference sequence (from which baits were designed for a particular locus) in the alignment written to file (default is FALSE).
#' @param ith.locus.start First locus in InputAlignmentFolder to process (default is 1).
#' @param ith.locus.end Last locus in InputAlignmentFolder to process (default is "all", when means process until no more loci to process).
#' @param locus.names.omit Names of loci to skip (default is NULL).
#' @param AA.pdist.drop.thresh Maximum mean pairwise p-distance for AA sequence of an individual (translated from CDS regions of input DNA sequences) to keep the individual in AA or CDS-containing DNA alignments.
#' @param trimto Optional character string with names of sequences in input alignments whose combined range should be used to trim the alignments. Default is NULL.
#' @return Writes up to nine different fasta-formatted sequence alignments for each input DNA alignment in InputAlignmentFolder, and partition files partitioned alignments. The nine output alignment types are:
#' (1) Partitioned version of the input DNA alignment. The alignment is partitioned into (if present): upstream noncoding DNA, each codon position (within CDS regions), and downstream noncoding DNA.
#' (2) Alignment containing only the CDS regions, which are partitioned by codon position.
#' (3) Alignment containing only the upstream noncoding DNA.
#' (4) Alignment containing only the downstream noncoding DNA.
#' (5) Alignment containing only noncoding DNA, partitioned into upstream and downstream regions.
#' (6) Alignment containing only the first codon positions of CDS regions.
#' (7) Alignment containing only the second codon positions of CDS regions.
#' (8) Alignment containing only the third codon positions of CDS regions.
#' (9) Alignment containing the amino acid sequence for the translated CDS region.
#'@export
make.partitioned.alignment  <- function(InputAlignmentFolder,output.dir,TargetCDS.path,bait.species.filename,ref.type="DNA",old.names=NA,new.names=NA,drop.reference=F,ith.locus.start=1,ith.locus.end="all",locus.names.omit=NULL,AA.pdist.drop.thresh=0.5,trimto=NULL){
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
	makeDirectories <- lapply(X=subdirectories,FUN=dir.check.create)  ### need to assign an object to avoid printing results of lapply 
	
	### Reads in files ###
	TargetDNA_CDS.regions       <- Biostrings::readDNAStringSet(TargetCDS.path,format="fasta")
	### reads the table specifying which species the probes were designed from for each locus 
	bait.species.table          <- data.table::fread(bait.species.filename)

	name.string.start           <- (str_locate_X(strings=names(TargetDNA_CDS.regions),pattern="_",X=3)+1)
	name.string.end             <- (str_locate_X(strings=names(TargetDNA_CDS.regions),pattern="_",X=4)-1)
	targetCDS.names             <- substring(names(TargetDNA_CDS.regions),first=name.string.start,last=name.string.end) ### target loci names of each sequence in TargetDNA_CDS.regions
	input.alignment.filenames   <- gtools::mixedsort(list.files(InputAlignmentFolder,full.names=T))
	

	### Next line assumes that files are named according to WeinelEntry names
	input.alignment.shortnames  <- gsub(".phy|.fa|.fasta","",basename(input.alignment.filenames))

	### Names of loci that have been aligned and that have a CDS region included in TargetDNA_CDS.regions
	shared.names                <- Reduce(intersect, list(input.alignment.shortnames,targetCDS.names))
	### names of loci to ignore
	if(!is.null(locus.names.omit)){
		shared.names  <- setdiff(shared.names,locus.names.omit)
	}
	
	if(file.exists(file.path(output.dir,"partitioned_alignments_made.txt"))){
		alignments.made <- utils::read.table(file=file.path(output.dir,"partitioned_alignments_made.txt"), header=T, colClasses="character")
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
	trimTo <- function(aln,nam){
		if(length(intersect(nam, names(aln)))>0){
			firstbase        <- as.character(subseq(Biostrings::DNAStringSet(x=gsub("-","",aln)),start=1,end=1))
			lastbase         <- as.character(subseq(Biostrings::reverse(Biostrings::DNAStringSet(x=gsub("-","",aln))),start=1,end=1))
			firstbases_pos   <- sapply(1:length(aln),function(x){stringr::str_locate(string=aln[x],pattern=firstbase[x])[1]})
			lastbases_pos    <- sapply(1:length(aln),function(x){REEs::str_locate_last(string=aln[x],pattern=lastbase[x])})
			aln_out        <- subseq(aln,start=min(firstbases_pos[names(aln) %in% nam]),end=max(lastbases_pos[names(aln) %in% nam]))
			
		} else {
			aln_out <- NA
		}
		aln_out
	}
	for(i in ith.locus.start:last.locus.process){
		#if(is.wholenumber(i/50)){
		#	print(i)
		#}
		print(i)
		locus.name.temp   <- shared.names[i]
		bait.species.temp <- bait.species.table$Species[which(bait.species.table$Bait==locus.name.temp)]
		### un-aligning sequences in the "novel" alignment
		novel             <- Biostrings::readDNAMultipleAlignment(input.alignment.filenames[which(input.alignment.shortnames==locus.name.temp)])
		novel2            <- Biostrings::DNAStringSet(x=gsub("-|\\?","",novel))
		if(ref.type=="DNA"){
			reference.cds        <- TargetDNA_CDS.regions[which(targetCDS.names==locus.name.temp)]
			names(reference.cds) <- paste(bait.species.temp,names(reference.cds),sep="_")
			final.locus          <- c(reference.cds,novel2)
		}
		# Skips locus if <4 sequences other than the reference CDS
		#if(FALSE) {
		#	if(length(final.locus)<5){
		#		next
		#	}
		#}
		alignment.names  <- names(final.locus)
		if(!all(is.na(old.names))){
			alignment.names  <- mgsub(old.names,new.names,alignment.names)
		}
		### Calling MAFFT
	#	alignment   <- REEs::mafft(final.locus,param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6")
		alignment   <- REEs::mafft(final.locus, param="--auto --adjustdirection --nwildcard --op 3 --ep 0.123 --quiet")
		alignment_A <- alignment
		### Trim alignment to the maximum width region that includes at least one individual in 'trimto' at the first and last base
		if(!is.null(trimto)){
			alignment   <- trimTo(aln=alignment,nam=trimto)
			alignment_B <- alignment
			#if(!names(reference.cds) %in% trimto){
			#	trimto <- c(names(reference.cds), trimto)
			#}
			#if(length(intersect(trimto, names(alignment)))>0){
			#	firstbase        <- as.character(subseq(Biostrings::DNAStringSet(x=gsub("-","",alignment)),start=1,end=1))
			#	lastbase         <- as.character(subseq(Biostrings::reverse(Biostrings::DNAStringSet(x=gsub("-","",alignment))),start=1,end=1))
			#	firstbases_pos   <- sapply(1:length(alignment),function(x){stringr::str_locate(string=alignment[x],pattern=firstbase[x])[1]})
			#	lastbases_pos    <- sapply(1:length(alignment),function(x){REEs::str_locate_last(string=alignment[x],pattern=lastbase[x])})
			#	alignment        <- subseq(alignment,start=min(firstbases_pos[names(alignment) %in% trimto]),end=max(lastbases_pos[names(alignment) %in% trimto]))
			#	alignment_B      <- alignment
			#}
		} else {
			alignment_B <- NULL
		}
		begin.exon       <- stringr::str_locate(string=alignment[1],pattern=as.character(subseq(reference.cds,start=1,end=1)))[1]
		end.exon         <- REEs::str_locate_last(string=alignment[1],pattern=as.character(subseq(reference.cds,start=(width(reference.cds)-1),end=width(reference.cds))))
		starts           <- c(1,begin.exon,begin.exon+1,begin.exon+2,end.exon+1)
		ends             <- c(begin.exon-1,end.exon,end.exon,end.exon,width(alignment[1]))
		widths           <- ends-(starts-1)
		extra.text       <- c("","\\3","\\3","\\3","")
		ranges           <- cbind(starts,ends,widths,extra.text)
		rownames(ranges) <- c("upstream.noncoding","CDS.1","CDS.2","CDS.3","downstream.noncoding")
		
		partition.groups        <- c(widths[1]>0,all(widths[2:4]>0),widths[5]>0)        ### logical vector indicating if upstream noncoding, CDS, and downstream noncoding data are present
		names(partition.groups) <- c("upstream.noncoding","CDS","downstream.noncoding")
		
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
				cds.alignment2        <- Biostrings::DNAStringSet(x=gsub("-","",cds.alignment))
				# cds.alignment2      <- REEs::mafft(cds.alignment2,param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6")
				cds.alignment2        <- REEs::mafft(cds.alignment2, param="--auto --adjustdirection --nwildcard --op 3 --ep 0.123 --quiet")
				if(!is.null(trimto)){
					cds.alignment2    <- trimTo(aln=cds.alignment2,nam=trimto)
				}
				alignment_C           <- cds.alignment2
			} else {
				partition.groups[2] <- F   ### updates decision from TRUE to FALSE on whether to write the CDS alignment, because too few individuals with CDS data
				alignment_C <- NULL
			}
		}
		#### Upstream noncoding alignment
		if(partition.groups[1]){
			if(begin.exon!=1){
				if(drop.reference==T){
					upstream.alignment    <- Biostrings::subseq(x=alignment,start=1,end=(begin.exon-1))[-1]
				} else {
					upstream.alignment    <- Biostrings::subseq(x=alignment,start=1,end=(begin.exon-1))
				}
				drop.nodata               <- stringr::str_count(upstream.alignment,"-") == width(upstream.alignment)
				if(length(which(!drop.nodata))>3){
					upstream.alignment         <- upstream.alignment[!drop.nodata]
					upstream.alignment2        <- Biostrings::DNAStringSet(x=gsub("-","",upstream.alignment))
					# upstream.alignment2      <- REEs::mafft(upstream.alignment2,param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6")
					upstream.alignment2        <- REEs::mafft(upstream.alignment2, param="--auto --adjustdirection --nwildcard --op 3 --ep 0.123 --quiet")
					if(!is.null(trimto)){
						upstream.alignment2    <- trimTo(aln=upstream.alignment2,nam=trimto)
					}
					alignment_D <- upstream.alignment2
				} else {
					partition.groups[1] <- F   ### updates decision from TRUE to FALSE on whether to write the upstream alignment, because too few individuals with upstream noncoding data
					alignment_D <- NULL
				}
			}
		}
		#### Downstream noncoding alignment
		if(partition.groups[3]){
			if(end.exon!=width(alignment[1])){
				if(drop.reference==T){
					downstream.alignment    <- Biostrings::subseq(x=alignment,start=(end.exon+1),end=width(alignment[1]))[-1]
				} else {
					downstream.alignment    <- Biostrings::subseq(x=alignment,start=(end.exon+1),end=width(alignment[1]))
				}				
				drop.nodata                 <- stringr::str_count(downstream.alignment,"-") == width(downstream.alignment)
				### Skips downstream noncoding alignment if 3 or less individuals
				if(length(which(!drop.nodata))>3){
					downstream.alignment    <- downstream.alignment[!drop.nodata]
					length(downstream.alignment)
					downstream.alignment2   <- Biostrings::DNAStringSet(x=gsub("-","",downstream.alignment))
					#downstream.alignment2  <- REEs::mafft(downstream.alignment2, param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6")
					downstream.alignment2   <- REEs::mafft(downstream.alignment2, param="--auto --adjustdirection --nwildcard --op 3 --ep 0.123 --quiet")
					if(!is.null(trimto)){
						downstream.alignment2    <- trimTo(aln=downstream.alignment2,nam=trimto)
					}
					alignment_E <- downstream.alignment2
				} else {
					### Updates decision from TRUE to FALSE on whether to write the downstream alignment, because too few individuals with downstream noncoding data.
					partition.groups[3] <- F
					alignment_E <- NULL
				}
			}
		}
		##########
		##### Remainder of Code below was added to include parts of the code in the R script "make_codon_alignments.R"
		##########

		if(partition.groups[2]){
			cds.temp2   <- cds.alignment2
			p1 <- which(is.wholenumber(c(0:(width(cds.temp2[1])-1))/3)) ### numerical vector of first codon position sites
			p2 <- intersect((p1+1),c(2:width(cds.temp2[1])))            ### numerical vector of second codon position sites
			p3 <- intersect((p2+1),c(3:width(cds.temp2[1])))            ### numerical vector of third codon position sites
			### CDS sequences aligned (extracted from whole alignment) but not yet re-aligned
			cds.temp3 <- Biostrings::subseq(cds.temp2,start=min(p1),end=max(p3))
			### CDS sequences not aligned
			cds.temp4 <- Biostrings::DNAStringSet(x=gsub("-","",cds.temp3))   ### unaligns sequences (after removing gaps) so that the sequences can be translated 
			
			### Translates CDS sequences, then does MAFFT multiple sequence alignment to make an alignment of amino acids
			aa.temp             <- suppressWarnings(Biostrings::translate(cds.temp4,if.fuzzy.codon="solve",no.init.codon=T))
			#aa.alignment.temp  <- REEs::mafft(aa.temp, param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6")
			aa.alignment.temp   <- REEs::mafft(aa.temp, param="--auto --adjustdirection --nwildcard --op 3 --ep 0.123 --quiet")
			# Nothing filtered by default
			aa.alignment.temp2  <- filter.alignment(aa.alignment.temp)
			# This alignment has no missing or ambiguous data. This alignment is used for calculating median pairwise p-distances, but is not the alignment written to file.
			aa.alignment.mdt0   <- filter.alignment(aa.alignment.temp, mdt=0)
			
			### Median pairwise p-distances for the AA alignment that can contain any amount of missing or ambiguous data
			distAA.mdt1         <- Biostrings::stringDist(aa.alignment.temp2,upper=T)
			#median.distAA.mdt1  <- apply(X=distAA.mdt1,MARGIN=1,FUN=median)/width(aa.alignment.temp2[1])
			### Median pairwise p-distances for the AA alignment containing no missing or ambiguous data
			if(!all(width(aa.alignment.mdt0)==0)){
				distAA.mdt0        <- Biostrings::stringDist(aa.alignment.mdt0)
				median.distAA.mdt0 <- apply(X=distAA.mdt0,MARGIN=1,FUN=median)/width(aa.alignment.mdt0[1])
			} else {
				distAA.mdt0        <- NA
				median.distAA.mdt0 <- NA
			}
			toDrop <- which(median.distAA.mdt0 > AA.pdist.drop.thresh) ### individuals with median p-distance greater than AA.pdist.drop.thresh should be dropped
			if(!!length(toDrop)){                                      ###| creates a character vector containing the names of individuals
				toDrop.names    <- names(median.distAA.mdt0)[toDrop]   ###| that should be dropped from the alignment because they are highly
			} else {                                                   ###| diverged from most other individuals
				toDrop.names <- NULL                                   ###|
			}                                                          ###|
			if(!!length(toDrop)){                                          ###| Drops individuals with unusually high amount of AA divergence unless all individuals dropped
				if(length(toDrop)==length(aa.alignment.temp2)){
					next
				}
				aa.alignment.temp3 <- aa.alignment.temp2[-toDrop]                    ###| sites-unfiltered AA alignment with highly diverged individuals removed
				aa.alignment.temp3 <- Biostrings::AAStringSet(x=gsub("-","",aa.alignment.temp3)) ###|
				names.aa.temp3     <- names(aa.alignment.temp3)                      ###|
				### rerun alignment algorithm
				# aa.alignment.temp3 <- REEs::mafft(aa.alignment.temp3,param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6")
				aa.alignment.temp3 <- REEs::mafft(aa.alignment.temp3,param="--auto --adjustdirection --nwildcard --op 3 --ep 0.123 --quiet")
			} else {                                               ###| if no individuals needed to be removed, then aa.alignment.temp3 is the same as aa.alignment.temp2
				aa.alignment.temp3 <- aa.alignment.temp2           ###|
			}                                                      ###|
			if(!length(aa.alignment.temp3)){  ##|Skips to next locus if no sites are retained
				next                          ##|
			}                                 ##|
			
			#### Writes target AA-only alignment
			alignments.made[i,6] <- "yes"
			Biostrings::writeXStringSet(x = aa.alignment.temp3, filepath=paste0(dir12,locus.name.temp,".fa"), append=FALSE,compress=FALSE, compression_level=NA, format="fasta")
			
			#### Updates CDS alignment to remove individuals with highly diverged AA sequences
			if(!!length(toDrop)){                          #|Drops individuals with unusually high amount of AA divergence
				cds.temp5 <- cds.temp2[-toDrop]            #|The difference between cds.temp5 and cds.temp6 is that cds.temp6
				cds.temp6 <- cds.temp3[-toDrop]            #|may have dropped one or two bases from the CDS region to ensure that
				cds.temp7 <- filter.alignment(cds.temp5)   #|the region sequence length is a multiple of three (i.e., a sequence
			} else {                                       #|of codons)
				cds.temp5 <- cds.temp2                     #|cds.temp7 removes columns with only missing data if they exist
				cds.temp6 <- cds.temp3                     #|
				cds.temp7 <- filter.alignment(cds.temp5)   #|
			}                                              #|
			if(length(cds.temp6)<2){
				next
			}
			cds.p1 <- filter.alignment(cds.temp5,keep.specific=p1) ### makes an alignment containing only first codon position sites (nodata sites dropped)
			cds.p2 <- filter.alignment(cds.temp5,keep.specific=p2) ### makes an alignment containing only second codon position sites (nodata sites dropped)
			cds.p3 <- filter.alignment(cds.temp5,keep.specific=p3) ### makes an alignment containing only third codon position sites (nodata sites dropped)
			
			### Writes the first, second, and third codon position alignments to file
			alignments.made[i,7:9] <- "yes"
	##		ape::write.dna(x = cds.p1, file=paste0(dir9,locus.name.temp,".phy"),format="sequential",append=F,nbcol=1,colw=100000000)
	##		ape::write.dna(x = cds.p2, file=paste0(dir10,locus.name.temp,".phy"),format="sequential",append=F,nbcol=1,colw=100000000)
	##		ape::write.dna(x = cds.p3, file=paste0(dir11,locus.name.temp,".phy"),format="sequential",append=F,nbcol=1,colw=100000000)
			
			Biostrings::writeXStringSet(x = cds.p1, file=paste0(dir9,locus.name.temp,".fa"))
			Biostrings::writeXStringSet(x = cds.p2, file=paste0(dir10,locus.name.temp,".fa"))
			Biostrings::writeXStringSet(x = cds.p3, file=paste0(dir11,locus.name.temp,".fa"))

			###########
			#### Checking if CDS frame preserved after dropping individuals with high AA divergence
			###########
			
			cds.p1.v2 <- filter.alignment(cds.temp7,keep.specific=p1)  ###|Like cds.p1,cds.p2,cds.p3, except that missing data columns
			cds.p2.v2 <- filter.alignment(cds.temp7,keep.specific=p2)  ###|were removed prior to extracting first, second, and third codon
			cds.p3.v2 <- filter.alignment(cds.temp7,keep.specific=p3)  ###|positions
			
			### The next steps write DNA alignments after dropping individuals with unusually high AA divergence,
			### but only if first, second, and third codon positions are exactly the same when no data columns are 
			### removed before or after extracting each codon position. This is a sanity check before writing
			### alignments.
			
			### Writing CDS-only alignment
			if(all(cds.p1==cds.p1.v2) & all(cds.p2==cds.p2.v2) & all(cds.p3==cds.p3.v2)){
				if(length(cds.temp7)>3){
					#ape::write.dna(x= cds.temp7,file=paste0(dir3,locus.name.temp,".phy"),format="sequential",append=F,nbcol=1,colw=100000000)
					Biostrings::writeXStringSet(x= cds.temp7,file=paste0(dir3,locus.name.temp,".fa"))
					alignments.made[i,2]  <- "yes"
					cds.starts            <- c(1:3)
					cds.ends              <- rep(width(cds.temp7[1]),3)
					cds.ranges            <- cbind(cds.starts,cds.ends)
					cds.partition.line    <- NULL
					for(j in 1:nrow(cds.ranges)){
						cds.partition.line[j] <- paste0("DNA, ", rownames(cds.ranges)[j]," = ",cds.ranges[j,1],"-",cds.ranges[j,2],"\\3")
					}
					write(cds.partition.line,file=paste0(dir4,locus.name.temp,"_parts.txt"))
				} else {
					partition.groups[2] <- F
				}
			} else {
				partition.groups[2] <- F
			}
		}
		### Writing upstream noncoding alignment
		if(partition.groups[1]){
			if(!!length(toDrop.names)){
				if(any(toDrop.names %in% names(upstream.alignment2))){
					toDrop.names.upstream <- toDrop.names[toDrop.names %in% names(upstream.alignment2)]
					upstream.alignment3   <- upstream.alignment2[-match(toDrop.names.upstream,names(upstream.alignment2))]
					upstream.alignment3   <- filter.alignment(upstream.alignment3)
				} else {
					### columns are required to have some non-gap, non-ambiguous data
					upstream.alignment3 <- filter.alignment(upstream.alignment2)
				}
			}	else {
				### columns are required to have some non-gap, non-ambiguous data
				upstream.alignment3 <- filter.alignment(upstream.alignment2)
			}
			if(length(upstream.alignment3)>3){
				#ape::write.dna(x = upstream.alignment3,file=paste0(dir5,locus.name.temp,".phy"),format="sequential",append=F,nbcol=1,colw=100000000)
				Biostrings::writeXStringSet(x = upstream.alignment3,file=paste0(dir5,locus.name.temp,".fa"))
				alignments.made[i,3]  <- "yes"
			} else {
				partition.groups[1] <- F
			}
		}
		### Writing downstream noncoding alignment
		if(partition.groups[3]){
			if(!!length(toDrop.names)){
				if(any(toDrop.names %in% names(downstream.alignment2))){
					toDrop.names.downstream <- toDrop.names[toDrop.names %in% names(downstream.alignment2)]
					downstream.alignment3   <- downstream.alignment2[-match(toDrop.names.downstream,names(downstream.alignment2))]
					downstream.alignment3   <- filter.alignment(downstream.alignment3)
				} else {
					downstream.alignment3   <- filter.alignment(downstream.alignment2)
				}
			}	else {
				downstream.alignment3 <- filter.alignment(downstream.alignment2)
			}
			if(length(downstream.alignment3)>3){
				#ape::write.dna(x = downstream.alignment3,file=paste0(dir6,locus.name.temp,".phy"),format="sequential",append=F,nbcol=1,colw=100000000)
				Biostrings::writeXStringSet(x = downstream.alignment3,file=paste0(dir6,locus.name.temp,".fa"))
				alignments.made[i,4]  <- "yes"
			} else {
				partition.groups[3] <- F
			}
		}
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
			Biostrings::writeXStringSet(x = noncoding.alignment3,file=paste0(dir7,locus.name.temp,".fa"))
			### updates the alignments.made matrix
			alignments.made[i,5]      <- "yes"
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
			write(noncoding.partition.line,file=paste0(dir8,locus.name.temp,"_parts.txt"))
		}
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
			
			ranges2      <- matrix(nrow=5,ncol=3)
			ranges2[,3]  <- c(1,2,2,2,3)
			ranges2[which(partition.groups)[1],1] <- 1
			ranges2[which(partition.groups)[1],2] <- ncol(dat)
			
			row.index             <- 1
			
			if(length(which(partition.groups)>1)){
				for(k in which(partition.groups)[-1]){
					row.index <- row.index+1
					dat       <- merge(dat, dat.list[[k]],by="rn", all=TRUE)
					ranges2[k,1] <- (ranges2[which(partition.groups)[row.index-1],2]+1)
					ranges2[k,2] <- ncol(dat)
					
				}
			}
			all.alignment2        <- REEs::na.replace(dat,"-")
			rn.temp               <- all.alignment2$rn
			all.alignment2        <- apply(all.alignment2[ ,!"rn"], 1, paste, collapse="")
			all.alignment2        <- Biostrings::DNAStringSet(all.alignment2)
			names(all.alignment2) <- rn.temp
			
			if(!all(is.na(ranges2[2,c(1,2)]))){
				ranges2[3,1]      <- ranges2[2,1]+1
				ranges2[4,1]      <- ranges2[2,1]+2
				ranges2[c(3:4),2] <- ranges2[2,2]
			}
			
			extra.text2       <- c("","\\3","\\3","\\3","")
			ranges3           <- cbind(ranges2,extra.text2)
			rownames(ranges3) <- c("upstream.noncoding","CDS.1","CDS.2","CDS.3","downstream.noncoding")
			ranges3           <- ranges3[c(1,2,2,2,3) %in% which(partition.groups),]
			#ape::write.dna(x= all.alignment2,file=paste0(dir1,locus.name.temp,".phy"),format="sequential",append=F,nbcol=1,colw=100000000)
			Biostrings::writeXStringSet(x= all.alignment2,file=paste0(dir1,locus.name.temp,".fa"))
			### updates the alignments.made matrix
			alignments.made[i,1]  <- "yes"
			### preparing to write the partition file
			partition.line        <- NULL
			for(j in 1:nrow(ranges3)){
				partition.line[j] <- paste("DNA, ", rownames(ranges3)[j]," = ",ranges3[j,1],"-",ranges3[j,2],ranges3[j,4],sep="")
			}
			write(partition.line,file=paste0(dir2,locus.name.temp,"_parts.txt"))
		}
		#### Writes table summarizing which alignments were written
		write.table(x=alignments.made,file=file.path(output.dir,"partitioned_alignments_made.txt"),sep="\t",quote=F)
	} # End i for loop
}
#' @examples
#' input.alignments.directory <- "~/Immune/unpartitioned/"
#' output.directory           <- "~/Immune/partitioned/"
#' target.loci                <- "~/Weinell_TargetLoci_Snakes_Final_targetCDS_v3.fa"
#' input.baits.table          <- "~/bait_species_table.txt"
#'
#' make.partitioned.alignment(InputAlignmentFolder=input.alignments.directory, output.dir=output.directory, TargetCDS.path=target.loci, bait.species.filename=input.baits.table)

