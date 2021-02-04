#' Get CDS Contained in Range
#'
#' Obtains CDS portion of a sequence from a genbank flatfile.
#' Also obtains AA sequences (target and full sequence).
#' For some loci, the start codon is poorly annotated in ncbi, but this function determines the start codon by comparing the translated CDS to the NCBI AA sequence; makes lists of the target loci with unusual or absent annotations.
#'
#'  @param  gb.object An object of class gbRecord
#'  @param  targetTable.object Input table containing target loci coordinates
#'  @param  output.filename Output filename
#'  @param  report.noCDS Logical indicating if the output object should also include two sets of amino acid sequences and information about missing/short/duplicate annotations. Default TRUE. If FALSE, the output object only includes the CDS regions contained within target loci.
#'  @return An object of class DNAStringSet if report.noCDS=FALSE, otherwise an object of class list (length = 10) that holds one DNAStringSet object, two AAStringSet objects, and seven character vectors.
#'    When report.noCDS = TRUE, the result is a list containing the following:
#'    A DNAStringSet object containing CDS sequences contained within (i.e. overlapping) the regions specified in targetTable.object,
#'    An AAStringSet object containing AA sequences coded by the CDS sequences in the DNAStringSet,
#'    An AAStringSet object containing the full protein sequences for proteins represented in the DNAStringSet,
#'    A vector of names of target loci for which another target has the same genomic coordinates,
#'    A vector of names of target loci that do not overlap with sequences in the Genbank flatfile,
#'    A vector of names of target loci without CDS annotation(s) in the Genbank record,
#'    A vector of names of target loci for which Genbank record includes a very short CDS annotation (< 12 bp). These short CDS and AA sequences were not written to file,
#'    A vector of names of target loci for which multiple translation frames of target CDS are found in full AA sequence,
#'    A vector of names of target loci for which no translation frames of target CDS are found in full AA sequence,
#'    A vector of short names for the 10 elements of the output: "CDS StringSet","Target AA StringSet","Full AA StringSet","duplicate target names","targets not in Genbank flatfile","targets without CDS annotation","targets with short CDS annotation","multiple translation frames","no matching translation frames","result names".
#'  @export 
CDS.from.gb            <- function(gb.object,targetTable.object,output.filename=NA,report.noCDS=TRUE){
	
	translate               <- Biostrings::translate ### needed to because multiple libraries have a translate function
	gbData                  <- gb.object
	loci.gb                 <- gsub("\\.\\.","-",names(gbData))  ### GenBank accession ID plus start and end region
	targetTable             <- targetTable.object
	accession.list          <- gsub("\\..","",targetTable$TargetName_NCBI)
	start.list              <- targetTable$TargetStart
	end.list                <- targetTable$TargetEnd
	negative.sense          <- which(targetTable$Sense_Target_ReferenceContig==2)
	targetTable.identifier1 <- paste(accession.list,":",start.list,"-",end.list, sep="")
	targetTable.identifier1[negative.sense] <- paste(accession.list[negative.sense],":complement(",start.list[negative.sense],"-",end.list[negative.sense],")", sep="")
	keep.loci1              <- match(unique(targetTable.identifier1),targetTable.identifier1)   ### a list of loci to keep (i.e., those that are not duplicates)
	duplicates              <- names(which(table(targetTable.identifier1)>1))
	duplicate.targets       <- targetTable$TargetName_ArborSci[setdiff(which(targetTable.identifier1 %in% duplicates),match(duplicates,targetTable.identifier1))]
	targetTable.reduced1    <- targetTable[keep.loci1,]
	targetTable.identifier2 <- targetTable.identifier1[keep.loci1]
	keep.loci2              <- which(targetTable.identifier2 %in% loci.gb)  ### only processes loci that exist in the flatfile
	targetTable.reduced     <- targetTable.reduced1[keep.loci2,]            ### loci found in genbank flatfile
	targetTable.remaining   <- targetTable.reduced1[-keep.loci2,]           ### loci not found in genbank flatfile
	targetNames.remaining   <- targetTable.remaining$TargetName_ArborSci    ### the unique names used for target loci (ie, WeinellEntry names) that were not represented in the GenBank flatfile
	targetTable.identifier  <- targetTable.identifier2[keep.loci2]          ### targetTable identifiers that are found in the gb flatfile
	targetNames             <- targetTable.reduced$TargetName_ArborSci      ### the unique names used for target loci (ie, WeinellEntry names)

	CDS.test                <- cbind(targetNames,rep("-",length(targetNames))) ### a 2-column matrix. first column = targetNames, and the second column is empty but will be filled with "yes" or "no" depending on whether or not a CDS feature (with AA sequence) exists
	
	for(i in 1:length(loci.gb)){
		## print every 50th i to keep progress
		if(is.wholenumber(i/50)){
			print(i)
		}
		cdsFeatures.temp <- gbData[[i]]["CDS"]  ### current CDS feature
		identifier.temp  <- which(targetTable.identifier==loci.gb[i])
		targetName.temp  <- targetNames[identifier.temp]
		locus.temp       <- paste(targetName.temp,"_",loci.gb[i],sep="")
		if(length(cdsFeatures.temp)!=0 & length(grep(pattern="translation",cdsFeatures.temp))!=0){
			CDS.test[identifier.temp,2] <- "yes"
			CDS_product.temp      <- biofiles::geneID(cdsFeatures.temp)                          #### Name of the CDS region
			if(all(is.na(CDS_product.temp))){
				CDS_product.temp  <- biofiles::proteinID(cdsFeatures.temp)
			}
			CDS_product.temp      <- gsub(", transcript variant.*","",CDS_product.temp)  #### Just removes this part of the name
			if(all(is.na(CDS_product.temp)) | length(CDS_product.temp)==0){
				CDS_product.temp  <- "MissingGeneName"
			}
			dna.temp          <- biofiles::getSequence(gbData[[i]])   ### entire target sequence
			clip              <- biofiles::ranges(cdsFeatures.temp)   ### GRanges object containing the ranges of the CDS sequence
			CDS.temp          <- dna.temp[clip]

			CDS.temp          <- unique(CDS.temp)                                            #### removes duplicate annotation if exists
			CDS.temp          <- CDS.temp[which(width(CDS.temp)==max(width(CDS.temp)))]
			CDS.temp          <- CDS.temp[1]
			CDS_product.temp  <- CDS_product.temp[which(width(CDS.temp)==max(width(CDS.temp)))][1]
			names(CDS.temp)   <- paste(CDS_product.temp,"_TargetCDS_of_",locus.temp,sep="")

			AA.temp           <- biofiles::translation(cdsFeatures.temp)
			names(AA.temp)    <- paste(CDS_product.temp,"_TargetAA_of_",locus.temp,sep="")
			AA.temp           <- unique(AA.temp)                                             #### removes duplicate annotation
			
			AA.temp           <- AA.temp[which(width(AA.temp)==max(width(AA.temp)))]
			AA.temp           <- AA.temp[1]
			
			cds.f1            <- CDS.temp
			cds.f2            <- XVector::subseq(x=CDS.temp,start=2)
			cds.f3            <- XVector::subseq(x=CDS.temp,start=3)
			cds.RC.f1         <- Biostrings::reverseComplement(CDS.temp)
			cds.RC.f2         <- XVector::subseq(x=reverseComplement(CDS.temp),start=2)
			cds.RC.f3         <- XVector::subseq(x=reverseComplement(CDS.temp),start=3)
			cds.frames        <- Biostrings::DNAStringSet(c(cds.f1,cds.f2,cds.f3,cds.RC.f1,cds.RC.f2,cds.RC.f3))
			
			aa.frames         <- suppressWarnings(translate(cds.frames,no.init.codon=T,if.fuzzy.codon="solve"))
			aa.frames         <- Biostrings::AAStringSet(gsub(pattern="\\*","X",aa.frames))
			### Removes potential stop codon from search
			aa.frames         <- XVector::subseq(aa.frames,start=1,end=(width(aa.frames)-1))
			if(any(width(aa.frames)==0) | (max(width(aa.frames))<3)){
				CDS.test[identifier.temp,2] <- "short"   ### changes this from "yes" to "short"
				next
			}
			
			frame.search      <- stringr::str_locate(AA.temp,as.character(aa.frames))[,1]
			frame.temp        <- unique(unlist(which(!is.na(frame.search))))
						
			if(length(frame.temp)==1){
				cds.target   <- cds.frames[frame.temp] 
				aa.target    <- aa.frames[frame.temp]
			}
			if(length(frame.temp)==2){
				CDS.test[identifier.temp,2]  <- "multipleFrames"
				next
			}
			if(length(frame.temp)==0){
				CDS.test[identifier.temp,2]  <- "noFramesMatch"
				next
			}
			
		} else {
			CDS.test[identifier.temp,2]      <- "no"
			next
		}
		if(i==1){
			CDS.all  <- cds.target
			AA.all   <- aa.target
			AA.full  <- AA.temp
		} else {
			CDS.all <- c(CDS.all,cds.target)
			AA.all  <- c(AA.all,aa.target)
			AA.full <- c(AA.full,AA.temp)
		}
	} # end of for loop i
	names(AA.all)        <- gsub("_TargetCDS_","_TargetAA_",names(AA.all))
	if(report.noCDS==T){
		noCDS            <- CDS.test[which(CDS.test[,2]=="no"),1]             ### list of target names for which Genbank record does not include a CDS annotation
		shortCDS         <- CDS.test[which(CDS.test[,2]=="short"),1]          ### list of target names for which Genbank record includes a very short CDS annotation (<12bp), such that the CDS target and AA target were not written to file
		multipleFrames   <- CDS.test[which(CDS.test[,2]=="multipleFrames"),1] ### list of target names for which multiple translation frames of target CDS are found in full AA sequence
		noFrames         <- CDS.test[which(CDS.test[,2]=="noFramesMatch"),1]  ### list of target names for which no translation frames of target CDS are found in full AA sequence
		notInGenbankFile <- targetNames.remaining                             ### list of target names for which the associated record is not yet included in GenBank flatfile
		duplicateTarget  <- duplicate.targets                                 ### list of target names for which another target has the same genomic coordinates
		result.names     <- c("CDS StringSet","Target AA StringSet","Full AA StringSet","duplicate target names","targets not in Genbank flatfile","targets without CDS annotation","targets with short CDS annotation","multiple translation frames","no matching translation frames","result names")
		result           <- list(CDS.all,AA.all,AA.full,duplicateTarget,notInGenbankFile,noCDS,shortCDS,multipleFrames,noFrames,result.names)
	}
	if(report.noCDS==F){
		result <- CDS.all
	}
	result
}
