#' @title Summarize a GenBank file of mitochondrial genome as a table
#'
#' Creates a table summizing the annotations held in the GenBank flatefile. This function is used in the function plotGBmtc.
#' 
#' @param pathGB Character string with path to the GenBank flatfile.
#' @return List of length five, which includes [[1]] organism name, [[2]], organelle name [[3]], length of the sequence [[4]] heavy and light strand labels, and [[5]] a table with features and their ranges.
#' @export summarizeGBmtc
summarizeGBmtc <- function(pathGB){
	mtcGB       <- read.gb(pathGB,progress=FALSE)
	features.df <- as.data.frame(do.call(rbind, lapply(X=1:length(mtcGB@features),function(x){ c(product=biofiles::product(mtcGB@features[x]),gene=biofiles::geneID(mtcGB@features[x]),start=as.numeric(start(mtcGB@features[x])),end=as.numeric(end(mtcGB@features[x])),sense=biofiles::strand(mtcGB@features[x]),"strand"=NA,key=mtcGB@features@.Data[[x]]@key)})))
	features.df$gene[is.na(features.df$gene)] <- gsub("tRNA-","",features.df$product[is.na(features.df$gene)])
	features.df2      <- features.df[which(!features.df$key %in% c("gene","source")),c("gene","start","end","sense","strand","key")]
	if(any(is.na(features.df2$gene))){
		features.df2$gene[which(is.na(features.df2$gene))] <- features.df2$key[which(is.na(features.df2$gene))]
	}
	features.df2$gene <- gsub("s-rRNA","12S",features.df2$gene)
	features.df2$gene <- gsub("l-rRNA","16S",features.df2$gene)
	features.df2$gene <- gsub(" ribosomal RNA","",features.df2$gene)

	dupGenes <- names(table(features.df2$gene)[which(table(features.df2$gene) > 1)])
	if(length(dupGenes)>0){
		for(x in 1:length(dupGenes)){
			features.df2$gene[which(features.df2$gene==dupGenes[x])] <- paste0(dupGenes[x],"_",1:length(which(features.df2$gene==dupGenes[x])))
		}
	}
	mode(features.df2$start) <- "numeric"
	mode(features.df2$end)   <- "numeric"
	if(nchar(gsub("[A,G]","",biofiles::getSequence(mtcGB))) < nchar(gsub("[A,G]","",Biostrings::reverseComplement(biofiles::getSequence(mtcGB))))) {
		strand1_weight="heavy"
		heavyStrand=1
		lightStrand=-1
	} else {
		strand1_weight="light"
		heavyStrand=-1
		lightStrand=1
	}
	features.df2$strand[features.df2$sense==heavyStrand] <- "heavy"
	features.df2$strand[features.df2$sense==lightStrand] <- "light"

	mtcSummary <- list(organism=mtcGB@features@.Data[[1]]@qualifiers["organism"],organelle=mtcGB@features@.Data[[1]]@qualifiers["organelle"],ntLength=biofiles::getLength(mtcGB),strands=c(heavy=heavyStrand,light=lightStrand),features=features.df2)
	mtcSummary
}


#' @title Summarize a GenBank file of mitochondrial genome as a table
#'
#' Creates a ggplot of the mitochondion that includes features and their locations on heavy and light strands
#' 
#' @param pathGB Character string with path to the GenBank flatfile.
#' @return List of length three, which includes [[1]] organism name, [[2]] ggplot of mitochondrion [[3]], data frame used to generate the ggplot
#' @export plotGBmtc
plotGBmtc <- function(pathGB, type="file", additionalDF=NULL, zoomout=1.1, radii=c(lightStrand=1, heavyStrand=0.85), linetypes=c(lightStrand=1,heavyStrand=1, CDS=1,stem_loop=1,misc_feature=1,tRNA=1,rRNA=1), widths=c(lightStrand=1,heavyStrand=2, CDS=4, stem_loop=4, misc_feature=4, tRNA=4, rRNA=4), colors=c(lightStrand="black",heavyStrand="black",CDS="green",stem_loop="orange",misc_feature="yellow",tRNA="purple",rRNA="brown"),textadj=c(heavy=0.95,light=1.04)){
	if(type=="accession"){
		accessionGB  <- pathGB
		subject.path <- tempfile()
		URL  <- sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=%s&rettype=gb&retmode=text",paste(accessionGB,collapse=","))
		conn <- utils::download.file(url=URL, destfile=subject.path,quiet=T)
		mtcSummary  <- summarizeGBmtc(subject.path)
	} else {
		mtcSummary  <- summarizeGBmtc(pathGB)
	}
	features.df <- mtcSummary$features
	strands.df  <- data.frame(arcname=factor(c("heavyStrand","lightStrand")),arcstart=c(0,0),arcend=c((2*pi),(2*pi)),strand=c("heavy","light"),type="source",radius=radii[c("heavyStrand","lightStrand")],linetype=linetypes[c("heavyStrand","lightStrand")],linecol=colors[c("heavyStrand","lightStrand")],linewidth=widths[c("heavyStrand","lightStrand")])
	genes.df    <- data.frame(arcname=mtcSummary$features$gene,arcstart=NA,arcend=NA,strand=features.df$strand,type=features.df$key,radius=NA,linetype=NA,linecol=NA,linewidth=NA)
	genes.df$arcstart <- ((features.df$start/mtcSummary$ntLength)*(2*pi))  # [features.df$strand==1]
	genes.df$arcend   <- ((features.df$end/mtcSummary$ntLength)*(2*pi))    # [features.df$strand==1]
	genes.df$radius[genes.df$strand=="heavy"] <- radii["heavyStrand"]
	genes.df$radius[genes.df$strand=="light"] <- radii["lightStrand"]
	

	for(i in 1:length(unique(genes.df$type))){
		tempType=unique(genes.df$type)[i]
		### track linetypes
		if(any(names(linetypes)==tempType)){
			genes.df$linetype[genes.df$type==tempType] <- linetypes[names(linetypes)==tempType]
		} else {
			genes.df$linetype[genes.df$type==tempType] <- 1
		}
		### track colors
		if(any(names(colors)==tempType)){
			genes.df$linecol[genes.df$type==tempType] <- colors[names(colors)==tempType]
		} else {
			genes.df$linecol[genes.df$type==tempType] <- "gray"
		}
		### track widths
		if(any(names(widths)==tempType)){
			genes.df$linewidth[genes.df$type==tempType] <- widths[names(widths)==tempType]
		} else {
			genes.df$linewidth[genes.df$type==tempType] <- widths["CDS"]
		}
	}
	if(is.null(additionalDF)){
		allArcs.df <- rbind(strands.df,genes.df)
	} else {
		allArcs.df <- rbind(strands.df,genes.df,additionalDF)
	}
	allArcs.df2 <- allArcs.df
	allArcs.df2$text <- gsub("^.+Strand$","",allArcs.df$arcname)
	Hadj <- textadj["heavy"]
	Ladj <- textadj["light"]
	allArcs.df2$Xtext <- NA
	allArcs.df2$Ytext <- NA
	allArcs.df2$Xtext[allArcs.df2$"strand"=="heavy"] <- (allArcs.df2$radius[allArcs.df2$"strand"=="heavy"]*Hadj)*cos((((allArcs.df2$arcstart[allArcs.df2$"strand"=="heavy"]+allArcs.df2$arcend[allArcs.df2$"strand"=="heavy"])/2)*(-1)+(pi/2)))
	allArcs.df2$Ytext[allArcs.df2$"strand"=="heavy"] <- (allArcs.df2$radius[allArcs.df2$"strand"=="heavy"]*Hadj)*sin((((allArcs.df2$arcstart[allArcs.df2$"strand"=="heavy"]+allArcs.df2$arcend[allArcs.df2$"strand"=="heavy"])/2)*(-1)+(pi/2)))
	allArcs.df2$Xtext[allArcs.df2$"strand"=="light"] <- (allArcs.df2$radius[allArcs.df2$"strand"=="light"]*Ladj)*cos((((allArcs.df2$arcstart[allArcs.df2$"strand"=="light"]+allArcs.df2$arcend[allArcs.df2$"strand"=="light"])/2)*(-1)+(pi/2)))
	allArcs.df2$Ytext[allArcs.df2$"strand"=="light"] <- (allArcs.df2$radius[allArcs.df2$"strand"=="light"]*Ladj)*sin((((allArcs.df2$arcstart[allArcs.df2$"strand"=="light"]+allArcs.df2$arcend[allArcs.df2$"strand"=="light"])/2)*(-1)+(pi/2)))
	
	allArcs.df2$quad1or4  <- allArcs.df2$Xtext > 0
	allArcs.df2$textangle <- ((((allArcs.df2$arcstart+allArcs.df2$arcend)/2)*(-1)+(pi/2)))*(180/pi)
	allArcs.df2$textangle[!allArcs.df2$quad1or4] <- allArcs.df2$textangle[!allArcs.df2$quad1or4]+180
	
	allArcs.df2$texthjust <- "right"
	allArcs.df2$texthjust[!allArcs.df2$quad1or4 & allArcs.df2$radius==allArcs.df2$radius[allArcs.df2$arcname=="heavyStrand"]] <- "left"
	allArcs.df2$texthjust[allArcs.df2$quad1or4 & allArcs.df2$radius==allArcs.df2$radius[allArcs.df2$arcname=="lightStrand"]]  <- "left"

	outerBoundary <- data.frame(arcname="outerBoundary", arcstart=0, arcend=2*pi, strand=NA, type=NA, radius=max(allArcs.df2$radius)*zoomout, linetype=1, linecol="white", linewidth=0, text="", Xtext=0, Ytext=0, quad1or4=FALSE, textangle=0, texthjust=0)
	organism.df   <- data.frame(arcname="organism", arcstart=0, arcend=2*pi, strand=NA, type=NA, radius=max(allArcs.df2$radius)*(zoomout/2), linetype=1, linecol="white", linewidth=0, text=paste0(mtcSummary$organism,"\nmitochondrion"), Xtext=0, Ytext=0, quad1or4=FALSE, textangle=0, texthjust=0)
	allArcs.df3 <- rbind(allArcs.df2, outerBoundary,organism.df)
	allArcs.df4 <- allArcs.df3[order(allArcs.df3$radius),]
	rownames(allArcs.df4) <- 1:nrow(allArcs.df4)
	allArcs.df4$plotOrder <- 1:nrow(allArcs.df4)
	allArcs.df4$arcname   <- factor(allArcs.df4$arcname, levels = allArcs.df4$arcname[order(allArcs.df4$plotOrder)])
	res1 <- ggplot2::ggplot(allArcs.df4) + ggforce::geom_arc(ggplot2::aes(x0 = 0, y0 = 0, r=radius, start = arcstart, end=arcend, size=arcname, color=arcname, linetype=arcname)) + ggplot2::scale_linetype_manual(values=allArcs.df4$linetype,guide="none") + ggplot2::scale_color_manual(values=allArcs.df4$linecol,guide="none") + ggplot2::scale_size_manual(values=allArcs.df4$linewidth, guide="none") + ggplot2::geom_text(ggplot2::aes(x=Xtext,y=Ytext,label=text, angle=textangle, hjust=texthjust), size=2) + ggplot2::theme_void()
	res1
	res2 <- allArcs.df4
	res  <- list(organism=mtcSummary$organism,plot=res1,data=res2)
	res
}

# library(REEs)
# library(biofiles)
# library(ggplot2)
# Trimeresurus.mtc.summary  <- REEs::summarizeGBmtc(pathGB="/users/jeff/Documents/Trimeresurus-albilabris_mitochondrial-genome_NC_022820.gb")
# Trimeresurus.mtc.plot     <- REEs::plotGBmtc(pathGB="/users/jeff/Documents/Trimeresurus-albilabris_mitochondrial-genome_NC_022820.gb")
# Trimeresurus.mtc.gb       <- Biofiles::read.gb("/users/jeff/Documents/Trimeresurus-albilabris_mitochondrial-genome_NC_022820.gb")
# Trimeresurus.mtc.sequence <- Biofiles::getSequence(Trimeresurus.mtc.gb)
# Biostrings::writeXStringSet(Trimeresurus.mtc.sequence,filepath="/users/jeff/Documents/reference.fa")
# genes <- Trimeresurus.mtc.summary$features
# for(i in 1:nrow(genes)){
# 	geneseq <- DNAStringSet(DNAString(as.character(Trimeresurus.mtc.sequence))[genes$start[i]:genes$end[i]])
# 	names(geneseq) <- genes$gene[i]
# 	if(i==1){
# 		allgenes <- geneseq
# 	} else {
# 		allgenes <- c(allgenes,geneseq)
# 	}
# }
# 
# writeXStringSet(allgenes,"/users/jeff/Documents/mtgenes.fa")

# '/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/Processed_Samples/Achalinus-spinalis_KU312258/Achalinus-spinalis_KU312258_singletons.fastq.gz'
# '/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/Processed_Samples/Achalinus-spinalis_KU312258/Achalinus-spinalis_KU312258_READ2_unmerged.fastq.gz'
# '/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/Processed_Samples/Achalinus-spinalis_KU312258/Achalinus-spinalis_KU312258_READ1_unmerged.fastq.gz'

# sampleName   <- "Cerberus-schneiderii_KU324534"
# read1        <- "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/Processed_Samples/Cerberus-schneiderii_KU324534/Cerberus-schneiderii_KU324534_READ1_unmerged.fastq.gz"
# read2        <- "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/Processed_Samples/Cerberus-schneiderii_KU324534/Cerberus-schneiderii_KU324534_READ2_unmerged.fastq.gz"
# reads.merged <- "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/Processed_Samples/Cerberus-schneiderii_KU324534/Cerberus-schneiderii_KU324534_singletons.fastq.gz"
# 
# referencePATH <- "/panfs/pfs.local/home/j926w878/work/mitogenomes/reference.fa"
# geneFilePATH  <- "/panfs/pfs.local/home/j926w878/work/mitogenomes/mtgenes.fa"
# out.dir       <- "/panfs/pfs.local/home/j926w878/work/mitogenomes/"
# 
# bbmapPATH  <- "/panfs/pfs.local/home/j926w878/programs/bbmap/bbmap.sh"
# cap3PATH   <- "/panfs/pfs.local/home/j926w878/programs/CAP3/cap3"
# spadesPATH <- "/panfs/pfs.local/home/j926w878/programs/SPAdes-3.12.0-Linux/bin/spades.py"

#' @title Get mitochondrial genome from reads
#'
#' 
#' 
#' @param sampleName Name of the sample
#' @param read1 Path to fasta file containing cleaned, unmerged (unmeargeable) read1 reads. This file has '_READ1_unmerged' in the filename if it was generated with the SnakeCap pipeline.
#' @param read2 Path to fasta file containing cleaned, unmerged (unmeargeable) read2 reads. This file has '_READ2_unmerged' in the filename if it was generated with the SnakeCap pipeline.
#' @param reads.merged Path to fasta file containing the cleaned, merged reads. This file has 'singletons' in the filename if it was generated with the SnakeCap pipeline.
#' @param referencePATH Path to fasta file containing mitogenomes of many snakes.
#' @param geneFilePATH Path to fasta file containing genes for one of the mitogenomes.
#' @param out.dir Path to directory where output should be saved.
#' @param bbmapPATH Path to bbmap bash script. Default "bbmap.sh", in which case bbmap. must be on your system path.
#' @param cap3PATH Path to cap3 executable. Default "cap3", in which case cap3 must be on your system path.
#' @param spadesPATH Path to spades python executable. Default "spades.py", in which case spades must be on your system path.
#' @param pblatPATH Path to pblat executable. Default "pblat", in which case pblat must be on your system path.
#' @return DNAStringSet; mitogenome is written to out.dir
#' @export get.mitogenome

get.mitogenome <- function(sampleName,read1,read2=NULL,reads.merged=NULL,referencePATH,geneFilePATH,out.dir,bbmapPATH="bbmap.sh'",cap3PATH="cap3",spadesPATH="spades.py",pblatPATH="pblat"){
	reference <- basename(referencePATH)
	#gene.file <- basename(geneFilePATH)
	#raw.dir   <- dirname(read1)
	dir.check.create(out.dir)
	#dir.check.create(file.path(out.dir,"Species_mtGenomes"))

	# Modify permissions to make sure that bbmap, cap3, and spades are executable
	#system(paste("chmod +x",bbmapPATH))
	#system(paste("chmod +x",cap3))
	#system(paste("chmod +x",spades.py))

	#Options
	options(stringsAsFactors = FALSE)
	options(warn=2)

	###################################################################
	### Step 1: Gather read data and assemble mitochondrial genomes ###
	###################################################################

	setwd(out.dir)
	#Pick out matching reads to mt Genomes. Output files will be in out.dir
	if(!is.null(read2)){
		system(sprintf("'%s'  -Xmx8g ref='%s' in1='%s' in2='%s' vslow k=12 minid=0.7 outm1=read1.fq outm2=read2.fq",bbmapPATH,referencePATH,read1,read2),ignore.stderr = T)
	} else {
		system(sprintf("'%s'  -Xmx8g ref='%s' in='%s' vslow k=12 minid=0.7 outm=read1.fq",bbmapPATH,referencePATH,read1),ignore.stderr = T)
	}
	if(!is.null(reads.merged)){
		system(sprintf("'%s' -Xmx8g ref='%s' in='%s' vslow k=12 minid=0.7 outm=singleton.fq",bbmapPATH,referencePATH,reads.merged), ignore.stderr = T)
	}
	k     <- c(9,13,21,33,55,77,99,127)
	k.val <- paste(k, collapse = ",")

	# creates an empty file called "current_seed.fasta" in out.dir
	system("touch current_seed.fasta")

	# Initial values for the while loop
	new.len        <- 0
	counter        <- 0
	repeat.counter <- 0
	seeding        <- T
	min.id         <- "0.7"

	# Running the while loop
	while (seeding == T){
		#Copy new reference to do recursively
		counter  <- counter+1
		prev.len <- new.len
		# Skips the first iteration since its already done.
		if (counter >= 2){
			#Pick out matching reads to mt Genomes
			if(!is.null(read2)){
				system(sprintf("'%s' -Xmx8g ref='current_seed.fasta' in1='%s' in2='%s' vslow k=12 minid=%s outm1=t_read1.fq outm2=t_read2.fq",bbmapPATH,referencePATH,read1,read2,min.id),ignore.stderr = T)
			} else {
				system(sprintf("'%s' -Xmx8g ref='current_seed.fasta' in='%s' vslow k=12 minid=%s outm=t_read1.fq",bbmapPATH,referencePATH,read1,min.id),ignore.stderr = T)
			}
			if(!is.null(reads.merged)){
				system(sprintf("'%s' -Xmx8g ref='current_seed.fasta' in='%s' vslow k=12 minid=%s outm=t_singleton.fq",bbmapPATH,referencePATH,reads.merged,min.id), ignore.stderr = T)
			}
			# concatenates 't_read1.fq' and 'o_read1.fq' and appends to 'read1.fq'
			system(sprintf("cat '%s/t_read1.fq' '%s/o_read1.fq' >> '%s/read1.fq'",out.dir,out.dir,out.dir))
			# removes current 't_read1.fq'
			invisible(file.remove(file.path(out.dir,"t_read1.fq")))
			# concatenates 't_read2.fq' and 'o_read2.fq' and appends to 'read2.fq'
			if(!is.null(read2)){
				system("cat t_read2.fq o_read2.fq >> read2.fq")
				# removes current 't_read2.fq'
				invisible(file.remove(file.path(out.dir,"t_read2.fq")))
			}
			if(!is.null(reads.merged)){
				# concatenates 't_singleton.fq' and 'o_singleton.fq' and appends to 'singleton.fq'
				system("cat t_singleton.fq o_singleton.fq >> singleton.fq")
				# removes current 't_singleton.fq'
				invisible(file.remove(file.path(out.dir,"t_singleton.fq")))
			}
		}
		#####################
		## Run SPADES on sample.
		if((!is.null(read2)) & (!is.null(reads.merged))){
			system(sprintf("'%s' --pe1-1 '%s/read1.fq' --pe1-2 '%s/read2.fq' --pe1-s '%s/singleton.fq' -o spades -k %s --careful -t 8 -m 8",spadesPATH,out.dir,out.dir,out.dir,k.val), ignore.stdout = T)
		}
		if((!is.null(read2)) & is.null(reads.merged)){
			system(sprintf("'%s' --pe1-1 '%s/read1.fq' --pe1-2 '%s/read2.fq' -o spades -k %s --careful -t 8 -m 8",spadesPATH,out.dir,out.dir,out.dir,k.val), ignore.stdout = T)
		}
		if(is.null(read2) & is.null(reads.merged)){
			system(sprintf("'%s' --s '%s/read1.fq' -o spades -k %s --careful -t 8 -m 8",spadesPATH,out.dir,out.dir,out.dir,k.val), ignore.stdout = T)
		}
		#Checks to see if one kmer failed or not
		while (file.exists(file.path(out.dir,"spades/contigs.fasta")) == F){
			#subtract Ks until it works
			#system("rm -r spades")
			invisible(unlink(file.path(out.dir,"spades"),recursive=T))
			k <- k[-length(k)]
			if(length(k) == 0) { break }
			k.val  <- paste(k, collapse = ",")
			min.id <- "0.6"
			if((!is.null(read2)) & (!is.null(reads.merged))){
				system(sprintf("'%s' --pe1-1 '%s/read1.fq' --pe1-2 '%s/read2.fq' --pe1-s '%s/singleton.fq' -o spades -k %s --careful -t 8 -m 8",spadesPATH,out.dir,out.dir,out.dir,k.val), ignore.stdout = T)
			}
			if((!is.null(read2)) & is.null(reads.merged)){
				system(sprintf("'%s' --pe1-1 '%s/read1.fq' --pe1-2 '%s/read2.fq' -o spades -k %s --careful -t 8 -m 8",spadesPATH,out.dir,out.dir,out.dir,k.val), ignore.stdout = T)
			}
			if(is.null(read2) & is.null(reads.merged)){
				system(sprintf("'%s' --s '%s/read1.fq' -o spades -k %s --careful -t 8 -m 8",spadesPATH,out.dir,out.dir,out.dir,k.val), ignore.stdout = T)
			}
		}#end while
		# If the k-mers are all run out, therefore nothing can be assembled
		if (length(k) == 0) { 
			paste("k-mer values all used up, cannot assemble!")
			# system("rm read1.fq t_read1.fq o_read1.fq")
			invisible(file.remove(file.path(out.dir,"read1.fq"),file.path(out.dir,"t_read1.fq"),file.path(out.dir,"o_read1.fq")))
			
			if(!is.null(read2)){
				# system("rm read2.fq t_read2.fq o_read2.fq")
				invisible(file.remove(file.path(out.dir,"read2.fq"),file.path(out.dir,"t_read2.fq"),file.path(out.dir,"o_read2.fq")))
			}
			if(!is.null(reads.merged)){
				#system("rm singleton.fq t_singleton.fq o_singleton.fq")
				invisible(file.remove(file.path(out.dir,"singleton.fq"),file.path(out.dir,"t_singleton.fq"),file.path(out.dir,"o_singleton.fq")))
			}
			#system("rm -r spades")
			invisible(unlink(file.path(out.dir,"spades"),recursive=T))
			seeding = F
		}# end if

		# renames read1.fq and read2.fq to read1.fq and read2.fq
		if (counter == 1){
			system(sprintf("mv '%s/read1.fq' '%s/o_read1.fq'",out.dir,out.dir))
			if(!is.null(read2)){
				#system("mv read2.fq o_read2.fq")
				system(sprintf("mv '%s/read2.fq' '%s/o_read2.fq'",out.dir,out.dir))
			}
			if(!is.null(reads.merged)){
				#system("mv singleton.fq o_singleton.fq")
				system(sprintf("mv '%s/singleton.fq' '%s/o_singleton.fq'",out.dir,out.dir))
			}
		}
		# copies contigs.fasta to current_seed.fasta
		system("cp spades/contigs.fasta current_seed.fasta")
		if (counter >= 2) {
			system("rm read1.fq")
			if(!is.null(read2)){
				system("rm read2.fq")
			}
			if(!is.null(reads.merged)){
				system("rm singleton.fq")
			}
		}
		# system("rm -r spades")
		invisible(unlink(file.path(out.dir,"spades"),recursive=T))
		reference <- "current_seed.fasta"

		#Check size
		temp.count <- scan(file = "current_seed.fasta", what = "character")
		new.len    <- sum(nchar(temp.count[-grep(">", temp.count)]))
		no.contigs <- length(temp.count[grep(">", temp.count)])
		print(paste0("iteration ", counter, " complete!"))
		print(paste0("new length: ", new.len, ". Old length: ", prev.len))
		if (new.len == prev.len || counter == 20){
			seeding<-F 
			system("rm o_read1.fq")
			if(!is.null(read2)){
				system("rm o_read2.fq")
			}
			if(!is.null(reads.merged)){
				system("rm o_singleton.fq")
			}
			print(paste("mitogenome complete after ", counter, " iterations!", sep = ""))
			min.id <- "0.7"
		}
		#If the file gets too large, its due to repeats
		if (new.len >= 23000){
			#runs cap3 to merge similar contigs (pull only clustered contigs out?)
			system(paste(cap3PATH," current_seed.fasta -z 1 -o 16 -e 11 -s 251", " > ","log.fasta.cap.txt", sep = "")) 
			#Reads in results files
			#temp.assembled <- scanFa(FaFile(paste("current_seed.fasta.cap.contigs", sep = "")))
			temp.assembled <- Biostrings::readDNAStringSet("current_seed.fasta.cap.contigs")
			#temp.singlets  <- scanFa(FaFile(paste("current_seed.fasta.cap.singlets", sep = "")))
			temp.singlets <- Biostrings::readDNAStringSet("current_seed.fasta.cap.singlets")
			keep.singlets <- temp.singlets[width(temp.singlets) >= 100]
			final.save    <- c(temp.assembled, keep.singlets)
			#Writes contigs for cap3
			#write.loci <- as.list(as.character(final.save))
			#write.fasta(sequences = write.loci, names = names(write.loci),"current_seed.fasta", nbchar = 1000000, as.string = T)
			Biostrings::writeXStringSet(final.save,filepath="current_seed.fasta")
			#Get cap3 files and deletes
			cap.files  <- list.files(pattern = "", full.names = F, recursive = F)
			cap.remove <- cap.files[grep(pattern = paste("fasta.cap*.", sep =""), x = cap.files)]
			system(paste("rm ", paste(cap.remove, collapse = " ") ))
			min.id <- "0.95"
			#makes sure this doesnt go on forever and ever
			repeat.counter<-repeat.counter+1
			if (repeat.counter >= 5){ 
				print(paste("repeat counter hit 5"))
				system("rm o_read1.fq")
				if(!is.null(read2)){
					system("rm o_read2.fq")
				}
				if(!is.null(reads.merged)){
					system("rm o_singleton.fq")
				}
				seeding <- F 
			}
		}#end length > 30,000 if
	}#end while
	
	### Save finished mitogenome
	# loads up fasta file
	#contigs <- scanFa(FaFile("current_seed.fasta"))
	contigs <- Biostrings::readDNAStringSet("current_seed.fasta")
	if(length(contigs) > 0){
		#Trys to merge contigs if there are more than 1
		if(length(contigs) >= 2){
			#runs cap3 to merge similar contigs (pull only clustered contigs out?)
			system(sprintf("'%s' current_seed.fasta -z 1 -o 16 -e 11 -s 251 > log.fasta.cap.txt",cap3PATH))
			
			#Reads in results files
			temp.assembled <- Biostrings::readDNAStringSet("current_seed.fasta.cap.contigs")
			temp.singlets  <- Biostrings::readDNAStringSet("current_seed.fasta.cap.singlets")

			keep.singlets  <- temp.singlets[width(temp.singlets) >= 100]
			contigs        <- c(temp.assembled, keep.singlets)
			
			#Get cap3 files and deletes
			cap.files  <- list.files(pattern = "", full.names = F, recursive = F)
			cap.remove <- cap.files[grep(pattern = paste("fasta.cap*.", sep =""), x = cap.files)]
			system(paste("rm ", paste(cap.remove, collapse = " ") ))
		}#end if
		
		if(sum(width(contigs)) <= 1000) { 
			print("less than 1000bp, not enough data to extract")
			next 
		}
		#Writes the mitochondrial genome file
		names(contigs)<- paste("sequence_", seq(1:length(contigs)), sep = "")
		#write.loci    <- as.list(as.character(contigs))
		#write.fasta(sequences = write.loci, names = names(write.loci),paste0("Species_mtGenomes/", sampleName, ".fa"), nbchar = 1000000, as.string = T)
		Biostrings::writeXStringSet(contigs,filepath=sprintf("%s/%s_mitocontigs.fa",out.dir,sampleName))
		system("rm current_seed.fasta")
	}
	system("rm -r ref")
	#### PBLAT search for genes in mitocontigs
	contigs <- Biostrings::readDNAStringSet(sprintf("%s/%s_mitocontigs.fa",out.dir,sampleName))
	system(sprintf("mpirun '%s' -threads=8 '%s' '%s' -tileSize=8 -minIdentity=60 -noHead -out=pslx '%s/mt_to_genes.pslx'" ,pblatPATH, sprintf("%s/%s_mitocontigs.fa",out.dir,sampleName), geneFilePATH, out.dir))

	setwd(out.dir)
	temp.count <- scan(file = "mt_to_genes.pslx", what = "character")
	if (length(temp.count) == 0){
		print("No matching mitochondrial genes were found.")
	} else {
		match.data <- data.table::fread(file.path(out.dir,"mt_to_genes.pslx"), sep = "\t", header = F, stringsAsFactors = FALSE)
		headers<-c("matches", "misMatches", "repMatches", "nCount", "qNumInsert", "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName", "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd", "blockCount", "blockSize", "qStarts", "tStarts", "qSeq", "tSeq")
		setnames(match.data, headers)
		loci.names <- unique(match.data$qName)
		sep.loci   <- DNAStringSet()
		for (j in 1:length(loci.names)){
			#pulls out data that matches to multiple contigs
			sub.data <- match.data[match.data$qName %in% loci.names[j],]
			sub.data <- sub.data[sub.data$matches == max(sub.data$matches),][1]
			if (sub.data$strand == "+"){
				#Cuts the node apart and saves separately
				sub.data$tStart<-sub.data$tStart-sub.data$qStart+1
				#Fixes ends
				sub.data$tEnd<-sub.data$tEnd+(sub.data$qSize-sub.data$qEnd)
			} else {
				sub.data$tStart<-sub.data$tStart-(sub.data$qSize-sub.data$qEnd)
				#Fixes ends
				sub.data$tEnd<-sub.data$tEnd+sub.data$qStart+1
			}
		
			#If it ends up with a negative start
			if (sub.data$tStart <= 0){
				sub.data$tStart <- 1
			}
			#Fixes if the contig is smaller than the full target locus
			if (sub.data$tEnd >= sub.data$tSize){
				sub.data$tEnd <- sub.data$tSize
			}
			#Gets start and end
			start.pos <- min(sub.data$tStart, sub.data$tEnd)
			end.pos   <- max(sub.data$tStart, sub.data$tEnd)
			  
			temp.contig    <- contigs[names(contigs) == sub.data$tName]
			new.seq        <- subseq(x = temp.contig, start = start.pos, end = end.pos)
			names(new.seq) <- sub.data$qName
			sep.loci       <- c(sep.loci, new.seq)
		} #end j loop
		#Writes the full mitochondrial genome file
		Biostrings::writeXStringSet(sep.loci, filepath = sprintf("%s/%s_mitogenome.fa",out.dir,sampleName))
		system(sprintf("rm '%s/mt_to_genes.pslx'",out.dir))
	}
}


