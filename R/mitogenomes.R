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

sampleName   <- "Cerberus-schneiderii_KU324534"
read1        <- "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/Processed_Samples/Cerberus-schneiderii_KU324534/Cerberus-schneiderii_KU324534_READ1_unmerged.fastq.gz"
read2        <- "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/Processed_Samples/Cerberus-schneiderii_KU324534/Cerberus-schneiderii_KU324534_READ2_unmerged.fastq.gz"
reads.merged <- "/panfs/pfs.local/home/j926w878/scratch/scratch_v3/SequenceCapture/SnakeCap_AllSamples/Processed_Samples/Cerberus-schneiderii_KU324534/Cerberus-schneiderii_KU324534_singletons.fastq.gz"

referencePATH <- "/panfs/pfs.local/home/j926w878/work/mitogenomes/reference.fa"
geneFilePATH  <- "/panfs/pfs.local/home/j926w878/work/mitogenomes/mtgenes.fa"
out.dir       <- "/panfs/pfs.local/home/j926w878/work/mitogenomes/"

bbmapPATH  <- "/panfs/pfs.local/home/j926w878/programs/bbmap/bbmap.sh"
cap3PATH   <- "/panfs/pfs.local/home/j926w878/programs/CAP3/cap3"
spadesPATH <- "/panfs/pfs.local/home/j926w878/programs/SPAdes-3.12.0-Linux/bin/spades.py"

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
#' @return DNAStringSet; mitogenome is written to out.dir
#' @export get.mitogenome
get.mitogenome <- function(sampleName,read1,read2,reads.merged,referencePATH,geneFilePATH,out.dir,bbmapPATH="bbmap.sh'",cap3PATH="cap3",spadesPATH="spades.py"){
	reference <- basename(referencePATH)
	gene.file <- basename(geneFilePATH)
	#raw.dir   <- dirname(read1)
	dir.check.create(out.dir)
	dir.check.create(file.path(out.dir,"Species_mtGenomes"))

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
	system(sprintf("'%s'  -Xmx8g ref='%s' in1='%s' in2='%s' vslow k=12 minid=0.7 outm1=read1.fq outm2=read2.fq",bbmapPATH,referencePATH,read1,read2),ignore.stderr = T)
	system(sprintf("'%s' -Xmx8g ref='%s' in='%s' vslow k=12 minid=0.7 outm=singleton.fq",bbmapPATH,referencePATH,reads.merged), ignore.stderr = T)

	k     <- c(9,13,21,33,55,77,99,127)
	k.val <- paste(k, collapse = ",")

	###################################################
	#### Not really sure what all of this is about: ###
	###################################################
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
		# Skips the first one [of???] since its already done.
		if (counter >= 2){
			#Pick out matching reads to mt Genomes
			system(sprintf("'%s' -Xmx8g ref='current_seed.fasta' in1='%s' in2='%s' vslow k=12 minid=%s outm1=t_read1.fq outm2=t_read2.fq",bbmapPATH,referencePATH,read1,read2,min.id),ignore.stderr = T)
			system(sprintf("'%s' -Xmx8g ref='current_seed.fasta' in='%s' vslow k=12 minid=%s outm=t_singleton.fq",bbmapPATH,referencePATH,reads.merged,min.id), ignore.stderr = T)
			# concatenates 't_read1.fq' and 'o_read1.fq' and appends to 'read1.fq'
			system("cat t_read1.fq o_read1.fq >> read1.fq")
			# concatenates 't_read2.fq' and 'o_read2.fq' and appends to 'read2.fq'
			system("cat t_read2.fq o_read2.fq >> read2.fq")
			# concatenates 't_singleton.fq' and 'o_singleton.fq' and appends to 'singleton.fq'
			system("cat t_singleton.fq o_singleton.fq >> singleton.fq")
			# removes current 't_read1.fq' 't_read2.fq' and 't_singleton.fq'
			system("rm t_read1.fq t_read2.fq t_singleton.fq")
		}
		# Run SPADES on sample.
		system(sprintf("'%s' --pe1-1 '%s/read1.fq' --pe1-2 '%s/read2.fq' --pe1-s '%s/singleton.fq' -o spades -k %s --careful -t 8 -m 8",spadesPATH,out.dir,out.dir,out.dir,k.val), ignore.stdout = T)
		#Checks to see if one kmer failed or not
		while (file.exists("spades/contigs.fasta") == F){
			#subtract Ks until it works
			system("rm -r spades")
			k <- k[-length(k)]
			if(length(k) == 0) { break }
			k.val  <- paste(k, collapse = ",")
			min.id <- "0.6"
			system(sprintf("'%s' --pe1-1 '%s/read1.fq' --pe1-2 '%s/read2.fq' --pe1-s '%s/singleton.fq' -o spades -k %s --careful -t 8 -m 8",spadesPATH,out.dir,out.dir,out.dir,k.val), ignore.stdout = T)
		}#end while
		# If the k-mers are all run out, therefore nothing can be assembled
		if (length(k) == 0) { 
			paste("k-mer values all used up, cannot assemble!")
			system("rm read1.fq read2.fq singleton.fq t_read1.fq t_read2.fq t_singleton.fq o_read1.fq o_read2.fq o_singleton.fq")
			system("rm -r spades")
			seeding = F 
		}# end if

		# renames read1.fq and read2.fq to read1.fq and read2.fq
		if (counter == 1){
			system("mv read1.fq o_read1.fq")
			system("mv read2.fq o_read2.fq")
			system("mv singleton.fq o_singleton.fq")
		}

		# copies contigs.fasta to current_seed.fasta
		system("cp spades/contigs.fasta current_seed.fasta")
		if (counter >= 2) {
			system("rm read1.fq read2.fq singleton.fq")
		}
		system("rm -r spades")
		reference <- "current_seed.fasta"

		#Check size
		temp.count <- scan(file = "current_seed.fasta", what = "character")
		new.len    <- sum(nchar(temp.count[-grep(">", temp.count)]))
		no.contigs <- length(temp.count[grep(">", temp.count)])
		print(paste0("iteration ", counter, " complete!"))
		print(paste0("new length: ", new.len, ". Old length: ", prev.len))
		if (new.len == prev.len || counter == 20){
			seeding<-F 
			system("rm o_read1.fq o_read2.fq o_singleton.fq")
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
				system("rm o_read1.fq o_read2.fq o_singleton.fq")
				seeding <- F 
			}
		}#end length > 30,000 if
	}#end while
	
	### Save finished genome
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
		#Writes the full mitochondrial genome file
		system("rm current_seed.fasta")
		names(contigs)<- paste("sequence_", seq(1:length(contigs)), sep = "")
		#write.loci    <- as.list(as.character(contigs))
		#write.fasta(sequences = write.loci, names = names(write.loci),paste0("Species_mtGenomes/", sampleName, ".fa"), nbchar = 1000000, as.string = T)
		Biostrings::writeXStringSet(contigs,filepath=paste0("Species_mtGenomes/", sampleName, ".fa"))
	}
	system("rm -r ref")
}



## ###### Make steps 2–4 separate functions
## ###########################################################################
## ### Step 2: Assess completeness of the mitochondrial genome and annotate ##
## ###########################################################################
## samples <- sampleName
## #PSLX headers
## headers<-c("matches", "misMatches", "repMatches", "nCount", "qNumInsert", "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName", 
##            "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd", "blockCount", "blockSize", "qStarts", "tStarts", "qSeq", "tSeq")
## 
## #Creates new directory and enters this working directory
## setwd(out.dir)
## dir.create("Species_Loci")
## spp.samples<-list.files("Species_mtGenomes/.")
## spp.samples<-gsub(".fa$", "", spp.samples)
## 
## for (i in 1:length(spp.samples)){
## 
## 	#Load in the data
## 	contigs<-readDNAStringSet(paste("Species_mtGenomes/", spp.samples[i], ".fa", sep = ""))   # loads up fasta file
## 
## 	#Matches samples to loci
## 	system(paste("mpirun pblat -threads=", threads, " Species_mtGenomes/", spp.samples[i], ".fa ",gene.file, " -tileSize=8 -minIdentity=60"," -noHead -out=pslx mt_to_genes.pslx", sep = ""), ignore.stdout = T)
## 	
## 	#Need to load in transcriptome for each species and take the matching transcripts to the database
## 	temp.count<-scan(file = "mt_to_genes.pslx", what = "character")
## 	if (length(temp.count) == 0){
## 		print("No matching mitochondrial genes were found.")
## 		next
## 	}
## 	match.data<-fread("mt_to_genes.pslx", sep = "\t", header = F, stringsAsFactors = FALSE)
## 	setnames(match.data, headers)
## 
## 	loci.names<-unique(match.data$qName)
## 	sep.loci<-DNAStringSet()
## 	for (j in 1:length(loci.names)){
## 		#pulls out data that matches to multiple contigs
## 		sub.data <- match.data[match.data$qName %in% loci.names[j],]
## 		sub.data <- sub.data[sub.data$matches == max(sub.data$matches),][1]
## 		if (sub.data$strand == "+"){
## 			#Cuts the node apart and saves separately
## 			sub.data$tStart<-sub.data$tStart-sub.data$qStart+1
## 			#Fixes ends
## 			sub.data$tEnd<-sub.data$tEnd+(sub.data$qSize-sub.data$qEnd)
## 		} else {
## 			sub.data$tStart<-sub.data$tStart-(sub.data$qSize-sub.data$qEnd)
## 			#Fixes ends
## 			sub.data$tEnd<-sub.data$tEnd+sub.data$qStart+1
## 		}
## 
## 		#If it ends up with a negative start
## 		if (sub.data$tStart <= 0){
## 			sub.data$tStart <- 1
## 		}
## 		#Fixes if the contig is smaller than the full target locus
## 		if (sub.data$tEnd >= sub.data$tSize){
## 			sub.data$tEnd <- sub.data$tSize
## 		}
## 		#Gets start and end
## 		start.pos <- min(sub.data$tStart, sub.data$tEnd)
## 		end.pos   <- max(sub.data$tStart, sub.data$tEnd)
## 		  
## 		temp.contig    <- contigs[names(contigs) == sub.data$tName]
## 		new.seq        <- subseq(x = temp.contig, start = start.pos, end = end.pos)
## 		names(new.seq) <- sub.data$qName
## 		sep.loci       <- append(sep.loci, new.seq)
## 	}#end j loop
## 
## 	#Writes the full mitochondrial genome file
## 	write.loci<-as.list(as.character(sep.loci))
## 	write.fasta(sequences = write.loci, names = names(write.loci),paste("Species_Loci/", spp.samples[i], "_mito_genes.fa", sep = ""), nbchar = 1000000, as.string = T)
## 
## 	system("rm mt_to_genes.pslx")
## 
## }#end i loop
## 
## ############################################
## ### Step 3: Create mitogenome alignments ###
## ############################################
## 
## # Alignment settings
## secondary.structure  <- TRUE                          ### If true, runs mafft-qinsi on mt regions that have secondary structure. Takes structure into acct.
## min.taxa             <- 3                             ### min number of individuals to keep an alignment
## min.prop             <- "0.25"                        ### min coverage per individual. e.g., if set to "0.25", for a 100bp gene, needs 25 bp to keep.
## min.len              <- "100"                         ### min length for trimming. Set to this value as you dont usually want to trim t-RNAs
## trim.cds             <- FALSE                         ### defaults to no trimming for coding sequence. Usually destroys mtGenes
## gblocks              <- FALSE                         ### If you want to use 
## trimal               <- TRUE                          ### If you want to use
## 
## setwd(out.dir)
## 
## #Sets up the loci to align
## ref.data      <- scanFa(FaFile(gene.file))
## species.names <- list.files("Species_Loci/.", full.names = F)
## species.names <- species.names[species.names != ""]
## dir.create("mtGenes_Fastas")
## dir.create("mtGenes_Aligned")
## 
## #Aligns each potential locus
## 
## for(i in 1:length(ref.data)){
## 	#######################################################
## 	### STEP 3.1: Gets the locus data from each species ###
## 	#######################################################
## 
## 	#Gets all species data
## 	final.gene<-DNAStringSet()
## 	for (j in 1:length(species.names)){
## 		#Looks for this gene in the species data
## 		spp.data<-scanFa(FaFile(paste("Species_Loci/", species.names[j], sep = "")))   # loads up fasta file
## 		spp.gene<-spp.data[names(spp.data) == names(ref.data)[i]]
## 		#Skips if none
## 		if (length(spp.gene) == 0){
## 			next
## 		}
## 		#Renames
## 		names(spp.gene)<-gsub("_mito_genes.fa", "", species.names[j])
## 		final.gene<-append(final.gene, spp.gene)
## 	}#end j loop
##   
##   ##############
##   ### STEP 3.2: Sets up for alignment
##   ##############
##   #Checks for a minimum length
##   final.gene<-final.gene[width(final.gene) >= width(ref.data)[i]*as.numeric(min.prop)]
##   
##   #Checks for minimum taxa number
##   if (length(names(final.gene)) <= min.taxa){
##     print(paste(names(ref.data)[i], " had too few taxa", sep = ""))
##     next
##   }
##   
##   #Adds reference locus
##   final.gene<-append(final.gene, ref.data[i])
##   names(final.gene)[length(final.gene)]<-"Nanorana_parkeri_genome"
##   final.loci<-as.list(as.character(final.gene))
##   
##   #Saves to folder to run with mafft
##   write.fasta(sequences = final.loci, names = names(final.loci), 
##               paste("mtGenes_Fastas/", names(ref.data)[i], ".fa", sep = ""), nbchar = 1000000, as.string = T)
##   
##   #####################################
##   ### STEP 3.3: Runs MAFFT to align ###
##   #####################################
##   
##   mafft.cmd<-"mafft"
##   if (names(ref.data)[i] == "12S_rRNA" || names(ref.data)[i] == "16S_rRNA"){
##     if (secondary.structure == TRUE){ mafft.cmd<-"mafft-qinsi" } else { mafft.cmd<-"mafft" }
##   }
##   
##   #Runs the mafft command 
##   system(paste(mafft.cmd, " --localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123"," --thread ", threads, " ", "mtGenes_Fastas/", names(ref.data)[i], ".fa"," > ", "mtGenes_Fastas/", names(ref.data)[i], "_align.fa", sep = ""))
##   
##   alignment<-scanFa(FaFile(paste("mtGenes_Fastas/", names(ref.data)[i], "_align.fa", sep = "")))   # loads up fasta file
##   
##   #Reverses alignment back to correct orientation
##   reversed<-names(alignment)[grep(pattern = "_R_", names(alignment))]
##   if (length(reversed[grep(pattern = "Nanorana_parkeri_genome", reversed)]) == 1){ alignment<-reverseComplement(alignment) }
##   
##   #Renames sequences to get rid of _R_
##   names(alignment)<-gsub(pattern = "_R_", replacement = "", x = names(alignment))
##   new.align<-strsplit(as.character(alignment), "")
##   mat.align<-lapply(new.align, tolower)
##   m.align<-as.matrix(as.DNAbin(mat.align))
##   
##   #Filters out weirdly divergent sequences
##   diff<-pairwise.inf.sites(as.character(m.align), "Nanorana_parkeri_genome")
##   bad.seqs<-names(diff)[which(diff >= 0.45)]
##   rem.align<-alignment[!names(alignment) %in% bad.seqs]
##   
##   # Moves onto next loop in there are no good sequences
##   if (length(rem.align) <= as.numeric(min.taxa)){ 
##     #Deletes old files
##     system(paste("rm ", "mtGenes_Fastas/", names(ref.data)[i], "_align.fa", sep = ""))
##     print(paste(names(ref.data)[i], " had too few taxa", sep = ""))
##     next }
##   
##   ### realign if bad seqs removed
##   if (length(bad.seqs) != 0  && width(ref.data)[i] >= 200){
##     #Aligns using mafft  
##     print(paste(names(ref.data)[i], " was realigned", sep = ""))
##     
##     #Saves to folder to run with mafft
##     final.loci<-as.list(as.character(rem.align))
##     
##     #Saves to folder to run with mafft
##     write.fasta(sequences = final.loci, names = names(final.loci),paste("mtGenes_Fastas/", names(ref.data)[i], ".fa", sep = ""), nbchar = 1000000, as.string = T)
## 
## 	system(paste(mafft.cmd, " --localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123"," --thread ", threads, " ", "mtGenes_Fastas/", names(ref.data)[i], ".fa"," > ", "mtGenes_Fastas/", names(ref.data)[i], "_align.fa", sep = ""))
## 	alignment <- scanFa(FaFile(paste(out.dir, "/mtGenes_Fastas/", names(ref.data)[i], "_align.fa", sep = "")))   # loads up fasta file
## 
## 	#Reverses alignment back to correction orientation
## 	reversed <- names(alignment)[grep(pattern = "_R_", names(alignment))]
## 	if (length(reversed[grep(pattern = "Nanorana_parkeri_genome", reversed)]) == 1){
## 		alignment<-reverseComplement(alignment)
## 	}
## 	
## 	#Renames sequences to get rid of _R_
## 	names(alignment) <- gsub(pattern = "_R_", replacement = "", x = names(alignment))
## 
##   } # end bad.seqs if
##   
##   #Removes the edge gaps
##   ref.aligned<-as.character(alignment['Nanorana_parkeri_genome'])
##   not.gaps<-str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
##   ref.start<-min(not.gaps)
##   ref.finish<-max(not.gaps)
##   trim.align<-subseq(alignment, ref.start, ref.finish)
## 
##   #readies for saving
##   write.temp<-strsplit(as.character(trim.align), "")
##   aligned.set<-as.matrix(as.DNAbin(write.temp) )
##   write.phy(aligned.set, file=paste("mtGenes_Aligned/", names(ref.data)[i], ".phy", sep = ""), interleave = F)
##   
## }#end i loop
## 
## ########################################################
## ### Step 4: Create alignments and partition by codon ###
## ########################################################
## 
## 
## taxa.remove          <- c("Nanorana_parkeri_genome")  ### If you dont want to keep the reference or other taxa (only needed for step 4)
## 
## #Create directory and loci to trim
## dir.create("mtGenes_Trimmed")
## locus.names<-list.files("mtGenes_Aligned/.")
## 
## #So it doesn't trim the cds
## if(trim.cds == FALSE){
## 	no.trim <- locus.names[grep("CDS", locus.names)]
## }
## 
## #Loops through each locus and does operations on them
## for (i in 1:length(locus.names)){
## 	
## 	#############################
## 	### STEP 4.1: Basic steps ###
## 	#############################
## 	
## 	align <- readAAMultipleAlignment(file = paste("mtGenes_Aligned/", locus.names[i], sep =""), format = "phylip")  #Reads in files
## 	
## 	tax.names <- rownames(align)
## 	if(length(!tax.names %in% taxa.remove)>0){
## 		tax.names <- tax.names[!tax.names %in% taxa.remove]  ### names of taxa to keep
## 	}
## 	new.align     <- strsplit(as.character(align), "")
## 	mat.align     <- lapply(new.align, tolower)
## 	m.align       <- as.matrix(as.DNAbin(mat.align))
## 	t.align       <- m.align[rownames(m.align) %in% tax.names,]
## 	save.rownames <- rownames(t.align)
## 	
## 	if (ncol(align) <= as.numeric(min.len)){                                                            #| Removes short loci
## 		write.phy(t.align, file= paste("mtGenes_Trimmed/", locus.names[i], sep = ""), interleave = F)   #| 
## 		next
## 	}
## 
## 	if (length(grep(locus.names[i], no.trim)) != 0) {                                                    #| So it doesnt trim the cds
## 		write.phy(t.align, file= paste("mtGenes_Trimmed/", locus.names[i], sep = ""), interleave = F)    #|
## 		next                                                                                             #| 
## 	}
## 	
## 	t.loci      <- as.character(as.list(t.align))                                                      #| makes alignments with introns removed
## 	w.loci      <- lapply(t.loci, toupper)                                                             #|
## 	write.align <- lapply(w.loci, c2s)                                                                 #|
## 	
## 	input.file  <- paste("mtGenes_Trimmed/", gsub(pattern = "\\..*", "", locus.names[i]), ".fa", sep = "")         # defines output filename
## 	write.fasta(sequences = write.align, names = names(write.align),input.file, nbchar = 1000000, as.string = T)   # writes no-intron alignment
## 	
## 	#########################
## 	### STEP 4.2: GBLOCKS ###
## 	#########################
## 	if (gblocks == TRUE){
## 		system(paste("Gblocks ", input.file, " -t=d -b1=50 -b2=50 -b5=h ", sep = ""))
## 		system(paste("rm ", input.file, " ", input.file, "-gb.htm", sep = ""))
## 		system(paste("mv ", input.file, "-gb ", input.file, sep = ""))
## 	}
## 	
## 	########################
## 	### STEP 4.3: TrimAI ###
## 	########################
## 	
## 	if(trimal == TRUE){
## 		#system(paste("trimal -in ", input.file, " -out ", input.file, "-tm ","-gt 0.75 -st 0.001 -cons 60 -resoverlap 0.75 -seqoverlap 50 -automated1", sep = ""))
## 		system(paste("trimal -in ", input.file, " -out ", input.file, "-tm -automated1", sep = ""))
## 		system(paste("rm ", input.file, sep = ""))
## 		system(paste("mv ", input.file, "-tm ", input.file, sep = ""))
## 	}
## 	
## 	###################################
## 	### STEP 4.4: Save as .phy file ###
## 	###################################
## 	
## 	locus.save.name<-gsub(pattern = ".fa", replacement = ".phy", x = input.file)
## 	alignment<-scanFa(FaFile(input.file))   # loads up fasta file
## 	
## 	temp<-names(alignment)[is.na(names(alignment)) == T]
## 	if (length(temp) > 0){ break }
## 	
## 	new.names<-c()
## 	for (j in 1:length(names(alignment))){ 
## 		new.names[j]<-save.rownames[grep(pattern = names(alignment)[j], x = save.rownames)]
## 	}
## 
## 	names(alignment) <- new.names
## 	
## 	#removes loci with too few taxa
## 	if (length(names(alignment)) <= as.numeric(min.taxa)){ 
## 		system(paste("rm ", input.file, sep = ""))
## 		print(paste(input.file, "deleted. Too few taxa after trimming."))
## 		write.phy(t.align, file= paste("mtGenes_Trimmed/", locus.names[i], sep = ""), interleave = F)
## 		next
## 	}
## 	
## 	write.temp<-strsplit(as.character(alignment), "")
## 	aligned.set<-as.matrix(as.DNAbin(write.temp) )
## 	
## 	#readies for saving
## 	write.phy(aligned.set, file= locus.save.name, interleave = F)
## 	system(paste("rm ", input.file, sep = ""))
## }

#####################
### END OF SCRIPT ###
#####################









