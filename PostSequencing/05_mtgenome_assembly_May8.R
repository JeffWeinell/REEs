#######################################
### Step 1: Settings for everything ###
#######################################

# Directory settings
work.dir      <- "/Users/Jeff/Documents/SnakeCap_Data/mitogenomes"                   ### Directory where you want to do the analysis; this should contain the "Processed_Samples" directory, which should itself contain the "assembly-reads" directory
out.dir       <- "/Users/Jeff/Documents/SnakeCap_Data/mitogenomes/mtGenomes"         ### Directory where you want the mitogenomes to be sent
raw.dir       <- "assembly-reads"                                                    ### Name of folder in each species Processed_Samples directory that contains assembly reads

referencePATH <- "/Users/Jeff/Documents/SnakeCap_Data/mitogenomes/reference.fa"  ### full path to file containing many mitogenomes
geneFilePATH  <- "/Users/Jeff/Documents/SnakeCap_Data/mitogenomes/mtGenes.fa"    ### full path to file containaining annotations (~37 genes) for each mitogenome

reference     <- "reference.fa"                                                  ### Name of the file containing many mitogenomes (this should be in work.dir)
gene.file     <- "mtGenes.fa"                                                    ### Name of the file containing many (~37) genes from mitogenomes (this should be in work.dir)

processed.dir <- "/Volumes/MyPassport/SequenceCapture/Results/Step2_contigs-spades/Processed_Samples/"  ### Processed_Samples director containing folders for each species, each of which contains the assembly-reads folder containing assembly reads

bbmapPATH     <- "/Applications/bbmap/bbmap.sh"                      ### path th bbmap.sh file
spadesPATH    <- "/Applications/SPAdes-3.13.0-Darwin/bin/spades.py"  ### path the spades.py file
cap3PATH      <- "/Applications/cap3.macosx.intel64/cap3"            ### path to cap3 program

# spadesOut   <- "/Users/Jeff/Documents/SnakeCap_Data/mitogenomes/spades"
# Run settings

threads              <- 8                             ### Number of threads
mem                  <- "8"                           ### GB of ram
min.id               <- "0.7"                         ### Initial value for matching raw reads to reference. Probably should leave alone, for problem taxa.

# Alignment settings
secondary.structure  <- TRUE                          ### If true, runs mafft-qinsi on mt regions that have secondary structure. Takes structure into acct.
min.taxa             <- 3                             ### min number of individuals to keep an alignment
min.prop             <- "0.25"                        ### min coverage per individual. e.g., if set to "0.25", for a 100bp gene, needs 25 bp to keep.
taxa.remove          <- c("Nanorana_parkeri_genome")  ### If you dont want to keep the reference or other taxa (only needed for step 4)
min.len              <- "100"                         ### min length for trimming. Set to this value as you dont usually want to trim t-RNAs
trim.cds             <- FALSE                         ### defaults to no trimming for coding sequence. Usually destroys mtGenes
gblocks              <- FALSE                         ### If you want to use 
trimal               <- TRUE                          ### If you want to use

###################################################
### Shouldnt need to change anything below here ###
###################################################

######################
### Load Libraries ###
######################

library(ape)
library(seqinr)
library(stringr)
library(data.table)
library(GenomicRanges)
library(Biostrings)
library(Rsamtools)

#Options
options(stringsAsFactors = FALSE)
options(warn=2) #for debugging warnings in loops

### load custom function
source("/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/SnakeCap_functions.R") ### loads in some necessary functions

###################################################################
### Step 1: Gather read data and assemble mitochondrial genomes ###
###################################################################

# Creates necessary directories if they dont exist

if(file.exists(out.dir)==F){
	dir.create(out.dir)
}
if(file.exists(paste(out.dir, "/", "Species_mtGenomes", sep = ""))==F){
	dir.create(paste(out.dir, "/", "Species_mtGenomes", sep = ""))
}

# Sets up the reads

setwd(processed.dir)
files   <- list.files(path = processed.dir, full.names = F, recursive = T)  ### all file names in all folders and subfolders of the "Processed_Samples" directory
samples <- list.dirs(path = processed.dir, full.names = F, recursive = F)   ### all folder names in the "Processed_Samples" directory (should be a folder for each sample)
reads   <- files[grep(pattern = raw.dir, x = files)]                        ### reads files that were used for assembly, these have "assembly-reads" as part of their name

setwd(out.dir)

for (i in 1:length(samples)){
#for (i in iVals){

	#Change to main directory
	
	setwd(out.dir)
		
	sample.reads <- reads[grep(samples[i], reads)]
	read1        <- paste(processed.dir, sample.reads[grep("READ1", sample.reads)], sep = "")
	read2        <- paste(processed.dir, sample.reads[grep("READ2", sample.reads)], sep = "")
	read3        <- paste(processed.dir, sample.reads[grep("singleton", sample.reads)], sep = "")
	
	#Pick out matching reads to mt Genomes
	system(paste(bbmapPATH," -Xmx8g ref=", referencePATH," in1=", read1, " in2=", read2, " vslow k=12 minid=",min.id," outm1=read1.fq outm2=read2.fq", sep = ""), ignore.stderr = T)
	system(paste(bbmapPATH," -Xmx8g ref=", referencePATH," in=", read3, " vslow k=12 minid=", min.id," outm=singleton.fq", sep = ""), ignore.stderr = T)
	
	system("touch current_seed.fasta")
	new.len        <- 0
	counter        <- 0
	repeat.counter <- 0
	seeding        <- T

	while (seeding == T){
		#Copy new reference to do recursively
		counter  <- counter+1
		prev.len <- new.len
	
		# Skips the first one since its already done
		if (counter >= 2){
			#Pick out matching reads to mt Genomes
			system(paste(bbmapPATH," -Xmx8g ref=current_seed.fasta", " in1=", read1, " in2=", read2, " vslow k=12 minid=",min.id," outm1=t_read1.fq outm2=t_read2.fq", sep = ""), ignore.stderr = T)
			system(paste(bbmapPATH," -Xmx8g ref=current_seed.fasta", " in=", read3, " vslow k=12 minid=", min.id," outm=t_singleton.fq", sep = ""), ignore.stderr = T)
			system(paste("cat t_read1.fq o_read1.fq >> read1.fq"))
			system(paste("cat t_read2.fq o_read2.fq >> read2.fq"))
			system(paste("cat t_singleton.fq o_singleton.fq >> singleton.fq"))
			system("rm t_read1.fq t_read2.fq t_singleton.fq")
		}

		# Run SPADES on sample
		k     <- c(9,13,21,33,55,77,99,127)
		k.val <- paste(k, collapse = ",")
		system(paste(spadesPATH," --pe1-1 /Users/Jeff/Documents/SnakeCap_Data/mitogenomes/mtGenomes/read1.fq --pe1-2 /Users/Jeff/Documents/SnakeCap_Data/mitogenomes/mtGenomes/read2.fq --pe1-s /Users/Jeff/Documents/SnakeCap_Data/mitogenomes/mtGenomes/singleton.fq"," -o /Users/Jeff/Documents/SnakeCap_Data/mitogenomes/mtGenomes/spades -k ",k.val," --careful -t ", threads, " -m ", mem, sep = ""), ignore.stdout = T) 
		
		#/Applications/SPAdes-3.13.0-Darwin/bin/spades.py --pe1-1 /Users/Jeff/Documents/SnakeCap_Data/mitogenomes/mtGenomes/read1.fq --pe1-2 /Users/Jeff/Documents/SnakeCap_Data/mitogenomes/mtGenomes/read2.fq --pe1-s /Users/Jeff/Documents/SnakeCap_Data/mitogenomes/mtGenomes/singleton.fq -o /Users/Jeff/Documents/SnakeCap_Data/mitogenomes/mtGenomes/spades -k 9 --careful -t 8 -m 8
		
		
		#Checks to see if one kmer failed or not
		while (file.exists("spades/contigs.fasta") == F){
			#subtract Ks until it works
			system("rm -r spades")
			k <- k[-length(k)]
			if(length(k) == 0) { break }
			k.val  <- paste(k, collapse = ",")
			min.id <- "0.6"
			system(paste(spadesPATH," --pe1-1 read1.fq --pe1-2 read2.fq --pe1-s singleton.fq"," -o spades -k ",k.val," --careful -t ", threads, " -m ", mem, sep = ""), ignore.stdout = T)
		}#end while

		# If the k-mers are all run out, therefore nothing can be assembled
		if (length(k) == 0) { 
			paste("k-mer values all used up, cannot assemble!")
			system("rm read1.fq read2.fq singleton.fq t_read1.fq t_read2.fq t_singleton.fq o_read1.fq o_read2.fq o_singleton.fq")
			system("rm -r spades")
			seeding = F 
		}# end if 
	
		if (counter == 1){
			system(paste("mv read1.fq o_read1.fq"))
			system(paste("mv read2.fq o_read2.fq"))
			system(paste("mv singleton.fq o_singleton.fq"))
		}
	
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

		print(paste("iteration ", counter, " complete!", sep = ""))
		print(paste("new length: ", new.len, ". Old length: ", prev.len, sep = ""))

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
			temp.assembled <- scanFa(FaFile(paste("current_seed.fasta.cap.contigs", sep = "")))
			temp.singlets  <- scanFa(FaFile(paste("current_seed.fasta.cap.singlets", sep = "")))
			keep.singlets  <- temp.singlets[width(temp.singlets) >= 100]
			final.save     <- append(temp.assembled, keep.singlets)
		
			#Writes contigs for cap3
			write.loci <- as.list(as.character(final.save))
			write.fasta(sequences = write.loci, names = names(write.loci),"current_seed.fasta", nbchar = 1000000, as.string = T)
		
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

	#Save finsihed genome
	contigs <- scanFa(FaFile("current_seed.fasta"))   # loads up fasta file
	
	#Skips if there are none
	if(length(contigs) == 0){
		next
	}

	#Trys to merge contigs if there are more than 1
	if(length(contigs) >= 2){
		#runs cap3 to merge similar contigs (pull only clustered contigs out?)
		system(paste(cap3PATH," current_seed.fasta -z 1 -o 16 -e 11 -s 251", " > ","log.fasta.cap.txt", sep = "")) 
		
		#Reads in results files
		temp.assembled<-scanFa(FaFile(paste("current_seed.fasta.cap.contigs", sep = "")))
		temp.singlets<-scanFa(FaFile(paste("current_seed.fasta.cap.singlets", sep = "")))
		keep.singlets<-temp.singlets[width(temp.singlets) >= 100]
		contigs<-append(temp.assembled, keep.singlets)
		
		#Get cap3 files and deletes
		cap.files<-list.files(pattern = "", full.names = F, recursive = F)
		cap.remove<-cap.files[grep(pattern = paste("fasta.cap*.", sep =""), x = cap.files)]
		system(paste("rm ", paste(cap.remove, collapse = " ") ))
	}#end if

	if(sum(width(contigs)) <= 1000) { 
		print("less than 1000bp, not enough data to extract")
		next 
	}

	#Writes the full mitochondrial genome file
	system("rm current_seed.fasta")
	names(contigs)<- paste("sequence_", seq(1:length(contigs)), sep = "")
	write.loci    <- as.list(as.character(contigs))
	write.fasta(sequences = write.loci, names = names(write.loci),paste("Species_mtGenomes/", samples[i], ".fa", sep = ""), nbchar = 1000000, as.string = T)

}#end i loop

system("rm -r ref")


###########################################################################
### Step 2: Assess completeness of the mitochondrial genome and annotate ##
###########################################################################

#PSLX headers
headers<-c("matches", "misMatches", "repMatches", "nCount", "qNumInsert", "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName", 
           "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd", "blockCount", "blockSize", "qStarts", "tStarts", "qSeq", "tSeq")

#Creates new directory and enters this working directory
setwd(out.dir)
dir.create("Species_Loci")
spp.samples<-list.files("Species_mtGenomes/.")
spp.samples<-gsub(".fa$", "", spp.samples)

for (i in 1:length(spp.samples)){

	#Load in the data
	contigs<-scanFa(FaFile(paste("Species_mtGenomes/", spp.samples[i], ".fa", sep = "")))   # loads up fasta file

	#Matches samples to loci
	system(paste("mpirun pblat -threads=", threads, " Species_mtGenomes/", spp.samples[i], ".fa ",gene.file, " -tileSize=8 -minIdentity=60"," -noHead -out=pslx mt_to_genes.pslx", sep = ""), ignore.stdout = T)
	
	#Need to load in transcriptome for each species and take the matching transcripts to the database
	temp.count<-scan(file = "mt_to_genes.pslx", what = "character")
	if (length(temp.count) == 0){
		print("No matching mitochondrial genes were found.")
		next
	}
	match.data<-fread("mt_to_genes.pslx", sep = "\t", header = F, stringsAsFactors = FALSE)
	setnames(match.data, headers)

	loci.names<-unique(match.data$qName)
	sep.loci<-DNAStringSet()
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
		sep.loci       <- append(sep.loci, new.seq)
	}#end j loop

	#Writes the full mitochondrial genome file
	write.loci<-as.list(as.character(sep.loci))
	write.fasta(sequences = write.loci, names = names(write.loci),paste("Species_Loci/", spp.samples[i], "_mito_genes.fa", sep = ""), nbchar = 1000000, as.string = T)

	system("rm mt_to_genes.pslx")

}#end i loop

############################################
### Step 3: Create mitogenome alignments ###
############################################

setwd(out.dir)

#Sets up the loci to align
ref.data      <- scanFa(FaFile(gene.file))
species.names <- list.files("Species_Loci/.", full.names = F)
species.names <- species.names[species.names != ""]
dir.create("mtGenes_Fastas")
dir.create("mtGenes_Aligned")

#Aligns each potential locus

for(i in 1:length(ref.data)){
	#######################################################
	### STEP 3.1: Gets the locus data from each species ###
	#######################################################

	#Gets all species data
	final.gene<-DNAStringSet()
	for (j in 1:length(species.names)){
		#Looks for this gene in the species data
		spp.data<-scanFa(FaFile(paste("Species_Loci/", species.names[j], sep = "")))   # loads up fasta file
		spp.gene<-spp.data[names(spp.data) == names(ref.data)[i]]
		#Skips if none
		if (length(spp.gene) == 0){
			next
		}
		#Renames
		names(spp.gene)<-gsub("_mito_genes.fa", "", species.names[j])
		final.gene<-append(final.gene, spp.gene)
	}#end j loop
  
  ##############
  ### STEP 3.2: Sets up for alignment
  ##############
  #Checks for a minimum length
  final.gene<-final.gene[width(final.gene) >= width(ref.data)[i]*as.numeric(min.prop)]
  
  #Checks for minimum taxa number
  if (length(names(final.gene)) <= min.taxa){
    print(paste(names(ref.data)[i], " had too few taxa", sep = ""))
    next
  }
  
  #Adds reference locus
  final.gene<-append(final.gene, ref.data[i])
  names(final.gene)[length(final.gene)]<-"Nanorana_parkeri_genome"
  final.loci<-as.list(as.character(final.gene))
  
  #Saves to folder to run with mafft
  write.fasta(sequences = final.loci, names = names(final.loci), 
              paste("mtGenes_Fastas/", names(ref.data)[i], ".fa", sep = ""), nbchar = 1000000, as.string = T)
  
  #####################################
  ### STEP 3.3: Runs MAFFT to align ###
  #####################################
  
  mafft.cmd<-"mafft"
  if (names(ref.data)[i] == "12S_rRNA" || names(ref.data)[i] == "16S_rRNA"){
    if (secondary.structure == TRUE){ mafft.cmd<-"mafft-qinsi" } else { mafft.cmd<-"mafft" }
  }
  
  #Runs the mafft command 
  system(paste(mafft.cmd, " --localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123"," --thread ", threads, " ", "mtGenes_Fastas/", names(ref.data)[i], ".fa"," > ", "mtGenes_Fastas/", names(ref.data)[i], "_align.fa", sep = ""))
  
  alignment<-scanFa(FaFile(paste("mtGenes_Fastas/", names(ref.data)[i], "_align.fa", sep = "")))   # loads up fasta file
  
  #Reverses alignment back to correct orientation
  reversed<-names(alignment)[grep(pattern = "_R_", names(alignment))]
  if (length(reversed[grep(pattern = "Nanorana_parkeri_genome", reversed)]) == 1){ alignment<-reverseComplement(alignment) }
  
  #Renames sequences to get rid of _R_
  names(alignment)<-gsub(pattern = "_R_", replacement = "", x = names(alignment))
  new.align<-strsplit(as.character(alignment), "")
  mat.align<-lapply(new.align, tolower)
  m.align<-as.matrix(as.DNAbin(mat.align))
  
  #Filters out weirdly divergent sequences
  diff<-pairwise.inf.sites(as.character(m.align), "Nanorana_parkeri_genome")
  bad.seqs<-names(diff)[which(diff >= 0.45)]
  rem.align<-alignment[!names(alignment) %in% bad.seqs]
  
  # Moves onto next loop in there are no good sequences
  if (length(rem.align) <= as.numeric(min.taxa)){ 
    #Deletes old files
    system(paste("rm ", "mtGenes_Fastas/", names(ref.data)[i], "_align.fa", sep = ""))
    print(paste(names(ref.data)[i], " had too few taxa", sep = ""))
    next }
  
  ### realign if bad seqs removed
  if (length(bad.seqs) != 0  && width(ref.data)[i] >= 200){
    #Aligns using mafft  
    print(paste(names(ref.data)[i], " was realigned", sep = ""))
    
    #Saves to folder to run with mafft
    final.loci<-as.list(as.character(rem.align))
    
    #Saves to folder to run with mafft
    write.fasta(sequences = final.loci, names = names(final.loci),paste("mtGenes_Fastas/", names(ref.data)[i], ".fa", sep = ""), nbchar = 1000000, as.string = T)

	system(paste(mafft.cmd, " --localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123"," --thread ", threads, " ", "mtGenes_Fastas/", names(ref.data)[i], ".fa"," > ", "mtGenes_Fastas/", names(ref.data)[i], "_align.fa", sep = ""))
	alignment <- scanFa(FaFile(paste(out.dir, "/mtGenes_Fastas/", names(ref.data)[i], "_align.fa", sep = "")))   # loads up fasta file

	#Reverses alignment back to correction orientation
	reversed <- names(alignment)[grep(pattern = "_R_", names(alignment))]
	if (length(reversed[grep(pattern = "Nanorana_parkeri_genome", reversed)]) == 1){
		alignment<-reverseComplement(alignment)
	}
	
	#Renames sequences to get rid of _R_
	names(alignment) <- gsub(pattern = "_R_", replacement = "", x = names(alignment))

  } # end bad.seqs if
  
  #Removes the edge gaps
  ref.aligned<-as.character(alignment['Nanorana_parkeri_genome'])
  not.gaps<-str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
  ref.start<-min(not.gaps)
  ref.finish<-max(not.gaps)
  trim.align<-subseq(alignment, ref.start, ref.finish)

  #readies for saving
  write.temp<-strsplit(as.character(trim.align), "")
  aligned.set<-as.matrix(as.DNAbin(write.temp) )
  write.phy(aligned.set, file=paste("mtGenes_Aligned/", names(ref.data)[i], ".phy", sep = ""), interleave = F)
  
}#end i loop

########################################################
### Step 4: Create alignments and partition by codon ###
########################################################

#Create directory and loci to trim
dir.create("mtGenes_Trimmed")
locus.names<-list.files("mtGenes_Aligned/.")

#So it doesn't trim the cds
if(trim.cds == FALSE){
	no.trim <- locus.names[grep("CDS", locus.names)]
}

#Loops through each locus and does operations on them
for (i in 1:length(locus.names)){
	
	#############################
	### STEP 4.1: Basic steps ###
	#############################
	
	align <- readAAMultipleAlignment(file = paste("mtGenes_Aligned/", locus.names[i], sep =""), format = "phylip")  #Reads in files
	
	tax.names <- rownames(align)
	if(length(!tax.names %in% taxa.remove)>0){
		tax.names <- tax.names[!tax.names %in% taxa.remove]  ### names of taxa to keep
	}
	new.align     <- strsplit(as.character(align), "")
	mat.align     <- lapply(new.align, tolower)
	m.align       <- as.matrix(as.DNAbin(mat.align))
	t.align       <- m.align[rownames(m.align) %in% tax.names,]
	save.rownames <- rownames(t.align)
	
	if (ncol(align) <= as.numeric(min.len)){                                                            #| Removes short loci
		write.phy(t.align, file= paste("mtGenes_Trimmed/", locus.names[i], sep = ""), interleave = F)   #| 
		next
	}

	if (length(grep(locus.names[i], no.trim)) != 0) {                                                    #| So it doesnt trim the cds
		write.phy(t.align, file= paste("mtGenes_Trimmed/", locus.names[i], sep = ""), interleave = F)    #|
		next                                                                                             #| 
	}
	
	t.loci      <- as.character(as.list(t.align))                                                      #| makes alignments with introns removed
	w.loci      <- lapply(t.loci, toupper)                                                             #|
	write.align <- lapply(w.loci, c2s)                                                                 #|
	
	input.file  <- paste("mtGenes_Trimmed/", gsub(pattern = "\\..*", "", locus.names[i]), ".fa", sep = "")         # defines output filename
	write.fasta(sequences = write.align, names = names(write.align),input.file, nbchar = 1000000, as.string = T)   # writes no-intron alignment
	
	#########################
	### STEP 4.2: GBLOCKS ###
	#########################
	if (gblocks == TRUE){
		system(paste("Gblocks ", input.file, " -t=d -b1=50 -b2=50 -b5=h ", sep = ""))
		system(paste("rm ", input.file, " ", input.file, "-gb.htm", sep = ""))
		system(paste("mv ", input.file, "-gb ", input.file, sep = ""))
	}
	
	########################
	### STEP 4.3: TrimAI ###
	########################
	
	if(trimal == TRUE){
		#system(paste("trimal -in ", input.file, " -out ", input.file, "-tm ","-gt 0.75 -st 0.001 -cons 60 -resoverlap 0.75 -seqoverlap 50 -automated1", sep = ""))
		system(paste("trimal -in ", input.file, " -out ", input.file, "-tm -automated1", sep = ""))
		system(paste("rm ", input.file, sep = ""))
		system(paste("mv ", input.file, "-tm ", input.file, sep = ""))
	}
	
	###################################
	### STEP 4.4: Save as .phy file ###
	###################################
	
	locus.save.name<-gsub(pattern = ".fa", replacement = ".phy", x = input.file)
	alignment<-scanFa(FaFile(input.file))   # loads up fasta file
	
	temp<-names(alignment)[is.na(names(alignment)) == T]
	if (length(temp) > 0){ break }
	
	new.names<-c()
	for (j in 1:length(names(alignment))){ 
		new.names[j]<-save.rownames[grep(pattern = names(alignment)[j], x = save.rownames)]
	}

	names(alignment) <- new.names
	
	#removes loci with too few taxa
	if (length(names(alignment)) <= as.numeric(min.taxa)){ 
		system(paste("rm ", input.file, sep = ""))
		print(paste(input.file, "deleted. Too few taxa after trimming."))
		write.phy(t.align, file= paste("mtGenes_Trimmed/", locus.names[i], sep = ""), interleave = F)
		next
	}
	
	write.temp<-strsplit(as.character(alignment), "")
	aligned.set<-as.matrix(as.DNAbin(write.temp) )
	
	#readies for saving
	write.phy(aligned.set, file= locus.save.name, interleave = F)
	system(paste("rm ", input.file, sep = ""))
}

#####################
### END OF SCRIPT ###
#####################
