library(ape)
library(seqinr)
library(stringr)
library(GenomicRanges)
library(Biostrings)
library(Rsamtools)
library(data.table)

options(stringsAsFactors = FALSE)
#options(warn=2) #for debugging warnings in loops

### load custom functions
source("/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/SnakeCap_functions.R") ### loads in some necessary functions

###################################
## Step 0 parameter and settings ##
###################################

threads <-8
resume  <- "yes" #If it crashes it will resume 

###### Parameter setup #########
min.taxa  <- 4    ### min number of taxa to write a trimAL alignment ### Carl used 3 here
min.aln   <- 80   ### minimum width of an alignment (also used in setup step as the minimum number of bases that an individual must have in an alignment to keep the individual)
min.len   <- 50   ### min size for an individual sample
edge.trim <- 0.50 ### proportion of taxa needed to keep column on ends
min.cov   <- 0.30 ### min percent of sample that must overlap with consensus

#If you want to remove any taxa (this can be empty)
taxa.remove  <- c()

### HOME VERSION
#work.dir     <- "/Users/chutter/Dropbox/Research/WIP/Anura_Phylogeny/Alignments"
#uce.file     <- "/Users/chutter/Dropbox/Research/WIP/FrogCap_Pipeline/Source_Files/Hutter_uce5k_loci.fa"
#probe.file   <- "/Users/chutter/Dropbox/Research/WIP/FrogCap_Pipeline/Source_Files/All_Combined_Loci_July24.fa"
#legacy.file  <- "/Users/chutter/Dropbox/Research/WIP/FrogCap_Pipeline/Source_Files/Final_Legacy_Consensus_Aug4.fa"

### CLUSTER VERSION
#work.dir     <- "/home/c111h652/scratch/Anura_Full/Alignments_all"
#uce.file     <- "/home/c111h652/scratch/Anura_Full/Hutter_uce5k_loci.fa"
#probe.file   <- "/home/c111h652/scratch/Anura_Full/All_Combined_Loci_July24.fa"
#legacy.file  <- "/home/c111h652/scratch/Anura_Full/Final_Legacy_Consensus_Aug4.fa"

work.dir                  <- "/Users/Jeff/Documents/SnakeCap_Data/Alignments/"            ### path to the directory containing the folder with input alignments. trimAL alignments will be placed in a folder in this directory.

### Upstream_noncoding_06-trimmed/alignmentFiles

input.alignments.dir      <- "All_data_with-Ophiophagus_someFlankingRegion_trimmed-somehow/alignmentFiles/"  ### name of folder containing input alignments
input.alignments.dir.path <- paste0(work.dir, input.alignments.dir)                                          ### path to folder containing input alignments "/Users/Jeff/Documents/SnakeCap_Data/Alignments/All-markers_withFlankingRegion_no-trimAL"
trimal.dir.path           <- "/Users/Jeff/miniconda3/pkgs/trimal-1.4.1-h04f5b5a_3/bin"                       ### path to the directory containing the trimal executable
trim.dir                  <- ""                                                                              ### name of folder containing output alignments
trim.dir.path             <- paste0(work.dir, "/", trim.dir)                                                 ### path to folder containing output alignments

#/Users/Jeff/Documents/SnakeCap_Data/Alignments/All_data_with-Ophiophagus_someFlankingRegion_trimmed-somehow/alignmentFiles/WeinellEntry1.phy
#/Users/Jeff/Documents/SnakeCap_Data/Alignments/Upstream_noncoding/alignmentFiles/WeinellEntry1.phy
###################################################
# Step 1. Trim the input alignments using trimAL ##
###################################################

dir.check.create(trim.dir.path)                                             ### creates the trimAL directory if it doesnt already exist
locus.names           <- list.files(input.alignments.dir.path)              ### filenames of the input alignments
input.alignment.paths <- list.files(input.alignments.dir.path,full.names=T) ### full paths to the input alignments

if (resume == "yes" & length(list.files(trim.dir.path))>0){
  done        <- list.files(trim.dir.path)                                  ### filenames of the output (trimAL) alignments already completed
  locus.names <- locus.names[!locus.names %in% done]                        ### input alignments not yet processed
}

#Loops through each locus
#for (i in 1:length(locus.names)){
for (i in 1:10){
  ##############
  #STEP 1: Setup steps
  ##############
  # Reads in files
  setwd(input.alignments.dir.path) ### sets the current directory as the directory containing the input alignments
  
  align <- DNAStringSet(readAAMultipleAlignment(file = input.alignment.paths[i], format = "phylip"))   ### reads in ith input alignment

  # Skips locus if the width of the input alignment is less than min.aln
  if(max(width(align)) <= as.numeric(min.aln)){
  	next
  }
  
  # Remove gap-only alignments
  gap.align <- strsplit(as.character(align), "")                             ### a list of character vectors, each containing the DNA sequence of an individual
  gap.count <- unlist(lapply(gap.align, function(x) length(x[x != "-"])))    ### number of non-gap characters per individual in alignment
  gap.rem   <- gap.count[gap.count <= as.numeric(min.aln)]                   ### number of non-gap characters for individuals with fewer non-gap characters than value of min.aln paramaeter
  
  # Remove the taxa that need to go
  rem.taxa  <- c(names(gap.rem), taxa.remove)                                ### names of individuals to be removed because to few non-gap characters
  rem.align <- align[!names(align) %in% rem.taxa]                            ### alignment after filtering individuals in rem.taxa

  # Skips locus if too few taxa
  if (length(rem.align) <= as.numeric(min.taxa)){
  	next
  }

  ############################################################################
  ## STEP 2: Runs trimming programs (3 different trimming methods are used) ##
  ############################################################################
  
  #A. Runs trimAL

  setwd(trim.dir.path)                                                 ### sets current directory as the directory where trimAL alignments will be saved
  trimal.align <- run.trimal(rem.align,trimal.exe.dir=trimal.dir.path) ### runs trimAL on rem.align using the -automated1 algorithm (which chooses between -gappyout and -strict methods)

##	This might be useful to do:
##	runs trimAL in a loop until input and output alignments are the same
##	current.align <- rem.align
##	while(max(width(trimal.align)) < max(width(current.align))){
##		current.align <- trimal.align
##		trimal.align  <- run.trimal(current.align,trimal.exe.dir=trimal.dir.path)
##	}
##	
  # skips locus if the width of the trimAL alignment is less than min.aln
  if(max(width(trimal.align)) <= as.numeric(min.aln)){
  	next
  }
  
  #B Replace poorly aligned subsequences within each window with an equal length string of gaps "-"
  red.align <- slice.trim(trimal.align, slice.size.bp = 80, threshold = 0.40)
  
  # Skips locus if too few taxa or alignment width too short
  if (length(red.align) <= as.numeric(min.taxa) | max(width(red.align)) <= as.numeric(min.aln)){
  	next
  }
  
  #C. Trims the alignment ends until the fraction of non-gap characters in each end column is greater than defined by edge.trim
  edge.align <- trim.ends(red.align, min.n.seq = ceiling(length(red.align) * edge.trim), codon.trim = F)
  
  #removes too short loci
  if (max(width(edge.align)) <= as.numeric(min.aln)){
  	next
  }
  #Gets consensus seq for trimming more
  con.seq <- make.consensus(edge.align, method = "majority")
  
  #Removes the edge gaps
  ok.seq <- gsub("\\+|-", "", as.character(con.seq))  ## deletes gaps in the consensus sequence
  
  # Skips locus if consensus sequence is less than min.aln
  # This means that the number of sites with more than 50% coverage among individuals must be > min.aln
  if ( nchar(ok.seq) <= as.numeric(min.aln)){
  	next
  }

  ##############
  #STEP 3: Cleanup and save
  ##############
  #Skips locus if too few taxa in alignment
  if (length(names(edge.align)) <= as.numeric(min.taxa)){ 
    print(paste0(locus.names[i], " deleted. Too few taxa after trimming.") )
    next
  }
  
  write.temp  <- strsplit(as.character(edge.align), "")  ### edge.align formatted as a list of character vectors
  aligned.set <- as.matrix(as.DNAbin(write.temp))        ### edge.align formatted as a DNAbin matrix
  
  #removes too short loci
  if (ncol(aligned.set) <= as.numeric(min.aln)){ 
    print(paste(locus.names[i], " deleted. Trimmed alignment length below threshold.", sep = "") )
    next 
  }
  
  #Removes samples that are too short individually
  len.temp <- as.character(as.list(aligned.set))                      ### aligned sequences as a list of character vectors
  len.loci <- lapply(len.temp, function (x) x[x != "-"])              ### for each individual, the set of non-gap bases
  spp.len  <- unlist(lapply(len.loci, function (x) length(x)))        ### number of non-gap bases per individual
  spp.rem  <- spp.len[spp.len <= max(spp.len) * as.numeric(min.cov)]  ### individuals with too low coverage
  spp.rem  <- append(spp.rem, spp.len[spp.len <= min.len])            ### individuals with too low of coverage or too few non-gap bases
  if (length(spp.rem) > 0){
  	aligned.set<-aligned.set[!rownames(aligned.set) %in% unique(names(spp.rem)),]
  }
  
  #Skips locus if too few taxa in alignment
  if (length(rownames(aligned.set)) <= as.numeric(min.taxa)){ 
    print(paste(locus.names[i], " deleted. Too few taxa after trimming.", sep = ""))
    next
  }
  
  #saves trimmed alignment of ith locus
  write.phy(aligned.set, file= paste0(gsub(pattern = "\\..*", "", locus.names[i]), ".phy"), interleave = F)

}


#### JLW verified everything works at least up to here.

#######################################################################
#######################################################################
### Step 2. Verify exons are present in loci expected to have exons ###
###         Trim to make CDS-only alignments.                       ###
###         Output exon table                                       ###
#######################################################################
#######################################################################


#Sets up new directory for this stuff
trim.dir    <-"exon-only_trimmed"
dir.create(paste0(work.dir, "/", trim.dir))
prot.dir    <-"exon-protein_trimmed"
dir.create(paste0(work.dir, "/", prot.dir))
setwd(paste0(work.dir, "/exon-only_untrimmed"))
locus.names <-list.files(".")

#Checks if to resume or not
if (resume == "yes"){
  done<-list.files(paste0(work.dir, "/", trim.dir))
  locus.names<-locus.names[!locus.names %in% done]
}

#Loops through each locus and does operations on them
for (i in 1:length(locus.names)){
  
  ##############
  #STEP 1: Load in data and basic filtering
  ##############
  #Skip if a UCE
  if (length(grep("uce", locus.names[i])) == 1){ next }
  
  #Reads in files
  setwd(paste0(work.dir, "/exon-only_untrimmed"))
  align<-DNAStringSet(readAAMultipleAlignment(file = locus.names[i], format = "phylip"))

  #Remove gap only sequences
  gap.align<-strsplit(as.character(align), "")
  gap.count<-unlist(lapply(gap.align, function(x) length(x[x != "-"]) ) )
  gap.rem<-gap.count[gap.count <= as.numeric(min.aln)]
  #Remove the taxa that need to go
  rem.taxa<-c(names(gap.rem), taxa.remove)
  red.align<-align[!names(align) %in% rem.taxa]
  
  #removes loci with too few taxa
  if ( length(red.align) <= as.numeric(min.taxa)){ next }
  #removes too short loci
  if (max(width(red.align)) <= as.numeric(min.aln)){ next }
  
  ##############
  #STEP 2: Trims the alignment 
  ##############  
  trimmed<-trim.ends(red.align, min.n.seq = ceiling(length(red.align) * edge.trim), codon.trim = T)
  save.names<-names(trimmed)
  
  #Gets consensus seq for trimming more
  con.seq<-make.consensus(trimmed, method = "majority")
  
  #Removes the edge gaps
  ref.aligned<-as.character(con.seq)
  not.gaps<-str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]

  #Finds weird gaps to fix
  temp.gaps<-as.numeric(1)
  for (k in 1:length(not.gaps)-1){ temp.gaps<-append(temp.gaps, not.gaps[k+1]-not.gaps[k]) }
  temp.gaps<-temp.gaps-1
  names(temp.gaps)<-not.gaps
  gap.spots<-temp.gaps[temp.gaps %% 3 != 0]
  
  del.col<-c()
  if (length(gap.spots) != 0){
    #Loop through each potential bad thing and fix
    for (k in 1:length(gap.spots)){
      del.col<-append(del.col, (as.numeric(names(gap.spots[k]))-gap.spots[k]):(as.numeric(names(gap.spots[k]))-1))
    }
  }#end if
  
  ##############
  #STEP 3: Fixes large gaps at ends of alignment
  ##############  
  #Looks for gaps that clearly wrong and not 3 BP 
  new.align<-strsplit(as.character(trimmed), "")
  x<-as.matrix(as.DNAbin(new.align))
  
  rem.n<-c()
  for (k in 1:ncol(x)){
    gaps<-table(as.character(x[,k]))
    per.gaps<-gaps[names(gaps) == "-"]/nrow(x)
    
    if (length(per.gaps) == 0){ next }
    
    #Records column when the gaps exceed this percentage
    if (per.gaps >= 0.75){ del.col<-append(del.col, k) }
  
    #Removes gap columns only consisting of Ns 
    n.gaps<-gaps[names(gaps) != "-"]
    if (length(n.gaps) == 1){
      if (names(n.gaps) == "n"){ rem.n<-append(rem.n, k)}
    }
    
  }#end k loop

  #combines columns to be deleted
  fin.del<-c(rem.n, del.col)
  if (length(fin.del) != 0){ x<-x[,-fin.del] }
  #Removes bad columsn and coverts alignment back to DNASTringSet
  char.align<-as.list(data.frame(t(as.character(x))))
  temp.align<-lapply(char.align, FUN = function(x) paste(x, collapse = ""))
  trimmed<-DNAStringSet(unlist(temp.align))
  names(trimmed)<-save.names

  ##############
  #STEP 4: Gathers table of best and longest stop codon free frames for each seq
  ##############  
  #Checks to make sure the codon position is correct
  save.frame<-data.frame()
  save.all<-data.frame()
  for (j in 1:length(trimmed)){
    #Finds open reading frames
    temp.codon<-find.orf(trimmed[j], codons = F, min.size = 80 )
    if(nrow(temp.codon) == 0){
      samp.frame<-cbind(Sample = names(trimmed[j]), FrameStart = 0, FrameEnd = 0, Size = 0, sppSize = 0, Frame = "0")
      save.frame<-rbind(save.frame, samp.frame)
      next
    }
    
    all.frame<-temp.codon[temp.codon$Size >= max(temp.codon$Size) * .70,]
    big.frame<-temp.codon[temp.codon$Size == max(temp.codon$Size),]

    if (nrow(big.frame) >= 2){
      #Picks the best from this order of things
      temp.stop<-big.frame[big.frame$Frame == "F1",]
      if (nrow(temp.stop) == 0){ temp.stop<-big.frame[big.frame$Frame == "R1",] }
      if (nrow(temp.stop) == 0){ temp.stop<-big.frame[big.frame$Frame == "F2",] }
      if (nrow(temp.stop) == 0){ temp.stop<-big.frame[big.frame$Frame == "R2",] }
      if (nrow(temp.stop) == 0){ temp.stop<-big.frame[big.frame$Frame == "F3",] }
      if (nrow(temp.stop) == 0){ temp.stop<-big.frame[big.frame$Frame == "R3",] }
      big.frame<-temp.stop
    }

    # #Saves teh data
    samp.frame<-cbind(Sample = names(trimmed[j]), big.frame)
    temp.size<-unlist(strsplit(as.character(trimmed[j]), ""), use.names = F)
    
    #Starts from the beginning and end to fill in end gaps
    sub.size<-0
    for (q in 1:length(temp.size)){ if (temp.size[q] == "-"){ sub.size<-sub.size+1 } else { break } }
    for (q in length(temp.size):1){ if (temp.size[q] == "-"){ sub.size<-sub.size+1 } else { break } }
    
    #Saves final data
    samp.frame<-cbind(samp.frame, sppSize = length(temp.size)-sub.size)
    save.frame<-rbind(save.frame, samp.frame)
    all.frame<-cbind(Sample = names(trimmed[j]), all.frame)
    all.frame<-cbind(all.frame, sppSize = length(temp.size)-sub.size)
    save.all<-rbind(save.all, all.frame)
  }#end j loop
  
  #Moves on if there are no frames found. Saves to Anon folder? 
  if (unique(save.frame$Frame)[1] == "0"){
    print(paste(locus.names[i], " sucked. No frames found.", sep = ""))
    next
  }
  
  ##############
  #STEP 5: Uses previous data to find a consistent frame
  ##############  
  #Looks at the overall data rather than the best indiv data to find a consistent frame
  temp.all<-save.all
  frame.names<-unique(temp.all$Frame)
  #Goes through the equally good frames and reduces to frames with the same range
  very.best<-data.frame()
  for (k in 1:length(frame.names)){
    temp.best<-temp.all[temp.all$Frame == frame.names[k],]
    starts<-table(temp.best$FrameStart)[table(temp.best$FrameStart) == max(table(temp.best$FrameStart))]
    ends<-table(temp.best$FrameEnd)[table(temp.best$FrameEnd) == max(table(temp.best$FrameEnd))]
    
    #Removes duplicates
    starts<-starts[as.numeric(names(starts)) == min(as.numeric(names(starts)))]
    ends<-ends[as.numeric(names(ends)) == min(as.numeric(names(ends)))]
    
    super.best<-temp.best[temp.best$FrameStart == as.numeric(names(starts)),]
    super.best<-super.best[super.best$FrameEnd == as.numeric(names(ends)),]
    very.best<-rbind(very.best, super.best)
  }#end k loop
  
  #Moves on if there are no frames found. Saves to Anon folder? 
  if (nrow(very.best) == 0){
    print(paste(locus.names[i], " sucked. No frames found.", sep = ""))
    next
  }

  ##############
  #STEP 6: Selects the best frame
  ##############  
  #Picks out the best frame
  best.frame<-table(very.best$Frame)[table(very.best$Frame) == max(table(very.best$Frame))]
  
  #If there are multiple good frames pick the biggest
  if (length(best.frame) != 1){ 
    temp.fix<-very.best[very.best$Frame %in% names(best.frame),]
    bigger<-temp.fix[temp.fix$Size == max(temp.fix$Size),]
    best.frame<-table(bigger$Frame)[table(bigger$Frame) == max(table(bigger$Frame))]
  }#end if
  
  #If they are same size just pick from this order
  if (length(best.frame) != 1){ 
    #Picks the best from this order of things
    temp.stop<-best.frame[names(best.frame) == "F1"]
    if (length(temp.stop) == 0){ temp.stop<-best.frame[names(best.frame) == "R1"] }
    if (length(temp.stop) == 0){ temp.stop<-best.frame[names(best.frame) == "F2"] }
    if (length(temp.stop) == 0){ temp.stop<-best.frame[names(best.frame) == "R2"] }
    if (length(temp.stop) == 0){ temp.stop<-best.frame[names(best.frame) == "F3"] }
    if (length(temp.stop) == 0){ temp.stop<-best.frame[names(best.frame) == "R3"] }
    best.frame<-temp.stop
    }
    
  #Checks the remaining ORF size
  temp.size<-save.all[save.all$Frame == names(best.frame),]
  samp.spp<-temp.size[temp.size$sppSize == max(temp.size$sppSize),]
  samp.seq<-trimmed[names(trimmed) == samp.spp$Sample[1]]
  samp.size<-nchar(gsub("-", "", as.character(samp.seq)))
  
  if (mean(temp.size$Size) <= samp.size*.5){ print(paste(locus.names[i], " was small.", sep = "")) }
  
  if (mean(temp.size$Size) <= samp.size*.25){ 
    print(paste(locus.names[i], " sucked. Too few sequence left.", sep = ""))
    next
  }
  
  if (best.frame <= length(trimmed) * .5){ 
    print(paste(locus.names[i], " sucked. No cosistent frame.", sep = ""))
    next
  }
  
  #Reverses if it needs to
  if (length(grep("R", names(best.frame))) != 0){
    new.align<-reverseComplement(trimmed)
  }else { new.align<-trimmed }
  
  ##############
  #STEP 7: Gets start and stop coordinates for each sequence and find best alignment
  ##############  
  #Gets trimming locations
  frame.ranges<-save.all[save.all$Frame == names(best.frame),]
  
  #Gets potential starts and ends
  starts<-table(frame.ranges$FrameStart)[table(frame.ranges$FrameStart) == max(table(frame.ranges$FrameStart))]
  ends<-table(frame.ranges$FrameEnd)[table(frame.ranges$FrameEnd) == max(table(frame.ranges$FrameEnd))]
  starts<-starts[starts == max(starts)]
  ends<-ends[ends == max(ends)]
  
  if (length(starts) != 1 || length(ends) != 1){ 
    frame.ranges<-save.all[save.all$Frame == names(best.frame),]
    starts<-table(frame.ranges$FrameStart)[table(frame.ranges$FrameStart) == max(table(frame.ranges$FrameStart))]
    ends<-table(frame.ranges$FrameEnd)[table(frame.ranges$FrameEnd) == max(table(frame.ranges$FrameEnd))]
    starts<-starts[starts == max(starts)]
    ends<-ends[ends == max(ends)]
  }
  
  if (length(starts) != 1 || length(ends) != 1){ 
    starts<-starts[as.numeric(names(starts)) == min(as.numeric(names(starts)))]
    ends<-ends[as.numeric(names(ends)) == max(as.numeric(names(ends)))]
  }
  
  ###################
  #STEP 8: Makes sure entire alignment is a multiple of 3
  ###################
  anu.start<-as.numeric(names(starts))
  new.end<-as.numeric(names(ends))
  new.len<-new.end-(anu.start-1)
  
  #Gets a new end to keep in multiples of 3 for proteins
  if (length(new.len[which(new.len %%3==0)]) == 0) {
    anu.end<-new.end-1
  } else { anu.end<-new.end }
  
  new.len<-anu.end-(anu.start-1)
  if (length(new.len[which(new.len %%3==0)]) == 0) {
    anu.end<-new.end-2
  } else { anu.end<-anu.end }
  
  #Trims sequence with new coords
  done.seq<-subseq(start = anu.start, end = anu.end, x = new.align)

  ###################
  #STEP 9: Trim out odd start/end bases
  ###################
  codon.seq<-DNAStringSet()
  for (k in 1:length(done.seq)){
    ref.aligned<-as.character(done.seq[k])
    
    #Chcecks at beginning of sequence
    not.gaps<-str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
    
    if (length(not.gaps) <= min.len){ next }
    
    #Checks if its odd, delete 1 base
    if ( (not.gaps[1]-1) %%3 == 2){ substr(ref.aligned, not.gaps[1], not.gaps[1])<-"-" }
    
    #Deletes 2 bases its off by
    if ( (not.gaps[1]-1) %%3 == 1){
      substr(ref.aligned, not.gaps[1], not.gaps[1])<-"-"
      substr(ref.aligned, not.gaps[1]+1, not.gaps[1]+1)<-"-"
    }#end if
    
    #checks for end of sequence
    not.gaps<-str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
    #counts characters
    char.len<-(not.gaps[length(not.gaps)]-not.gaps[1])+1
    end.pos<-not.gaps[length(not.gaps)]
    #removes odd characters at ends
    if ( char.len %%3 == 1){ substr(ref.aligned, end.pos, end.pos)<-"-" }
    
    if ( char.len %%3 == 2){ 
      substr(ref.aligned, end.pos-1, end.pos-1)<-"-"
      substr(ref.aligned, end.pos, end.pos)<-"-"
    } #end if 
    
    #Saves final seqs
    save.seq<-DNAStringSet(ref.aligned)
    names(save.seq)<-names(done.seq[k])
    codon.seq<-append(codon.seq, save.seq)
    
  } #END K  
  
  ###################
  #STEP 10: Change stop codons to N
  ###################
  #Finds stop codons to replace
  n.seq<-DNAStringSet(gsub("-", "N", as.character(codon.seq)))
  stop.seq<-DNAStringSet()
  for (k in 1:length(n.seq)){
    stop.data<-find.orf(n.seq[k], codon = T, min.size = 80)
    stop.data<-stop.data[stop.data$Frame == "F1",]
    
    #Skips if there are more than 3 stop codons
    if (nrow(stop.data) >= 3){ next }
    
    if (stop.data$Start[1] == 0){ stop.seq<-append(stop.seq, n.seq[k]) } else {
      #Goes through each codon
      ref.aligned<-as.character(n.seq[k])
      for (y in 1:nrow(stop.data)){      
        #Saves final seqs
        substr(ref.aligned, stop.data$Start[y], stop.data$Start[y])<-"N"
        substr(ref.aligned, stop.data$Start[y]+1, stop.data$Start[y]+1)<-"N"
        substr(ref.aligned, stop.data$Start[y]+2, stop.data$Start[y]+2)<-"N"
      }#end Y LOOP
      
      #Saves final data
      save.seq<-DNAStringSet(ref.aligned)
      names(save.seq)<-names(n.seq[k])
      stop.seq<-append(stop.seq, save.seq)
    }#end if state
   }# END K loop
  
  if (length(stop.seq) <= min.taxa){ next }
  
  ###################
  #FINAL STEP: Save everything after some final spp and length filtering
  ###################
  #Removes sequences that are less than a certain coverage
  t.align<-strsplit(as.character(stop.seq), "")
  len.loci<-lapply(t.align, function (x) x[x != "N"])
  spp.len<-unlist(lapply(len.loci, function (x) length(x)))
  spp.rem<-spp.len[spp.len <= max(spp.len) * as.numeric(min.cov)]  
  spp.rem<-append(spp.rem, spp.len[spp.len <= min.len]  )
  
  if (length(spp.rem) > 0){ 
    red.align<-t.align[!names(t.align) %in% unique(names(spp.rem))]
  } else { red.align<-t.align }
  
  #writes alignment
  mat.align<-lapply(red.align, tolower)
  write.align<-as.matrix(as.DNAbin(mat.align))
  
  #readies for saving
  write.phy(write.align, file=paste0(work.dir, "/", trim.dir, "/", locus.names[i]), interleave = F)
  
  #Saves protein sequence version
  trans.prot<-translate(stop.seq, genetic.code=GENETIC_CODE, if.fuzzy.codon="solve")
  mat.align<-lapply(trans.prot, tolower)
  write.align<-as.matrix(as.AAbin(trans.prot))
  write.phy(write.align, file=paste0(work.dir, "/", prot.dir, "/", locus.names[i]), interleave = F)
  
  #CHECKS FOR STOP CODNS
  stop.align<-strsplit(as.character(trans.prot), "")
  stop.loci<-lapply(stop.align, function (x) x[x == "*"])
  stop.len<-unlist(lapply(stop.loci, function (x) length(x)))
  stop.pres<-stop.len[stop.len >= 1]  
  
  if (length(stop.pres) != 0){ stop("STOP CODON") }
  #IQTREE TEST
 # system(paste0("iqtree -s ", work.dir, "/", trim.dir, "/", locus.names[i], 
  #             " -nt ", threads, " -m MFP -st CODON -rcluster 10 -msub nuclear"))
  
  
}#end i loop
  

####################################################################
####################################################################
##### Step 3. Save UCE and bad exon loci separately        #########
####################################################################
####################################################################

#Blast two probe set files together to get uce loci to keep
setwd(work.dir)

#Make blast database for the probe loci
system(paste("makeblastdb -in ", probe.file, " -parse_seqids -dbtype nucl ",
             " -out probe_blast_db", sep = ""))

#Matches samples to loci
system(paste("blastn -task dc-megablast -db probe_blast_db",
             " -query ", uce.file, " -out uce_match.txt", 
             " -outfmt 6 -num_threads ", threads, sep = ""))

#headers for the blast db
headers<-c("qName", "tName", "pident", "matches", "misMatches", "gapopen", 
           "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore")

#Load in matches
match.data<-fread("uce_match.txt", sep = "\t", header = F, stringsAsFactors = FALSE)
setnames(match.data, headers)
match.data<-match.data[match.data$evalue <= 0.05,]
uce.names<-unique(match.data$tName)
uce.names<-paste(uce.names, ".phy", sep = "")
system("rm probe_blast* uce_match.txt")

#Sets up new directory for this stuff
trim.dir<-"uce_trimmed"
dir.create(paste(work.dir, "/", trim.dir, sep = ""))
setwd(paste(work.dir, "/all-markers_untrimmed", sep = ""))
file.names<-list.files(".")
locus.names  <- uce.names[uce.names %in% file.names]
#locus.names <- file.names

#Loops through each locus and does operations on them
for (i in 1:length(locus.names)){
  
  ##############
  #STEP 1: Basic steps
  ##############
  #Reads in files
  setwd(paste(work.dir, "/all-markers_untrimmed", sep = ""))
  align<-readAAMultipleAlignment(file = paste(work.dir, "/all-markers_untrimmed/", locus.names[i], sep =""), format = "phylip")
  align<-DNAStringSet(align)
  
  #Remove gap only alignments
  gap.align<-strsplit(as.character(align), "")
  gap.count<-unlist(lapply(gap.align, function(x) length(x[x != "-"]) ) )
  gap.rem<-gap.count[gap.count <= as.numeric(min.aln)]
  
  #Remove the taxa that need to go
  rem.taxa<-c(names(gap.rem), taxa.remove)
  rem.align<-align[!names(align) %in% rem.taxa]
  
  #removes loci with too few taxa
  if ( length(rem.align) <= as.numeric(min.taxa)){ next }
  
  #removes too short loci
  if (max(width(rem.align)) <= as.numeric(min.aln)){ next }
  
  ##############
  #STEP 2: Runs trimming programs
  ##############
  #A. Runs Trimal
  trimal.align<-run.trimal(rem.align)
  
  #B Slice up alignment
  red.align<-slice.trim(trimal.align, slice.size.bp = 100, threshold = 0.40)
  
  #C. Trims the alignment ends
  edge.align<-trim.ends(red.align, min.n.seq = ceiling(length(red.align) * edge.trim), codon.trim = F)
  #Gets consensus seq for trimming more
  con.seq<-make.consensus(edge.align, method = "majority")
  #Removes the edge gaps
  ok.seq<-gsub("\\+|-", "", as.character(con.seq))
  
  #removes loci with too few taxa
  if ( nchar(ok.seq) <= as.numeric(min.aln)){ next  }
  #removes too short loci
  if (max(width(edge.align)) <= as.numeric(min.aln)){ next }
  
  ##############
  #STEP 3: Cleanup and save
  ##############
  #removes loci with too few taxa
  if (length(names(edge.align)) <= as.numeric(min.taxa)){ 
    print(paste(locus.names[i], " deleted. Too few taxa after trimming.", sep = "") )
    next
  }
  #string splitting
  write.temp<-strsplit(as.character(edge.align), "")
  aligned.set<-as.matrix(as.DNAbin(write.temp) )
  
  #removes too short loci
  if (ncol(aligned.set) <= as.numeric(min.aln)){ 
    print(paste(locus.names[i], " deleted. Trimmed alignment length below threshold.", sep = "") )
    next 
  }
  
  #Removes samples that too short individually
  len.temp<-as.character(as.list(aligned.set))
  len.loci<-lapply(len.temp, function (x) x[x != "-"])
  spp.len<-unlist(lapply(len.loci, function (x) length(x)))
  spp.rem<-spp.len[spp.len <= max(spp.len) * as.numeric(min.cov)]  
  spp.rem<-append(spp.rem, spp.len[spp.len <= min.len]  )
  if (length(spp.rem) > 0){  aligned.set<-aligned.set[!rownames(aligned.set) %in% unique(names(spp.rem)),] }
  
  #removes loci with too few taxa
  if (length(rownames(aligned.set)) <= as.numeric(min.taxa)){ 
    print(paste(locus.names[i], " deleted. Too few taxa after trimming.", sep = ""))
    next
  }
  
  #readies for saving
  write.phy(aligned.set, file= paste0(work.dir, "/", trim.dir, "/", gsub(pattern = "\\..*", "", locus.names[i]), ".phy"), interleave = F)
  
}


###########################################################
###########################################################
## Step 4. Save the intron regions along separately     ###
###########################################################
###########################################################


### GET EXON_INTRON FILES NOT USED FOR EXON 

#Sets up new directory for this stuff
setwd(work.dir)
trim.dir<-"intron-only_trimmed"
dir.create(paste(work.dir, "/", trim.dir, sep = ""))
notrim.dir<-"intron-only_untrimmed"
dir.create(paste(work.dir, "/", notrim.dir, sep = ""))
setwd(paste(work.dir, "/exon-only_trimmed", sep = ""))
exon.names<-list.files(".")

#Checks nad removes some missing from others
setwd(paste(work.dir, "/all-markers_untrimmed", sep = ""))
all.names<-list.files(".")
locus.names<-all.names[all.names %in% exon.names]

#Checks if to resume or not
if (resume == "yes"){
  done<-list.files(paste0(work.dir, "/", trim.dir))
  locus.names<-locus.names[!locus.names %in% done]
}

#Loops through each locus and does operations on them
for (i in 1:length(locus.names)){
  
  ##############
  #STEP 1: Set up and match loci
  ##############
  #Reads in files
  setwd(paste(work.dir, "/exon-only_trimmed", sep = ""))
  align<-readAAMultipleAlignment(file = locus.names[i], format = "phylip")
  align<-DNAStringSet(align)
  
  #Gets consensus seq for trimming more
  con.seq<-make.consensus(align, method = "majority")
  #Removes the edge gaps
  ok.seq<-DNAStringSet(gsub("\\+|-", "", as.character(con.seq)))
  names(ok.seq)<-paste("Reference_Locus")
  
  setwd(paste(work.dir, "/all-markers_untrimmed", sep = ""))
  intron.align<-DNAStringSet(readAAMultipleAlignment(file = locus.names[i], format = "phylip"))

  ##############
  #STEP 2: Runs MAFFT to add
  ##############
  setwd(paste(work.dir, "/", trim.dir, sep = ""))
  alignment<-run.mafft(unaligned.contigs = intron.align, add.contigs = ok.seq, rev.dir = T,
                       algorithm = "add", save.name = gsub(".phy", "",locus.names[i]), delete.files = T)
  
  #Aligns and then reverses back to correction orientation
  reversed<-names(alignment)[grep(pattern = "_R_", names(alignment))]
  if (length(reversed[grep(pattern = "Reference_Locus", reversed)]) == 1){ alignment<-reverseComplement(alignment) }
  names(alignment)<-gsub(pattern = "_R_", replacement = "", x = names(alignment))
  
  #Gets the divergence to make sure not crazy
  diff<-pairwise.inf.sites(alignment, "Reference_Locus")
  bad.seqs<-names(diff)[which(diff >= 0.5)]
  rem.align<-alignment[!names(alignment) %in% bad.seqs]
  
  # Moves onto next loop in there are no good sequences
  if (length(rem.align) <= as.numeric(min.taxa)){ 
    print(paste(locus.names[i], " had too few taxa", sep = ""))
    next }
  
  #removes too short loci
  #if (max(width(rem.align)) <= as.numeric(min.aln)){ next }
  
  ##############
  #STEP 3: Removes exon from the intron part
  ##############
  #Removes the edge gaps
  ref.aligned<-as.character(alignment['Reference_Locus'])
  not.gaps<-str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
  ref.start<-min(not.gaps)
  ref.finish<-max(not.gaps)
  
  #Finds weird gaps to fix
  temp.gaps<-as.numeric(1)
  for (k in 1:length(not.gaps)-1){ temp.gaps<-append(temp.gaps, not.gaps[k+1]-not.gaps[k]) }
  temp.gaps<-temp.gaps-1
  names(temp.gaps)<-not.gaps
  bad.gaps<-which(temp.gaps >= 30)
  front.gaps<-bad.gaps[bad.gaps <= length(not.gaps) *0.10]
  end.gaps<-bad.gaps[bad.gaps >= length(not.gaps) *0.90]

  #Fix big gaps if there are any
  if (length(front.gaps) != 0){ 
    temp.change<-(max(as.numeric(names(front.gaps))-ref.start))-(max(front.gaps)-1)
    ref.start<-ref.start+temp.change
  }#end gap if
  
  #Fix big gaps if there are any
  if (length(end.gaps) != 0){ 
    add.bp<-length(temp.gaps)-min(end.gaps)
    #add.bp<-(ref.finish-min(as.numeric(names(end.gaps))))
    min.gaps<-temp.gaps[min(end.gaps)]
    temp.change<-as.numeric(names(min.gaps))-as.numeric(min.gaps)
    ref.finish<-temp.change+add.bp
  }#end gap if
  
  #Cuts out the intron pieces
  intron.left<-subseq(alignment, 1, ref.start-1)
  intron.right<-subseq(alignment, ref.finish+1, width(alignment))
  save.names<-names(alignment)
  
  #Merges the alignments
  intron.align<-DNAStringSet(paste0(as.character(intron.left), as.character(intron.right)))
  names(intron.align)<-save.names
  intron.align<-intron.align[names(intron.align) != "Reference_Locus"]
  
  #Remove gap only alignments
  gap.align<-strsplit(as.character(intron.align), "")
  gap.count<-unlist(lapply(gap.align, function(x) length(x[x != "-"]) ) )
  gap.rem<-gap.count[gap.count <= as.numeric(min.aln)]
  
  #Remove the taxa that need to go
  rem.taxa<-c(names(gap.rem), taxa.remove)
  rem.align<-intron.align[!names(intron.align) %in% rem.taxa]
  
  #removes loci with too few taxa
  if ( length(rem.align) <= as.numeric(min.taxa)){ next }
  
  #removes too short loci
  if (max(width(rem.align)) <= as.numeric(min.aln)){ next }
  
  #Saves a copy of not modified intron only
  intron.notrim<-strsplit(as.character(rem.align), "")
  save.intron<-as.matrix(as.DNAbin(intron.notrim) )
  write.phy(save.intron, file= paste0(work.dir, "/", notrim.dir, "/", locus.names[i]), interleave = F)
  
  ##############
  #STEP 4: Runs trimming programs
  ##############
  #A. Runs Trimal
  trimal.align<-run.trimal(rem.align)
  
  if (length(trimal.align) == 0){ next }
  
  #B Slice up alignment
  red.align<-slice.trim(trimal.align, slice.size.bp = 100, threshold = 0.5)
  
  if (length(red.align) == 0){ next }
  
  #C. Trims the alignment ends
  edge.align<-trim.ends(red.align, min.n.seq = ceiling(length(red.align) * edge.trim), codon.trim = F)
  #Gets consensus seq for trimming more
  con.seq<-make.consensus(edge.align, method = "majority")
  #Removes the edge gaps
  ok.seq<-gsub("\\+|-", "", as.character(con.seq))
  
  #removes loci with too few taxa
  if ( nchar(ok.seq) <= as.numeric(min.aln)){ next  }
  #removes too short loci
  if (max(width(edge.align)) <= as.numeric(min.aln)){ next }
  
  ##############
  #STEP 5: Cleanup and save
  ##############
  #removes loci with too few taxa
  if (length(names(edge.align)) <= as.numeric(min.taxa)){ 
    print(paste(locus.names[i], " deleted. Too few taxa after trimming.", sep = "") )
    next
  }
  #string splitting
  write.temp<-strsplit(as.character(edge.align), "")
  aligned.set<-as.matrix(as.DNAbin(write.temp) )
  
  #removes too short loci
  if (ncol(aligned.set) <= as.numeric(min.aln)){ 
    print(paste(locus.names[i], " deleted. Trimmed alignment length below threshold.", sep = "") )
    next 
  }
  
  #Removes samples that too short individually
  len.temp<-as.character(as.list(aligned.set))
  len.loci<-lapply(len.temp, function (x) x[x != "-"])
  spp.len<-unlist(lapply(len.loci, function (x) length(x)))
  spp.rem<-spp.len[spp.len <= max(spp.len) * as.numeric(min.cov)]  
  spp.rem<-append(spp.rem, spp.len[spp.len <= min.len]  )
  if (length(spp.rem) > 0){  aligned.set<-aligned.set[!rownames(aligned.set) %in% unique(names(spp.rem)),] }
  
  #removes loci with too few taxa
  if (length(rownames(aligned.set)) <= as.numeric(min.taxa)){ 
    print(paste(locus.names[i], " deleted. Too few taxa after trimming.", sep = ""))
    next
  }
  
  #readies for saving
  write.phy(aligned.set, file= paste0(work.dir, "/", trim.dir, "/", gsub(pattern = "\\..*", "", locus.names[i]), ".phy"), interleave = F)
  
}



################################################################
#### Step 5. Group linked loci and exons from same gene   ######
################################################################


#Sets up new directory for this stuff
trim.dir<-"locus-combined"
exon.dir<-"exon-only_trimmed"
dir.create(paste(work.dir, "/", trim.dir, sep = ""))
setwd(paste(work.dir, "/", exon.dir, sep = ""))

#Finds loci from the file names
file.names<-list.files(".")
exon.names<-file.names[grep("-ex", file.names)]
temp.names<-gsub(".*_", "", exon.names)
temp.names<-gsub("-.*", "", temp.names)
temp.names<-gsub(".phy", "", temp.names)
locus.names<-unique(temp.names[duplicated(temp.names)])

#Loops through each locus and does operations on them
for (i in 1:length(locus.names)){
  
  ##############
  #STEP 1: Get the taxa in the alignments for these loci
  ##############
  setwd(paste(work.dir, "/", exon.dir, sep = ""))
  locus.files<-exon.names[grep(paste0("_",locus.names[i],"-"), exon.names)]
  locus.files<-locus.files[order(gsub(".*-ex", "", locus.files))]
  
  if (length(grep("-ex2.phy", locus.files)) > 1){ stop("grep FAIL") }
  
  taxa.names<-c()
  for (j in 1:length(locus.files)){
    #Reads in files
    align<-readAAMultipleAlignment(file = paste(work.dir, "/", exon.dir, "/", locus.files[j], sep =""), format = "phylip")
    taxa.names<-append(taxa.names, rownames(align))
  }
  
  taxa.names<-unique(taxa.names)
  
  ##############
  #STEP 2: Get different loci
  ##############
  
  combined.align<-DNAStringSet()
  for (j in 1:length(locus.files)){
    #Reads in files
    align<-readAAMultipleAlignment(file = paste(work.dir, "/", exon.dir, "/", locus.files[j], sep =""), format = "phylip")
    align<-DNAStringSet(align)

    add.taxa<-taxa.names[!taxa.names %in% names(align)]
    
    blank.align<-DNAStringSet()
    if (length(add.taxa) != 0){
      for (y in 1:length(add.taxa)){
        blank.align<-append(blank.align, DNAStringSet(paste0(rep("-", max(width(align))), collapse = "")) )
      }
      names(blank.align)<-add.taxa
    }#end rem seqs if
    
    #Saves the slices and cats
    new.align<-append(align, blank.align)
    new.align<-new.align[order(names(new.align))]
    save.names<-names(new.align)
    combined.align<-DNAStringSet(paste0(as.character(combined.align), as.character(new.align)))
    names(combined.align)<-save.names
  }#end j loop
  
  ##############
  #STEP 3: Cleanup and save
  ##############

  #string splitting
  write.temp<-strsplit(as.character(combined.align), "")
  aligned.set<-as.matrix(as.DNAbin(write.temp) )
  
  #Removes samples that too short individually
  len.temp<-as.character(as.list(aligned.set))
  len.loci<-lapply(len.temp, function (x) x[x != "-"])
  spp.len<-unlist(lapply(len.loci, function (x) length(x)))
  spp.rem<-spp.len[spp.len <= max(spp.len) * as.numeric(min.cov)]  
  spp.rem<-append(spp.rem, spp.len[spp.len <= min.len]  )
  if (length(spp.rem) > 0){  aligned.set<-aligned.set[!rownames(aligned.set) %in% unique(names(spp.rem)),] }
  
  #removes loci with too few taxa
  if (length(rownames(aligned.set)) <= as.numeric(min.taxa)){ 
    print(paste(locus.names[i], " deleted. Too few taxa after trimming.", sep = ""))
    next
  }
  
  #readies for saving
  write.phy(aligned.set, file= paste0(work.dir, "/", trim.dir, "/", gsub(pattern = "\\..*", "", locus.names[i]), ".phy"), interleave = F)
  
}#end i loop
  
####################################################
####################################################
####### Step 6. Save Legacy loci separately     ####
####################################################
####################################################
# # 
# # 
# # #Blast two probe set files together to get uce loci to keep
# # setwd(work.dir)
# # 
# # #Make blast database for the probe loci
# # system(paste("makeblastdb -in ", probe.file, " -parse_seqids -dbtype nucl ",
# #              " -out probe_blast_db", sep = ""))
# # 
# # #Matches samples to loci
# # system(paste("blastn -task dc-megablast -db probe_blast_db",
# #              " -query ", legacy.file, " -out leg_match.txt", 
# #              " -outfmt 6 -num_threads ", threads, sep = ""))
# # 
# # #headers for the blast db
# # headers<-c("qName", "tName", "pident", "matches", "misMatches", "gapopen", 
# #            "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore")
# # 
# # #Load in matches
# # match.data<-fread("leg_match.txt", sep = "\t", header = F, stringsAsFactors = FALSE)
# # setnames(match.data, headers)
# # match.data<-match.data[match.data$evalue <= 0.05,]
# # 
# # #Gets rid of very poor matches
# # filt.data<-match.data[match.data$matches > 40,]
# # filt.data<-filt.data[filt.data$evalue <= 0.05,]
# # 
# # #Fixes direction and adds into data
# # filt.data[,qDir:= as.character("0")]
# # #Finds out if they are overlapping
# # for (k in 1:nrow(filt.data)){
# #   if (filt.data$tStart[k] > filt.data$tEnd[k]){
# #     filt.data$qDir[k]<-"-"
# #     new.start<-min(filt.data$tStart[k], filt.data$tEnd[k])
# #     new.end<-max(filt.data$tStart[k], filt.data$tEnd[k])
# #     filt.data$tStart[k]<-new.start
# #     filt.data$tEnd[k]<-new.end
# #   } else { filt.data$qDir[k]<-"+" }
# # }#end k loop
# # 
# # #Get the sizes from the contig name
# # probe.loci<-scanFa(FaFile(probe.file))
# # contigs<-scanFa(FaFile(legacy.file))
# # new.qsize<-width(contigs)[pmatch(filt.data$qName, names(contigs), duplicates.ok = T)]
# # filt.data[,qSize:=as.numeric(new.qsize)]
# # #Gets the sizes from the probes
# # new.tsize<-width(probe.loci)[pmatch(filt.data$tName, names(probe.loci), duplicates.ok = T)]
# # filt.data[,tSize:=as.numeric(new.tsize)]
# # #Removes matches that barely match to the full contig
# # filt.data<-filt.data[filt.data$matches >= filt.data$tSize *0.20,]
# # 
# # #Looks to remove duplicates
# # dup.match<-unique(filt.data[duplicated(filt.data$qName),]$qName)
# # 
# # keep.match<-c()
# # for (i in 1:length(dup.match)){
# #   
# #   sub.match<-filt.data[filt.data$qName %in% dup.match[i],]
# #   
# #   #Skips if they all match to each other
# #   if (length(unique(sub.match$qName)) == 1 && length(unique(sub.match$tName)) == 1){ next }
# #   
# #   #Looks for a pretty high bit score
# #   b.match<-sub.match[sub.match$bitscore/max(sub.match$bitscore) >= 0.6,]
# #   if (nrow(b.match) == 1){ 
# #     keep.match<-rbind(keep.match, b.match)
# #     next
# #   }
# #   
# #   c.match<-sub.match[sub.match$qSize/sub.match$tSize >= 0.6,]
# #   if (nrow(c.match) == 1){ 
# #     keep.match<-rbind(keep.match, c.match)
# #     next
# #   }
# # }
# # 
# # red.data<-filt.data[!filt.data$qName %in% dup.match,]
# # new.data<-rbind(red.data, keep.match)
# # 
# # #Gets the files to save separately
# # leg.names<-unique(new.data$tName)
# # leg.names<-paste(leg.names, ".phy", sep = "")
# # system("rm probe_blast* leg_match.txt")
# # 
# # #Sets up new directory for this stuff
# # trim.dir<-"legacy-markers_trimmed"
# # dir.create(paste0(work.dir, "/", trim.dir))
# # setwd(paste0(work.dir, "/exon-only_trimmed"))
# # file.names<-list.files(".")
# # locus.names<-leg.names[leg.names %in% file.names]
# # 
# # #Saves the new markers with new names 
# # for (i in 1:length(locus.names)){
# #   
# #   temp.data<-new.data[new.data$tName %in% gsub(".phy", "", locus.names[i]),]
# #   temp.data$qName<-gsub("_", "-", temp.data$qName)
# #   new.name<-paste0(temp.data$tName, "_", temp.data$qName)
# #   system(paste0("cp ", locus.names[i]," ", work.dir, "/", trim.dir, "/", new.name, ".phy"))
# # 
# # }
# # 
# # # END SCRIPT
# # 
