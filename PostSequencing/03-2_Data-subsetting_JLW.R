library(ape)
library(seqinr)
library(stringr)
library(GenomicRanges)
library(Biostrings)
library(Rsamtools)
library(rMSA)
library(data.table)
library(DECIPHER)


#########
## This script does the following:
## –- Extract UCE or anchored phylogenomics data (or any non-SNP-based next gen data) that might overlap with SnakeCap data
## ----- Extract individuals from alignment file and save to its own file with gaps removed
## ----- BLASTn is then conducted with target loci as the "query", and the extracted data from previous step as the "subject" (this step performed independently of this R script, typically via terminal)
## ----- BLASTn Hit table is used to extract the loci that match SnakeCap for the individual of interest, and these data are appended to a file 
#########


input.file       <-  "/Users/Jeff/Google Drive/Herpetology/Pyron2014_Colubroidea_333loci_DNA-Data_SequentialFormat.phy"
alignment        <-  read.dna(file=input.file,format="sequential")
partition.matrix <-  fread("/Users/Jeff/Google Drive/Herpetology/Pyron2014_333Loci_PartitionMatrix.txt")
outfile          <-  "/Users/Jeff/Documents/SnakeCap_Data/Pyron2014_Anchored-subset/Storeria_PyronData.fa"
taxa             <-  attributes(alignment)$dimnames[[1]]

for(i in 1:nrow(partition.matrix)){
	alignment.temp <- alignment[3,partition.matrix$start[i]:partition.matrix$end[i]]
	#alignment.temp <- as.list(as.character(alignment.temp))
	temp.names <- paste("Locus",i,"-",taxa[18],sep="")
	attributes(alignment.temp)$dimnames[[1]] <- temp.names
	if(i==1){
		append=F
	} else{
		append=T
	}
	#write.fasta(sequences=alignment.temp,file.out="/Users/Jeff/Documents/SnakeCap_Data/PyronData.fa",names<-temp.names,nbchar=1000000,as.string=F,open=open)
	write.dna(x=alignment.temp,file=outfile,format="fasta",append=append,nbcol=1,colw=1000000000)
}

#######
## After running Blast (between Storeria and Thamnophis), now extracting the matches for the species of interest
#######

input.file <- "/Users/Jeff/Google Drive/Herpetology/Pyron2014_Colubroidea_333loci_DNA-Data_SequentialFormat.phy"
alignment <- read.dna(file=input.file,format="sequential")
partition.matrix <- fread("/Users/Jeff/Google Drive/Herpetology/Pyron2014_333Loci_PartitionMatrix.txt") 
outfile <- "/Users/Jeff/Documents/SnakeCap_Data/Pyron2014_Anchored-subset/PyronData_AnchoredLoci_Keep.fa"
taxa <- attributes(alignment)$dimnames[[1]]


keep.loci    <- c(4,11,122,146,173,176,218,300)
new.names    <- c("WeinellEntry2636","WeinellEntry852","WeinellEntry1726","WeinellEntry5316","WeinellEntry1399","WeinellEntry2815","WeinellEntry1132","WeinellEntry1431")
keep.species <- c(3,6,7,8,11,12,15,16,18,19,24,28,32) ##

#Acrochorus:7, Prosymna:16, Xenodermus:19, Gonionotophis:28, Pseudaspis:8, Atractaspis:11, Calliophis:32, Lapemis:3, Aparallactus:6, Pareas:12, Psammophis:15, Hormonotus:24, Storeria:18

for(i in 1:length(keep.loci)){
	#locus.temp <- keep.loci
	start.temp <- partition.matrix$start[keep.loci[i]]
	end.temp <- partition.matrix$end[keep.loci[i]]
	alignment.temp <- alignment[keep.species,start.temp:end.temp]
	#alignment.temp <- as.list(as.character(alignment.temp))
	temp.names <- paste(new.names[i],"-","AnchoredLocus",keep.loci[i],"_",taxa[keep.species],sep="")
	attributes(alignment.temp)$dimnames[[1]] <- temp.names
	if(i==1){
		append=F
	} else{
		append=T
	}
	#write.fasta(sequences=alignment.temp,file.out="/Users/Jeff/Documents/SnakeCap_Data/PyronData.fa",names<-temp.names,nbchar=1000000,as.string=F,open=open)
	write.dna(x=alignment.temp,file=outfile,format="fasta",append=append,nbcol=1,colw=1000000000)
}


##########

#/Users/Jeff/Documents/SnakeCap_Data/Streicher2016_UCE-subset/Lampropeltis_Streicher2016_allData.fa

input.file <- "/Volumes/MyPassport/SequenceCapture/Results/Streicher2016_data_snakes_only_RAxML.phy"
alignment  <- RemoveGaps(unmasked(readDNAMultipleAlignment(file = input.file, format = "phylip")))

#partition.matrix <- fread("/Users/Jeff/Documents/SnakeCap_Data/BlastQueries/Lampropeltis-getula_UCEs_vs_target-loci.txt") 

out.dir <- "/Users/Jeff/Documents/SnakeCap_Data/Streicher2016_UCE-subset/"
outfiles <- paste(out.dir,names(alignment),"_Streicher2016_allData.fa",sep="") #"/Users/Jeff/Documents/SnakeCap_Data/Streicher2016_UCE-subset/Streicher2016_SnakeCap1_subset.fa"

#taxa <- attributes(alignment)$dimnames[[1]]

#keep.loci    <- c(1:nrow(partition.matrix))
#new.names <- as.character(unlist(partition.matrix[,2]))
keep.species <- c(1,2,18,19,25,29) ## corresponds to the order of the taxa in the alignment
outfiles <- outfiles[keep.species]

for(i in 1:length(keep.species)){
	alignment.temp <- alignment[keep.species[i]]	
	temp.names <- names(alignment.temp)
	writeXStringSet(x=alignment.temp, filepath=outfiles[i], append=F, format="fasta",width=20000)
}


##################
#######
# Subsetting the Streicher data by locus
# After running Blast (between species UCE data and Thamnophis SnakeCap loci), now extracting the matches for the species of interest

#input.file <- "/Users/Jeff/Documents/SnakeCap_Data/Streicher2016_UCE-subset/Lampropeltis_Streicher2016_allData.fa"
input.file <- outfiles[6]
alignment  <- RemoveGaps(unmasked(readDNAMultipleAlignment(file = input.file, format = "fasta")))   ### may switch to this soon

partition.matrix <- fread(paste("/Users/Jeff/Documents/SnakeCap_Data/BlastQueries/",names(alignment),"_UCEs_vs_target-loci.txt",sep=""))
outfile <- "/Users/Jeff/Documents/SnakeCap_Data/Streicher2016_UCE-subset/Streicher2016_SnakeCap1_subset.fa"


keep.loci    <- c(1:nrow(partition.matrix))
new.names <- as.character(unlist(partition.matrix[,2]))

for(i in 1:length(keep.loci)){
	start.temp       <- as.numeric(unlist(partition.matrix[keep.loci[i],7]))
	end.temp         <- as.numeric(unlist(partition.matrix[keep.loci[i],8]))
	#alignment.temp  <- alignment[keep.species,start.temp:end.temp]
	alignment.temp <- DNAMultipleAlignment(x=alignment,start=start.temp,end=end.temp)  ### may switch to this soon
	
	temp.names <- paste(new.names[i],"-", rownames(alignment.temp), sep="")
	#attributes(alignment.temp)$dimnames[[1]] <- temp.names
	rownames(alignment.temp) <- temp.names
	#if(i==1){
	#	append=F
	#} else{
	#	append=T
	#}
	writeXStringSet(x=unmasked(alignment.temp), filepath=outfile, append=T, format="fasta",width=20000)
}