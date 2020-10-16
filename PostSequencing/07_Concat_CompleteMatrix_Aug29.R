#The two R packages needed
library(Biostrings)
library(data.table)

#This is always good to have
options(stringsAsFactors = FALSE)

#################
### Functions ###
#################
### charTest function checks that a sequence has at least one of each nucleotide types A,C,G,T
### ls = a list containing one or more character vectors
### this function is mostly useful after first using the uniqueLetters function

charTest <- function(ls){
		all(c("A","C","G","T") %in% ls)
}

###############################################################################
###############################################################################
######################           INFO                 #########################
###############################################################################
###############################################################################

# AMAS is the only outside program needed for this. Its python 3 based.

# https://github.com/marekborowiec/AMAS

# Just move the AMAS.py file to some place your path can access

###############################################################################
###############################################################################
######################    1.  Parameter setup         #########################
###############################################################################
###############################################################################

###### PARAMETER SETUP ####
#Note that you only need 1 thread for this, takes 15 min

# TO DO PARAMETER LIST
# min.per<-0.25 #min percent of taxa to keep an alignment instead of or alongside min number taxa (to do)
# min.pis<-0.05 #Filters out the bottom percentage of loci with the lowest parsimoney informative sites (to do)
# taxon.comp<-0.10 #Filters out taxa that are below a certain threshold. To do. 
# loci.sets<-"both" #can also choose exon_only or exon_intron (only both works now)
# gblock<-"no" #gblock entire alignment (not yet implemented, might be better in alignment script)
# Get missing data percentage per taxa on a bp by bp basis, maybe keep based on this per gene

#Percent completeness matrices:
#Means that for a locus to be kept, it must have X percent total samples in the alignment for that locus

per.comp    <- c("50", "70", "75", "90", "95") #Can include as many as you like, with whatever percentages (in whole numbers)
folder.keep <- "no" #saves or deletes folder of individual loci created for each data matrix generated

#Add taxa you would like to remove from alignment before concatenation
#Note that removing taxa can lead to gappy alignments, and should be realigned if its too many

taxa.remove<-c() ## need to include this even though its usually empty!

#working directory where the loci are located

work.dir <- "/Volumes/MyPassport/SequenceCapture/Results/Alignments_includeOphiophagus-Anchored-UCEs" #Folder with target_only and all-loci folders
#tar.dir  <- "target-only_untrimmed"
tar.dir  <- "all-markers_untrimmed"
out.dir  <- "Concatenated_Matrices" #Whatever you want to call it

min.taxa   <-4    #min number of taxa to keep an alignment
min.length <-100  #Filters out loci with less than this number of bp

#Sets working directory and creates output directory if it doesnt already exist
setwd(work.dir)
if (file.exists(out.dir) == F){
	dir.create(out.dir)
}

#Create log file with useful data
log.con<-file("concat_log.txt", open="a")
cat("Concatenation script run log", file = log.con)
cat("\n", file = log.con)
cat("\n", file = log.con)

###############################################################################
###############################################################################
######################      1.. Concat                 ########################
###############################################################################
###############################################################################

locus.names<-unique(list.files(path = paste0(tar.dir, "/.")))
bad.loci<-c()
taxa.count<-data.table(Locus = locus.names, Count = as.numeric(0))

for (i in 1:length(locus.names)){
	#Reads in files
	# align <- readAAMultipleAlignment(file = paste0(work.dir, "/", tar.dir, "/", locus.names[i]), format = "phylip")

	align  <- unmasked(readDNAMultipleAlignment(file = paste0(work.dir, "/", tar.dir, "/", locus.names[i]), format = "phylip"))
	
	#########################################
	# A. Excludes specific taxa and/or loci based on selected parameters
	#########################################
	#removes too short loci
	if (width(align)[1] < min.length){
		bad.loci <- append(bad.loci, locus.names[i])
	}
	
	#removes loci with too few taxa
	if ( length(tax.names) < min.taxa){
		bad.loci <- append(bad.loci, locus.names[i])
	}
	
	#Remove taxa that lack at least one each of A,C,G,and T bases (usually these are individuals with only gaps)
	uniqueChar  <- lapply(X=align,FUN=uniqueLetters)         ### gets the list of unique character states for each individual
	taxa.test   <- unlist(lapply(X=uniqueChar,FUN=charTest)) ### gets a list of logicals indicating if the individuals has A,C,G, and T nucleotides all present
	align       <- DNAMultipleAlignment(align[c(taxa.test)]) ### removes individuals that lack an A,G,C, or T

	#Use taxa remove
	tax.names <- rownames(align)
	tax.names <- tax.names[!tax.names %in% taxa.remove]
	
	
	
	
	
		
	#Add other removal parameters later

	####
	####
	
	#########################################
	# B. Loci number vs. Parsimony Informative Sites
	#########################################
	
	#ALL
	#aligned.meta<-fread(file = "Final_Datasets_2017/aligned_loci_stats_all.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
	
	#pis.data<-aligned.meta[order(aligned.meta$NI_PISper, decreasing = F),]
	#plot(pis.data$NI_PISper, rep(1:nrow(pis.data)))
	
	#########################################
	# C. Generates loci stats for percent completeness matrices 
	#########################################
	
	taxa.count[i,2]<-length(tax.names)
	
}
 
#########################################
# Create loci lists for each data matrix
######################################### 
 
fin.tax <- taxa.count[!taxa.count$Locus %in% bad.loci,]  
 
datasets<-vector("list")
max.tax<-max(fin.tax$Count)-length(taxa.remove)
for (k in 1:length(per.comp)){
	datasets[[k]]<-fin.tax[(fin.tax$Count/max.tax)*100 > per.comp[k]]$Locus
}
 
names(datasets)<-paste(per.comp, "matrix", sep="_")
 
#Adds some log output
cat("Exon only dataset loci counts for each matrix", file = log.con)
cat("\n", file = log.con)
matrix.count<-data.table(Matrix = names(unlist(lapply(datasets, length))), Loci = unlist(lapply(datasets, length)))

for (k in 1:nrow(matrix.count)){
	cat(paste(matrix.count[k,1]), file = log.con)
	cat("\n", file = log.con)
	cat(paste(matrix.count[k,2]), file = log.con)
	cat("\n", file = log.con)
}

#Run AMAS python concatenation for each of the data matrices

per.comp <- per.comp[which(unlist(lapply(datasets,FUN=length))>0)] ### updates the per.comp vector so that the loop only tries to create datasets that can actually exist

for (j in 1:length(per.comp)){
	#Creates folders of loci files
	setwd(paste0(work.dir, "/", out.dir))
	concat.files<-paste(tar.dir, per.comp[j], sep = "_")
	if(file.exists(concat.files)==F){
		dir.create(concat.files)
	}
	
	#copies files from prior folder into new folder
	start<-1
	end<-1000
	for (k in 1:ceiling(length(datasets[[j]])/1000) ) {
		data.files<-datasets[[j]][start:end]
		data.files<-data.files[is.na(data.files) !=T]
		system(paste("cp ", paste(work.dir, "/", tar.dir, "/", data.files, collapse = " ", sep ="")," ", work.dir, "/", out.dir, "/", concat.files, sep = ""))
		start<-start+1000
		end<-end+1000
	}

	#Runs the AMAS command to concatenate
	setwd(concat.files)

	#Runs differently if taxa are to be deleted
	#system(paste0("AMAS.py concat -f phylip -d dna -i *phy -u phylip --part-format raxml"))
	system(paste0("python /Applications/AMAS-master/amas/AMAS.py concat -f phylip -d dna -i *phy -u phylip --part-format raxml"))

	#Copies and pastes the concatenated files to output directory, deletes and renames
	system(paste0("cp concatenated.out partitions.txt ", work.dir, "/", out.dir))
	system(paste0("rm concatenated.out partitions.txt"))

	setwd(paste0(work.dir, "/", out.dir))
	system(paste0("mv concatenated.out ", concat.files, ".phy"))
	system(paste0("mv partitions.txt ", concat.files, "_parts.txt"))

	#Command to remove taxa from large alignment
	if (length(taxa.remove) != 0){
		system(paste0("python /Applications/AMAS-master/amas/AMAS.py remove -x ", paste(taxa.remove, collapse = " ")," -d dna -f phylip -i ", concat.files, ".phy", " -u phylip --part-format raxml"))
		#Delete previous file and keep new file
		system(paste0("rm ", concat.files, ".phy"))
		system(paste0("mv reduced_", concat.files, ".phy-out.phy ",concat.files, ".phy"))
	}

	#Deletes folder if desired
	if (folder.keep == "no"){
		system(paste0("rm -r ", work.dir, "/",out.dir, "/", concat.files))
	}
}

#moves log file
system(paste0("cp ../concat_log.txt concat_log.txt"))
system(paste0("rm ../concat_log.txt"))

### END SCRIPT
