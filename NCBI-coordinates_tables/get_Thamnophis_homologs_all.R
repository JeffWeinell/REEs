packages.to.load <- c("R.methodsS3","R.oo","assertthat","Rcpp","tibble","magrittr","lazyeval","DBI","BH","dplyr","R.utils","data.table","utils","BiocGenerics","bitops","S4Vectors","IRanges","RCurl","XVector","zlibbioc","GenomeInfoDb","GenomeInfoDbData","GenomicRanges","Biostrings","lambda.r","futile.options","snow","futile.logger","BiocParallel","Rsamtools","ape","rentrez","rMSA","stringr","stringi","biofiles")
invisible(lapply(packages.to.load, FUN=library, character.only = TRUE))

###### load functions
source("/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/SnakeCap_functions.R") ### loads in some necessary functions

### I blasted (blastn) all of the ddRAD loci from Weinell_TargetLoci_Snakes_Final_18April2019.fa to Thamnophis sirtalis (database = refseq_genomes); default search parameters.
### Then I download the hit tables as a single CSV file.
### Then, I used the following code:

############
## Step 1 ##
############

Thermophis.Thamnophis_hit.table           <- read.table(file="/Users/Jeff/Downloads/Thermophis_vs_Thamnophis_Alignment-HitTable_26March2020.csv",colClasses="character",sep=",")
colnames(Thermophis.Thamnophis_hit.table) <- c("TargetName","Thamnophis.accession","percent.identity","alignment.length","num.mismatches","num.gap.opens","query.start","query.end","subject.start","subject.end","e.value","bit.score")
Thermophis.Thamnophis.best.table          <- Thermophis.Thamnophis_hit.table[match(unique(Thermophis.Thamnophis_hit.table$TargetName),Thermophis.Thamnophis_hit.table$TargetName),]
targetTable                               <- read.table(file="/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/TargetTable1_25March2020.txt",header=T,colClasses="character")

### Removes hits with alignment length < 100
Thermophis.Thamnophis.best.table          <- Thermophis.Thamnophis.best.table[which(as.numeric(Thermophis.Thamnophis.best.table$alignment.length)>=100),]

### get lengths of target loci (query lengths)

target.loci   <- readDNAStringSet(filepath="/Users/Jeff/Documents/SnakeCap_Data/Weinell_TargetLoci_Snakes_Final_18April2019.fa")
matched.loci  <- target.loci[Thermophis.Thamnophis.best.table$TargetName]
query.lengths <- width(matched.loci)

query.five.prime.not.hit.length  <- (as.numeric(Thermophis.Thamnophis.best.table$query.start)-1)
query.three.prime.not.hit.length <- (query.lengths-as.numeric(Thermophis.Thamnophis.best.table$query.end))

negative.subject.sense <- which(as.numeric(Thermophis.Thamnophis.best.table$subject.start) > as.numeric(Thermophis.Thamnophis.best.table$subject.end))

subject.start <- as.numeric(Thermophis.Thamnophis.best.table$subject.start)
subject.end   <- as.numeric(Thermophis.Thamnophis.best.table$subject.end)
query.start   <- as.numeric(Thermophis.Thamnophis.best.table$query.start)
query.end     <- as.numeric(Thermophis.Thamnophis.best.table$query.end)

subject.start.new <- vector(length=length(subject.start))
subject.end.new   <- vector(length=length(subject.end))

subject.start.new.with.buffer <- vector(length=length(subject.start))
subject.end.new.with.buffer   <- vector(length=length(subject.end))
edge.buffer.length <- 100
### also adds additional buffer zone to ends
for(i in 1:length(subject.start)){
	if(subject.start[i] > subject.end[i]){
		subject.start.new[i] <- subject.start[i] + query.five.prime.not.hit.length[i]
		subject.end.new[i]   <- subject.end[i] - query.three.prime.not.hit.length[i]
		subject.start.new.with.buffer[i] <- subject.start.new[i] + edge.buffer.length
		subject.end.new.with.buffer[i]   <- subject.end.new[i] - edge.buffer.length
	}
	if(subject.start[i] < subject.end[i]){
		subject.start.new[i]             <- subject.start[i] - query.five.prime.not.hit.length[i]
		subject.end.new[i]               <- subject.end[i] + query.three.prime.not.hit.length[i]
		subject.start.new.with.buffer[i] <- subject.start.new[i] - edge.buffer.length
		subject.end.new.with.buffer[i]   <- subject.end.new[i] + edge.buffer.length
	}
}

subject.contig.sense <- rep(1,length(subject.start.new))
subject.contig.sense[negative.subject.sense] <- 2
retrieve.subject.seqs.table <- cbind("target.loci" = as.character(Thermophis.Thamnophis.best.table$TargetName),"Thamnophis.accession" = as.character(Thermophis.Thamnophis.best.table$Thamnophis.accession),subject.start.new,subject.end.new,subject.start.new.with.buffer,subject.end.new.with.buffer,subject.contig.sense)

## COMMENTED TO AVOID OVERWRITING ## write.table(x=retrieve.subject.seqs.table,file="/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/retrieve_ddRAD-loci_Thamnophis_putative-homologs-of-Thermophis_with-100bp-buffer.txt",quote=F,sep="\t",row.names=F)
### this should download the putative homologs plus 100bp buffer region on each end. File is a genbank flatfile.
## COMMENTED TO AVOID OVERWRITING ## get_ncbi_sequences(outfile="/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/ddRAD-loci_Thamnophis_putative-homologs-of-Thermophis_with-100bp-buffer_flatfile.gb",accessionList=retrieve.subject.seqs.table[,2],startList=retrieve.subject.seqs.table[,5],endList=retrieve.subject.seqs.table[,6],strandList=retrieve.subject.seqs.table[,7],rettype="gb")

#################
#### Step 2 #####
#################
### Using functions in the Biofiles package, read the genbank flatfile created in step 1, and make a table indicating whether or not each locus has an annotated gene, mRNA, and/or CDS feature.
### For loci with gene, mRNA, or CDS features, use MAFFT to align putative homolog +100bp buffer region to the target sequence
### Get genomic coordinates for the region of the homolog+buffer sequence that aligns to the target sequence (and save these in a table as text file)
### use get_ncbi_sequences function to download a genbank flatfile containing the aligned-homolog of the target sequence

### reading the gb flatfile created in Step 1
gbData   <- read.gb("/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/ddRAD-loci_Thamnophis_putative-homologs-of-Thermophis_with-100bp-buffer_flatfile.gb",progress=T)

### re-importing or redefining some objects that were created in Step 1, and are also needed here
retrieve.subject.seqs.table               <- read.table(file="/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/retrieve_ddRAD-loci_Thamnophis_putative-homologs-of-Thermophis_with-100bp-buffer.txt",header=T,sep="\t",colClasses="character")
Thermophis.Thamnophis_hit.table           <- read.table(file="/Users/Jeff/Downloads/Thermophis_vs_Thamnophis_Alignment-HitTable_26March2020.csv",colClasses="character",sep=",")
colnames(Thermophis.Thamnophis_hit.table) <- c("TargetName","Thamnophis.accession","percent.identity","alignment.length","num.mismatches","num.gap.opens","query.start","query.end","subject.start","subject.end","e.value","bit.score")
Thermophis.Thamnophis.best.table          <- Thermophis.Thamnophis_hit.table[match(unique(Thermophis.Thamnophis_hit.table$TargetName),Thermophis.Thamnophis_hit.table$TargetName),]
targetTable                               <- read.table(file="/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/TargetTable1_25March2020.txt",header=T,colClasses="character")
Thermophis.Thamnophis.best.table          <- Thermophis.Thamnophis.best.table[which(as.numeric(Thermophis.Thamnophis.best.table$alignment.length)>=100),]
target.loci                               <- readDNAStringSet(filepath="/Users/Jeff/Documents/SnakeCap_Data/Weinell_TargetLoci_Snakes_Final_18April2019.fa")
matched.loci                              <- target.loci[Thermophis.Thamnophis.best.table$TargetName]
query.lengths                             <- width(matched.loci)

### identifying sequences with gene, mRNA, and/or CDS regions
features.test.table           <- cbind(names(matched.loci),rep("-",length(matched.loci)),rep("-",length(matched.loci)),rep("-",length(matched.loci)))
colnames(features.test.table) <- c("target.locus","gene","mRNA","CDS")
for(i in 1:nrow(features.test.table)){
	gene.temp <- gbData[[i]]["gene"]
	mRNA.temp <- gbData[[i]]["mRNA"]
	cds.temp  <- gbData[[i]]["CDS"]
	if(length(gene.temp)==0){
		features.test.table[i,"gene"] = "no"
	} else {
		features.test.table[i,"gene"] = "yes"
	}
	if(length(mRNA.temp)==0){
		features.test.table[i,"mRNA"] = "no"
	} else {
		features.test.table[i,"mRNA"] = "yes"
	}
	if(length(cds.temp)==0){
		features.test.table[i,"CDS"] = "no"
	} else {
		features.test.table[i,"CDS"] = "yes"
	}
}

#####
accession.homolog     <- retrieve.subject.seqs.table[,2]
subject.start.homolog <- vector(length=nrow(features.test.table))
subject.end.homolog   <- vector(length=nrow(features.test.table))
for(i in 1:nrow(features.test.table)){
#for(i in 1:10){
	dna.temp              <- biofiles::getSequence(gbData[[i]])      ### change 1 to i
	target.and.subject    <- c(matched.loci[i],dna.temp)             ### change 1 to i
	#DNAStringSet(rMSA::mafft(target.and.subject,param="--localpair --maxiterate 1000 --quiet --lop 3 --lep 0.123 --lexp 0.3 --thread 6"))
	## MAFFT alignment. Relatively high gap open and gap extension penalties used, whereas no penalty for gap offsets (i.e., no penalty for a string of gaps at the ends of the sequences)
	## These parameters were used because we expect that gaps will be present at the end of the target sequence, because we included a buffer (i.e., putatively nonhomologous region) on each end of the putative homolog
	alignment.temp        <- DNAStringSet(rMSA::mafft(target.and.subject,param="--localpair --maxiterate 1000 --quiet --lop -3 --lep 2 --lexp -3 --thread 6"))
	names(alignment.temp) <- names(target.and.subject)

	### finds first and last columns of alignment in which the target sequence is not a gap, and then trims the alignment to that range
	first.char.target <- str_locate(string=alignment.temp[1],pattern= as.character(subseq(matched.loci[i],start=1,end=1)))
	last.char.target  <- str_locate_last(string=alignment.temp[1],pattern= as.character(subseq(matched.loci[i],start=width(matched.loci[i]))))
	alignment.temp.trimmed.to.target <- subseq(alignment.temp,start=first.char.target,end=last.char.target)

	### make vectors containing either nongap or gap positions for target and homolog sequences of trimmed-to-target alignments
	nongap.positions.target  <- str_locate_all(string=alignment.temp.trimmed.to.target[1],pattern="[A-Z]")[[1]][,1]
	gap.positions.target     <- str_locate_all(string=alignment.temp.trimmed.to.target[1],pattern="-")[[1]][,1]
	nongap.positions.homolog <- str_locate_all(string=alignment.temp.trimmed.to.target[2],pattern="[A-Z]")[[1]][,1]
	gap.positions.homolog    <- str_locate_all(string=alignment.temp.trimmed.to.target[2],pattern="-")[[1]][,1]

	## Make a 4x2 character matrix containing:
	## [Row1,Column1]: "target.locus.name target.column.numbers"    ### A character string containing the name of the locus and "target.column.numbers", separated by a space.
	## [Row2,Column1]: "target.locus.name aligned.target.sequence"  ### A character string containing the name of the locus and "aligned.target.sequence", separated by a space.
	## [Row3,Column1]: "target.locus.name homolog.column.numbers"   ### A character string containing the name of the locus and "homolog.column.numbers", separated by a space.
	## [Row4,Column1]: "target.locus.name aligned.homolog.sequence" ### A character string containing the name of the locus and "aligned.homolog.sequence", separated by a space.
	## [Row1,Column2]: space-delimited character string of nucleotide-containing column numbers for the target sequence. A gap "-" is used for columns without characters
	## [Row2,Column2]: space-delimited character string of the aligned nucleotide sequence for the target sequence (i.e., includes with gaps "-";  alignment.temp.trimmed.to.target[1] collapsed with sep=" ").
	## [Row3,Column2]: space-separated character string of character column numbers for the Thamnophis putative homolog sequence. A gap "-" is used for columns without characters
	## [Row4,Column2]: space-delimited character string of the aligned nucleotide sequence of the Thamnophis putative homolog sequence. (i.e., alignment.temp.trimmed.to.target[2] collapsed with sep=" ").
	#####
	# The above matrix is generated for each locus and Thamnophis homolog, and then appended using rbind to create a matrix of matrices that can be saved to a text file and parsed later
	# The matrix will be useful to transform the Thamnophis feature table annotations (in Step 3) onto the target locus (when plotting using "graph_target_and_features.R"), and later when partitioning data according to type (e.g., CDS vs. non CDS, codon positions, ect.).
	
	target.column.numbers                            <- vector(length=width(alignment.temp.trimmed.to.target[1]),mode="character")
	target.column.numbers[nongap.positions.target]   <- c(1:length(nongap.positions.target))
	target.column.numbers[gap.positions.target]      <- "-"	
	homolog.column.numbers                           <- vector(length=width(alignment.temp.trimmed.to.target)[2],mode="character")
	homolog.column.numbers[nongap.positions.homolog] <- c(1:length(nongap.positions.homolog))
	homolog.column.numbers[gap.positions.homolog]    <- "-"
	
	row1.mat.temp <- cbind(paste(names(matched.loci[i]),"target.column.numbers"),paste(target.column.numbers,collapse=" "))
	row2.mat.temp <- cbind(paste(names(matched.loci[i]),"aligned.target.sequence"),paste(strsplit(as.character(alignment.temp.trimmed.to.target[1]),split="")[[1]],collapse=" "))
	row3.mat.temp <- cbind(paste(names(matched.loci[i]),"homolog.column.numbers"),paste(homolog.column.numbers,collapse=" "))
	row4.mat.temp <- cbind(paste(names(matched.loci[i]),"aligned.homolog.sequence"),paste(strsplit(as.character(alignment.temp.trimmed.to.target[2]),split="")[[1]],collapse=" "))
	target.vs.homolog.character.range.matrix         <- rbind(row1.mat.temp,row2.mat.temp,row3.mat.temp,row4.mat.temp) ### table that will be written
	if(i==1){
		append.if <- F
	} else {
		append.if <- T
	}
	## COMMENTED TO AVOID OVERWRITING ## write.table(target.vs.homolog.character.range.matrix,file="/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/Thermophis.vs.Thamnophis.transform.range.matrix.txt",quote=T,row.names=F,col.names=F,append=append.if,sep="\t")
	
	####### now back to making the table of Thamnophis homolog accession information (contig name and contig range).
	five.prime.subject.to.trim  <- gsub(pattern="-",replacement="",x= subseq(alignment.temp[2],start=1,end=(first.char.target-1)[1]))
	three.prime.subject.to.trim <- gsub(pattern="-",replacement="",x= subseq(alignment.temp[2],start=(last.char.target+1)[1]))

	if(retrieve.subject.seqs.table[i,7]==1){
		subject.start.homolog[i] <- as.numeric(retrieve.subject.seqs.table[i,5])+(width(five.prime.subject.to.trim))
		subject.end.homolog[i]   <- as.numeric(retrieve.subject.seqs.table[i,6])-(width(three.prime.subject.to.trim))
		subject.name.temp        <- paste0(accession.homolog[i],":",subject.start.homolog[i],"-",subject.end.homolog[i])
	}
	if(retrieve.subject.seqs.table[i,7]==2){
		subject.start.homolog[i] <- as.numeric(retrieve.subject.seqs.table[i,5])-(width(five.prime.subject.to.trim))
		subject.end.homolog[i]   <- as.numeric(retrieve.subject.seqs.table[i,6])+(width(three.prime.subject.to.trim))
		subject.name.temp        <- paste0(accession.homolog[i],":complement(",subject.start.homolog[i],"-",subject.end.homolog[i],")")
	}

	names(alignment.temp.trimmed.to.target)[2] <- paste("Thamnophis_sirtalis",subject.name.temp,"putative homolog of Thermophis baileyi",names(alignment.temp.trimmed.to.target)[1])
	subject.sequence.to.write <- DNAStringSet(gsub(pattern="-",replacement="",alignment.temp.trimmed.to.target[2]))
	if(i==1){
		append.choice <- F
	} else {
		append.choice <- T
	}
 
	## COMMENTED TO AVOID OVERWRITING ## writeXStringSet(x=subject.sequence.to.write,filepath="/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/ddRAD-loci_Thamnophis_putative-homologs-of-Thermophis_v2.fasta",append=append.choice,format="fasta")
}

subject.coordinates.table.temp    <- cbind(accession.homolog,subject.start.homolog,subject.end.homolog,subject.contig.sense)
subject.coordinates.table.new     <- subject.coordinates.table.temp
subject.coordinates.table.new[,2] <- apply(X=subject.coordinates.table.temp[,2:3],MARGIN=1,FUN=min)
subject.coordinates.table.new[,3] <- apply(X=subject.coordinates.table.temp[,2:3],MARGIN=1,FUN=max)
subject.NCBI.coordinates          <- paste0(subject.coordinates.table.new[,1],":",subject.coordinates.table.new[,2],"-",subject.coordinates.table.new[,3])
subject.table                     <- cbind(features.test.table,subject.NCBI.coordinates,subject.coordinates.table.new)

ddrad.target.names                <- targetTable$TargetName_ArborSci[which(targetTable$Locus.Type=="ddRAD-like")]
loci.missing                      <- setdiff(ddrad.target.names,subject.table[,"target.locus"])                    ### target loci not found with blastn
missing.loci.subject.table        <- cbind(loci.missing,matrix(data="—",ncol=8,nrow=length(loci.missing)))
subject.table.semifinal           <- rbind(subject.table,missing.loci.subject.table)
subject.table.final               <- subject.table.semifinal[match(ddrad.target.names,subject.table.semifinal[,1]),]
## COMMENTED TO AVOID OVERWRITING ## write.table(x=subject.table.final,file="/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/ddRAD-loci_Thamnophis_putative-homologs-of-Thermophis_NCBI-coordinates_table.txt",quote=F,sep="\t",row.names=F)

### This should download the putative homologs. File is a genbank flatfile.
## COMMENTED TO AVOID OVERWRITING ## get_ncbi_sequences(outfile="/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/ddRAD-loci_Thamnophis_putative-homologs-of-Thermophis_flatfile_v2_29March2020.gb",accessionList=subject.table[,"accession.homolog"],startList=subject.table[,"subject.start.homolog"],endList=subject.table[,"subject.end.homolog"],strandList=subject.table[,"subject.contig.sense"],rettype="gb")

############
## Step 3 ##
############

### reading the gb flatfile
# gbData.noBuffer   <- read.gb("/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/ddRAD-loci_Thamnophis_putative-homologs-of-Thermophis_flatfile.gb")
# # gbData.noBuffer     <- biofiles::gbRecord("/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/ddRAD-loci_Thamnophis_putative-homologs-of-Thermophis_flatfile.gb",progress=T)
## ^^^ produces a warning message "In mclapply(seq_len(n), do_one, mc.preschedule = mc.preschedule,  :scheduled cores 4, 3 encountered errors in user code, all values of the jobs will be affected"
## vvvv Re-loaded the gb flatfile after removing the gb record of "WeinellEntry5594", which was causing a glitch
#gbData.noBuffer_v2  <- biofiles::gbRecord("/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/ddRAD-loci_Thamnophis_putative-homologs-of-Thermophis_flatfile.gb",progress=T)

#### vvvv Re-loaded yet again the gb flatfile after re-downloading gb data and removing the gb record of "WeinellEntry5594", which was still causing a glitch
## UNCOMMENT WHEN NEEDED ## gbData.noBuffer_v2  <- biofiles::gbRecord("/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/ddRAD-loci_Thamnophis_putative-homologs-of-Thermophis_flatfile_v2_29March2020.gb",progress=T)

### re-importing or redefining some objects that were created in Step 1, and are also needed here
retrieve.subject.seqs.table               <- read.table(file="/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/retrieve_ddRAD-loci_Thamnophis_putative-homologs-of-Thermophis_with-100bp-buffer.txt",header=T,sep="\t",colClasses="character")
Thermophis.Thamnophis_hit.table           <- read.table(file="/Users/Jeff/Downloads/Thermophis_vs_Thamnophis_Alignment-HitTable_26March2020.csv",colClasses="character",sep=",")
colnames(Thermophis.Thamnophis_hit.table) <- c("TargetName","Thamnophis.accession","percent.identity","alignment.length","num.mismatches","num.gap.opens","query.start","query.end","subject.start","subject.end","e.value","bit.score")
Thermophis.Thamnophis.best.table          <- Thermophis.Thamnophis_hit.table[match(unique(Thermophis.Thamnophis_hit.table$TargetName),Thermophis.Thamnophis_hit.table$TargetName),]
targetTable                               <- read.table(file="/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/TargetTable1_25March2020.txt",header=T,colClasses="character")
namesTable                                <- read.table(file="/Users/Jeff/Documents/SnakeCap_Data/NamesTable_All-Loci_3March2020.txt",header=T,colClasses="character")
Thermophis.Thamnophis.best.table          <- Thermophis.Thamnophis.best.table[which(as.numeric(Thermophis.Thamnophis.best.table$alignment.length)>=100),]
target.loci                               <- readDNAStringSet(filepath="/Users/Jeff/Documents/SnakeCap_Data/Weinell_TargetLoci_Snakes_Final_18April2019.fa")
matched.loci                              <- target.loci[Thermophis.Thamnophis.best.table$TargetName]
query.lengths                             <- width(matched.loci)

### re-importing and generating some tables because closed R
subject.table.final <- read.table("/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/ddRAD-loci_Thamnophis_putative-homologs-of-Thermophis_NCBI-coordinates_table.txt",header=T,colClasses="character")
subject.table       <- subject.table.final[match(names(matched.loci),subject.table.final$target.locus),]
#subject.table.temp <- subject.table[-184,] ### removed the row with "WeinellEntry5594" because this was removed from the gb file
subject.table.temp  <- subject.table[-grep("WeinellEntry5594",subject.table[,"target.locus"]),] ### removed the row with "WeinellEntry5594" because this was removed from the gb file

### Re-checking the features for the putative homologs of the target sequences (without the buffer regions).
matched.loci.names                     <- subject.table.temp[,"target.locus"]
features.test.table.noBuffer           <- cbind(matched.loci.names,rep("-",length(matched.loci.names)),rep("-",length(matched.loci.names)),rep("-",length(matched.loci.names)))
colnames(features.test.table.noBuffer) <- c("target.locus","gene","mRNA","CDS")

for(i in 1:nrow(features.test.table.noBuffer)){
	gene.temp <- gbData.noBuffer_v2[[i]]["gene"]
	mRNA.temp <- gbData.noBuffer_v2[[i]]["mRNA"]
	cds.temp  <- gbData.noBuffer_v2[[i]]["CDS"]
	if(length(gene.temp)==0){
		features.test.table.noBuffer[i,"gene"] = "no"
	} else {
		features.test.table.noBuffer[i,"gene"] = "yes"
	}	
	if(length(mRNA.temp)==0){
		features.test.table.noBuffer[i,"mRNA"] = "no"
	} else {
		features.test.table.noBuffer[i,"mRNA"] = "yes"
	}	
	if(length(cds.temp)==0){
		features.test.table.noBuffer[i,"CDS"] = "no"
	} else {
		features.test.table.noBuffer[i,"CDS"] = "yes"
	}	
}

subject.table.temp[,1:4]   <- features.test.table.noBuffer
gbData.noBuffer.with.genes <- gbData.noBuffer_v2[which(subject.table.temp[,"gene"]=="yes")]  ### gbData with only the loci with a gene feature
subject.table.with.genes   <- subject.table.temp[which(subject.table.temp[,"gene"]=="yes"),] ### subject table with only the loci with a gene feature

feature.location.table.with.genes           <- matrix(data="–",ncol=10,nrow=nrow(subject.table.with.genes))
colnames(feature.location.table.with.genes) <- c("Gene.ID","Gene.Name","Gene.Range","Gene.target.Sense","mRNA.Range","mRNA.target.Sense","CDS.Range","CDS.target.Sense","CDS_codon_start","target_codon_start")

translate <- Biostrings::translate ### needed to because multiple libraries have a translate function

for(i in 1:nrow(feature.location.table.with.genes)){
	source.features.temp   <- gbData.noBuffer.with.genes[[i]]["source"]
	gene.features.temp     <- gbData.noBuffer.with.genes[[i]]["gene"]
	mRNA.features.temp     <- gbData.noBuffer.with.genes[[i]]["mRNA"]
	CDS.features.temp      <- gbData.noBuffer.with.genes[[i]]["CDS"]
	
	#if(nrow(dbxref(gene.features.temp))>1){
	#	feature.location.table.with.genes[i,"Gene.ID"]           <- "multiple gene IDs; locus skipped; come back to this later"
	#	next
	#}
	
	#gene.ID.temp             <- as.character(dbxref(gene.features.temp))
	#gene.name.temp           <- as.character(geneID(gene.features.temp))	
	gene.ID.temp              <- paste(as.matrix(dbxref(gene.features.temp)),collapse=";")
	gene.name.temp            <- paste(as.matrix(geneID(gene.features.temp)),collapse=";")
	gene.location.temp        <- cbind(as.character(fuzzy(gene.features.temp)[,1]),as.character(start(gene.features.temp)),as.character(fuzzy(gene.features.temp)[,2]),as.character(end(gene.features.temp)))
	gene.location.temp[,1]    <- mgsub(x=c("TRUE","FALSE"),y=c("<",""),z=gene.location.temp[,1])
	gene.location.temp[,3]    <- mgsub(x=c("TRUE","FALSE"),y=c(">",""),z=gene.location.temp[,3])
	gene.location.temp.string <- paste(paste0(gene.location.temp[,1],gene.location.temp[,2],"..",gene.location.temp[,3],gene.location.temp[,4]),collapse=";")
	gene.sense                <- paste(strand(gene.features.temp),collapse=";")

	feature.location.table.with.genes[i,"Gene.ID"]           <- gene.ID.temp
	feature.location.table.with.genes[i,"Gene.Name"]         <- gene.name.temp
	feature.location.table.with.genes[i,"Gene.Range"]        <- paste(gene.location.temp.string,collapse=";")
	feature.location.table.with.genes[i,"Gene.target.Sense"] <- gene.sense

	if(length(mRNA.features.temp)!=0){
		mRNA.location.temp        <- cbind(as.character(fuzzy(mRNA.features.temp)[,1]),as.character(unlist(start(mRNA.features.temp))),as.character(fuzzy(mRNA.features.temp)[,2]),as.character(unlist(end(mRNA.features.temp))))
		mRNA.location.temp[,1]    <- mgsub(x=c("TRUE","FALSE"),y=c("<",""),z=mRNA.location.temp[,1])
		mRNA.location.temp[,3]    <- mgsub(x=c("TRUE","FALSE"),y=c(">",""),z=mRNA.location.temp[,3])
		mRNA.location.temp.string <- paste(paste0(mRNA.location.temp[,1],mRNA.location.temp[,2],"..",mRNA.location.temp[,3],mRNA.location.temp[,4]),collapse=";")

		mRNA.sense         <- strand(mRNA.features.temp)
		feature.location.table.with.genes[i,"mRNA.Range"]        <- mRNA.location.temp.string
		mRNA.sense.string  <- paste(as.character(unlist(mRNA.sense)),collapse=";")
		feature.location.table.with.genes[i,"mRNA.target.Sense"] <- mRNA.sense.string
	}
	if(length(CDS.features.temp)!=0){
		#full.translation.temp    <- translation(CDS.features.temp)
		if("translation" %in% names(qualif(CDS.features.temp))){
			full.translation.temp    <- translation(CDS.features.temp)
		}
		CDS.location.temp        <- cbind(as.character(fuzzy(CDS.features.temp)[,1]),as.character(unlist(start(CDS.features.temp))),as.character(fuzzy(CDS.features.temp)[,2]),as.character(unlist(end(CDS.features.temp))))
		CDS.location.temp[,1]    <- mgsub(x=c("TRUE","FALSE"),y=c("<",""),z=CDS.location.temp[,1])
		CDS.location.temp[,3]    <- mgsub(x=c("TRUE","FALSE"),y=c(">",""),z=CDS.location.temp[,3])
		CDS.location.temp.string <- paste(paste0(CDS.location.temp[,1],CDS.location.temp[,2],"..",CDS.location.temp[,3],CDS.location.temp[,4]),collapse=";")
		CDS.sense                <- unlist(strand(CDS.features.temp))
		CDS.codon.start          <- vector(length=length(CDS.sense),mode="character")
		target.codon.start       <- vector(length=length(CDS.sense),mode="character")
		for(j in 1:length(CDS.sense)){
			if(CDS.sense[j]==1){
				full.target_CDS.plus.sense  <- getSequence(source.features.temp)  ### full target region and same strand as the CDS strand
				#CDS.plus.sense             <- getSequence(CDS.features.temp)    ### only the CDS region
				CDS.plus.sense <- subseq(full.target_CDS.plus.sense,start=as.numeric(CDS.location.temp[j,2]),end=as.numeric(CDS.location.temp[j,4]))
			} else{
				full.target_CDS.plus.sense  <- reverseComplement(getSequence(source.features.temp))
				#CDS.plus.sense             <- reverseComplement(getSequence(CDS.features.temp))
				CDS.plus.sense              <- reverseComplement(subseq(getSequence(source.features.temp),start=as.numeric(CDS.location.temp[j,2]),end=as.numeric(CDS.location.temp[j,4])))
			}
			if(width(CDS.plus.sense)<15){
				CDS.codon.start[j]    <- "skipped because CDS < 15bp"
				target.codon.start[j] <- "skipped because CDS < 15bp"
				next
			}
			cds1 <- CDS.plus.sense
			cds2 <- subseq(CDS.plus.sense,start=2)
			cds3 <- subseq(CDS.plus.sense,start=3)
			cds.frames.list <- DNAStringSet(c(cds1,cds2,cds3))
	
			cds.translated      <- translate(subseq(cds.frames.list,start=rep(1,3),end=c(mround(width(cds1),base=3,direction="down"),mround(width(cds2),base=3,direction="down"),mround(width(cds3),base=3,direction="down"))),if.fuzzy="solve")
			AA.seq.frames.list  <- gsub(pattern="^\\*+|\\*+$",replacement="",cds.translated)     ### removes beginning or terminal string of stop codons

			internal.stop.check <- unlist(lapply(X=str_locate_all(string=AA.seq.frames.list,pattern="\\*"),FUN=nrow)) ### number of internal stop codons if translating for each possible frame

			if(any(internal.stop.check==0) & "translation" %in% names(qualif(CDS.features.temp))){
				possible.reading.frames <- which(internal.stop.check==0)
				targetAA.in.fullAA      <- matrix(nrow=3,ncol=length(full.translation.temp)) ## empty matrix to be filled
				## now checking which translated reading frames are in the NCBI translations
				for(k in 1:3){
					if(internal.stop.check[k]>0){
						targetAA.in.fullAA[k,] <- 0
					} else {
						search.pattern.temp         <- as.character(AA.seq.frames.list[k])
						targetAA.in.fullAA[k,]      <- unlist(lapply(X=str_locate_all(string=full.translation.temp,pattern=search.pattern.temp),FUN=nrow))
					}
				}
					if(any(targetAA.in.fullAA!=0)){
					possible.reading.frames2        <- rbind(internal.stop.check==0,apply(X=targetAA.in.fullAA,MARGIN=1,FUN=function(y){any(y>0)}))
					CDS.codon.start.temp            <- which(apply(X=possible.reading.frames2,MARGIN=2,FUN=all))
					CDS.codon.start[j]              <- paste(CDS.codon.start.temp,collapse=",")
					target.cds1                     <- full.target_CDS.plus.sense
					target.cds2                     <- subseq(full.target_CDS.plus.sense,start=2)
					target.cds3                     <- subseq(full.target_CDS.plus.sense,start=3)
					target.frames.list.cds.plus     <- DNAStringSet(c(target.cds1,target.cds2,target.cds3))
					full.target.translated          <- translate(subseq(target.frames.list.cds.plus,start=rep(1,3),end=c(mround(width(target.cds1),base=3,direction="down"),mround(width(target.cds2),base=3,direction="down"),mround(width(target.cds3),base=3,direction="down"))),if.fuzzy="solve")
					target.codon.start.temp         <- vector(length=length(CDS.codon.start.temp),mode="character")
					for(z in 1:length(CDS.codon.start.temp)){
						target.codon.start.temp[z]  <- which(unlist(lapply(X=str_locate_all(string=full.target.translated,pattern=as.character(AA.seq.frames.list[CDS.codon.start.temp[z]])),FUN=nrow)) !=0)
					}
					target.codon.start[j]           <- paste(target.codon.start.temp,collapse=",")				
				} else {
					CDS.codon.start[j]    <- "—"
					target.codon.start[j] <- "–"
				}
			} else {
				CDS.codon.start[j]    <- "—"
				target.codon.start[j] <- "–"
				next
			}
		}
		feature.location.table.with.genes[i,"CDS.Range"]          <- CDS.location.temp.string
		feature.location.table.with.genes[i,"CDS.target.Sense"]   <- paste(CDS.sense,collapse=";")
		feature.location.table.with.genes[i,"CDS_codon_start"]    <- paste(CDS.codon.start,collapse=";")
		feature.location.table.with.genes[i,"target_codon_start"] <- paste(target.codon.start,collapse=";")
	}
}

subject.table.new                    <- cbind(subject.table.with.genes,feature.location.table.with.genes)
ddrad.target.names                   <- targetTable$TargetName_ArborSci[which(targetTable$Locus.Type=="ddRAD-like")]
loci.missing                         <- setdiff(ddrad.target.names,subject.table.new[,"target.locus"])                    ### target loci not found with blastn
missing.loci.subject.table           <- cbind(loci.missing,matrix(data="—",ncol=(ncol(subject.table.new)-1),nrow=length(loci.missing)))
colnames(missing.loci.subject.table) <- colnames(subject.table.new)
subject.table.semifinal.new          <- rbind(subject.table.new,missing.loci.subject.table)
subject.table.final                  <- subject.table.semifinal.new[match(ddrad.target.names,subject.table.semifinal.new[,1]),]

## COMMENTED TO AVOID OVERWRITING ## write.table(x=subject.table.final,file="/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/ddRAD-loci_Thamnophis_putative-homologs-of-Thermophis_NCBI-coordinates_table_MoreInfo_v2_30March2020.txt",quote=F,sep="\t",row.names=F)
## COMMENTED TO AVOID OVERWRITING ## write.table(x=subject.table.final,file="/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/ddRAD-loci_Thamnophis_putative-homologs-of-Thermophis_NCBI-coordinates_table_MoreInfo_v3_30March2020.txt",quote=F,sep="\t",row.names=F)

