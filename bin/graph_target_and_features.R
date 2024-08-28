#### load necessary functions and libraries
source("~/SnakeCap_functions.R")
packages.to.load <- c("R.methodsS3","R.oo","assertthat","Rcpp","tibble","magrittr","lazyeval","DBI","BH","dplyr","R.utils","data.table","utils","BiocGenerics","bitops","S4Vectors","IRanges","RCurl","XVector","zlibbioc","GenomeInfoDb","GenomeInfoDbData","GenomicRanges","Biostrings","lambda.r","futile.options","snow","futile.logger","BiocParallel","Rsamtools","ape","rentrez","rMSA","stringr","stringi","biofiles","tidyr")
invisible(lapply(packages.to.load, FUN=library, character.only = TRUE))

####
## loads fasta-formatted DNA sequences of target.loci and probes used to capture target loci
####
target.loci             <- readDNAStringSet(file="Weinell_TargetLoci_Snakes_Final_18April2019.fa")
probes                  <- readDNAStringSet(file="Weinell_FinalProbeSet_20020Probes_7-Oct-2018.fasta")

####
## loads start/stop/sense coordinates (location qualifiers) of features (gene, mRNA, and CDS regions) for each target locus. A different table is used for each type of target locus.
####

Exon.Immune.annotations <- read.table("https://github.com/JeffWeinell/SnakeCap/raw/main/NCBI-coordinates_tables/Exon.and.Immune-Loci_Thamnophis_NCBI-coordinates_table_MoreInfo_12April2020_v3.txt",colClasses="character",header=T,sep="\t")                           ### table made from NCBI annotations using code in /Users/Jeff/Documents/SnakeCap_Data/Notes_12April2020.txt
UCE.annotations         <- read.table("https://github.com/JeffWeinell/SnakeCap/raw/main/NCBI-coordinates_tables/UCE-Loci_Thamnophis_NCBI-coordinates_table_MoreInfo_12April2020_v2.txt",colClasses="character",header=T,sep="\t")                                       ### table made from NCBI annotations using code in /Users/Jeff/Documents/SnakeCap_Data/Notes_12April2020.txt
ddRAD.annotations       <- read.table("https://github.com/JeffWeinell/SnakeCap/raw/main/NCBI-coordinates_tables/ddRAD-loci_Thamnophis_putative-homologs-of-Thermophis_NCBI-coordinates_table_MoreInfo_v3_30March2020.txt",colClasses="character",header=T,sep="\t")     ### table made from NCBI annotations using code in /Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/get_Thamnophis_homologs_all.R
vision.annotations      <- read.table("https://github.com/JeffWeinell/SnakeCap/raw/main/NCBI-coordinates_tables/Vision-Loci_Thamnophis_NCBI-coordinates_table_MoreInfo_12April2020_v2.txt",colClasses="character",header=T,sep="\t")                                    ### table made from NCBI annotations using code in /Users/Jeff/Documents/SnakeCap_Data/Notes_12April2020.txt
scalation.annotations   <- read.table("https://github.com/JeffWeinell/SnakeCap/raw/main/NCBI-coordinates_tables/Scalation-Loci_NCBI-coordinates_table_MoreInfo_12April2020.txt",colClasses="character",header=T,sep="\t")                              ### table made from annotations in paper rather than from NCBI feature tables

####
## merge the annotation tables into a single annotation table
####
features.table          <- rbind(Exon.Immune.annotations,UCE.annotations,scalation.annotations,vision.annotations,ddRAD.annotations)

# Write/read the combined feature table.
# write.table(x=features.table,file="All-Loci_NCBI-coordinates_table_MoreInfo_13April2020.txt",quote=F,sep="\t",row.names=F)
# features.table        <- read.table("https://github.com/JeffWeinell/SnakeCap/raw/main/NCBI-coordinates_tables/All-Loci_NCBI-coordinates_table_MoreInfo_13April2020.txt",colClasses="character",header=T,sep="\t")

### This next table is used to transform gene/mRNA/CDS feature ranges that are annotated for Thamnophis sirtalis onto the homologous region for loci not designed from T. sirtalis.
### So far, this table only includes range transormation info for Thamnophis vs. Thermophis (ddRAD-like loci). The table was created with the script get_Thamnophis_homologs_all.R
transform.ranges.table <- read.table("https://github.com/JeffWeinell/SnakeCap/raw/main/NCBI-coordinates_tables/Thermophis.vs.Thamnophis.transform.range.matrix.txt",colClasses="character",header=F,sep="\t")

loci.names  <- names(target.loci) ### Names of target loci. These have the form "WeinellEntryXXX", where XXX is a unique number identifier (one to 4 digits)
probe.names <- names(probes)      ### Names of probes. These have the form "WeinellEntryXXX_StartPosition_EndPosition" (first base of target locus = position "0")

loci.per.plot         <- 20               ### Number loci to include on each page of the pdf that will be generated. 20 seems to look nice.
rev.loci.names        <- rev(loci.names)  ### Reverse orders the list of loci.names. This is done because loci will be plotted from bottom to top, so we want to plot in reverse order so that numerical order is plotted from top to bottom of each pdf page
rev.order.target.loci <- rev(target.loci) ### Reverse orders the DNAStringSet of target loci (doesnt reverse each sequence, just the order of the set of sequences)

locus.of.probe   <- gsub("_.+","",probe.names)
origin.of.probe  <- as.numeric(gsub(".+_","",probe.names))                                 ### extracts the first position -1 of each probe relative to its target locus
probe.start      <- origin.of.probe+1                                                      ### start position of each probe relative to its target locus
probe.end        <- origin.of.probe+120                                                    ### end position of each probe relative to its target locus
n.full.plots     <- mround(length(loci.names)/loci.per.plot,base=1,direction="down")       ### number of plots with number of loci plotted equal to loci.per.plot
n.loci.last.plot <- ((length(loci.names)/loci.per.plot)-n.full.plots)*loci.per.plot        ### number of loci on the last plot (if non-zero)
k.i.mat          <- matrix(data=1:(n.full.plots*loci.per.plot),ncol=loci.per.plot,byrow=T) ### matrix of "i" values. Number of columns determines loci/plot. Number of rows determines number of plots.

#####
### number of plots (calculated from value of loci.per.plot, or pre-specified to only plot as subset of loci)
#####
if(!n.loci.last.plot==0){
	num.plots     <- nrow(k.i.mat)+1
	i.values.last <- (max(k.i.mat)+1):length(loci.names)
} else {
	num.plots    <- nrow(k.i.mat)
}
## num.plots <- 1

########
my.plots         <- list(); length(my.plots)= num.plots                                    ### empty list that will hold plot
#for(k in 1:num.plots){ ### test for "WeinellEntry5717" (k=156, j=11) because lots of features. See if some can be dropped.
for(k in 130:135){ ### plots scalation loci
	#### Define dimensions of current plot by calculating number of loci per plot and maximum length of loci on current plot
	if(k>nrow(k.i.mat)){
		x.axis.temp    <- seq(from=0,to=max(width(target.loci[i.values.last])),length=10)*1.2
	} else {
		x.axis.temp    <- seq(from=0,to=max(width(target.loci[k.i.mat[k,]])),length=10)*1.2
	}
	x.axis.temp[1] <- x.axis.temp[3]*(-1)
	y.axis.temp    <- seq(from=0,to=(loci.per.plot*10),length=10)
	par(mai=c(0.5,0,0.5,0),omi=c(0.5,0,0,0.25),xpd=TRUE)
	plot(x=x.axis.temp,y=y.axis.temp,col="white",axes=F,xlab="",ylab="")                          ### blank plotting area to plot on
	if(k>nrow(k.i.mat)){
		x.max    <- mround(x=ceiling(max(width(target.loci[i.values.last]))),base=100,direction="up")   ### longest locus among the kth set of loci, rounded up to the nearest 100  
	} else {
		x.max    <- mround(x=ceiling(max(width(target.loci[k.i.mat[k,]]))),base=100,direction="up")     ### longest locus among the kth set of loci, rounded up to the nearest 100  
	}
	x.median <- median(c(0,x.max))
	
	#### Draws x axis, tick marks, and x-axis labels
	segments(x0=0,y0=3,x1=x.max,lwd=1,lend=2)                                                     ### Draws the x-axis
	x.small.ticks  <- seq(from=0,to=(x.max-100),by=100)                                           ### Calculates positions on x-axis where small tick marks will be drawn (at frequent intervals) 
	x.medium.ticks <- seq(from=0,to=x.max,by=1000)                                                ### Calculates positions on x-axis where medium tick marks will be drawn (at intermediate intervals)
	x.big.ticks    <- c(0,x.max)                                                                  ### Calculates positions on x-axis where large tick marks will be drawn (at beginning and end of x-axis)
	segments(x0=x.small.ticks,x1=x.small.ticks,y0=rep(3,length(x.small.ticks)),y1=rep(3.5,length(x.small.ticks)))       ### Draws the small x-axis tick marks
	segments(x0=x.medium.ticks,x1=x.medium.ticks,y0=rep(3,length(x.medium.ticks)),y1=rep(4,length(x.medium.ticks)))     ### Draws the medium x-axis tick marks
	segments(x0=x.big.ticks,x1=x.big.ticks,y0=rep(3,length(x.big.ticks)),y1=rep(4.5,length(x.big.ticks)))               ### Draws the large x-axis tick marks
	text(x=0,y=1,"0",cex=0.5)                                                                     ### Draws min x-axis label
	text(x=x.max,y=1,x.max,cex=0.5)                                                               ### Draws max x-axis label
	text(x=x.medium.ticks,y=rep(1,length(x.medium.ticks)),x.medium.ticks,cex=0.5)                 ### Draws 1kb x-axis labels (at intermediate tick marks)
	text(x=x.median,y=-3,"sequence length [nt]",cex=0.7)                                          ### Draws main x-axis label
	x.text           <- median(c(x.axis.temp[1],0))                                               ### Calculates horizontal position where loci names will be written
	
	#### Calculates which loci will be plotted on the current plot (and in what order)
	if(k>nrow(k.i.mat)){
		i.values         <- rev(i.values.last)
	} else {
		i.values         <- rev(k.i.mat[k,])
	}
	
	#### Draws a legend for kth plot (below the x-axis)
	legend.mat.colnames           <- c("start.plus","end.plus","sense","feature.type","y.coords","segment.colors","lwd","arrow.direction","arrow.color","arrow.size","arrow.transparency","cap.transparency","line.type")
	feature.types                 <- c("target.sequence","probes","probe.coverage","gene.region","mRNA.region","CDS.region")
	legend.mat                    <- matrix(data="-",ncol=length(legend.mat.colnames),nrow=length(feature.types))
	colnames(legend.mat)          <- legend.mat.colnames
	legend.mat[,"feature.type"]   <- feature.types
	legend.mat[,"y.coords"]       <- "-15"
	legend.mat[,"segment.colors"] <- c("black","blue","darkgreen","green","purple","green")
	legend.mat[,"lwd"]            <- c(1.25,1,1,1.25,2.5,6)
	legend.mat[,"arrow.size"]     <- c(0,0,0,0.25,0.4,0.4)
	legend.mat[,"line.type"]      <- c(1,1,1,2,1,1)
	width.temp                    <- ((x.max*0.75)/6)
	x.end.legend                  <- rollSum(rep(width.temp,6))
	x.start.legend                <- x.end.legend-(0.9*width.temp)
	legend.mat[,"start.plus"]     <- x.start.legend
	legend.mat[,"end.plus"]       <- x.end.legend
	lim.legend.mat                <- legend.mat[,c("start.plus","end.plus")]
	mode(lim.legend.mat)          <- "numeric"
	x.legend.labels               <- apply(X=lim.legend.mat,MARGIN=1,FUN=mean)
	y.legend.labels               <- -19
	legend.labels                 <- c("target sequence","probes","probe coverage","gene region","mRNA region","CDS region")
	segments(x0=as.numeric(legend.mat[,"start.plus"]),y0=as.numeric(legend.mat[,"y.coords"]),x1=as.numeric(legend.mat[,"end.plus"]),lwd=as.numeric(legend.mat[,"lwd"]),lend=2,col=legend.mat[,"segment.colors"],lty=as.numeric(legend.mat[,"line.type"]))
	text(legend.labels, x=x.legend.labels,y=y.legend.labels,cex=0.5)
	
	##### Calculates y-position where text and tracks will be drawn
	y.temp.all      <- ((1:length(i.values))*10)
	y.text.all      <- (y.temp.all+1.5)
	y.target.all    <- (y.temp.all+0)
	y.probes1.all   <- (y.temp.all+5)
	y.probes2.all   <- (y.temp.all+6)
	y.coverage.all  <- (y.temp.all+7)
	y.gene.all      <- (y.temp.all+1.5)
	y.mRNA.all      <- (y.temp.all+1.5)
	y.CDS.all       <- (y.temp.all+1.5)
	
	### Adds loci names to the current plot
	text(x=x.text,y=y.text.all,labels=loci.names[i.values],cex=0.75)
	
	#### This loop draws the annotation tracks for each locus of the kth plot.
	#### For a particular target locus, if annotation tracks and the target locus sequence are derived from different individuals, then a transformation table is used to calculate the transformed locations of annotations.
	for(j in 1:length(i.values)){
		i                  <- i.values[j]
		locus.name.temp    <- loci.names[i]
		probe.matches.temp <- which(locus.of.probe==locus.name.temp)
		probes.temp.start  <- probe.start[probe.matches.temp]
		probes.temp.end    <- probe.end[probe.matches.temp]
		
		y.temp             <- j*10        ### minimum y value of tracks for ith locus
		y.text             <- y.temp+1.5  ### vertical position where the ith locus name will be written
		y.target           <- y.temp+0    ## vertical position where the target locus will be drawn
		y.probes1          <- y.temp+5    ## vertical position where the first set of non-overlapping probes will be drawn
		y.probes2          <- y.temp+6    ## vertical position where the second set of non-overlapping probes will be drawn
		y.coverage         <- y.temp+7    ## vertical position where the coverage track will be drawn
		y.gene             <- y.temp+1.5  ## = (y.coverage+1)
		y.mRNA             <- y.temp+1.5  ## = (y.gene+1)
		y.CDS              <- y.temp+1.5  ## = (y.mRNA+1)

		probes2.unique.start.end <- unique(c(probes.temp.start[c(T,F)],probes.temp.end[c(T,F)]))
		probes1.unique.start.end <- unique(c(probes.temp.start[c(F,T)],probes.temp.end[c(F,T)]))
		coverage                 <- as.matrix(reduce(IRanges(start=probes.temp.start,probes.temp.end)))
		
		#### plotting gene, mRNA, and CDS features
		features.row.temp  <- which(features.table[,"target.locus"]==locus.name.temp)
		
		if(features.table[features.row.temp,"gene"]=="yes"){
			gene.ranges.string      <- features.table[features.row.temp,"Gene.Range"]
			gene.sense.string       <- features.table[features.row.temp,"Gene.target.Sense"]
			gene.ranges.mat         <- do.call(rbind,apply(X=mat.strsplit(x=gene.ranges.string,split=";|,",byrow=F,ncol=1),MARGIN=1,FUN=function(y){mat.strsplit(y,ncol=2,byrow=T,split="\\.\\.")}))
			gene.fuzzy.mat          <- matrix(c(1:length(gene.ranges.mat)) %in% grep(pattern="<|>",x=gene.ranges.mat),ncol=2)
			gene.sense.mat          <- mat.strsplit(gene.sense.string,split=";|,",byrow=F,ncol=1)
			gene.data.mat.tmp1      <- cbind(gene.ranges.mat,gene.fuzzy.mat,gene.sense.mat)
			gene.data.mat.tmp2      <- matrix(gsub("<|>","",gene.data.mat.tmp1),nrow=nrow(gene.data.mat.tmp1))
			### removes duplicate rows because no need to plot multiple genes on top of each other. Alternatively, keep all genes but filter exons with same range for different isoforms
			if(nrow(gene.data.mat.tmp2)>1){
				gene.data.mat.tmp3      <- unique(gene.data.mat.tmp2)
			} else {
				gene.data.mat.tmp3      <- gene.data.mat.tmp2
			}
			gene.data.mat           <- cbind(gene.data.mat.tmp3,rep("gene",nrow(gene.data.mat.tmp3)))
			colnames(gene.data.mat) <- c("start.plus","end.plus","start.fuzzy","end.fuzzy","sense","feature.type")
		}
		if(features.table[features.row.temp,"mRNA"]=="yes"){
			mRNA.ranges.string      <- features.table[features.row.temp,"mRNA.Range"]
			mRNA.sense.string       <- features.table[features.row.temp,"mRNA.target.Sense"]
			mRNA.ranges.mat         <- do.call(rbind,apply(X=mat.strsplit(x=mRNA.ranges.string,split=";|,",byrow=F,ncol=1),MARGIN=1,FUN=function(y){mat.strsplit(y,ncol=2,byrow=T,split="\\.\\.")}))
			mRNA.fuzzy.mat          <- matrix(c(1:length(mRNA.ranges.mat)) %in% grep(pattern="<|>",x=mRNA.ranges.mat),ncol=2)
			mRNA.sense.mat          <- mat.strsplit(mRNA.sense.string,split=";|,",byrow=F,ncol=1)
			mRNA.data.mat.tmp1      <- cbind(mRNA.ranges.mat,mRNA.fuzzy.mat,mRNA.sense.mat)
			mRNA.data.mat.tmp2      <- matrix(gsub("<|>","",mRNA.data.mat.tmp1),nrow=nrow(mRNA.data.mat.tmp1))
			mRNA.ranges.and.senses  <- paste(mRNA.data.mat.tmp2[,1],mRNA.data.mat.tmp2[,2],mRNA.data.mat.tmp2[,5])
			mRNAs.keep              <- match(mRNA.ranges.and.senses,mRNA.ranges.and.senses)
			mRNA.data.mat.tmp3      <- matrix(mRNA.data.mat.tmp2[mRNAs.keep,],ncol=ncol(mRNA.data.mat.tmp2))         ### removes duplicate rows because no need to plot mRNAs with same range but different transcript variants.
			mRNA.data.mat           <- cbind(mRNA.data.mat.tmp3,rep("mRNA",nrow(mRNA.data.mat.tmp3)))
			colnames(mRNA.data.mat) <- c("start.plus","end.plus","start.fuzzy","end.fuzzy","sense","feature.type")
		}
		if(features.table[features.row.temp,"CDS"]=="yes"){
			CDS.ranges.string      <- features.table[features.row.temp,"CDS.Range"]
			CDS.sense.string       <- features.table[features.row.temp,"CDS.target.Sense"]
			CDS.ranges.mat         <- do.call(rbind,apply(X=mat.strsplit(x=CDS.ranges.string,split=";|,",byrow=F,ncol=1),MARGIN=1,FUN=function(y){mat.strsplit(y,ncol=2,byrow=T,split="\\.\\.")}))
			CDS.fuzzy.mat          <- matrix(c(1:length(CDS.ranges.mat)) %in% grep(pattern="<|>",x=CDS.ranges.mat),ncol=2)
			CDS.sense.mat          <- mat.strsplit(CDS.sense.string,split=";|,",byrow=F,ncol=1)
			CDS.data.mat.tmp1      <- cbind(CDS.ranges.mat,CDS.fuzzy.mat,CDS.sense.mat)
			CDS.data.mat.tmp2      <- matrix(gsub("<|>","",CDS.data.mat.tmp1),nrow=nrow(CDS.data.mat.tmp1))
			if(nrow(CDS.data.mat.tmp2)>1){
				CDS.data.mat.tmp3 <- unique(CDS.data.mat.tmp2)
			} else {
				CDS.data.mat.tmp3 <- CDS.data.mat.tmp2
			}
			CDS.data.mat           <- cbind(CDS.data.mat.tmp3,rep("CDS",nrow(CDS.data.mat.tmp3)))
			colnames(CDS.data.mat) <- c("start.plus","end.plus","start.fuzzy","end.fuzzy","sense","feature.type")
		}
		
		target.data.mat    <- cbind(start.plus=1,end.plus=width(target.loci[i]),start.fuzzy="FALSE",end.fuzzy="FALSE",sense="1",feature.type="target.sequence")
		probes1.data.mat   <- cbind(start.plus=probes.temp.start[c(F,T)],end.plus=probes.temp.end[c(F,T)],start.fuzzy=rep("FALSE"),end.fuzzy=rep("FALSE"),sense="1",feature.type="probes.layer1")
		probes2.data.mat   <- cbind(start.plus=probes.temp.start[c(T,F)],end.plus=probes.temp.end[c(T,F)],start.fuzzy=rep("FALSE"),end.fuzzy=rep("FALSE"),sense="1",feature.type="probes.layer2")
		coverage.data.mat  <- cbind(start.plus=coverage[,1],end.plus=((coverage[,1]-1)+coverage[,2]),start.fuzzy=rep("FALSE"),end.fuzzy=rep("FALSE"),sense="1",feature.type="probe.coverage")
		feature.names.temp <- c("gene.data.mat","mRNA.data.mat","CDS.data.mat")[which(features.table[features.row.temp,c("gene","mRNA","CDS")]=="yes")]
		features.mat.temp  <- do.call(rbind,lapply(X=feature.names.temp,FUN=get))
		all.data.mat       <- rbind(features.mat.temp,target.data.mat,probes1.data.mat,probes2.data.mat,coverage.data.mat)
		all.data.mat       <- unique(all.data.mat)
		
		##########
		## splitting transform.ranges.table into a submatrix for ith locus (if ith target locus was designed from species other than Thamnophis)
		##########
		
		if(length(grep(paste0(locus.name.temp," "),transform.ranges.table[,1]))>0){
			need.to.transform          <- T
			mat.locus.temp             <- transform.ranges.table[grep(paste0(features.table[features.row.temp,1]," "),transform.ranges.table[,1]),]
			mat.locus.temp.split       <- cbind(column.split(mat.locus.temp[,1]),column.split(mat.locus.temp[,2]))
			mode(mat.locus.temp.split) <- "character"
		} else {
			need.to.transform          <- F
			seq.temp    <- do.call(rbind,strsplit(as.character(target.loci[i]),split=""))
			column1     <- rep(locus.name.temp,4)
			column2     <- c("target.column.numbers","aligned.target.sequence","homolog.column.numbers","aligned.homolog.sequence")
			columns3etc <- rbind(c(1:length(seq.temp)),seq.temp,c(1:length(seq.temp)),seq.temp)
			mat.locus.temp.split       <- cbind(column1,column2,columns3etc)
			mode(mat.locus.temp.split) <- "character"
			rownames(mat.locus.temp.split) <- colnames(mat.locus.temp.split) <-NULL
		}
		
		### Defining some empty vectors to be filled
		xmin.homolog.plus <- vector(nrow(all.data.mat),mode="numeric") ## empty vector that will contain xmin for the untransformed ranges
		xmax.homolog.plus <- vector(nrow(all.data.mat),mode="numeric") ## empty vector that will contain xmax for the untransformed ranges
		xmin.target.plus  <- vector(nrow(all.data.mat),mode="numeric") ## empty vector that will contain xmin for the transformed ranges
		xmax.target.plus  <- vector(nrow(all.data.mat),mode="numeric") ## empty vector that will contain xmax for the transformed ranges
		
		##########
		unaligned.target.column.numbers.temp  <- c(1:width(target.loci[i]))
		unaligned.homolog.column.numbers.temp <- c(1:length(which(mat.locus.temp.split[3,3:ncol(mat.locus.temp.split)]!="-")))
		aligned.target.column.numbers.temp    <- mat.locus.temp.split[1,3:ncol(mat.locus.temp.split)]
		aligned.homolog.column.numbers.temp   <- mat.locus.temp.split[3,3:ncol(mat.locus.temp.split)]
		plus.sense.rows                       <- which(all.data.mat[,"sense"]=="1")
		minus.sense.rows                      <- which(all.data.mat[,"sense"]=="-1")
		
		if(length(plus.sense.rows)>0){
			for(z in 1:length(plus.sense.rows)){
				if(all.data.mat[plus.sense.rows[z],"feature.type"] %in% c("target.sequence","probes.layer1","probes.layer2","probe.coverage")){
					xmin.target.plus[plus.sense.rows[z]]  <- as.numeric(all.data.mat[plus.sense.rows[z],"start.plus"])
					xmax.target.plus[plus.sense.rows[z]]  <- as.numeric(all.data.mat[plus.sense.rows[z],"end.plus"])
				}
				if(all.data.mat[plus.sense.rows[z],"feature.type"] %in% c("gene","mRNA","CDS")){
					#if(need.to.transform == T){
					#	xmin.target.plus[plus.sense.rows[z]]  <- as.numeric(all.data.mat[plus.sense.rows[z],"start.plus"])
					#	xmax.target.plus[plus.sense.rows[z]]  <- as.numeric(all.data.mat[plus.sense.rows[z],"end.plus"])
					#} else {
						start.feature.plus.temp               <- as.numeric(all.data.mat[plus.sense.rows[z],"start.plus"])
						end.feature.plus.temp                 <- as.numeric(all.data.mat[plus.sense.rows[z],"end.plus"])
						xmin.homolog.plus[plus.sense.rows[z]] <- start.feature.plus.temp
						xmax.homolog.plus[plus.sense.rows[z]] <- end.feature.plus.temp
						#start.temp <- mat.locus.temp.split[1,match(as.character(xmin.homolog.plus[plus.sense.rows[z]]),mat.locus.temp.split[3,])]
						#end.temp   <- mat.locus.temp.split[1,match(as.character(xmax.homolog.plus[plus.sense.rows[z]]),mat.locus.temp.split[3,])]
						target.col.range <- c(mat.locus.temp.split[1,match(as.character(xmin.homolog.plus[plus.sense.rows[z]]),mat.locus.temp.split[3,]):match(as.character(xmax.homolog.plus[plus.sense.rows[z]]),mat.locus.temp.split[3,])])
						start.temp       <- min(as.numeric(target.col.range[which(target.col.range!="-")]))
						end.temp         <- max(as.numeric(target.col.range[which(target.col.range!="-")]))
						#end.temp   <- max(mat.locus.temp.split[1,match(as.character(xmin.homolog.plus[plus.sense.rows[z]]),mat.locus.temp.split[3,]):match(as.character(xmax.homolog.plus[plus.sense.rows[z]]),mat.locus.temp.split[3,])])
						xmin.target.plus[plus.sense.rows[z]]  <- start.temp
						xmax.target.plus[plus.sense.rows[z]]  <- end.temp
					#}
				}
			}
		}
		
		###
		if(length(minus.sense.rows)>0){
			for(z in 1:length(minus.sense.rows)){
				start.feature.plus.temp                 <- as.numeric(all.data.mat[minus.sense.rows[z],"start.plus"])
				end.feature.plus.temp                   <- as.numeric(all.data.mat[minus.sense.rows[z],"end.plus"])
				xmin.homolog.plus[minus.sense.rows[z]]  <- rev(unaligned.homolog.column.numbers.temp)[end.feature.plus.temp]
				xmax.homolog.plus[minus.sense.rows[z]]  <- rev(unaligned.homolog.column.numbers.temp)[start.feature.plus.temp]
				target.col.range <- c(mat.locus.temp.split[1,match(as.character(xmin.homolog.plus[minus.sense.rows[z]]),mat.locus.temp.split[3,]):match(as.character(xmax.homolog.plus[minus.sense.rows[z]]),mat.locus.temp.split[3,])])
				start.temp       <- min(as.numeric(target.col.range[which(target.col.range!="-")]))
				end.temp         <- max(as.numeric(target.col.range[which(target.col.range!="-")]))
				#start.temp <- mat.locus.temp.split[1,match(as.character(xmin.homolog.plus[minus.sense.rows[z]]),mat.locus.temp.split[3,])]
				#end.temp   <- mat.locus.temp.split[1,match(as.character(xmax.homolog.plus[minus.sense.rows[z]]),mat.locus.temp.split[3,])]
				xmin.target.plus[minus.sense.rows[z]] <- rev(unaligned.target.column.numbers.temp)[as.numeric(end.temp)]
				xmax.target.plus[minus.sense.rows[z]] <- rev(unaligned.target.column.numbers.temp)[as.numeric(start.temp)]
			}
		}
		
		#### Makes a table containing parameters for plotting tracks
		track.order         <- c("target.sequence","gene","mRNA","CDS","probes.layer1","probes.layer2","probe.coverage")
		segment.colors      <- mgsub(c("gene","mRNA","CDS","probes.layer1","probes.layer2","probe.coverage","target.sequence"),c("green","purple","green","blue","blue","darkgreen","black"),all.data.mat[,"feature.type"])
		lwd                 <- as.numeric(mgsub(c("gene","mRNA","CDS","probes.layer1","probes.layer2","probe.coverage","target.sequence"),c(1.25,2.5,6,1,1,1,1.25),all.data.mat[,"feature.type"]))
		arrow.direction     <- mgsub(c("-1","1"),c("\U25C4" ,"\U25BA"),all.data.mat[,"sense"]) ### need to update this so that it only applies if feature type = "gene","mRNA",or "CDS"		
		arrow.color         <- mgsub(c("gene","mRNA","CDS","probes.layer1","probes.layer2","probe.coverage","target.sequence"),c("black","black","black","black","black","black","black"),all.data.mat[,"feature.type"])
		arrow.size          <- mgsub(c("gene","mRNA","CDS","probes.layer1","probes.layer2","probe.coverage","target.sequence"),c("0.25","0.4","0.4","0","0","0","0"),all.data.mat[,"feature.type"])
		arrow.transparency  <- mgsub(c("gene","mRNA","CDS","probes.layer1","probes.layer2","probe.coverage","target.sequence"),c("0","0","0","100","100","100","100"),all.data.mat[,"feature.type"])
		cap.transparency    <- mgsub(c("gene","mRNA","CDS","probes.layer1","probes.layer2","probe.coverage","target.sequence"),c("100","100","100","0","0","100","100"),all.data.mat[,"feature.type"])
		line.type           <- mgsub(c("gene","mRNA","CDS","probes.layer1","probes.layer2","probe.coverage","target.sequence"),c("2","1","1","1","1","1","1"),all.data.mat[,"feature.type"])
		if(any(all.data.mat[,"feature.type"] %in% c("gene","mRNA","CDS"))){
			row.order       <- mgsub(track.order,c(1,3,3,3,5,6,7),all.data.mat[,"feature.type"])
		} else {
			row.order       <- mgsub(track.order,c(1,0,0,0,2,3,4),all.data.mat[,"feature.type"])
		}
		y.coords            <- ((as.numeric(row.order)-1)+y.temp)
		track.width         <- as.numeric(mgsub(c("gene","mRNA","CDS","probes.layer1","probes.layer2","probe.coverage","target.sequence"),c(1.25,2.5,6,1,1,1,1.25),all.data.mat[,"feature.type"]))
		target.plus.mat     <- cbind(all.data.mat,xmin.target.plus,xmax.target.plus,y.coords,segment.colors,lwd,arrow.direction,arrow.color,arrow.size,arrow.transparency,cap.transparency,line.type,row.order)
		target.plus.mat     <- target.plus.mat[order(target.plus.mat[,"row.order"]),]
		
		###################
		##### plots tracks, direction pointers, and probe end-caps
		segments(x0=as.numeric(target.plus.mat[,"xmin.target.plus"]),y0=as.numeric(target.plus.mat[,"y.coords"]),x1=as.numeric(target.plus.mat[,"xmax.target.plus"]),lwd=as.numeric(target.plus.mat[,"lwd"]),lend=2,col=target.plus.mat[,"segment.colors"],lty=as.numeric(target.plus.mat[,"line.type"]),lend=1)
		### plots reading direction pointer of gene/mRNA/CDS features
		points(x=c(as.numeric(target.plus.mat[,"xmin.target.plus"]),as.numeric(target.plus.mat[,"xmax.target.plus"])),y=rep(as.numeric(target.plus.mat[,"y.coords"]),2),pch=rep(target.plus.mat[,"arrow.direction"],2),cex=rep(as.numeric(target.plus.mat[,"arrow.size"]),2),col= "black")
		### plots white end-caps at ends of probes (allows adjacent probes to be seen as distinct)
		text(x = c(as.numeric(target.plus.mat[,"xmin.target.plus"]),as.numeric(target.plus.mat[,"xmax.target.plus"])), y = rep(as.numeric(target.plus.mat[,"y.coords"]),2), labels= "*", cex=0.3, col = rep(t_col("white", perc = as.numeric(target.plus.mat[,"cap.transparency"])),2))
	} #end j loop
	
	### adds page number to lower right corner of plot
	mtext(side = 1, text = paste0(k,"/",num.plots), outer = TRUE,cex=0.75,adj=1)
	
	### saves kth (current) plot to a list that can be called later to export all plots as a multipage pdf
#	my.plot <- recordPlot()
#	dev.off()
#	my.plots[[k]] <- my.plot
}

### On OsX: the quartz function is used to save the plot
### On Windows: the cairo_pdf functon is used to save the plot

pdf.pathname <- "Target-loci_Coverage_graph_14April2020_scalation.pdf"
if(Sys.info()["sysname"]=="Darwin"){
	quartz(type = "pdf", file = pdf.pathname)
} else {
	cairo_pdf(pdf.pathname, onefile=TRUE, family="DejaVu Sans")
}
#for(k in 1:num.plots){
for(k in 130:135){
	replayPlot(my.plots[[k]])
}
dev.off()

