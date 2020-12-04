#####
## Dependencies
###
## ape, Biostrings, data.table, stringr, ...
#####


############
### str_locate_last function
#####
### Dependencies:
## libraries: stringr
############
## parameters
### strings = one or more character strings to search within
### pattern = a charcter string to search for
####
# returns the character position of the last occurrence of pattern
# if the pattern is >1 character long, the function returns the position
# of the last character of the last instance of the pattern

str_locate_last <- function(strings,pattern){
 	all.locs <- str_locate_all(strings,pattern)
 	last.locs <- list()
 	length(last.locs) <- length(strings)
 	for(i in 1:length(strings)){
 		last.locs[i] <- max(all.locs[[i]][,2])
 	}
 	last.locs <- unlist(last.locs)
 	last.locs
}

############
### str_locate_X function
############
### Dependencies:
## libraries: stringr
############
## parameters
### strings = one or more character strings to search within
### pattern = a charcter string to search for
### X       = which to locate (first, second, third,... etc.) Must supply an integer
####
# returns the character position of the Xth occurrence of pattern
# if the pattern is >1 character long, the function returns the 
# first and last position of the Xth instance of the pattern

str_locate_X <- function(strings,pattern,X){
 	all.locs <- str_locate_all(strings,pattern)
 	X.locs <- matrix(ncol=2,nrow=length(strings))
 	for(i in 1:length(strings)){
 		X.locs[i,1]      <- all.locs[[i]][X,1]
 		X.locs[i,2]      <- all.locs[[i]][X,2]
 	} 	
 	if(all(X.locs[,1]==X.locs[,2])){
 		X.locs <- c(X.locs[,1])
 	}
 	X.locs
}

###########
### nameFromPath function
###########
### Dependencies:
## libraries: stringr
############
## get the lowest level directory name or filename from the full pathname
#########

nameFromPath <- function(fullpath){
	pathLength <- nchar(fullpath)
	all.locs   <- str_locate_all(fullpath,"/")[[1]][,1]
	if(all(max(all.locs)==pathLength)){
		start.loc <- (all.locs[(length(all.locs)-1)]+1)
		end.loc   <- (pathLength-1)
		name <- substring(fullpath,first=start.loc,last=end.loc)
	} else {
		start.loc <- (max(all.locs)+1)
		name      <- substring(fullpath,first=start.loc)
	}
	name
}

############
### countChars function
############
## parameters
### x = character string
### pattern = a string to count in x

countChars <- function(x,pattern) {
	counter <- gsub(pattern,"",x)
	nchar(x)-nchar(counter)
}

########
## is.overlap
########
## test if two integer ranges overlap
########
## r1 = vector containing the limits of the first range (intermediate values optional)
## r2 = vector containing the limits of the second range (intermediate values optional)
########

is.overlap <- function(r1,r2){
	 min(r1) <= max(r2) && min(r2) <= max(r1)
}

########
## read.tree.multipleFiles function
########
### Dependencies
##    packages: ape
########
### Given a list of newick tree filenames or filepaths, returns a multiPhylo object containing all of the trees
########

read.tree.multipleFiles <- function(filepaths){
	treenames.local <- paste0("tree",seq(1:length(filepaths)))
	for(i in 1:length(filepaths)){
		assign(treenames.local[i],read.tree(filepaths[i]))
		if(i==1){
			all.trees <- get(treenames.local[i])
		} else {
			all.trees <- c(all.trees,get(treenames.local[i]))
		}
	}
	all.trees
}

############
### intersect.all function
############
### returns a vector of elements shared among a list of vectors
####
## parameters
### x = a list of vectors

intersect.all <- function(x){
	result <- Reduce(intersect,x)
	result
}

############
## mround function
####
## parameters
####
## x         = number to be rounded
## base      = what number to round by (e.g., base = 5 if "round to the nearest 5")
## direction = "nearest" (default), alternatively "up" or "down"

mround <- function(x,base=1,direction="nearest"){
	res <- base*round(x/base)
	if(direction=="up" & res < x){
		res <- res+base
	}
	if(direction=="down" & res > x){
		res <- res-base
	}
	res
}

############
### rollSum function
############
## parameters
#### x = a numerical vector
## value is a numerical vector, a rolling summation of the input vector x
####

rollSum <- function(x){
	rollingSum <- list()
	length(rollingSum) <- length(x)
	for(i in 1:length(x)){
		rollingSum[i] <- sum(x[1:i])
	}
	rollingSum <- unlist(rollingSum)
	rollingSum
}

######## string.list function
string.list <- function(y){
	paste.collapse <- function(z){
		z <- paste(z,collapse="")
		z
	}
	if(class(y)[1]!="matrix"){
		#y <- as.character(y)
		y <- as.matrix(y)
	}
	y <- apply(X=y,MARGIN=1,FUN=paste.collapse)
	y
}

######## mgsub function
### find and replace value in a vector with corresponding values of a second vector, in an object
mgsub <- function(x,y,z){
	for(i in 1:length(x)){
		z <- gsub(pattern=x[i],replacement=y[i],z)
	}
	z
}

#mgsub <- function(patt,repl,subj,complete.matches.only=FALSE){
#	if(complete.matches.only){
#		patt <- paste0("^",patt,"$") # "^" indicates beginning of string and "$" indicates end of string
#	}
#	for(i in 1:length(patt)){
#		subj <- gsub(pattern=patt[i],replacement=repl[i],subj)
#		
#	}
#	subj
#}

#################
######## mgrep function
### find location of each value in a vector query within a vector subject
mgrep <- function(query,subject){
	z <- NULL
	for(i in 1:length(query)){
		z[i] <- grep(pattern=query[i],x=subject)
	}
	z
}

######## is.wholenumber function
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
	abs(x - round(x)) < tol
}

######## mat.strsplit function
## split a string into a matrix containing the splitted strings
####
## x = character string to be split
## split = character delimeter to split by # default is to split at spaces
## ncol = "auto" {or numerical value} ## if "auto", the number of columns will be equal to the number of split string components (i.e., nrow=1)
## byrow=T {or F to fill by columns}

mat.strsplit <- function(x,split=" ",ncol="auto",byrow=T){
	
	split.pattern <- split
	step1         <- strsplit(x,split=split.pattern) ### a list of vectors
	
	if(byrow==T){
		result <- do.call(rbind,step1)
	}
	if(byrow==F){
		result <- do.call(cbind,step1)
	}
	
	#if(ncol=="auto" & byrow==T){
	#	result <- matrix(step1,ncol=length(step1),byrow=T)
	#}
	#if(ncol!="auto" & byrow==T) {
	#	n=ncol
	#	result <- matrix(step1,ncol=n,byrow=T)
	#}
	#if(ncol=="auto" & byrow==F){
	#	result <- matrix(step1,ncol=length(step1),byrow=F)
	#}
	#if(ncol!="auto" & byrow==F) {
	#	n = ncol
	#	result <- matrix(step1,ncol=n,byrow=F)
	#}
	result
}

#### column.split function
## split a matrix column (of class character) into multiple columns by a delimiter
####
## x = matrix columnd to be split
## split = " " character delimeter to split by  ### default is to split at spaces
## ncol = "auto" {or numerical value} ## if "auto", the number of columns will be equal to the number of split string components (i.e., nrow=1)
## byrow=T {or F to fill by columns}

column.split <- function(x,split=" ",ncol="auto",byrow=T){
	if(is.vector(x,mode="character")){
		x <- matrix(x,ncol=1)
	}
	result.temp <- apply(X=x,MARGIN=1,FUN=mat.strsplit)
	result      <- do.call(rbind, result.temp)
	result
}

### t_col function
#### makes a transparent version of a named color
##   color = color name
## percent = % transparency                 ### default is 50
##    name = an optional name for the color ### default is NULL
######
### source of function:
## Transparent colors
## Mark Gardener 2015
## www.dataanalytics.org.uk
######
t_col <- function(color, percent = 50, name = NULL) {
	rgb.val <- col2rgb(color)
	t.col   <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],max = 255,alpha = (100 - percent) * 255 / 100,names = name)
	invisible(t.col)
}

######## arrowLine function
## draw a line with arrows in it; e.g. -->-->-->
## this was taken from Stack Overflow (user "Sacha Epskamp")
arrowLine <- function(x0,y0,x1,y1=y0,nArrow=1,...){
	lines(c(x0,x1),c(y0,y1),...)
	Ax=seq(x0,x1,length=nArrow+1)
	Ay=seq(y0,y1,length=nArrow+1)
	for (i in 1:nArrow){
		arrows(Ax[i],Ay[i],Ax[i+1],Ay[i+1],...)
	}
}

######## dir.check.create function
## checks if a directory exists, and if not, creates it
## parent directories are also created if they do not already exist
dir.check.create <- function(directory){
	dir.parts <- unlist(strsplit(directory,"/"))
	for(i in 1:(length(dir.parts)-1)){
		dir.temp <- paste(dir.parts[1:(i+1)],collapse="/")
		if (!file.exists(dir.temp)){
			dir.create(dir.temp)
		}
	}
}

#####
## percent.missing.data function
#####
### calculates the percent of missing data in an alignment
### by default, missing characters include any of the following: "-","n","b","h","d","v","k","s","r","w","y","?"
### Therefore, ambiguous characters are treated as missing characters by default
#####
percent.missing.data <- function(alignment){
	store.class <- class(alignment)
	if(store.class %in% c("DNAMultipleAlignment","DNAStringSet")){
		alignment <- as.DNAbin(DNAMultipleAlignment(alignment))
		data.type <- "DNA"
	}
	if(store.class %in% c("AAMultipleAlignment","AAStringSet")){
		alignment <- as.AAbin(AAMultipleAlignment(alignment))
		data.type <- "AA"
	}
	x      <- as.character(alignment)
	y      <- table(x)
	n      <- c("-", "n", "b", "h", "d", "v", "k", "s", "r", "w","y","?")
	y.n    <- y[(names(y) %in% n)]
	result <- round((sum(y.n)/sum(y)*100),digits=2)
	result
}

######
## filter.alignment function
######
## filters individuals or sites (columns) from an alignment based on percent missing or ambiguous data per site, variation per site, frequency of common alleles per site, frequency of rare alleles per site, or specific sites (e.g., every second, third,
## Running with default values returns the input alignment.
######
## alignment               # input alignment
## mdt             = 1     # missing data threshold: percent of non-gap, non-ambiguous data required among individuals to keep an alignment column
## min.alleles     =       ### Not yet implemented, but hoping to use this to remove invariant sites
## min.var         = 1     # minimum number of individuals required to have the "rarest allele"; if min.freq.rare = 1 and min.freq.common = 1, then invariant sites are maintained.
## min.freq.common = 1     # a site is removed if the most common allele occurs in fewer than this many individuals
## min.freq.rare   = 1     # a site is removed if the rarest allele occurs in fewer than this many individuals
## keep.specific   = NA    # providing a numeric range or set allows you to extract specific sites
## drop.nodata.samples = T # drops an individual if it doesnt have any data for the sites that meet the other criteria
######
### examples:
## filter.alignment(alignment,min.var=2,min.freq.rare=1) ### returns the alignment containing only variable sites
## filter.alignment(alignment,min.var=2,min.freq.rare=2) ### returns the alignment containing only parsimony informative sites

filter.alignment <- function (alignment,mdt=1,min.var=1,min.freq.common=1,min.freq.rare=1,keep.specific=NA,drop.nodata.samples = T){
	######
	## Defining internal functions needed
	## for the whole function to work

	filter.md <- function(alignment.site) {                       ### alignment.site is one column of a DNA alignment
		y    <- table(alignment.site)                         ### This is the number of each type of character in a given column
		md.y <- (sum(y[which(names(y) %in% c("-"))])/sum(y))
		if(md.y < mdt|md.y==0){
			z = TRUE
		} else {
			z = FALSE
		}
		z
	}

    ## & length(y)>=min.hap

	#### pars.inf.dna function
	pars.inf.dna <- function(alignment.site) {	                          # alignment.site is one column of a DNA alignment
		y <- table(alignment.site)                                        # This is the number of each type of character in a given column
		n <- c("-", "n", "b", "h", "d", "v", "k", "s", "r", "w","y","?")  # A list of possible ambiguous or missing characters to ignore
		y <- y[!(names(y) %in% n)]
		# if(any(y>1) & length(y)>1){
		if(any(y>=min.freq.common) & length(y[y>=min.freq.rare])>=min.var){
			y <- TRUE
		} else {
			y <- FALSE
		}
		y
	}
	
	#### pars.inf.aa function
	pars.inf.aa <- function(alignment.site) {	   #|alignment.site is a column of an AA alignment
		y   <- table(alignment.site)               #|this is the number of each type of character in a given column
		n   <- c("-","X")                          #|a list of possible ambiguous or missing characters to ignore
		#AAs <- c(unique(GENETIC_CODE))      
		y <- y[!(names(y) %in% n)]
		#if(any(y>1) & length(y)>1){                #|tests if the column has at least two different characters, and if at least on of the characters if shared by more than one individul
		if(any(y>=min.freq.common) & length(y>=min.freq.rare)>=min.var){
			y <- TRUE                              
		} else {
			y <- FALSE
		}
		y
	}
	
	#### string.list function
	string.list <- function(y){
		paste.collapse <- function(z){
			z <- paste(z,collapse="")
			z
		}
		if(class(y)[1]!="matrix"){
			#y <- as.character(y)
			y <- as.matrix(y)
		}
		y <- apply(X=y,MARGIN=1,FUN=paste.collapse)
		y
	}
		
	#### charTest2.dna function
	#### tests whether at least one a, g, c, or t is present in a sequence
	charTest2.dna <- function(ls){
		ls <- tolower(ls)
		any(c("a","c","g","t") %in% ls)
	}
	
	#### charTest.aa function
	#### tests whether at least one non-missing and non-ambiguous amino acid is present in a sequence
	charTest.aa <- function(ls){
		#ls <- tolower(ls)
		#all(c("a","c","g","t") %in% ls)
		any(GENETIC_CODE %in% ls)
	}

	#### Start of actual function  ####

	store.class <- class(alignment)
	
	if(store.class %in% c("DNAMultipleAlignment","DNAStringSet")){
		alignment <- as.DNAbin(DNAMultipleAlignment(alignment))
	}
	
	if(store.class %in% c("AAMultipleAlignment","AAStringSet")){
		alignment <- as.AAbin(AAMultipleAlignment(alignment))
	}
	
	if(store.class == "DNAbin"){
		data.type <- "DNA"
	}
	
	if(store.class == "AAbin"){
		data.type <- "AA"
	}
	
	x <- as.character(alignment)
	
	if(is.na(all(keep.specific))){
		keep.cols3 <- rep(TRUE,ncol(x))
	} else {
		keep.cols3 <- c(1:ncol(x)) %in% keep.specific
	}

	if(data.type == "DNA"){
		keep.cols1 <- unlist(apply(x, 2, pars.inf.dna))        ### 2 is the columns margin to apply the function over
		keep.cols2 <- unlist(apply(x, 2, filter.md))           ### 
		keep.cols  <- keep.cols1 & keep.cols2 & keep.cols3     ### TRUE if keep.cols1, keep.cols2, and keep.cols3 are all TRUE for a particular site
		out        <- as.matrix(x[,keep.cols])
		keep.rows  <- apply(X=out,MARGIN=1,FUN=charTest2.dna)  ### A list of logicals indicating whether each individual has at least one non-missing and non-ambiguous character
	}
	
	if(data.type == "AA"){
		keep.cols1 <- unlist(apply(x, 2, pars.inf.aa))         ### 2 is the columns margin to apply the function over
		keep.cols2 <- unlist(apply(x, 2, filter.md))           ### 
		keep.cols  <- keep.cols1 & keep.cols2 & keep.cols3     ### TRUE if keep.cols1, keep.cols2, and keep.cols3 are all TRUE for a particular site
		out        <- as.matrix(x[,keep.cols])
		keep.rows  <- apply(X=out,MARGIN=1,FUN=charTest.aa)    ### A list of logicals indicating whether each individual has at least one non-missing and non-ambiguous character
	}

	if(drop.nodata.samples == T){
		out        <- out[keep.rows,]                          ### Only keep individuals that have some data present (i.e., at least one A, C, G, or T).
	}
	
	if(store.class=="DNAbin"){
		out <- as.DNAbin(out)
	}
	if(store.class=="DNAStringSet"){
		out <- DNAStringSet(string.list(out))
	}
	if(store.class=="DNAMultipleAlignment"){
		out <- DNAMultipleAlignment(string.list(out))
	}
	if(store.class=="AAbin"){
		out <- as.AAbin(out)
	}
	if(store.class=="AAStringSet"){
		out <- AAStringSet(string.list(out))
	}
	if(store.class=="AAMultipleAlignment"){
		out <- AAMultipleAlignment(string.list(out))
	}
	out
}

#######
## bait.coverage function
####
# targets       # set of target sequences
# baits         # set of baits, each associated with a target sequence. Multiple baits can be associated with one target sequence.
#               # baits must be names using the following format: targetName_startPosistionInTarget
# distance = 0  # distance between target site and bait (probe) to be considered covered. Default is 0.
# bait.length = 120 # length of baits (probes), 60 or 90 nt are also common
bait.coverage <- function(targets,baits,distance=0,bait.length=120,fraction=T) {
	bait.matrix           <- do.call(rbind,strsplit(names(baits),split="_"))
	bait.matrix[,2]       <- as.numeric(bait.matrix[,2])+1
	bait.matrix           <- cbind(bait.matrix,as.numeric(bait.matrix[,2])+(bait.length-1))
	rownames(bait.matrix) <- names(baits)
	colnames(bait.matrix) <- c("target","probe.start.in.target","probe.end.in.target")
	result                <- vector(mode="numeric",length=length(targets))
	buffer                <- distance
	for(i in 1:length(targets)){
		bait.matrix.temp       <- bait.matrix[which(bait.matrix[,"target"]==names(targets[i])),,drop=F]
		ranges.matrix.temp     <- cbind(as.numeric(bait.matrix.temp[,"probe.start.in.target"]),as.numeric(bait.matrix.temp[,"probe.end.in.target"]))
		buffered.ranges.matrix <- cbind((ranges.matrix.temp[,1]-buffer),(ranges.matrix.temp[,2]+buffer))
		covered.temp           <- unique(as.numeric(apply(X=buffered.ranges.matrix,MARGIN=1,FUN=function(X){seq(from=X[1],to=X[2])})))
		target.temp.sites      <- 1:width(targets[i])
		covered                <- intersect(target.temp.sites,covered.temp)
		if(fraction){
			bait.coverage          <- length(covered)/width(targets[i])
		} else {
			bait.coverage          <- length(covered)
		}
		result[i]              <- bait.coverage
	}
	result
}

##########
### na.replace function
#######
### used to replace NA values of data.table with a desired value (usually "-")
#######
na.replace <- function(v,value=x){
		v[is.na(v)] <- value
		v
}

##############
## simple.concatenate.DNA function
######
## Used to concatenate multiple DNAStringSet or DNAMultipleAlignment files
## Concatenation order is the same as the input argument order
######
simple.concatenate.DNA <- function(...,make.partitions.table=F){
	list.of.alignments <- list(...)
	for(i in seq(length(list.of.alignments))){
		names(list.of.alignments)[i]  <-  as.character(sys.call()[[i+1]])
	}
	list.of.loci       <- names(list.of.alignments)
	
	# alignments.temp <- lapply(X=list.of.alignments,FUN= function(input){data.table(as.matrix(input),keep.rownames = TRUE)})
	# cols.temp       <- lapply(X=alignments.temp,FUN=function(input){c(1:(ncol(input)-1))})
	
	for(i in 1:length(list.of.alignments)){
		alignment.temp           <- data.table(as.matrix(list.of.alignments[[i]]),keep.rownames = TRUE)
		locus.name.temp          <- list.of.loci[i]
		col.temp                 <- c(1:(ncol(alignment.temp)-1))
		colnames(alignment.temp) <- c("rn",paste(locus.name.temp,col.temp,sep="."))                          ### redefining column names of data.table so that duplicate colnames wont exist when concatenating

		if(i==1){
			alignment.current   <- alignment.temp
			start.current       <- 1
			end.current         <- (ncol(alignment.current)-1)
			app                 <- FALSE
			rn                  <- alignment.current$rn
			alignment.current   <- apply(alignment.current[ ,!"rn"], 1, paste, collapse="")
			alignment.current   <- data.table(cbind(rn,alignment.current))
		
		} else {
			start.current      <- (end.current+1)
			width.current      <- max(nchar(alignment.current$alignment.current))
			alignment.current1 <- merge(alignment.current, alignment.temp, by="rn", all=TRUE)  ### concatenated data frame
			alignment.current1$alignment.current <- na.replace(alignment.current1$alignment.current,paste(rep("-",width.current),collapse=""))
			#alignment.current <- na.replace(merge(alignment.current, alignment.temp2, by="rn", all=TRUE),"-")  ### concatenated data frame
			alignment.current2 <- na.replace(alignment.current1,"-")
			alignment.current <- alignment.current2
			end.current       <- (start.current+(ncol(alignment.temp2)-2))
			rn                <- alignment.current$rn
			alignment.current <- apply(alignment.current[ ,!"rn"], 1, paste, collapse="")
			alignment.current <- data.table(cbind(rn,alignment.current))
			app               <- TRUE
		}
	}
	alignment.write    <- DNAStringSet(apply(alignment.current[ ,!"rn"], 1, paste, collapse=""))
	names(alignment.write) <- alignment.current$rn
	alignment.write
}

##############
### concatenate.alignments function
######
### used to concatenate multiple DNA or AA alignments stored in
### separate files of a folder; by default, creates a partition file
######
### if DNA, input alignments must be in sequential phylip format
### if AA, input alignment must be in interleaved fasta format
######
## Dependencies (and their dependencies):
##    custom functions: filter.alignment, is.wholenumber, na.replace
##    packages: Biostrings, ape, data.table
###############
## parameters
###
## alignment.folder         ### a directory containing gene alignents in separate files
## type="DNA"               ### type of sequence data. "AA" is also possible
## seqs.format="phylip"     ### can also be "fasta"
## alignment.output=NA
## partition.output=NA
## show.progress=T
## istart=NA
## iend=NA
## only.variable.sites=F
concatenate.alignments    <- function(alignment.folder,type="DNA",seqs.format="phylip",alignment.output=NA,partition.output=NA,show.progress=T,istart=NA,iend=NA,only.variable.sites=F){
	alignment.filenames   <- list.files(alignment.folder,full.names=T)
	alignment.shortnames  <- gsub("\\..+","",list.files(alignment.folder))

	if(is.na(istart)){
		istart <- 1
	}
	if(is.na(iend)){
		iend <- length(alignment.filenames)
	}
	
	for(i in istart:iend){
		locus.name.temp <- alignment.shortnames[i]
		if(show.progress==T & is.wholenumber(i/50)){
			print(paste(i," of ",length(alignment.filenames),sep=""))
		}
		if(type=="DNA"){
			alignment.temp   <- readDNAMultipleAlignment(alignment.filenames[i],format=seqs.format)
		}
		if(type=="AA"){
			alignment.temp   <- readAAMultipleAlignment(alignment.filenames[i],format="fasta")
		}
		if(all(names(unmasked(alignment.temp)) %in% c(1:nrow(alignment.temp)))){
			next
		} else {
			if(only.variable.sites==T){
				alignment.temp2      <- filter.alignment(alignment.temp,min.var=2)
			}
			if(only.variable.sites==F){
				alignment.temp2      <- filter.alignment(alignment.temp)
			}
			
		}
		if(ncol(alignment.temp2)==0){
			next
		} else {
			alignment.temp2  <- data.table(as.matrix(alignment.temp2),keep.rownames = TRUE)                ### coerces alignment into data.table
			col.temp  <- c(1:(ncol(alignment.temp2)-1))                                                    ### subtract 1 because first column of data.table is sample name
			### Next line is used for SnakeCap loci but wont harm other data such as Sanger data
			key.temp  <- gsub("WeinellEntry","WE",locus.name.temp)                                         ### shorter version of SnakeCap locus name but wont be generalizable across datasets
			colnames(alignment.temp2) <- c("rn",paste(key.temp,col.temp,sep="."))                          ### redefining column names of data.table so that duplicate colnames wont exist when concatenating
		}
		if(i==istart){
			alignment.current   <- alignment.temp2
			start.current       <- 1
			end.current         <- (ncol(alignment.current)-1) 
			app                 <- FALSE
			rn                  <- alignment.current$rn
			alignment.current   <- apply(alignment.current[ ,!"rn"], 1, paste, collapse="")
			alignment.current   <- data.table(cbind(rn,alignment.current))
		
		} else {
			start.current     <- (end.current+1)
			width.current     <- max(nchar(alignment.current$alignment.current))
			alignment.current1 <- merge(alignment.current, alignment.temp2, by="rn", all=TRUE)  ### concatenated data frame
			#alignment.current2  <- na.replace(alignment.current1$alignment.current,rep("-",width.current))
			alignment.current1$alignment.current <- na.replace(alignment.current1$alignment.current,paste(rep("-",width.current),collapse=""))
			#alignment.current <- na.replace(merge(alignment.current, alignment.temp2, by="rn", all=TRUE),"-")  ### concatenated data frame
			alignment.current2 <- na.replace(alignment.current1,"-")
			alignment.current <- alignment.current2
			end.current       <- (start.current+(ncol(alignment.temp2)-2))
			rn                <- alignment.current$rn
			alignment.current <- apply(alignment.current[ ,!"rn"], 1, paste, collapse="")
			alignment.current <- data.table(cbind(rn,alignment.current))
			app               <- TRUE
		}
		if(type=="DNA"){
			test.rows <- DNAStringSet(apply(alignment.current[ ,!"rn"], 1, paste, collapse=""))
		}
		if(type=="AA"){
			test.rows <- AAStringSet(apply(alignment.current[ ,!"rn"], 1, paste, collapse=""))
		}
		if(!all(width(test.rows)==width(test.rows[1]))){
			print(paste(i, "error",sep=" "))
		}
		text.current = paste(type,", ",locus.name.temp,"_VarSites = ",start.current,"-",end.current,sep="")
		if(!is.na(partition.output)){
			write(text.current,file=partition.output,append=app)
		}
	}
	alignment.temp3        <- apply(alignment.current[ ,!"rn"], 1, paste, collapse="")
	if(type=="DNA"){
		alignment.write    <- DNAStringSet(alignment.temp3)
	}
	if(type=="AA"){
		alignment.write    <- AAStringSet(alignment.temp3)
	}
	names(alignment.write) <- alignment.current$rn
	if(!is.na(alignment.output)){
		writeXStringSet(x=alignment.write, filepath=alignment.output, append = F)
	}
	alignment.write
}


###############################
### concatenateDNAStringSets ###
###############################
### parameters:
## ... = multiple DNAStringSet objects separated by commas
#####
concatenateDNAStringSets <- function(...){
	x              <- list(...)
	x.names        <- lapply(X=x,FUN=names)        ### list of character vectors, each containing the names of the sequences in the corresponding DNAStringSet
	x.names2       <- lapply(X=x.names,FUN=sort)   ### each character vector of names is sorted alphabetically
	names.matrix   <- do.call(rbind,x.names2)      ### character matrix. each row is a vector in x.names2
	unique.names   <- unique(names.matrix)         ### 
	if(nrow(unique.names)>1){
		result = print("one or more sequence names not in all input alignments")
	} else {
		x.characters   <- lapply(X=x,FUN=as.character)
		x.mats         <- do.call(cbind,x.characters)
		concat.list    <- apply(X=x.mats,MARGIN=1,FUN=function(y){paste(y,collapse="")})
		result         <- DNAStringSet(concat.list)
	}
	result
}

##################################
### writeDNAStringSet function ###
##################################
## This is like writeXStringSet (Biostrings package) except data can be written in sequential fasta format.
## Requires the "ape" package.
### parameters:
## x                      ### an object of class DNAStringSet
## filepath               ### where to save the file
## append = F             ### whether or not to add data the end of a file
## format = "fasta"       ### no options currently possible...
## charsPerLine = 100000  ### maximum number of characters to write on each line
######
writeDNAStringSet <- function(x,filepath,append=F,format="fasta",charsPerLine=100000){
	data         <- as.DNAbin(x)
	param.append <- append
	if(format=="phylip"){
		param.format = "sequential"
	}
	if(format=="fasta"){
		param.format = "fasta"
	}
	write.dna(x=data,file=filepath,format=param.format,append=param.append,nbcol=1,colsep="",colw=charsPerLine)
}

##########################
### filter.annotationTable function
#####################
## parameters # values last used are shown
### input.gff   <- "ref_Thamnophis_sirtalis-6.0_top_level_JLW.gff3"                 ### path to unfiltered genome annotation file
### output.gff  <- "CDS_ref_Thamnophis_sirtalis-6.0_top_level_JLW_longer120bp.gff3" ### path and filename where to write filtered annotation table
### region.type <- "CDS"
### min.length  <- 120
filter.annotationTable <- function(input.gff,output.gff,region.type="CDS",min.length=120) {
	unfiltered.gff   <- fread(input.gff)
	filtered.gff1A   <- unfiltered.gff[which(unfiltered.gff$region==region.type),]
	widths1A         <- abs(as.numeric(filtered.gff1A$start)-as.numeric(filtered.gff1A$end))	
	filtered.gff1B   <- filtered.gff1A[which(widths1A>=min.length),]
	write.table(x=filtered.gff1B,file=output.gff,sep="\t",row.names=F)
}

##########################
### get.exome.from.annotationTable function ###
#####################
## parameters # values last used are shown
### species.name      <- "Thamnophis_sirtalis"
### genome.filepath   <- "Thamnophis_sirtalis_GCF_001077635.1_genome_renamed_sequential.fas"   ### directory to genome
### input.gff         <- "CDS_ref_Thamnophis_sirtalis-6.0_top_level_JLW_longer120bp.gff3"      ### path to filtered annotation table
### output.dir.exome  <- "Exomes_TempFolder_3Nov2019/"                                         ### directory to put exomes into
### additional.ID     <- "Scaffold-Name-Key.txt"                                               ### filename of a table that cross-references contig names to names in GFF file
get.exome.from.annotationTable <- function(species.name,genome.filepath,input.gff,output.dir.exome,additional.ID) {
	filtered.gff1B    <- fread(input.gff)
	refseq.names     <- unlist(filtered.gff1B[,1])
	ScaffoldKey      <- fread(additional.ID)                                                                 ### this is needed because the names in the gff file are not exactly the same as the ones in the genome file; they both link to a common Genbank Accession # though
	alt.scaff.names  <- ScaffoldKey$ScaffoldName[match(refseq.names,ScaffoldKey$RefSeq.ScaffoldAccession)]   ### links the names in the gff file to the names of the CDS matches

	subject.id        <- alt.scaff.names
	subject.start     <- as.numeric(filtered.gff1B$start)
	subject.end       <- as.numeric(filtered.gff1B$end)
	new.names         <- paste(subject.id,"_",refseq.names,":",subject.start,"-",subject.end,sep="")	

	indexFa(genome.filepath)                                           ### create an index of file 'foo.fasta'; this avoids having to actually copy or move the file to a new directory
	fa <- FaFile(genome.filepath)
	gr <- as(seqinfo(fa), "GRanges")

	scaff.matches.all   <- match(subject.id, names(gr))
	start.all           <- apply(X=cbind(subject.start,subject.end),MARGIN=1,FUN=min)
	end.all             <- apply(X=cbind(subject.start,subject.end),MARGIN=1,FUN=max)
	subranges           <- IRanges(start=start.all ,end=end.all,names=subject.id)
	gsubranges          <- GRanges(seqnames=subject.id,ranges=subranges)
	genome.scaff        <- getSeq(fa, gsubranges)
	names(genome.scaff) <- new.names
	writeXStringSet(x = genome.scaff, filepath=paste(output.dir.exome,species.name,"_exome_longer120bp.fas",sep=""), append=F, format="fasta")
}

##########################
### get.loci.from.annotationTable function ###
#####################
### identical to the function get.exome.from.annotationTable, but used more generally.
### I used this function to extract the MHC loci from the Thamnophis sirtalis genome
#######
## parameters # values last used are shown
### species.name      <- "Thamnophis_sirtalis"
### genome.filepath   <- "Thamnophis_sirtalis_GCF_001077635.1_genome_renamed_sequential.fas"   ### directory to genome
### input.gff         <- "ref_Thamnophis_sirtalis-6.0_top_level_MHC.gff3"                      ### path to filtered annotation table containing only the loci of interest
### output.file       <- "MHC-loci.fas"                                                        ### directory to put exomes into
### additional.ID     <- "Scaffold-Name-Key.txt"                                               ### filename of a table that cross-references contig names to names in GFF file

get.loci.from.annotationTable <- function(species.name,genome.filepath,input.gff,output.file,additional.ID) {
	filtered.gff1B   <- fread(input.gff)
	refseq.names     <- unlist(filtered.gff1B[,1])
	ScaffoldKey      <- fread(additional.ID)                                                                 ### this is needed because the names in the gff file are not exactly the same as the ones in the genome file; they both link to a common Genbank Accession # though
	alt.scaff.names  <- ScaffoldKey$ScaffoldName[match(refseq.names,ScaffoldKey$RefSeq.ScaffoldAccession)]   ### links the names in the gff file to the names of the CDS matches

	subject.id        <- alt.scaff.names
	subject.start     <- as.numeric(filtered.gff1B$start)
	subject.end       <- as.numeric(filtered.gff1B$end)
	new.names         <- paste(subject.id,"_",refseq.names,":",subject.start,"-",subject.end,sep="")	

	indexFa(genome.filepath)                                           ### create an index of file 'foo.fasta'; this avoids having to actually copy or move the file to a new directory
	fa <- FaFile(genome.filepath)
	gr <- as(seqinfo(fa), "GRanges")

	scaff.matches.all   <- match(subject.id, names(gr))
	start.all           <- apply(X=cbind(subject.start,subject.end),MARGIN=1,FUN=min)
	end.all             <- apply(X=cbind(subject.start,subject.end),MARGIN=1,FUN=max)
	subranges           <- IRanges(start=start.all ,end=end.all,names=subject.id)
	gsubranges          <- GRanges(seqnames=subject.id,ranges=subranges)
	genome.scaff        <- getSeq(fa, gsubranges)
	names(genome.scaff) <- new.names
	writeXStringSet(x = genome.scaff, filepath=output.file, append=F, format="fasta")
}

############
## reportBestMatches function
############
## parameters
#### input.dir      <- "TBlastXResults/" ### must have that last slash mark
#### output.dir     <- input.dir
#### species.names  <- c("Anolis_carolinensis","Crotalus_mitchellii","Pogona_vitticeps","Gekko_japonicus","Ophiophagus_hannah","Vipera_berus","Crotalus_horridus","Thamnophis_sirtalis","Python_bivittatus","Protobothrops_mucrosquamatus","Pantherophis_guttatus")
#### blastMethod    <- "tblastx"   ### change to "blastn" if type = "UCEs"
#### locusType      <- "exons"     ### other options include "UCEs"...
reportBestMatches <- function(input.dir,species.names,output.dir=NA,blastMethod="tblastx",locusType="exons") {
	if(is.na(output.dir)){
		output.dir <- input.dir
	}
	if(substring(input.dir,first=nchar(input.dir))!="/"){
		input.dir <- paste(input.dir,"/",sep="")
	}
	if(substring(output.dir,first=nchar(output.dir))!="/"){
		output.dir <- paste(output.dir,"/",sep="")
	}
	for(i in 1:length(species.names)){
		species.temp      <- species.names[i]
		data.50           <- fread(input=paste(input.dir,species.temp,".",blastMethod,".",locusType,".50hits.txt",sep=""),sep="\t",header=F)
		colnames(data.50) <- c("gseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
		data.50.ordered   <- data.50[order(data.50$gseqid, data.50$bitscore)]
		best.matches      <- match(unique(data.50.ordered$gseqid), data.50.ordered$gseqid)
		best.data         <- data.50[best.matches,]
		write.table(best.data,file=paste(output.dir,species.temp,".",blastMethod,".",locusType,".best.txt",sep=""),sep=",",row.names=F,col.names=T)
	}
}

##########################
### get.exome.from.blastTable function ###
#####################
## parameters   # values last used are shown
#### species           <- "Vipera_berus"   ### species names
#### genome.filepath   <- "Vipera_berus_GCA_000800605.1_Vber.be_1.0_genomic.fna"  ### genome filepath
#### input.blastTable  <- "Vipera_berus.tblastx.exons.best.txt"                   ### BLAST hit table filepath
#### output.dir.exome  <- "Exomes_TempFolder_3Nov2019/"                           ### directory to put exomes into
get.exome.from.blastTable <- function(species,genome.filepath,input.blastTable,output.dir.exome) {

	data.best           <- fread(file=input.blastTable,sep=",",header=T)
	
	query.subject.id    <- paste(as.character(data.best$gseqid),"_Subject=",as.character(data.best$sseqid),"_",as.character(data.best$sstart),"_",as.character(data.best$send),sep="")
	subject.id          <- gsub(".+_","",as.character(data.best$sseqid))
	subject.start       <- data.best$sstart
	subject.end         <- data.best$send

	indexFa(genome.filepath)                                           ### create an index of file 'foo.fasta'; this avoids having to actually copy or move the file to a new directory
	fa <- FaFile(genome.filepath)
	gr <- as(seqinfo(fa), "GRanges")

	scaff.matches.all   <- match(subject.id, names(gr))
	start.all           <- apply(X=cbind(subject.start,subject.end),MARGIN=1,FUN=min)
	end.all             <- apply(X=cbind(subject.start,subject.end),MARGIN=1,FUN=max)
	subranges           <- IRanges(start=start.all ,end=end.all,names=subject.id)
	gsubranges          <- GRanges(seqnames=subject.id,ranges=subranges)
	exome.scaff         <- getSeq(fa, gsubranges)
	names(exome.scaff)  <- query.subject.id
	writeXStringSet(x = exome.scaff, filepath=paste(output.dir.exome,species,"_exome.fas",sep=""), append=F, format="fasta")
}

###########
### get.UCEs.from.blastTable function
###########
### (This function is identical to get.exome.from.blastTable function)
## parameters   # values last used are shown
#### species           <- "Vipera_berus"                                               ### species names
#### genome.filepath   <- "Vipera_berus_GCA_000800605.1_Vber.be_1.0_genomic.fna"       ### genome filepath
#### input.blastTable  <- ""                                                           ### BLAST hit table filepath
#### output.dir        <- ""                                                           ### directory to save UCE-containing files
get.UCEs.from.blastTable <- function(species,genome.filepath,input.blastTable,output.dir) {
	data.best           <- fread(file=input.blastTable,sep=",",header=T)
	query.subject.id    <- paste(as.character(data.best$gseqid),"_Subject=",as.character(data.best$sseqid),"_",as.character(data.best$sstart),"_",as.character(data.best$send),sep="")
	subject.id          <- gsub(".+_","",as.character(data.best$sseqid))
	subject.start       <- data.best$sstart
	subject.end         <- data.best$send
	indexFa(genome.filepath)                                           ### create an index of file 'foo.fasta'; this avoids having to actually copy or move the file to a new directory
	fa <- FaFile(genome.filepath)
	gr <- as(seqinfo(fa), "GRanges")
	scaff.matches.all   <- match(subject.id, names(gr))
	start.all           <- apply(X=cbind(subject.start,subject.end),MARGIN=1,FUN=min)
	end.all             <- apply(X=cbind(subject.start,subject.end),MARGIN=1,FUN=max)
	subranges           <- IRanges(start=start.all ,end=end.all,names=subject.id)
	gsubranges          <- GRanges(seqnames=subject.id,ranges=subranges)
	UCE.scaff           <- getSeq(fa, gsubranges)
	names(UCE.scaff)    <- query.subject.id
	writeXStringSet(x = UCE.scaff, filepath=paste(output.dir,species,"_UCEs.fas",sep=""), append=F, format="fasta")
}

##############
## makeExomeStatsTable function
##############
## parameters
#### exomes.filepaths     # paramA <- list.files(path="~/Exomes/",full.names=T)
#### is.primary.exome     # paramB <- 1  ### indicates which of the species names correspond to the "subject" species during the tblastn step. This currently only works if set to 1.
#### annotationTable.path # paramC <- "CDS_ref_Thamnophis_sirtalis-6.0_top_level_JLW_withGenBankAcc_longer120bp.gff3"
#### output.dir           # paramD <- "Exomes_TempFolder_3Nov2019/"  ### location to write Exome Stats Table, currently set as folder for testing
#### species              # paramE <- c("Thamnophis_sirtalis","Ophiophagus_hannah","Crotalus_mitchellii","Python_bivittatus","Vipera_berus","Crotalus_horridus","Protobothrops_mucrosquamatus","Pantherophis_guttatus","Anolis_carolinensis","Pogona_vitticeps","Gekko_japonicus")
#### subgroup <-          # paramF <- c("Thamnophis_sirtalis","Ophiophagus_hannah","Crotalus_mitchellii","Python_bivittatus","Vipera_berus","Crotalus_horridus","Protobothrops_mucrosquamatus","Pantherophis_guttatus")  ### a subset of the full list of species to estimate stats for (here, subgroup includes only the snakes)
## annotationTable.path   # paramC <- "CDS_ref_Thamnophis_sirtalis-6.0_top_level_JLW_longer120bp.gff3"
makeExomeStatsTable <- function(exomes.filepaths,annotationTable.path,species,subgroup,output.dir,is.primary.exome=1,i.start=1,i.stop=NA){
	species              <- c(species[is.primary.exome],species[-is.primary.exome])         ### reorders the list of species such that the primary.species is first in the list
	is.subgroup          <- mgrep(query=subgroup,subject=species)                           ### a list of numbers indicating which of the species are also in the subgroup
	exomes.filepaths     <- exomes.filepaths[mgrep(species,exomes.filepaths)]               ### puts the exomes.filepaths in the same order as species
	exome.names          <- paste("exome",c(1:length(exomes.filepaths)),sep="")
	matches.names        <- paste("matches.exome",c(1:length(exomes.filepaths)),sep="")
	
	for(i in 1:length(exome.names)){                                                       #|reads in exomes
		assign(x=exome.names[i],value=readDNAStringSet(filepath=exomes.filepaths[i]))  #|and assigns them
	}                                                                                      #|to object names

	annotationTable  <- fread(input=annotationTable.path)   ### reads in annotation table, which is associated with the primary exome
	colnames(annotationTable)[9] <- "ID"

	CDS.names        <- names(exome1)
	unique.CDS.names <- unique(CDS.names)
	to.keep          <- match(unique.CDS.names, CDS.names)
	CDS.names        <- unique.CDS.names
	exome1           <- exome1[to.keep]   ### updates exome1 to remove duplicate sequences (e.g., because the region has multiple annotations as a result of multiple isomers). This may have already been done in updated filter.annotationTable function
	
	for(j in c(1:length(matches.names))){
		assign(x=matches.names[j],value=mgsub(paste("_Subject=",species,".*",sep=""),rep("",length(species)), names(get(exome.names[j]))))
	}

	CDS.names <- intersect.all(lapply(X=matches.names,FUN=get))  ### This is the set of shared exons

	annotationTable.identifier  <- paste(unlist(annotationTable[,1]),"_",unlist(annotationTable$start),"_", unlist(annotationTable$end),sep="")  #|For each shared exon, extracts the 
	CDS.locus.identifier        <- gsub(pattern="_Subject.*\\.1_","_",x=CDS.names)                                                               #|gene name from the last column of the 
	match.identifier            <- match(CDS.locus.identifier,annotationTable.identifier)                                                        #|filtered annotation table. This part of the code
	gene.names.temp             <- gsub(".*;gene=","gene=",unlist(annotationTable$ID))                                                           #|may fail if other annotation table formats were used.
	gene.names.temp2            <- gsub(";.*","",gene.names.temp)                                                                                #|
	gene.names                  <- gene.names.temp2[match.identifier]                                                                            #|

	CDS.ranges  <- gsub(".*\\.1_","",CDS.names)                                                          ### character string location of CDS within contig, with format "start_end"
	CDS.lengths <- (abs(as.numeric(gsub(".*_","",CDS.ranges))-as.numeric(gsub("_.*","",CDS.ranges)))+1)  ### length of each exon (in primary exome, T. sirtalis for SnakeCap) for shared exons
	
	if(is.na(i.start) | i.start > length(CDS.names)){
		i.start <- 1                  ### default i.start, ie which exons to start the loop at
	}
	if(is.na(i.stop) | i.stop > length(CDS.names)){
		i.stop  <- length(CDS.names)   ### default i.stop, ie which exons to stop the loop at
	}
	
	for(i in i.start:i.stop){
		#tempMatrix     <- matrix(data=0, nrow=1, ncol=19)
		tempMatrix      <- matrix(data=0, nrow=1, ncol=12+length(species))
		tempMatrix[1]   <- CDS.names[i]    ### stats.data.exome[,1]  <- CDS.names
		tempMatrix[19]  <- gene.names[i]   ### stats.data.exome[,19] <- gene.names
		
		temp.exon.names <- paste("temp.exon.exome",c(1:length(species)),sep="")
		
		for(k in 1:length(temp.exon.names)){
			assign(x=temp.exon.names[k],value = get(exome.names[k])[grep(CDS.names[i],names(get(exome.names[k])))])
		}
		
		temp.seqs <- DNAStringSet(get(temp.exon.names[1]))
		for(z in 2:length(temp.exon.names)){
			temp.seqs <- c(temp.seqs,DNAStringSet(get(temp.exon.names[z])))
		}
		num.species   <- length(temp.seqs)     ### number of species with the ith exon
		tempMatrix[2] <- num.species           ### stats.data.exome[i,2] = num.species
		
		if (num.species>1){
			alignment           <- rMSA::mafft(temp.seqs,param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6")  ### aligns the ith exon of each species
			rownames(alignment) <- names(temp.seqs)  ### needed because the mafft function truncates sequence names
			alignment.width     <- ncol(alignment)   ### number of sites in alignment
			distances           <- 100*as.matrix(dist.dna(as.DNAbin(alignment),model="raw",pairwise.deletion=T))     ### pairwise matrix of percent genetic distance
			pident              <- 100-distances                                                                     ### pairwise matrix of percent genetic identity
			mean.pident         <- mean(pident[1,]) ### stats.data.exome[i,6] = mean(pident[1,])                     ### calculates mean percent identity from primary species for ith exon
			tempMatrix[6]       <- mean.pident
			tempMatrix[7] = pident[1,1]                                                                              ### percent identity of primary species to itself...this better be 100
			
			for(y in 2:length(species)){
				if(length(grep(pattern=species[y],rownames(pident)))!=0) {
					tempMatrix[y+6] = pident[1,grep(pattern=species[y],rownames(pident))]
				} else {
					tempMatrix[y+6] = as.numeric(-1)
				}
			}		
			
		}
		if (num.species > 3){
			count.cover <- 0                                                                       #| This counts the number of sites
			for(w in 1:ncol(alignment)){                                                           #| of exon i with at least
				temp.cover        <- subseq(DNAStringSet(alignment),start=w,end=w)                 #| 4 individuals with non-missing data
				samples.at.site.w <- length(temp.cover)-length(grep(pattern="-",temp.cover))       #|
				if(samples.at.site.w>3){                                                           #|
					count.cover<-count.cover+1                                                     #|
				}                                                                                  #|
			}	                                                                                   #|
			count.cover <- as.numeric(count.cover)                                                 #|

			if(count.cover == 0){                                                   #|If no sites have ≥4 individuals with non-missing data,
				absolutePIS <- 0                                                    #|then the number of parsimony informative sites is zero for exon i.
			}  else {                                                               #|Otherwise, calculates the number of parsimony informative
				absolutePIS <- ips::pis(as.DNAbin(alignment),what="absolute")       #|sites (using pis function of ips package) for exon i
			}                                                                       #|

			if(absolutePIS == 0){                           #|If no parsinomy informative sites, then percent of sites that are parsimony informative
				percentPIS <- 0                             #|is set to zero, otherwise divides number of parsimony informative sites
			} else {                                        #|by the total number of sites to get the percent of sites that are partsimony
				percentPIS <-  absolutePIS/count.cover      #|informative for exon i
			}                                               #|
			
			tempMatrix[3]  = count.cover ### number of characters at locus with at least 4 individuals represented
			tempMatrix[4]  = absolutePIS ### number of parsimony informative sites at locus
			tempMatrix[5]  = percentPIS  ### percent of characters with at least 4 individuals represented that are parsimony informative
			tempMatrix[18] = alignment.width
		}
		### Stats added from the first part of Step 5
		mean.var.sites      <- ((100-as.numeric(mean.pident))/100)*as.numeric(CDS.lengths[i])  ### Calculates the average number of variable sites; this value isnt used later
		min.pident.all      <- min(as.numeric(tempMatrix[7:(length(species)+6)]))              ### minimum pident to exon of primary species (i.e., minimum pident to T. sirtalis for SnakeCap)
		if(all(is.na(subgroup))){
			min.pident.subgroup <- "-1"
		} else {
			min.pident.subgroup <- min(as.numeric(tempMatrix[(is.subgroup+6)]))                      ### minimum pident to exon of primary species among a subset of the species (for SnakeCap, subgroup = snakes)
		}
		
		tempMatrix[20] <- CDS.lengths[i]
		tempMatrix[21] <- mean.var.sites
		tempMatrix[22] <- min.pident.all
		tempMatrix[23] <- min.pident.subgroup
		
		if(i==1){
			colnames.part1  <- paste(species[is.primary.exome],".CDS.names",sep="")                      #|Determines what the column names of
			colnames.part2  <- c("num.Species","CountCover","absolutePIS","percentPIS","mean.pident")    #|the stats matrix should be
			colnames.part3  <- paste("pident.",species,sep="")                                           #|
			colnames.part4 <- c("alignment.width","gene.name",paste("CDS.length.",species[1],sep=""),"mean.variable.sites","min.pident.all","min.pident.subgroup")   #|
			colnames(tempMatrix) <- c(colnames.part1,colnames.part2,colnames.part3,colnames.part4)       # sets the column names
			write.table(tempMatrix,file=paste(output.dir,"stats_exome_data_TBLASTX.txt",sep=""),sep=",",col.names=T,row.names=F,append=F) ### writes the column names and tempMatrix for exon i = 1
		}
		if(i!=1){
			write.table(tempMatrix,file=paste(output.dir,"stats_exome_data_TBLASTX.txt",sep=""),sep=",",col.names=F,row.names=F,append=T) ### writes tempMatrix for exon i > 1
		}
	}
}


######################
## pick.loci function
######################
## parameters
####
### min.pident.keep         <- c(65,100)  ### minimum and maximum values for min.pident necessary to keep locus
### max.capture.coverage    <- 1200000    ### the maximum number of nucleotides that you expect to be able to capture given constrains of your particular probe Kit
### write.stats.tables      <- FALSE      ### You will need to change this to TRUE; I have set to FALSE to avoid overwriting my earlier results
### plot.results            <- FALSE      ### Change to TRUE if you want to plot your result; haven't plotted this in a while so I'm not sure if this works
### statsTable.path         <- "stats_exome_data_TBLASTX.txt"
### output.dir              <- input.dir
### primary.species         <- "Thamnophis_sirtalis"
### species.subgroup        <- c(7:14) ### either a character vector of species names that are used in column names, or an integer vector specifying which columns are for species included in the subgroup
### use.min.pident.subgroup <- T
### fast.stat               <- "pident"   ### alternatively, "percentPIS". Whether to use mean percent identity or mean percent of sites parsimony informative when sorting loci before choosing those those with summed length < max.capture.coverage
####### Next like runs the function under the values that were used for SnakeCap
### result.pident <- pick.loci(statsTable.path = "stats_exome_data_TBLASTX.txt", output.dir = NULL, primary.species = "Thamnophis_sirtalis", species.subgroup = c(7:14), use.min.pident.subgroup = T, min.pident.keep = c(65,100), max.capture.coverage = 1200000, write.stats.tables = T, plot.results = T, fast.stat = "pident")
pick.loci <- function(statsTable.path, output.dir, primary.species,species.subgroup=NULL, min.pident.keep=c(65,100), max.capture.coverage=1200000, write.stats.tables=F, plot.results=T,use.min.pident.subgroup=F,fast.stat="pident"){
	stats.data.exome          <- fread(input=statsTable.path,sep=",",header=T)                                           ### reads in the full stats table generated in step 4
	stats.data.exome.ordered  <- stats.data.exome[with(stats.data.exome, order(gene.name, mean.pident, decreasing=F)),]  ### Sorts loci by gene name, then by increasing percent identity to Thamnophis sirtalis
	
	###### Eventually you can comment out these next two if statements #####
	all.pident.columns           <- grep("pident.",colnames(stats.data.exome))[which(grep("pident.",colnames(stats.data.exome)) < grep("gene.name",colnames(stats.data.exome)))]
	if(length(which(colnames(stats.data.exome.ordered)=="min.pident.all"))==0){
		min.pident.all           <- apply(X=stats.data.exome.ordered[,..all.pident.columns],MARGIN=1,FUN=min) #| [NOTE: Must have the ".." before the object "all.pident.columns" to specify that it is an integer index.] This line Calculates minimum percent identity (among all species to primary species)
		stats.data.exome.ordered <- cbind(stats.data.exome.ordered,min.pident.all)            #| Adds the minimum percent identity (all species) stat as a column if it doesnt already exist
	}
	if(length(which(colnames(stats.data.exome.ordered)=="min.pident.subgroup"))==0 & !is.null(species.subgroup)){
		if(is.integer(species.subgroup)){
			subgroup.pident.columns  <- species.subgroup
		}
		if(is.character(species.subgroup)){
			expected.colnames <- paste("pident.",species.subgroup,sep="")
			subgroup.pident.columns  <- match(expected.colnames,colnames(stats.data.exome))
		}
		min.pident.subgroup      <- apply(X=stats.data.exome.ordered[,..subgroup.pident.columns],MARGIN=1,FUN=min)    #| Calculates minimum percent identity (among subgroup species to primary species)
		stats.data.exome.ordered <- cbind(stats.data.exome.ordered,min.pident.subgroup)                               #| Adds the minimum percent identity (subgroup species) stat as a column if it doesnt already exist
	}
	
	if(use.min.pident.subgroup==T){
		min.pident.filter <- stats.data.exome.ordered$min.pident.subgroup
	} else {
		min.pident.filter <- stats.data.exome.ordered$min.pident.all
	}
	to.keep1                      <- which(min.pident.filter >= min.pident.keep[1] & min.pident.filter < min.pident.keep[2])  ####| Filters out loci if minimum percent identity to Thamnophis sirtalis is < 65%, or = 100%
	stats.data.exome.ordered      <- stats.data.exome.ordered[to.keep1,]                                                                                                        ####|
	
	to.keep2                      <- match(unique(stats.data.exome.ordered$gene.name),stats.data.exome.ordered$gene.name)  #### List of rows containing the fastest-evolving exon per gene
	stats.data.FastestExonPerGene <- stats.data.exome.ordered[to.keep2,]                                                   #### Creates a table containing only the fastests-evolving exons per gene
	
	if(write.stats.tables == T){
		write.table(stats.data.FastestExonPerGene,file=paste(output.dir,"stats_data_FastestExonPerGene.txt",sep=""),sep=",",col.names=T,row.names=F,append=F)
	}
	
	if(fast.stat == "pident"){
		mean.pident.fast                      <- stats.data.FastestExonPerGene$mean.pident                             ###| Orders exon stats table by increasing mean percent identity to T. sirtalis
		stats.data.FastestExonPerGene.ordered <- stats.data.FastestExonPerGene[order(mean.pident.fast),]               ###|
	}
	
	if(fast.stat == "percentPIS"){
		mean.percentPIS                      <- stats.data.FastestExonPerGene$percentPIS                               ###| Orders exon stats table by increasing mean percent identity to T. sirtalis
		stats.data.FastestExonPerGene.ordered <- stats.data.FastestExonPerGene[order(percentPIS),]                     ###|
	}
	
	###### Eventually you can comment out this next if statement ######
	if(length(which(colnames(stats.data.FastestExonPerGene.ordered)=="CDS.lengths"))==0){
		CDS.ranges  <- gsub(".*\\.1_","",stats.data.FastestExonPerGene.ordered$Thamnophis.CDS.names)          ###| Only need these to lines if
		CDS.lengths <- (abs(as.numeric(gsub(".*_","",CDS.ranges))-as.numeric(gsub("_.*","",CDS.ranges)))+1)   ###| exon stats table doesnt already
		stats.data.FastestExonPerGene.ordered <- cbind(stats.data.FastestExonPerGene.ordered,CDS.lengths)     ###| include a column for CDS.lengths
	}
	
	roll.sum                           <- rollSum(stats.data.FastestExonPerGene.ordered$CDS.lengths)                       ### Calculates the rolling sum of exon lengths
	stats.data.FastestExonPerGene.best <- stats.data.FastestExonPerGene.ordered[which(roll.sum <= max.capture.coverage),]  ### Keeps the loci in the first N rows of the table, which have roll.sum <= max.capture.coverage
	
	if(write.stats.tables == T){
		write.table(stats.data.FastestExonPerGene.best,file=paste(output.dir,"stats_data_FastestExonPerGene_best.txt",sep=""),sep=",",col.names=T,row.names=F,append=F)
	}
	
	mean.pident.hist               <- hist(stats.data.exome$mean.pident,breaks=round(max(stats.data.exome$mean.pident)-min(stats.data.exome$mean.pident)),plot=F)                                      ###| Histogram data to be plotted later
	mean.pident.hist.keep1         <- hist(stats.data.exome.ordered$mean.pident,breaks=round(max(stats.data.exome.ordered$mean.pident)-min(stats.data.exome.ordered$mean.pident)),plot=F)                      ###| After filtering to keep only loci within the range of min.pident.keep
	mean.pident.hist.keep2         <- hist(stats.data.FastestExonPerGene$mean.pident,breaks=round(max(stats.data.FastestExonPerGene$mean.pident)-min(stats.data.FastestExonPerGene$mean.pident)),plot=F)            ###| After filtering to keep only fastest exon per gene
	mean.pident.hist.best          <- hist(stats.data.FastestExonPerGene.best$mean.pident,breaks=round(max(stats.data.FastestExonPerGene.best$mean.pident)-min(stats.data.FastestExonPerGene.best$mean.pident)),plot=F)  ###| After filtering to keep only the fastest exons with concatenated length < max.capture.coverage
	
	#mean.pident.dens               <- density(stats.data.exome$mean.pident)                    ###| Density data to be plotted later
	#mean.pident.dens.keep1         <- density(stats.data.exome.ordered$mean.pident)            ###| After filtering to keep only loci within the range of min.pident.keep
	#mean.pident.dens.keep2         <- density(stats.data.FastestExonPerGene$mean.pident)       ###| After filtering to keep only fastest exon per gene
	#mean.pident.dens.best          <- density(stats.data.FastestExonPerGene.best$mean.pident)  ###| After filtering to keep only the fastest exons with concatenated length < max.capture.coverage
		
	if (plot.results==T){
		par(mfrow=c(4,1), mar=c(1,1,1,1),oma=c(4,4,0,0))  ### 4 rows and 1 column of plots
		cexSize   <- 1              # text size
		x.lab2   <- "% identical"   # xlabel 2
		y.lab2   <- "# loci"        # y label 2
		x.lab    <- NULL            #| x label, y label, and plot title
		y.lab    <- NULL            #|
		main.lab <- NULL            #|
		
		x.lim.min <- mround(x=min(stats.data.exome$mean.pident),base=5,direction="down")  ### determines the mininum x-axis value
		
		plot(mean.pident.hist, xlim=c(x.lim.min,100), xlab=NULL, ylab=y.lab2, main = main.lab,cex.main=cexSize,col="grey")  ### Histogram
		Axis(side=2)
		text(x=c(((par("usr")[2]-par("usr")[1])/10)+par("usr")[1]),y= c(((par("usr")[4]-par("usr")[3])/2)+par("usr")[3]),labels=c("Exon Set 1 (all annotated exons)"),cex=1,col=1,adj=c(0,NA))
		lines(x=rep(mean(stats.data.exome$mean.pident),2),y=c(0,par("usr")[4]),col="green")
		
		plot(mean.pident.hist.keep1, xlim=c(x.lim.min,100), xlab=NULL, ylab=y.lab2, main = main.lab,cex.main=cexSize,col="grey")
		Axis(side=2)
		text(x=c(((par("usr")[2]-par("usr")[1])/10)+par("usr")[1]),y= c(((par("usr")[4]-par("usr")[3])/2)+par("usr")[3]),labels= c("Exon Set 2: subset of Exon Set 1 in which:"),cex=1,col=1,adj=c(0,NA))
		text(x=c(((par("usr")[2]-par("usr")[1])/10)+par("usr")[1]),y= c(((par("usr")[4]-par("usr")[3])/2.5)+par("usr")[3]),labels= paste("100% > minimum similarity to ",gsub("_"," ",primary.species)," ≥ ",min.pident.keep[1],"%",sep=""),cex=1,col=1,adj=c(0,NA))
		lines(x=rep(mean(stats.data.exome.ordered$mean.pident),2),y=c(0,par("usr")[4]),col="green")
		
		plot(mean.pident.hist.keep2, xlim=c(x.lim.min,100), xlab=NULL, ylab=y.lab2, main = main.lab,cex.main=cexSize,col="grey")
		Axis(side=2)
		text(x=c(((par("usr")[2]-par("usr")[1])/10)+par("usr")[1]),y= c(((par("usr")[4]-par("usr")[3])/2)+par("usr")[3]),labels=c("Exon Set 3: fastest exon/gene of Exon Set 2"),cex=1,col=1,adj=c(0,NA))
		lines(x=rep(mean(stats.data.FastestExonPerGene$mean.pident),2),y=c(0,par("usr")[4]),col="green")
		
		plot(mean.pident.hist.best, xlim=c(x.lim.min,100), xlab=NULL, ylab=y.lab2, main = main.lab,cex.main=cexSize,col="grey")
		Axis(side=2)
		text(x=c(((par("usr")[2]-par("usr")[1])/10)+par("usr")[1]),y= c(((par("usr")[4]-par("usr")[3])/2)+par("usr")[3]),labels= c("Exon Set 4: fastest exons of Exon Set 3 in which:"),cex=1,col=1,adj=c(0,NA))
		text(x=c(((par("usr")[2]-par("usr")[1])/10)+par("usr")[1]),y= c(((par("usr")[4]-par("usr")[3])/2.5)+par("usr")[3]),labels= paste("sum(exon lengths) < ",max.capture.coverage,"nt",sep=""),cex=1,col=1,adj=c(0,NA))
		lines(x=rep(mean(stats.data.FastestExonPerGene.best$mean.pident),2),y=c(0,par("usr")[4]),col="green")
		
		mtext(paste("% of sites identical to ", gsub("_"," ",primary.species)),c(SOUTH=1),line=2, cex=1, outer=TRUE)
		mtext(c("number of exons"),c(WEST=2),line=1.5, cex=1, outer=TRUE)
	}
	result.table.names <- c("Table1. Unique CDS regions above min length","Table2. Subset of Table1 with pident above cuttoff","Table3. Subset of Table2 with only fastest exons per gene","Table4. Fastest exon subset of Table3 with summed length below 120Mbp")
	result <- list(stats.data.exome,stats.data.exome.ordered,stats.data.FastestExonPerGene,stats.data.FastestExonPerGene.best,result.table.names)
	result
}

######################
## align.and.concatenate.best.exons function
######################
## parameters
### exomes.filepaths     # paramA <- list.files(path="~/Exomes/",full.names=T)
### is.primary.exome     # paramB <- 1  ### indicates which of the species names correspond to the "subject" species during the tblastn step. This currently only works if set to 1.
### statsTable.path      # paramC <- "stats_data_FastestExonPerGene_best.txt"
### output.dir           # paramD <- "Exomes_TempFolder_3Nov2019/" ### location to write Exome Stats Table, currently set as folder for testing
### species              # paramE <- c("Thamnophis_sirtalis","Ophiophagus_hannah","Crotalus_mitchellii","Python_bivittatus","Vipera_berus","Crotalus_horridus","Protobothrops_mucrosquamatus","Pantherophis_guttatus","Anolis_carolinensis","Pogona_vitticeps","Gekko_japonicus")
### i.start <- 1
### i.stop  <- NA
align.and.concatenate.best.exons <- function(exomes.filepaths,is.primary.exome,statsTable.path,output.dir,species,i.start = 1,i.stop = NA){
	statsTable           <- fread(statsTable.path,sep=",",header=T)                         ### reads in the info on fastest exon per gene dataset with total length <1.2Mbp
	species              <- c(species[is.primary.exome],species[-is.primary.exome])         ### reorders the list of species such that the primary.species is first in the list
	exomes.filepaths     <- exomes.filepaths[mgrep(species,exomes.filepaths)]               ### puts the exomes.filepaths in the same order as species
	exome.names          <- paste("exome",c(1:length(exomes.filepaths)),sep="")
	#matches.names        <- paste("matches.exome",c(1:length(exomes.filepaths)),sep="")
	for(i in 1:length(exome.names)){                                                       #|reads in exomes
		assign(x=exome.names[i],value=readDNAStringSet(filepath=exomes.filepaths[i]))      #|and assigns them
	}                                                                                      #|to object names

	CDS.names        <- names(exome1)
	unique.CDS.names <- unique(CDS.names)
	to.keep          <- match(unique.CDS.names, CDS.names)
	CDS.names        <- unique.CDS.names
	exome1           <- exome1[to.keep]

	locus.names <- paste(unlist(statsTable[,1]),unlist(statsTable$gene.name),sep="_")
	CDS.names   <- unlist(statsTable[,1])

	if(is.na(i.start) | i.start > length(CDS.names)){
		i.start <- 1                  ### default i.start, ie which exons to start the loop at
	}
	if(is.na(i.stop) | i.stop > length(CDS.names)){
		i.stop  <- length(CDS.names)   ### default i.stop, ie which exons to stop the loop at
	}
	
	dir.check.create(paste(output.dir,"exonAlignments",sep=""))
	
	for(i in i.start:i.stop){
	
		temp.exon.names <- paste("temp.exon.exome",c(1:length(species)),sep="")
		
		for(k in 1:length(temp.exon.names)){
			assign(x=temp.exon.names[k],value = get(exome.names[k])[grep(CDS.names[i],names(get(exome.names[k])))])
		}
		
		temp.seqs <- DNAStringSet(get(temp.exon.names[1]))
		for(z in 2:length(temp.exon.names)){
			temp.seqs <- c(temp.seqs,DNAStringSet(get(temp.exon.names[z])))
		}
		alignment           <- rMSA::mafft(temp.seqs,param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6")  ### aligns the ith exon of each species
		rownames(alignment) <- names(temp.seqs)  ### needed because the mafft function truncates sequence names
		
		writeXStringSet(x=alignment, filepath=paste(output.dir,"exonAlignments/",locus.names[i],".fas",sep=""), append = F)
		
		if(i == 1){
			alignment2 <- alignment          ### concatenated alignment for i = 1
			start.temp <- 1                  ### site where ith partition starts for i = 1
			end.temp <- width(alignment)[1]  ### site where ith partition ends for i = 1
			app=FALSE                        ### states that we do not want to the append partition info to next line of partition file for the first loop
		}
		if(i != 1){
			start.temp <- width(alignment2)[1]+1      ### width of the (i-1)th alignment2 + 1
			alignment2 <- xscat(alignment,alignment2) ### concatenates the ith alignment with the (i-1)th alignment2
			end.temp   <- width(alignment2)[1]        ### width of the ith alignment2
			app=TRUE                                  ### states that we will want to the append partition info to next line of partition file
		}
		temp.write <- paste("DNA, ",locus.names[i]," = ",start.temp,"-",end.temp,sep="")
		temp.write <- gsub("Subject=","Subject_",temp.write)
		temp.write <- gsub("gene=","gene_",temp.write)
		write(temp.write, paste(output.dir,"fastestExonPerGene_BestExons_partitionFile.txt",sep="") ,append=app)
		if(i%%10==0){print(paste(i,"of",length(CDS.names),"complete",Sys.time(),sep=" "))}
	}
	names(alignment2) <- species
	writeXStringSet(alignment2, paste(output.dir,"FastestExonPerGene_BestExons_aligned_concatenated.txt",sep=""), append = F) #writes the aligned and concatenated target loci to file
}

##########
### makeBlastDB function
#########
### R wrapper for the ncbi makeblastdb command
### This function makes an ncbi blast database (which is needed to run blast against)
## parameters
###
## makeblastdb.path <- "~/ncbi-blast-2.5.0+/bin/makeblastdb"                     ### path to the script called "makeblastdb"
## subject.path     <- "~/GCA_000737285.1_CrotMitch1.0_genomic.fna"              ### path to the subject sequences that will be made into a database
makeBlastDB    <- function(makeblastdb.path,subject.path){
	command.part1     <- makeblastdb.path              #|character strings that will be 
	command.part2     <- "-in"                         #|pasted together into a "command" string
	command.part3     <- subject.path                  #|that can be executed by calling terminal
	command.part4     <- "-parse_seqids -dbtype nucl"  #|
	command           <- paste(command.part1,command.part2,command.part3,command.part4)  ### pastes the parts into a character string that can be executed in terminal
	system(command)   ### calls terminal to execute the character string "command"
}

##########
### blastnR function
#########
### R wrapper for the ncbi blastn command
### This function uses a set of query sequences and blasts them against a set of subject sequences
### and produces a hit table with up to 50 hits per query
## parameters
### blastn.path   <- "~/ncbi-blast-2.5.0+/bin/blastn"
### subject.path  <- "GCA_000737285.1_CrotMitch1.0_genomic.fna"
### query.path    <- "micrurus_UCEs.fa"
### output.path   <- "Crotalus-mitchellii.blastn.exons.50hits.txt"
# These next three lines are the full set of subject paths and the query path used for UCEs for SnakeCap
### subject.path <- genome.subjects.paths  <- c("~/GCA_000737285.1_CrotMitch1.0_genomic.fna","~/Protobothrops_mucrosquamatus_GCF_001527695.2_P.Mucros_1.0_genomic.fna","~/Pantherophus_guttatus_GCA_001185365.1_PanGut1.0_genomic.fna","~/Python_bivittatus_GCF_000186305.1_Python_molurus_bivittatus-5.0.2_genomic.fna","~/Thamnophis_sirtalis_GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.fna","~/ContigFiles/Ophiophagus_hannah_GCA_000516915.1_OphHan1.0_genomic.fna","~/Vipera_berus_GCA_000800605.1_Vber.be_1.0_genomic.fna","~/Crotalus_horridus_GCA_001625485.1_ASM162548v1_genomic.fna")
### query.path   <- micrurus.UCEs.path     <- "~/StreicherWiens2017_UCEs_fasta/micrurus_UCEs.fa"
### output.path  <- output.50hits.tables.paths <- c("~/Crotalus-mitchellii.blastn.exons.50hits.txt","~/Protobothrops.blastn.exons.50hits.txt","~/Pantherophis.blastn.exons.50hits.txt","~/Python.blastn.exons.50hits.txt","~/Thamnophis.blastn.exons.50hits.txt","~/Ophiophagus.blastn.exons.50hits.txt","~/Vipera.blastn.exons.50hits.txt","~/Crotalus-horridus.blastn.exons.50hits.txt")
blastnR <- function(blastn.path,subject.path,query.path,output.path){
	DB.extensions        <- c(".nog",".nsd",".nsi",".nsq",".nhr",".nin",".fai",".gz")  ### double-check that these are the files produced by makeBlastDB
	makeblastdb.filepath <- gsub("/blastn","/makeblastdb",blastn.path)                 ### unless these files have been moved around
	subject.db.files     <- paste(subject.path,DB.extensions,sep="")                   ### a list of files that should be present if local NCBI database exists (has been made) for the subject sequences
	
	if(!all(file.exists(subject.db.files))){             #| Makes NCBI database if one doesnt already exist
		makeBlastDB(makeblastdb.filepath,subject.path)   #| 
	}                                                    #|
	
	command.part1 <- blastn.path
	command.part2 <- "-db"
	command.part3 <- subject.path
	command.part4 <- "-query"
	command.part5 <- "-out"
	command.part6 <- output.path
	command.part7 <- "-evalue 1e-5 -outfmt 6 -max_target_seqs 50 -num_threads 16"
	command       <- paste(command.part1,command.part2,command.part3,command.part4,command.part5,command.part6,command.part7)
	system(command)   ### calls terminal to execute the character string "command"
}

######################
## align.bestHit.UCEs function
######################
## parameters
### species.UCEs.filepaths   # paramA <- list.files(path="~/UCEs.In.Snake.Genomes/",full.names=T)
### is.primary.species       # paramB <- 1                                                                                ### indicates which of the species names correspond to the "subject" species during the tblastn step. This currently only works if set to 1.
### output.dir               # paramD <- "~/UCEs_TempFolder_7Nov2019" ### location to write Exome Stats Table, currently set as folder for testing
### species                  # paramE <- c("Thamnophis_sirtalis","Ophiophagus_hannah","Crotalus_mitchellii","Python_bivittatus","Vipera_berus","Crotalus_horridus","Protobothrops_mucrosquamatus","Pantherophis_guttatus")
### i.start <- 1
### i.stop  <- NA
align.bestHit.UCEs <- function(species.UCEs.filepaths,is.primary.species,output.dir,species,i.start = 1,i.stop = NA){
	species                <- c(species[is.primary.species],species[-is.primary.species])         ### reorders the list of species such that the primary.species is first in the list
	species.UCEs.filepaths <- species.UCEs.filepaths[mgrep(species,species.UCEs.filepaths)]       ### puts the species.UCEs.filepaths in the same order as species
	UCEs.species           <- paste("UCEs",c(1:length(species.UCEs.filepaths)),sep="")
	UCE.names              <- paste("UCE.names",c(1:length(species.UCEs.filepaths)),sep="")
	
	for(i in 1:length(UCEs.species)){                                                        #|reads in exomes
		assign(x=UCEs.species[i],value=readDNAStringSet(filepath=species.UCEs.filepaths[i])) #|and assigns them
	}                                                                                        #|to object names
	
	for(i in 1:length(UCEs.species)){                                                        
		assign(x=UCE.names[i],value=  gsub(".*_uce","UCE",names(get(UCEs.species[i])))) #| gets set of UCE names for each species and assigned the set to an object name
	}
	
	UCEs.shared <- intersect.all(lapply(UCE.names,get))  ### list of UCEs found in all genomes
	
	UCEs.species
	test <- lapply(x=lapply(UCE.names,get),match,table=UCEs.shared)
	
	UCEs.species.shared <- paste(UCEs.species,".shared",sep="")

	for(i in 1:length(species)){
		temp.UCEs <- get(UCE.names[i])
		matches <- match(UCEs.shared,temp.UCEs)
		assign(UCEs.species.shared[i],get(UCEs.species[i])[matches])
	}
	
	geti <- function(x,index=i){  #| Gets the ith element of an object with name "x"
		get(x)[index]             #| This is a specifial version of get() that is
	}                             #| sometimes useful when combined with lapply() (only useful in i loops though)
	
	
	if(is.na(i.stop)){
		i.stop  <- length(UCEs.shared)
	}
	
	for(i in i.start:i.stop){
		temp.seqs  <- lapply(UCEs.species.shared,geti)
		
		temp.seqs2 <- DNAStringSet(temp.seqs[[1]])
		for(z in 2:length(temp.seqs)){
			temp.seqs2 <- c(temp.seqs2,DNAStringSet(temp.seqs[[z]]))
		}
		alignment           <- rMSA::mafft(temp.seqs2,param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6")  ### aligns the ith exon of each species
		rownames(alignment) <- names(temp.seqs2)                  ### needed because the mafft function truncates sequence names
		alignment           <- DNAStringSet(alignment)
		writeXStringSet(alignment, paste(output.dir,UCEs.shared[i],".fasta",sep=""), append = F)
	}
}

###### xprod.combn function
### input: list.of.vects = a list of vectors (aka, sets)
### value: a character vector of each possible combinations that can be created by drawing one sample from each input set,
###        each element has the form "set1sample set2sample". Can then use strsplit(value[i],split=" ") to get each seperately
###        
xprod.combn <- function(list.of.vects) {
	nloops  <- length(list.of.vects)-1
	
	for(i in 1:nloops){
		if(i==1){
			res.temp <- c(list.of.vects[[1]])
		}
		res.temp <- c(outer(res.temp, c(list.of.vects[[i+1]]), FUN=paste))
	}
	res.temp
}

###### xprod.combn.mat function
### input: list.of.vects = a list of vectors (aka, sets; numerical or character)
###### NOTE: the previous version of this function required that the input be a list of vectors, whereas this version takes any number of vectors as input and creates the list
#### value: a two column matrix, with each row containing a unique combination from from set1 (column 1 element) and set2 (column 2 element)
###
xprod.combn.mat <- function(...) {
	list.of.vects <- list(...)
	nloops        <- (length(list.of.vects)-1)
	for(i in 1:nloops){
		if(i==1){
			res.temp <- c(list.of.vects[[1]])
		}
		res.temp <- c(outer(res.temp, c(list.of.vects[[i+1]]), FUN=paste))
	}
	res.temp <- strsplit(res.temp,split=" ")
	res.mat <- matrix(unlist(res.temp), ncol=2, byrow=T)
	res.mat
}

###########
## run.mafft function
###########
run.mafft<-function(unaligned.contigs, add.contigs = NULL, algorithm = "localpair", rev.dir = T, save.name = NULL, threads = 6, delete.files = T){
  save.contigs<-as.list(as.character(unaligned.contigs))
  if (is.null(save.name) == T) { save.name<-paste(sample(LETTERS, 5, replace = T), collapse = "")}
  if (rev.dir == T){
  		adjust.direction <- "--adjustdirection"
  	} else {
  		adjust.direction <- ""
  	}
  
  #Adds a sequence into the alignment. Saves much computation.
  if (algorithm == "add"){
    #Saves to folder to run with mafft
    write.fasta(sequences = save.contigs, names = names(save.contigs), paste(save.name, ".fa", sep = ""), nbchar = 1000000, as.string = T)
    
    #Saves to folder to run with mafft
    add.save<-as.list(as.character(add.contigs))
    write.fasta(sequences = add.save, names = names(add.save), "add_sequences.fa", nbchar = 1000000, as.string = T)
    
    #Runs MAFFT to align
    system(paste("mafft --",algorithm, " add_sequences.fa ", adjust.direction, " --maxiterate 1000 ", save.name, ".fa > ", save.name, "_align.fa", sep = ""), ignore.stderr = T)
    
    alignment<-scanFa(FaFile(paste(save.name, "_align.fa", sep = "")))   # loads up fasta file
    unlink(paste(save.name, ".fa", sep = ""))
    unlink("add_sequences.fa")    
  }#end -add
  
  #Does Regular MAFFT Local Pair
  if (algorithm == "localpair"){
    #Saves to folder to run with mafft
    write.fasta(sequences = save.contigs, names = names(save.contigs), paste(save.name, ".fa", sep = ""), nbchar = 1000000, as.string = T)
    
    #Runs MAFFT to align
    system(paste("mafft --",algorithm, " --maxiterate 1000 ", adjust.direction, " --quiet --op 3 --ep 0.123"," --thread ", threads, " ", save.name, ".fa > ", save.name, "_align.fa", sep = ""))
    
    alignment <- scanFa(FaFile(paste(save.name, "_align.fa", sep = "")))   # loads up fasta file
    unlink(paste(save.name, ".fa", sep = ""))
  }#end local pair
  
  if (delete.files == T){
    unlink(paste(save.name, "_align.fa", sep = ""))
    return(alignment)
  } else { return(alignment) }
}#function end

### pairwise.inf.sites function
#####
### 
pairwise.inf.sites <- function(x, y) {
  temp.align <- strsplit(as.character(x), "")
  mat.align  <- lapply(temp.align, tolower)
  m.align    <- as.matrix(as.DNAbin(mat.align))
  
  #Filters out weirdly divergent sequences
  new.align  <- as.character(m.align)
  new.align[new.align == "n"]<-"-"
  new.align[is.na(new.align) == T]<-"-"
  ref<-new.align[rownames(new.align) == y,]
  summary.data<-c()
  all.pars<-c()
  all.over<-c()
  for (z in 1:nrow(new.align)) {
    #Site counter
    pars         <- 0
    overlap      <- 0
    tar          <- new.align[z,]
    combined     <- matrix(NA_character_, ncol = max(length(ref), length(tar)), nrow =2)
    combined[1,] <- ref
    combined[2,] <- tar
    for (k in 1:ncol(combined)) {
      #Pulls out column of data
      seq.col<-vector("character", length = nrow(combined))
      seq.col<-combined[,k]
      #not equal to -
      f.char<-seq.col[seq.col != '-'] 
      #don't count missing seq
      if (length(f.char) <= 1) {
      	next
      }
      if (length(f.char) >= 2){
        overlap <- overlap+1
        if (f.char[1] != f.char [2]) {
        	pars<-pars+1
        }
      }#end if
    }#ends informative sites loop
    all.pars<-append(all.pars, pars)
    all.over<-append(all.over, overlap)
  }# ends seq loop
  #Summarizes and returns data
  summary.data<-all.pars/all.over
  summary.data[is.nan(summary.data)]<-0
  names(summary.data)<-rownames(new.align)
  return(summary.data)
}


### write.phy function
write.phy <- function (x, file = "", interleave = FALSE, strict = FALSE){
  str2cha <- function(x) {
    unlist(strsplit(x, ""))
  }
  datatype <- ifelse(is.numeric(x[1, 1]), "continuous", "nc")
  ntax <- nrow(x)
  nchar <- ncol(x)
  taxnames <- rownames(x)
  if (strict) {
    taxnames <- substring(taxnames, 1, truncate)
    missing <- 10 - unlist(lapply(strsplit(taxnames, ""), 
                                  length))
    for (i in seq(along = taxnames)) taxnames[i] <- paste(taxnames[i], 
                                                          paste(rep("*", missing[i]), collapse = ""), sep = "")
    if (any(duplicated(taxnames))) 
      cat("WARNING: Truncation of taxon names created", 
          "identical strings.")
  }
  else {
    xx <- nchar(taxnames)
    diff <- max(xx) - xx + 3
    for (i in 1:ntax) taxnames[i] <- paste(taxnames[i], paste(rep(" ", 
                                                                  diff[i]), collapse = ""), sep = "")
  }
  if (!interleave) 
    interleave <- nchar
  nbpart <- ceiling(nchar/interleave)
  pt <- matrix(nrow = nbpart, ncol = 2)
  pt[1, ] <- c(1, interleave)
  if (nbpart > 1) 
    for (i in 2:(dim(pt)[1])) {
      pt[i, ] <- c(pt[i - 1, 2] + 1, pt[i - 1, 2] + interleave)
      pt[nbpart, 2] <- nchar
    }
  phy <- paste(ntax, nchar)
  for (i in seq(along = pt[, 1])) {
    sm <- as.character(x[, pt[i, 1]:pt[i, 2]])
    if (is.null(dim(sm))) 
      sm <- as.matrix(sm, ncol = 1)
    sm <- apply(sm, 1, paste, collapse = "")
    if (i == 1) 
      sm <- paste(taxnames, sm)
    if (i < max(seq(along = pt[, 1]))) 
      sm <- c(sm, "")
    phy <- c(phy, sm)
  }
  if (file == "") {
    cat(phy, sep = "\n")
  }
  else {
    write(phy, file = file)
  }
}

##################
### charTest.dna function
### checks that an individual has at least one of each these nucleotide types A,C,G,T in its sequence
### ls = a list containing one or more character vectors
### this function is mostly useful after first using the uniqueLetters function
charTest.dna <- function(ls){
		ls <- tolower(ls)
		all(c("a","c","g","t") %in% ls)
}

##################
### charTest.aa function
### checks that an individual has at least one of each these amino acid types A,C,G,T in its sequence
### ls = a list containing one or more character vectors
charTest.aa <- function(ls){
		#ls <- tolower(ls)
		#all(c("a","c","g","t") %in% ls)
		any(GENETIC_CODE %in% ls)
}

##################
### trim.alignment function
### trims the ends of an alignment containing less than a threshold number of individuals (default = 1)
### alignment = DNA sequence alignment of class DNAbin, DNAMultipleAlignment, or DNAStringSet
######
### value = a trimmed alignment with the same class as the input (untrimmed) alignment, which may be DNAbin, DNAMultipleAlignment, or DNAStringSet
trim.alignment <- function (alignment, threshold=1){
	
	store.class <- class(alignment)
	
	if(store.class %in% c("DNAMultipleAlignment","DNAStringSet")){
		alignment <- as.DNAbin(DNAMultipleAlignment(alignment))
	}
	
	meets.threshold <- function(x) {
		x <- table(x)                                                     # This is the number of each type of character in a given column
		n <- c("-", "n", "b", "h", "d", "v", "k", "s", "r", "w","y","?")  # A list of possible characters to check against
		if (sum(x[!names(x) %in% n]) > threshold){                        #| These lines are specifying whether the 
			x <- TRUE                                                     #| the ith column has at least theshold characters that are not ambiguous or gaps
		} else {                                                          #|
			x <- FALSE                                                    #| the ith column has fewer than the theshold number of non-ambiguous/non-gap characters
		}                                                                 #|
	}
		
	string.list <- function(y){
		paste.collapse <- function(z){
			z <- paste(z,collapse="")
			z
		}
		y <- as.character(y)
		y <- apply(X=y,MARGIN=1,FUN=paste.collapse)
		y
	}
	
	x         <- as.character(alignment)
	logicals  <- unlist(apply(x, 2, meets.threshold))
	new.start <- min(which(logicals))
	new.end   <- max(which(logicals))
	out       <- alignment[,new.start:new.end]
	
	keep.rows  <- apply(X=out,MARGIN=1,FUN=charTest.dna)    ### A list of logicals indicating whether each individual has at least one non-missing character
	out        <- out[keep.rows,]                       ### Only keeps individual that have some data present
	
	if(store.class=="DNAStringSet"){
		out <- DNAStringSet(string.list(out))
	}
	
	if(store.class=="DNAMultipleAlignment"){
		out <- DNAMultipleAlignment(string.list(out))
	}
	
	out
}

###################################
### codonsDNAStringSet function ###
###################################
## x                     ### A DNAStringSet object (representing an alignment)
## which.codon.positions ### An integer in {1,2,3}, or an integer vector containing the values 1, 2, and/or 3.
## predict.frame = F     ### Whether or not to try to guess the reading frame, based on the relative number of snps in each frame.
##                       ### The third frame is assumed to be the frame with the most snps. This method works best for longer loci, faster-evolving loci,
##                       ### and when the alignment has many individuals, and/or moderate to high genetic distances among individuals.
## frameAtFirstchar=NULL ### the frame (i.e., codon position) of the first character of the sequence (if known). If unknown, leave as NULL, in which case the first base treated as if it is the first codon position.
codonsDNAStringSet <- function(x,which.codon.positions,predict.reading.frame=F,frameAtFirstchar=NULL){
	x.mat               <- do.call(rbind,strsplit(as.character(x),split=""))
	if(predict.reading.frame==T){
		x.frame1.mat    <- x.mat[,c(T,F,F)]
		x.frame2.mat    <- x.mat[,c(F,T,F)]
		x.frame3.mat    <- x.mat[,c(F,F,T)]
		x.frame1.set    <- DNAStringSet(apply(X=x.frame1.mat,MARGIN=1,FUN = function(y){paste(y,collapse="")}))
		x.frame2.set    <- DNAStringSet(apply(X=x.frame2.mat,MARGIN=1,FUN = function(y){paste(y,collapse="")}))
		x.frame3.set    <- DNAStringSet(apply(X=x.frame3.mat,MARGIN=1,FUN = function(y){paste(y,collapse="")}))
		snps.per.frame  <- c(width(snpsDNAStringSet(x.frame1.set)[1]),width(snpsDNAStringSet(x.frame2.set)[1]),width(snpsDNAStringSet(x.frame3.set)[1]))
		frame3          <- which(snps.per.frame==max(snps.per.frame))
		frame1          <- c(frame3-2,frame3+1)[c(frame3-2,frame3+1) %in% c(1:3)]
		frame2          <- c(frame3-1,frame3+2)[c(frame3-1,frame3+2) %in% c(1:3)]
		frames          <- c(frame1,frame2,frame3)
		codon.positions <- frames %in% which.codon.positions
	}
	if(is.null(frameAtFirstchar)==F){
		char1.frame     <- frameAtFirstchar
		char2.frame     <- c(char1.frame-2,char1.frame+1)[c(char1.frame-2,char1.frame+1) %in% c(1:3)]
		char3.frame     <- c(char2.frame-2,char2.frame+1)[c(char2.frame-2,char2.frame+1) %in% c(1:3)]
		frames          <- c(char1.frame,char2.frame,char3.frame)
		codon.positions <- frames %in% which.codon.positions
	}
	if(predict.reading.frame==F & is.null(frameAtFirstchar)==T) {
		codon.positions <- c(1:3) %in% which.codon.positions
	}
	positions.mat <- x.mat[,codon.positions]
	positions.set <- DNAStringSet(apply(X=positions.mat,MARGIN=1,FUN = function(y){paste(y,collapse="")}))
	positions.set
}

##################
### snp.alignment function
### removes non-SNP columns of a DNA alignment
### alignment = either DNA sequence alignment of class DNAbin, DNAMultipleAlignment, or DNAStringSet
###             or AA sequene alignment of class AAbin, AAMultipleAlignment, or AAStringSet
##################
### value = an alignment of snps with the same class as the input (untrimmed) alignment
### Remember that snps (herein) = sites with at least two alleles, each of which occurs at least twice
### Therefore, snps = parsimony informative sites...
snp.alignment <- function (alignment){
	######
	## Defining internal functions needed
	## for the whole function to work
	#### pars.inf.dna function ### tests if an alignment site is parsimony informative
	pars.inf.dna <- function(alignment.site) {	   # alignment.site is one column of a DNA alignment
		y <- table(alignment.site)                 # This is the number of each type of character in a given column
		n <- c("-", "n", "b", "h", "d", "v", "k", "s", "r", "w","y","?")  # A list of possible ambiguous or missing characters to ignore
		y <- y[!(names(y) %in% n)]
		if(any(y>1) & length(y)>1){
			y <- TRUE
		} else {
			y <- FALSE
		}
		y
	}	
	#### pars.inf.aa function
	pars.inf.aa <- function(alignment.site) {	   #|alignment.site is a column of an AA alignment
		y   <- table(alignment.site)               #|this is the number of each type of character in a given column
		n   <- c("-","X")                          #|a list of possible ambiguous or missing characters to ignore
		AAs <- c(unique(GENETIC_CODE))
		y   <- y[!(names(y) %in% n)]
		if(any(y>1) & length(y)>1){                #|tests if the column has at least two different characters,
			y <- TRUE                              #|and if at least on of the characters if shared by more than one individual
		} else {
			y <- FALSE
		}
		y
	}
	#### string.list function
	string.list <- function(y){
		paste.collapse <- function(z){
			z <- paste(z,collapse="")
			z
		}
		if(class(y)!="matrix"){
			y <- as.character(y)
		}
		y <- apply(X=y,MARGIN=1,FUN=paste.collapse)
		y
	}
	#### charTest.dna function ### tests if an individual has at least one of each nucleotide
	charTest.dna <- function(ls){
		ls <- tolower(ls)
		any(c("a","c","g","t") %in% ls)  ### changed "all" to "any"
	}
	#### charTest.aa function
	charTest.aa <- function(ls){
		#ls <- tolower(ls)
		#all(c("a","c","g","t") %in% ls)
		any(GENETIC_CODE %in% ls)
	}
	#### Start of actual function  ####
	store.class <- class(alignment)
	if(store.class %in% c("DNAMultipleAlignment","DNAStringSet")){
		alignment <- as.DNAbin(DNAMultipleAlignment(alignment))
		data.type <- "DNA"
	}
	if(store.class %in% c("AAMultipleAlignment","AAStringSet")){
		alignment <- as.AAbin(AAMultipleAlignment(alignment))
		data.type <- "AA"
	}
	x <- as.character(alignment)
	if(data.type == "DNA"){
		keep.cols  <- unlist(apply(x, MARGIN=2, pars.inf.dna))
		out        <- x[,keep.cols]
		keep.rows  <- apply(X=out,MARGIN=1,FUN=charTest.dna)    ### A list of logicals indicating whether each individual has at least one non-missing character

	}
	if(data.type == "AA"){
		keep.cols  <- unlist(apply(x, 2, pars.inf.aa))          ### 2 is the margin to apply the function over
		out        <- x[,keep.cols]
		keep.rows  <- apply(X=out,MARGIN=1,FUN=charTest.aa)     ### A list of logicals indicating whether each individual has at least one non-missing character
	}
	out        <- out[keep.rows,]                               ### Only keeps individuals that have some data present
	if(store.class=="DNAbin"){
		out <- as.DNAbin(out)
	}
	if(store.class=="DNAStringSet"){
		out <- DNAStringSet(string.list(out))
	}
	if(store.class=="DNAMultipleAlignment"){
		out <- DNAMultipleAlignment(string.list(out))
	}
	if(store.class=="AAbin"){
		out <- as.AAbin(out)
	}
	if(store.class=="AAStringSet"){
		out <- AAStringSet(string.list(out))
	}
	if(store.class=="AAMultipleAlignment"){
		out <- AAMultipleAlignment(string.list(out))
	}
	out
}

##################
### snpsDNAStringSet function
### removes non-SNP columns of a DNA alignment of class DNAStringSet
### Same as snp.alignment function, but more efficient although only works on DNAStringSet objects
#####
### alignment = DNA sequence alignment of class DNAStringSet
### remove.NoDataSeqs = T  #### filter the SNPs alignment to only include individuals with at least one A, C, G, or T 
##################
### value = an DNAStringSet alignment of snps
### Remember that snps = sites with at least two alleles, each of which occurs at least twice
### 
snpsDNAStringSet <- function(alignment,remove.NoDataSeqs = T){
	keep.cols <- rep(F,width(alignment[1]))
	for(i in 1:width(alignment[1])){
		chars.temp <- table(subseq(alignment,start=i,width=1))
		nucs.temp  <- chars.temp[names(chars.temp) %in% c("A","C","G","T") & chars.temp > 1]
		if(length(nucs.temp)>1){
			keep.cols[i] <- T
		}
	}
	alignment.mat       <- do.call(rbind,strsplit(as.character(alignment),split=""))
	alignment.mat.snps  <- alignment.mat[,keep.cols]
	alignment.list.snps <- apply(X=alignment.mat.snps,MARGIN=1,FUN=function(y){paste(y,collapse="")})
	alignment.set.snps.temp  <- DNAStringSet(alignment.list.snps)
	if(remove.NoDataSeqs == T){
		keep.rows <- which(countChars(alignment.set.snps.temp,"A|C|G|T")>0)
		alignment.set.snps <- alignment.set.snps.temp[keep.rows]
	} else {
		alignment.set.snps <- alignment.set.snps.temp
	}
	alignment.set.snps
}

#########
### find.orf function
#########
find.orf <- function(input.seq, codons = F, min.size = 80){
  #Sets up data
  # input.seq<-trimmed[j]
  codon.table <- data.frame(Start = rep(0,6), End = rep(0,6), Frame = c("F1", "F2", "F3", "R1", "R2", "R3"))
  for.seq     <- as.character(input.seq)
  
  #Gets codon stuff
  TAA <- matchPattern("TAA", for.seq)
  TGA <- matchPattern("TGA", for.seq)
  TAG <- matchPattern("TAG", for.seq)
  
  #Forward Frame 1
  result1 <- TAA[(TAA@ranges@start+2) %% 3 == 0]   
  result2 <- TGA[(TGA@ranges@start+2) %% 3 == 0]    
  result3 <- TAG[(TAG@ranges@start+2) %% 3 == 0]    
  
  starts <- c(result1@ranges@start, result2@ranges@start, result3@ranges@start)
  ends   <- c(result1@ranges@start+2, result2@ranges@start+2, result3@ranges@start+2)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "F1",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "F1")
    codon.table<-rbind(codon.table, temp.table)
  } 
  
  #Forward Frame 2
  result1 <- TAA[(TAA@ranges@start+1) %% 3 == 0]   
  result2 <- TGA[(TGA@ranges@start+1) %% 3 == 0]    
  result3 <- TAG[(TAG@ranges@start+1) %% 3 == 0]    
  
  starts<-c(result1@ranges@start-1, result2@ranges@start-1, result3@ranges@start-1)
  ends<-c(result1@ranges@start+1, result2@ranges@start+1, result3@ranges@start+1)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "F2",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "F2")
    codon.table<-rbind(codon.table, temp.table)
  }
  
  #Forward Frame 3
  result1<-TAA[(TAA@ranges@start) %% 3 == 0]   
  result2<-TGA[(TGA@ranges@start) %% 3 == 0]    
  result3<-TAG[(TAG@ranges@start) %% 3 == 0]    
  
  starts<-c(result1@ranges@start-2, result2@ranges@start-2, result3@ranges@start-2)
  ends<-c(result1@ranges@start, result2@ranges@start, result3@ranges@start)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "F3",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "F3")
    codon.table<-rbind(codon.table, temp.table)
  }
  
  #Sets up data
  rev.seq<-as.character(reverseComplement(input.seq))
  
  #Gets codon stuff
  TAA<-matchPattern("TAA", rev.seq)
  TGA<-matchPattern("TGA", rev.seq)
  TAG<-matchPattern("TAG", rev.seq)
  
  #Rev Frame 1
  result1<-TAA[(TAA@ranges@start+2) %% 3 == 0]   
  result2<-TGA[(TGA@ranges@start+2) %% 3 == 0]    
  result3<-TAG[(TAG@ranges@start+2) %% 3 == 0]    
  
  starts<-c(result1@ranges@start, result2@ranges@start, result3@ranges@start)
  ends<-c(result1@ranges@start+2, result2@ranges@start+2, result3@ranges@start+2)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "R1",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "R1")
    codon.table<-rbind(codon.table, temp.table)
  }
  
  #Rev Frame 2
  result1<-TAA[(TAA@ranges@start+1) %% 3 == 0]   
  result2<-TGA[(TGA@ranges@start+1) %% 3 == 0]    
  result3<-TAG[(TAG@ranges@start+1) %% 3 == 0]    
  
  starts<-c(result1@ranges@start-1, result2@ranges@start-1, result3@ranges@start-1)
  ends<-c(result1@ranges@start+1, result2@ranges@start+1, result3@ranges@start+1)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "R2",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "R2")
    codon.table<-rbind(codon.table, temp.table)
  }
  
  #Rev Frame 3
  result1<-TAA[(TAA@ranges@start) %% 3 == 0]   
  result2<-TGA[(TGA@ranges@start) %% 3 == 0]    
  result3<-TAG[(TAG@ranges@start) %% 3 == 0]    
  
  starts<-c(result1@ranges@start-2, result2@ranges@start-2, result3@ranges@start-2)
  ends<-c(result1@ranges@start, result2@ranges@start, result3@ranges@start)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "R3",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "R3")
    codon.table<-rbind(codon.table, temp.table)
  } #end if 
  
  if (codons == T) { return(codon.table) }
  
  if (codons == F) {  
    frames<-unique(codon.table$Frame)
    orf.frame<-data.frame()
    for (x in 1:length(frames)){
      temp.codon<-codon.table[codon.table$Frame == frames[x],]
      temp.codon<-temp.codon[order(temp.codon$Start),]
      
      if (temp.codon$Start[1] == 0){
        temp.start<-as.numeric(gsub("F|R", "", temp.codon$Frame))
        add.frame<-data.frame(FrameStart = temp.start, FrameEnd = width(input.seq), 
                              Size = (width(input.seq)-temp.start)+1, Frame = frames[x])
        orf.frame<-rbind(orf.frame, add.frame)
        next
      }
      #Goes through each of the given directions codons and converts to frame ranges
      temp.frame<-data.frame()
      for (y in 1:(nrow(temp.codon)+1)){
        #First y the start is 1, otherwise take from previous end
        if (y == 1){ frame.start<-as.numeric(gsub("F|R", "", temp.codon$Frame[y])) } else { frame.start<-temp.frame$FrameEnd[y-1]+4 }
        
        #Gets end by subtracting from the codon start
        frame.end<-temp.codon$Start[y]-1
        temp.frame<-rbind(temp.frame, data.frame(FrameStart = frame.start, FrameEnd = frame.end))
      } # end y loop
      
      temp.frame$FrameEnd[nrow(temp.frame)]<-width(input.seq)
      
      #Adds all the data together
      add.frame<-cbind(temp.frame, Size = (temp.frame$FrameEnd-temp.frame$FrameStart)+1, Frame = frames[x])
      orf.frame<-rbind(orf.frame, add.frame)
      
    } #end x loop
    
    orf.frame<-orf.frame[orf.frame$Size >= min.size,]
    return(orf.frame)
  } # end else
}# END FUNCTION

##########
## trim.ends function
######
trim.ends <- function (x, min.n.seq = 4, codon.trim = T){
  #Converts DNAStringSet to something usable
  #x<-trimmed
  #min.n.seq <- 4
  new.align  <- strsplit(as.character(x), "") ### input alignment stored as a list of character vectors
  mat.align  <- lapply(new.align, tolower)    ### same as new.align but with lowercase DNA characters
  x <- as.matrix(as.DNAbin(mat.align))        ### DNAbin format of mat.align

  if (!inherits(x, "DNAbin")) {
    stop("'x' is not of class 'DNAbin'")
  }
  if (!is.matrix(x)) {
    stop("'x' must be a matrix")
  }
  
  replaceWithN <- function(x) {
    id <- (x == as.raw(4))    ### a matrix of logicals, having the same size as x. TRUE for each gap character, FALSE for each non-gap character.
    ## tests if the first or last alignment column of has any gaps
    if (length(id) > 0 & any(id[c(1, length(id))])) {
      id <- which(id)   ### a vector containing the element index/numbers of gap characters of the alignment
      
      getIndex <- function(x) {
        for (i in (seq_along(id) - 1)) {
          if (any(id[1:(i + 1)] != (1:(i + 1)))) 
            break
        }
        id <- rev(id)
        jj <- head(id, 1)
        j <- jj - 1
        for (k in seq_along(id)[-1]) {
          if (any(id[1:k] != (jj:j))) 
            break
          j <- j - 1
        }
        j <- j + 1
        id <- c(0:i, j:jj)
        id[id != 0]
      }
      id    <- getIndex(id)
      x[id] <- as.raw(240)   ### 2-digit byte format of a character
    }
    return(x)
  }
  
  #Does stuff
  x        <- t(apply(x, 1, replaceWithN))
  class(x) <- "DNAbin"
  b        <- as.raw(c(136, 40, 72, 24))
  percentInformation <- function(x, b) {
    length(which(x %in% b))
  }
  m <- apply(x, 2, percentInformation, b)
  if (max(m) < min.n.seq) {
  	stop("alignment contains less sequences then required")
  }
  m <- range(which(m >= min.n.seq))
  
  #Forward Frame 2
  if (codon.trim == T){
    if ((m[1]-1) %% 3 == 0){
    	m[1]<-m[1]
    }
    if ((m[1]-1) %% 3 == 1){
    	m[1]<-m[1]+2
    }
    if ((m[1]-1) %% 3 == 2){
    	m[1]<-m[1]+1
    }
  }
  
  m  <- seq(from = m[1], to = m[2])
  x2 <- as.matrix(x[, m])             ### alignment after trimming ends
  #Converts back
  save.names <- rownames(x2)
  
  #Removes N end gaps
  x3 <- as.list(data.frame(t(as.character(x2))))
  for (y in 1:length(x3)){
   #Starts from the beginning and end to fill in end gaps
    for (q in 1:length(x3[[y]])){ if (x3[[y]][q] == "n"){ x3[[y]][q]<-"-" } else { break } }
    for (q in length(x3[[y]]):1){ if (x3[[y]][q] == "n"){ x3[[y]][q]<-"-" } else { break } }
  }#end x loop
  #Saves final stuff
  temp.align <- lapply(x3, FUN = function(x) paste(x, collapse = ""))
  align.out  <- DNAStringSet(unlist(temp.align))
  names(align.out) <- save.names
  return(align.out)
}

#######
## slice.trim function
#######
## uses a sliding window method to replace poorly aligned subsections of sequences of an alignment (these are replaced with a string of "-")
## slice.size.bp = window width
## threshold = maximum fraction of sites that can be parsimony informative sites (mean of pairwise scores calculated for each individual) for the sequence to be considered aligned "good"
## individuals with fewer than 20bp in the slice.trim alignment are removed
slice.trim <- function(input.align, slice.size.bp = 100, threshold = 0.45){  
  #makes consensus sequence for comparison
  #input.align<-trimal.align
  input.con <- make.consensus(input.align, method = "majority")   ### majority-rule consensus sequence of input.align
  names(input.con)<-"Reference_Locus"                             ### name of the sequence
  
  comb.align  <- append(input.align, input.con)                   ### alignment that includes input.align and the consensus reference locus
  
  #Gets slice information ready
  slice.no    <- ceiling(max(width(input.align))/slice.size.bp)   ### number of slices
  slice.start <- 1
  slice.end   <- slice.size.bp
  
  #checks to see if its out of bounds
  if (slice.end > max(width(input.align))){ 
    slice.end   <-max(width(input.align))
  }#end if check
  output.align <- DNAStringSet()
  for (x in 1:slice.no){
    #Slice alignment into number of slices 
    sliced.align <- subseq(comb.align, start = slice.start, end = slice.end)
    #Checks for badly aligned sequences 
    bad.align    <- pairwise.inf.sites(sliced.align, "Reference_Locus")     ### mean pairwise percent of sites that are parsimony informative for each individual in the slice
    #Remove bad sequence chunks
    rem.seqs     <- bad.align[bad.align >= threshold]                       ### badly aligned sequences to remove (that have mean pairwise percent parsimony informative greater than the threshold)
    good.align   <- sliced.align[!names(sliced.align) %in% names(rem.seqs)] ### the sequences not considered badly aligned 
    #Makes replacement gap seqs for the bad ones
    blank.align <- DNAStringSet()
    if (length(rem.seqs) != 0){
      for (y in 1:length(rem.seqs)){
        blank.align<-append(blank.align, DNAStringSet(paste0(rep("-", slice.end-slice.start+1), collapse = "")) )
      }
      names(blank.align)<-names(rem.seqs)
    }#end rem seqs if
    
    #Saves the slices and cats
    save.slice          <- append(good.align, blank.align)
    save.slice          <- save.slice[order(names(save.slice))]
    save.names          <- names(save.slice)
    output.align        <- DNAStringSet(paste0(as.character(output.align), as.character(save.slice)))
    names(output.align) <- save.names
    
    #Gets new start and stop
##  slice.start <- slice.start+100  #|These two lines were in the original function, but would cause some bases to be overlooked if
##  slice.end   <- slice.end+100    #|slice.size.bp is set to less than 100. Updated versions of these lines are below.
    slice.start <- slice.start+slice.size.bp
    slice.end   <- slice.end+slice.size.bp

    #checks to see if the next slice would be out of bounds
    if (slice.end > max(width(input.align))){ 
      slice.end<-max(width(input.align))
      if (slice.end-slice.start <= 25){
      	break
      } else {
        save.slice          <- subseq(comb.align, start = slice.start, end = slice.end)
        save.slice          <- save.slice[order(names(save.slice))]
        save.names          <- names(save.slice)
        output.align        <- DNAStringSet(paste0(as.character(output.align), as.character(save.slice))) ### concatenates previous slice with current slice
        names(output.align) <- save.names
		break
      }
    }#end if
  }#end x loop
  
  #Removes reference
  output.align <- output.align[names(output.align) != "Reference_Locus"]
  #removes gap only taxa
  str.splitted <- strsplit(as.character(output.align), "")
  x.align      <- as.matrix(as.DNAbin(str.splitted))
  len.temp     <- as.character(as.list(x.align))
  len.loci     <- lapply(len.temp, function (x) x[x != "-"])
  spp.len      <- unlist(lapply(len.loci, function (x) length(x)))         ### number of non-gap characters in each sequence
  spp.rem      <- spp.len[spp.len <= 20]                                   ### individuals to remove from the alignment because they have fewer than 20 non-gap characters in the alignment
  return.align <- output.align[!names(output.align) %in% names(spp.rem)]   ### slice.trim alignment without spp.rem individuals
  return(return.align)
}#end FUNCTION


#######
## run.trimal function
#######
## reads a phylip DNA alignment (input.align)
## saves a copy of input.align in fasta format
## runs trimAL on the fasta formatted alignment, and saves the trimAL alignment in the current directory
## deletes the fasta version of input.align
## reads trimAL alignment into R
## deletes the trimAL alignment from the current directory
## value returned by this function is the trimAL alignment
run.trimal <- function(input.align, method = "auto",locus.name=locus.names[i],trimal.exe.dir=getwd(),trimal.exe.name="trimal"){
  # input.align  <- rem.align
  # Finds probes that match to two or more contigs
  save.rownames <- names(input.align)                                       ### names of taxa in input alignment
  write.align   <- as.list(as.character(input.align))                       ### a list of character strings (each character string is the alignment DNA sequence of an individual)
  ## Writes input alignment as a fasta file in current directory
  locus.name.fa <- paste0(gsub(pattern = "\\..*", "", locus.name), ".fa")   ### name of temportary fasta alignment that is the same as input.align
  write.fasta(sequences = write.align, names = names(write.align),file.out=locus.name.fa, nbchar = 1000000, as.string = T)
  ## name of write.align (other than format, this is the same alignment as input.align)
  input.file      <- locus.name.fa  
  output.file     <- paste0(input.file,"-tm")
  trimAL.exe.path <- paste0(trimal.exe.dir,"/",trimal.exe.name)                             ### full path to the trimal executable
  system(paste0(trimAL.exe.path," -in ", input.file, " -out ", output.file," -automated1")) ### executes trimAL
  system(paste0("rm ", input.file))                                                         ### removes the fasta version of input.align
  
  ## if trimAL alignment does not exist, delete the fasta version of input.align and print a warning
  if (file.exists(output.file) == F) { 
    system(paste0("rm ", input.file))
    print(paste0("deleted. Not enough overlapping data in alignment.") )
    return(DNAStringSet())
  } else {
  	system(paste("mv",output.file,input.file))  ### renames the trimAL alignment filename
  }
  
  out.align <- scanFa(FaFile(input.file))       ### reads the trimAL alignment
  
  #Fixes any terrible NA names introduced by trimal
  new.names<-c()
  for (j in 1:length(names(out.align))){ 
    new.names[j] <- save.rownames[grep(pattern = names(out.align)[j], x = save.rownames)]
  }
  
  temp <- names(out.align)[is.na(names(out.align)) == T]
  if (length(temp) > 0){
  	stop("there are NAs in the names")
  }
  names(out.align) <- new.names
  unlink(input.file) ### deletes the trimAL alignment (because its not in its final form yet)
  return(out.align)  ### the value of the function is out.align
}#end function


#####
## make.consensus function
#####
make.consensus <- function (input.alignment, method = c("majority","threshold","IUPAC","profile"), threshold = 0.6, warn.non.IUPAC = FALSE, type = c("DNA", "RNA")) {
  
  #input.alignment<-trimmed
  #Converts alignment to matrix of characters to be used
  new.align <- strsplit(as.character(input.alignment), "")                           ### alignment as a list of vectors
  align.in  <- matrix(unlist(new.align), ncol = length(new.align[[1]]), byrow = T)   ### alignment as a matrix 
  
  #Does based on method
  method <- match.arg(method) ### this only works when excuting the function
  
  if (method == "IUPAC") {
    type <- match.arg(type)
    res <- apply(align.in, 2, bma, warn.non.IUPAC = warn.non.IUPAC, 
                 type = type)
    names(res) <- NULL
  }
  if (method == "majority") {
    majority <- function(x) names(which.max(table(x)))
    res <- apply(align.in, 2, majority)
    names(res) <- NULL
  }
  if (method == "profile") {
    obsvalue <- levels(factor(align.in))
    nrow     <- length(obsvalue)
    row.names(align.in) <- NULL
    res <- apply(matali, 2, function(x) table(factor(x, levels = obsvalue)))
  }
  if (method == "threshold") {
    profile <- consensus(align.in, method = "profile")
    profile.rf <- apply(profile, 2, function(x) x/sum(x))
    res <- rownames(profile.rf)[apply(profile.rf, 2, which.max)]
    res <- ifelse(apply(profile.rf, 2, max) >= threshold, 
                  res, NA)
    names(res) <- NULL
  }
  
  out.consensus<-DNAStringSet(paste0(res, collapse = ""))
  names(out.consensus)<-"Consensus_Sequence"
  return(out.consensus)
}

### charTest function checks that a sequence has at least one of each nucleotide types A,C,G,T
### ls = a list containing one or more character vectors
### this function is mostly useful after first using the uniqueLetters function
charTest <- function(ls){
		ls <- tolower(ls)
		all(c("a","c","g","t") %in% ls)
}

##################
## percent.gaps function
##################	
## alignment must be an Xstringset
## returns the percentage of the alignment that is a gap "-"
percent.gaps <- function(alignment){
	area.alignment <- ((unique(width(alignment)))*length(alignment))
	noGaps.seqs <- RemoveGaps(alignment)
	result <- (area.alignment-sum(width(noGaps.seqs)))/area.alignment
	result
}

##################
### pis function
##################
### alternative parsimony informative sites calculation
### this is the pis function in phyloch package
#####
### x = DNAbin alignment
###
pis.new <- function (x, abs = TRUE, use.ambiguities = FALSE,as.percent=F){
	if(class(x)=="DNAStringSet"){
		x <- as.DNAbin(DNAMultipleAlignment(x))
	}
	pars.inf <- function(x) {
		x <- table(x)
		x <- x[x > 1]                                                     # This is the number of each type of character in a given column
		n <- c("-", "n", "b", "h", "d", "v", "k", "s", "r", "w","y","?")  # A list of possible characters to check against
		if (length(x[!names(x) %in% n]) > 1){                             #| These lines are specifying whether the
			x <- TRUE                                                 #| the ith column has at least two DIFFERENT but not ambiguous characters,
		} else {                                                          #|
			x <- FALSE                                                #|
		}                                                                 #|
	}
	nbchar <- dim(x)[2]
	x <- as.character(x)
	out <- apply(x, 2, pars.inf)
	out <- length(out[out])
	if (!abs){
		out <- round(out/nbchar * 100, digits = 2) #### If abs = FALSE, the output is the percentage of sites with at least two different and not ambiguous character states
	}
	if(abs ==F & as.percent==F){
		out <- out/100
	}
	out
}

##########
## alignment.stats function
## make alignment info table containing info about each locus alignment prior to concatenating them
## align.dir is the directory containing the alignments
##########
## Dependencies
##  functions: percent.gaps, pis.new
##  packages: ape, Biostrings, data.table
alignment.stats <- function(align.dir,outfile,seqs.format="phylip"){
	ext         <- c("phy","fa"); names(ext) <- c("phylip","fasta")
	ext.temp    <- as.character(ext[seqs.format])
	#files       <- list.files(path=align.dir, pattern="*.phy",full.names=T)  ### list of the names of alignment files
	#short.names <- list.files(path=align.dir, pattern="*.phy")         ###list of the names of alignment files
	
	files       <- list.files(path=align.dir, pattern=paste("*.",ext.temp,sep=""),full.names=T)  ### list of the names of alignment files
	short.names <- list.files(path=align.dir, pattern=paste("*.",ext.temp,sep=""))         ###list of the names of alignment files

	stats.table <- matrix(data="0",nrow=length(files),ncol=10)
	colnames(stats.table) <- c("locus","n_individuals","n_characters","percent_gaps","fraction_sites_pars_inform","mean_pDistance","all_some_overlap","stat7","stat8","stat9","stat10")[1:ncol(stats.table)]
	
	for(i in 1:length(files)){
		
		multipleAlignment.temp <- readDNAMultipleAlignment(file = files[i], format = seqs.format)
		DNAbin.temp    <- as.DNAbin(multipleAlignment.temp)   ### DNAbin version of the alignment
		alignment.temp <- unmasked(multipleAlignment.temp)
		
		distances        <- dist.dna(DNAbin.temp,model="raw",pairwise.deletion=T)
		mean.distance    <- round(mean(distances[which(distances!="NaN")]),digits=3)
		all.some.overlap <- all(distances!="NaN")
		name.temp   <- gsub(pattern=".phy","",x=short.names[i])

		stats.table[i,1] = name.temp
		stats.table[i,2] = as.numeric(length(alignment.temp))
		stats.table[i,3] = as.numeric(width(alignment.temp)[1])
#		stats.table[i,4] = as.numeric(round(percent.gaps(alignment.temp),digits=5))
		stats.table[i,5] = pis.new(alignment.temp,abs=F) ### percent of sites that are informative
		stats.table[i,6] = mean.distance
		stats.table[i,7]= all.some.overlap
		# stats.table[i,6] = round(mean(dist.dna(DNAbin.temp,model="raw",pairwise.deletion=T)),digits=3)   ### mean pairwise distance among individuals
		# stats.table[i,6]= ((1-as.numeric(round(percent.gaps(alignment.temp),digits=5)))*as.numeric(width(alignment.temp)[1])*as.numeric(length(alignment.temp)))
		# stats.table[i,8]=
		# stats.table[i,9]=
		
		if(i==1){
			write(x=colnames(stats.table),ncolumns=ncol(stats.table),file=outfile,append=T)
		}
		write(x=stats.table[i,],ncolumns=ncol(stats.table),file=outfile,append=T)
	}
	stats.table
}

#####
## align.alignment.columns function
#####
## for asthetics only
## uses the linux commands "column" and "sed" to ensure that columns of the alignment are lined up in the text file
## even if sequence names are of different lengths
####
## phylip.alignment = filename or path to sequence alignment in sequential phylip format
## nspaces = minimum number of spaces separating columns (this is the number of spaces between the longest sequence name and the first alignment column)
####
# this doesnt quite work yet because there is a limit on the character length of columns, and DNA seqs often too long
align.alignment.columns <- function(phylip.alignment,npaces=6){
	### define temporary filename
	tmp.name        <- paste0(phylip.alignment,"_tmp")
	
	### align the columns of the sequences in the alignment
	command.string1 <- paste0("column -t -s '      ' '",phylip.alignment,"' > '",tmp.name,"'")
	system(command.string1)

#	command.string1 <- paste0("awk 'NR==FNR{for(i=1;i<=NF;i++) max[i] = length($i) > max[i] ? length($i) : max[i]; next} { for(i=1;i<=NF;i++) printf '%-'max[i]'s  ', $i; printf '\n'}' ",phylip.alignment," ",phylip.alignment," > ",tmp.name)
#	system(command.string1)

	### replace multiple spaces with a single space on first line (header line of phylip alignment)
	command.string2 <- paste0("sed -i '' '1 s/ * / /' '",tmp.name,"'")
	system(command.string2)
	
	### deletes strings of whitespace at ends of lines
	command.string3 <- paste0("sed -i '' 's/ *$//' '",tmp.name,"'")
	system(command.string3)
	
	### renames the temporary file as the name of the original file
	command.string4 <- paste("mv",tmp.name,phylip.alignment)
	system(command.string4)
}

###########################
#### scanFa.multiple function
# used to pull in fasta formatted data stored in multiple files
# directory = path to folder holding the fasta formatted files
scanFa.multiple <- function(directory){
	scan.fa <- function(filename){
		res <- scanFa(FaFile(filename))
		res
	}
	
	files <- list.files(directory,full.names=T)
	for(i in 1:length(files)){
		temp.data <- scan.fa(files[i])
		if(i==1){
			all.data <- temp.data
		}
		if(i>1){
			all.data <- append(all.data,temp.data)
		}
		#detach(temp.data)
	}
	all.data
}


			
###########################
## make.partitioned.alignment function
###########################
## Dependencies:
##    R packages: Biostrings, stringr, ape, data.table
#############
## parameters (probably need to update this to make it more generalizable)
###########################
## output.dir                     # where to save alignments (alignments for each type of data will be in different subdirectory)
## InputAlignmentFolder           # folder containing input alignments, which are not partitioned by coding region
## TargetCDS.path                 # full path to the fasta file containing only the CDS sequence of target loci (from which probes were designed). Sequence names must have the following format: "<GeneName>_TargetCDS_of_<TargetLocusName>_<AnyAdditionalIformation>", where <GeneName> and <TargetLocusName> are replaced with the actual names, and <AnyAdditionalIformation> can be a string of any characters
## bait.species.filename          # full path to the file containing a two column table in which the first column contains name of locus and the second column contains name of species that the probes for that locus were designed from
## ref.type = "DNA"               # indicates that the input alignments are DNA sequences
## old.names=NA                   # character vector of sample names that should be changed to those defined in new.names parameter
## new.names=NA                   # character vector of sample names to change from those defined in old.names parameter
## drop.reference=F               # whether or not to include the reference sequence in the alignment written to file
## ith.locus.start=1              # first locus to process
## ith.locus.end = "all"          # last locus to process. "all" means process until no more loci to process.
## locus.names.omit=NULL          # set to "WeinellEntry5179" for SnakeCap analyses, because the wrong target sequence was used when designing probes for this locus. I intended to target the Beta-keratin 2 gene, but instead targeted a non-coding region 80,000bp away.
## AA.pdist.drop.thresh=0.5       # maximum mean pairwise p-distance for AA sequence of an individual to keep the individual in AA alignments or alignments containing CDS regions

#########
# Example Usage:
#     source("~/SnakeCap_Functions.R") # load in the functions in this file
#     library(Biostrings)
#     library(stringr)
#     library(ape)
#     library(data.table)
#     make.partitioned.alignment(InputAlignmentFolder="~/Immune/unpartitioned/", output.dir="~/Immune/partitioned/", TargetCDS.path="~/Weinell_TargetLoci_Snakes_Final_targetCDS_v3.fa", bait.species.filename="~/bait_species_table.txt")
#########

make.partitioned.alignment       <- function(InputAlignmentFolder,output.dir,TargetCDS.path,bait.species.filename,ref.type="DNA",old.names=NA,new.names=NA,drop.reference=F,ith.locus.start=1,ith.locus.end="all",locus.names.omit=NULL,AA.pdist.drop.thresh=0.5){
	#source(libs)            ### loads required libraries from the libraries file
	#source(functions)       ### loads required functions from the functions file
	
	### makes necessary output subdirectories if they dont exist ###
	dir1  <- paste(output.dir,"All_parts/alignmentFiles/",sep="")
	dir2  <- paste(output.dir,"All_parts/partitionFiles/",sep="")
	dir3  <- paste(output.dir,"CDS_only/alignmentFiles/",sep="")
	dir4  <- paste(output.dir,"CDS_only/partitionFiles/",sep="")
	dir5  <- paste(output.dir,"Upstream_noncoding/alignmentFiles/",sep="")
	dir6  <- paste(output.dir,"Downstream_noncoding/alignmentFiles/",sep="")
	dir7  <- paste(output.dir,"All_noncoding/alignmentFiles/",sep="")
	dir8  <- paste(output.dir,"All_noncoding/partitionFiles/",sep="")
	dir9  <- paste(output.dir,"CDS_FirstCodonPosition/alignmentFiles/",sep="")
	dir10 <- paste(output.dir,"CDS_SecondCodonPosition/alignmentFiles/",sep="")
	dir11 <- paste(output.dir,"CDS_ThirdCodonPosition/alignmentFiles/",sep="")
	dir12 <- paste(output.dir,"AminoAcids/alignmentFiles/",sep="")

	subdirectories  <- c(dir1,dir2,dir3,dir4,dir5,dir6,dir7,dir8,dir9,dir10,dir11,dir12)
	makeDirectories <- lapply(X=subdirectories,FUN=dir.check.create)  ### need to assign an object to avoid printing results of lapply 
	
	### Reads in files ###
	TargetDNA_CDS.regions       <- readDNAStringSet(TargetCDS.path,format="fasta")
	bait.species.table          <- fread(bait.species.filename) ### reads the table specifying which species the probes were designed from for each locus 

	#targetCDS.names             <- gsub(pattern=".*_Weinell",replacement="Weinell",x=names(TargetDNA_CDS.regions))  ###| used this if working with beta-keratin CDS
	#targetCDS.names             <- gsub(pattern="_.+",replacement="",x=targetCDS.names)                             ###| 
	#test1 <- gsub(pattern=".*_TargetCDS_of_",replacement="",x=names(TargetDNA_CDS.regions))
	#test2 <- gsub(pattern="_.+",replacement="",x=test1)
	#test3 <- strsplit(names(TargetDNA_CDS.regions))
	#test4 <- stri_list2matrix(test3,byrow=T)

	name.string.start           <- (str_locate_X(strings=names(TargetDNA_CDS.regions),pattern="_",X=3)+1)
	name.string.end             <- (str_locate_X(strings=names(TargetDNA_CDS.regions),pattern="_",X=4)-1)
	targetCDS.names             <- substring(names(TargetDNA_CDS.regions),first=name.string.start,last=name.string.end) ### target loci names of each sequence in TargetDNA_CDS.regions
	input.alignment.filenames   <- list.files(InputAlignmentFolder,full.names=T)
	input.alignment.shortnames  <- gsub(".phy","",list.files(InputAlignmentFolder))

	## 27/27 immune genes are in TargetDNA_CDS.regions
	## 1652/1652 WholeExon genes are in TargetDNA_CDS.regions  ### potentially re-add the ones that were missing to include the potentially one or two bases not included (for those with frame not 0)
	## 54/95 scalation genes missing from TargetDNA_CDS.regions
	## 29/119 vision genes missing from TargetDNA_CDS.regions
	## 85/??? exon-containing UCEs present in TargetDNA_CDS.regions
	## 0/??? exon-containing ddRAD-like loci present in TargetDNA_CDS.regions

	shared.names                <- Reduce(intersect, list(input.alignment.shortnames,targetCDS.names)) ### updates shared.names so that they only processes loci that have been aligned and that have a known CDS region (which is included in TargetDNA_CDS.regions)
	if(!is.null(locus.names.omit)){
		shared.names  <- setdiff(shared.names,locus.names.omit)   ### doesnt process loci in locus.names.omit
	}
	
	if(file.exists(paste(output.dir,"partitioned_alignments_made.txt",sep=""))){
		alignments.made <- read.table(file=paste(output.dir,"partitioned_alignments_made.txt",sep=""),header=T,colClasses="character")
	} else {
		alignments.made <- matrix(data="no",nrow=length(shared.names),ncol=9)
		rownames(alignments.made)   <- shared.names
		colnames(alignments.made)   <- c("all.data","CDS","5'NC","3'NC","Noncoding","AA","FirstCodon","SecondCodon","ThirdCodon")
	}
	if(ith.locus.end=="all"){
		last.locus.process          <- length(shared.names)
	} else {
		last.locus.process          <- ith.locus.end
	}
	for(i in ith.locus.start:last.locus.process){
		#if(is.wholenumber(i/50)){
		#	print(i)
		#}
		print(i)
		locus.name.temp   <- shared.names[i]
		bait.species.temp <- bait.species.table$Species[which(bait.species.table$Bait==locus.name.temp)]
	
		novel             <- readDNAMultipleAlignment(input.alignment.filenames[which(input.alignment.shortnames==locus.name.temp)],format="phylip")
		novel2            <- DNAStringSet(x=gsub("-","",novel))    ### un-aligning sequences in the "novel" alignment

		if(ref.type=="DNA"){
			reference.cds        <- TargetDNA_CDS.regions[which(targetCDS.names==locus.name.temp)]
			names(reference.cds) <- paste(bait.species.temp,names(reference.cds),sep="_")
			final.locus             <- c(reference.cds,novel2)
		}
		
		if(length(final.locus)<5){        ##| Skips locus if fewer than 5 sequences in final.locus
			next                          ##| (i.e. need at least 4 sequences other than the reference CDS)
		}                                 ##|
	
		alignment.names  <- names(final.locus)
		if(!all(is.na(old.names))){
			alignment.names  <- mgsub(old.names,new.names,alignment.names)
		}
		
		alignment        <- rMSA::mafft(final.locus,param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6") ### alternative way of calling MAFFT using rMSA package
		alignment        <- as(alignment, "DNAStringSet")
		names(alignment) <- alignment.names
	
		begin.exon       <- str_locate(string=alignment[1],pattern=as.character(subseq(reference.cds,start=1,end=1)))[1]
		end.exon         <- str_locate_last(string=alignment[1],pattern=as.character(subseq(reference.cds,start=(width(reference.cds)-1),end=width(reference.cds))))
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
					cds.alignment     <- subseq(x=alignment,start=begin.exon,end=end.exon)[-1]
				} else {
					cds.alignment     <- subseq(x=alignment,start=begin.exon,end=end.exon)
			}
			drop.nodata               <- str_count(cds.alignment,"-") == width(cds.alignment)    ### specifies individuals with no data
			if(length(which(!drop.nodata))>3){
				cds.alignment         <- cds.alignment[!drop.nodata]
				cds.alignment2        <- DNAStringSet(x=gsub("-","",cds.alignment))
				cds.alignment2        <- rMSA::mafft(cds.alignment2,param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6")
				cds.alignment2        <- as(cds.alignment2, "DNAStringSet")
				names(cds.alignment2) <- names(cds.alignment)
			} else {
				partition.groups[2] <- F   ### updates decision from TRUE to FALSE on whether to write the CDS alignment, because too few individuals with CDS data
			}
		}
		
		#### Upstream noncoding alignment
		if(partition.groups[1]){
			if(begin.exon!=1){
				if(drop.reference==T){
					upstream.alignment    <- subseq(x=alignment,start=1,end=(begin.exon-1))[-1]
				} else {
					upstream.alignment    <- subseq(x=alignment,start=1,end=(begin.exon-1))
				}
				drop.nodata                <- str_count(upstream.alignment,"-") == width(upstream.alignment)
				if(length(which(!drop.nodata))>3){
					upstream.alignment         <- upstream.alignment[!drop.nodata]
					upstream.alignment2        <- DNAStringSet(x=gsub("-","",upstream.alignment))
					upstream.alignment2        <- rMSA::mafft(upstream.alignment2,param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6")
					upstream.alignment2        <- as(upstream.alignment2, "DNAStringSet")
					names(upstream.alignment2) <- names(upstream.alignment)
				} else {
					partition.groups[1] <- F   ### updates decision from TRUE to FALSE on whether to write the upstream alignment, because too few individuals with upstream noncoding data
				}
			}
		}
		
		#### Downstream noncoding alignment
		if(partition.groups[3]){
			if(end.exon!=width(alignment[1])){
				if(drop.reference==T){
					downstream.alignment    <- subseq(x=alignment,start=(end.exon+1),end=width(alignment[1]))[-1]
				} else {
					downstream.alignment    <- subseq(x=alignment,start=(end.exon+1),end=width(alignment[1]))
				}				
				drop.nodata                  <- str_count(downstream.alignment,"-") == width(downstream.alignment)
				### skips downstream noncoding alignment if 3 or less individuals
				if(length(which(!drop.nodata))>3){
					downstream.alignment         <- downstream.alignment[!drop.nodata]
					length(downstream.alignment)
					downstream.alignment2        <- DNAStringSet(x=gsub("-","",downstream.alignment))
					downstream.alignment2        <- rMSA::mafft(downstream.alignment2,param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6")
					downstream.alignment2        <- as(downstream.alignment2, "DNAStringSet")
					names(downstream.alignment2) <- names(downstream.alignment)
				} else {
					partition.groups[3] <- F   ### updates decision from TRUE to FALSE on whether to write the downstream alignment, because too few individuals with downstream noncoding data
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
			cds.temp3 <- subseq(cds.temp2,start=min(p1),end=max(p3))
			### CDS sequences not aligned
			cds.temp4 <- DNAStringSet(x=gsub("-","",cds.temp3))   ### unaligns sequences (after removing gaps) so that the sequences can be translated 
			
			### Translates CDS sequences, then does MSA to make an alignment of amino acids
			aa.temp                  <- suppressWarnings(translate(cds.temp4,if.fuzzy.codon="solve",no.init.codon=T))
			aa.alignment.temp        <- rMSA::mafft(aa.temp,param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6")
			aa.alignment.temp        <- AAStringSet(aa.alignment.temp)
			names(aa.alignment.temp) <- names(aa.temp)

			aa.alignment.temp2  <- filter.alignment(aa.alignment.temp)       # Nothing filtered by default
			aa.alignment.mdt0   <- filter.alignment(aa.alignment.temp,mdt=0) # This alignment has no missing or ambiguous data. This alignment is used for calculating median pairwise p-distances, but is not the alignment written to file.
			
			distAA.mdt1          <- stringDist(aa.alignment.temp2)
			distAA.mdt0          <- stringDist(aa.alignment.mdt0)
			median.distAA.mdt1   <- apply(X=distAA.mdt1,MARGIN=1,FUN=median)/width(aa.alignment.temp2[1])  ### median pairwise p-distances for the AA alignment that can contain any amount of missing or ambiguous data
			median.distAA.mdt0   <- apply(X=distAA.mdt0,MARGIN=1,FUN=median)/width(aa.alignment.mdt0[1])   ### median pairwise p-distances for the AA alignment containing no missing or ambiguous data
			
			toDrop <- which(median.distAA.mdt0 > AA.pdist.drop.thresh) ### individuals with median p-distance greater than AA.pdist.drop.thresh should be dropped
			
			if(!!length(toDrop)){                                    ###| creates a character vector containing the names of individuals
				toDrop.names    <- names(median.distAA.mdt0)[toDrop] ###| that should be dropped from the alignment because they are highly
			} else {                                                 ###| diverged from most other individuals
				toDrop.names <- NULL                                 ###|
			}                                                        ###|
			
			if(!!length(toDrop)){                                          ###| Drops individuals with unusually high amount of AA divergence unless all individuals dropped
				if(length(toDrop)==length(aa.alignment.temp2)){
					next
				}
				aa.alignment.temp3 <- aa.alignment.temp2[-toDrop]                    ###| sites-unfiltered AA alignment with highly diverged individuals removed
				aa.alignment.temp3 <- AAStringSet(x=gsub("-","",aa.alignment.temp3)) ###|
				names.aa.temp3     <- names(aa.alignment.temp3)                      ###|
				aa.alignment.temp3 <- rMSA::mafft(aa.alignment.temp3,param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6") ### rerun alignment algorithm
				aa.alignment.temp3 <- AAStringSet(aa.alignment.temp3)
				names(aa.alignment.temp3) <- names.aa.temp3
			} else {                                               ###| if no individuals needed to be removed, then aa.alignment.temp3 is the same as aa.alignment.temp2
				aa.alignment.temp3 <- aa.alignment.temp2           ###|
			}                                                      ###|
			if(!length(aa.alignment.temp3)){  ##|Skips to next locus if no sites are retained
				next                          ##|
			}                                 ##|
			
			#### Writes target AA-only alignment
			alignments.made[i,6] <- "yes"
			writeXStringSet(x = aa.alignment.temp3, filepath=paste0(dir12,locus.name.temp,".fa"), append=FALSE,compress=FALSE, compression_level=NA, format="fasta")
			
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
			write.dna(x = cds.p1,file=paste0(dir9,locus.name.temp,".phy"),format="sequential",append=F,nbcol=1,colw=100000000)
			write.dna(x = cds.p2,file=paste0(dir10,locus.name.temp,".phy"),format="sequential",append=F,nbcol=1,colw=100000000)
			write.dna(x = cds.p3,file=paste0(dir11,locus.name.temp,".phy"),format="sequential",append=F,nbcol=1,colw=100000000)
			
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
					ape::write.dna(x= cds.temp7,file=paste0(dir3,locus.name.temp,".phy"),format="sequential",append=F,nbcol=1,colw=100000000)
					alignments.made[i,2]  <- "yes"
					cds.starts            <- c(1:3)
					cds.ends              <- rep(width(cds.temp7[1]),3)
					cds.ranges            <- cbind(cds.starts,cds.ends)
					cds.partition.line    <- NULL
					for(j in 1:nrow(cds.ranges)){
						cds.partition.line[j] <- paste("DNA, ", rownames(cds.ranges)[j]," = ",cds.ranges[j,1],"-",cds.ranges[j,2],"\\3",sep="")
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
					upstream.alignment3 <- filter.alignment(upstream.alignment2)  ### columns are required to have some non-gap, non-ambiguous data
				}
			}	else {
				upstream.alignment3 <- filter.alignment(upstream.alignment2)      ### columns are required to have some non-gap, non-ambiguous data
			}
			if(length(upstream.alignment3)>3){
				ape::write.dna(x = upstream.alignment3,file=paste0(dir5,locus.name.temp,".phy"),format="sequential",append=F,nbcol=1,colw=100000000)
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
				ape::write.dna(x = downstream.alignment3,file=paste0(dir6,locus.name.temp,".phy"),format="sequential",append=F,nbcol=1,colw=100000000)
				alignments.made[i,4]  <- "yes"
			} else {
				partition.groups[3] <- F
			}
		}
		### Writing all-noncoding alignment (upstream and downstream non-coding regions concatenated)
		if(partition.groups[1] & partition.groups[3]){
			upstream.dataTable          <- data.table(as.matrix(upstream.alignment3),keep.rownames = TRUE)
			downstream.dataTable        <- data.table(as.matrix(downstream.alignment3),keep.rownames = TRUE)
			noncoding.alignment3        <- merge(upstream.dataTable, downstream.dataTable, by="rn", all=TRUE)
			noncoding.alignment3        <- na.replace(noncoding.alignment3,"-")
			rn.temp                     <- noncoding.alignment3$rn
			noncoding.alignment3        <- apply(noncoding.alignment3[ ,!"rn"], 1, paste, collapse="")
			noncoding.alignment3        <- DNAStringSet(noncoding.alignment3)
			names(noncoding.alignment3) <- rn.temp
			ape::write.dna(x = noncoding.alignment3,file=paste0(dir7,locus.name.temp,".phy"),format="sequential",append=F,nbcol=1,colw=100000000)
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
				noncoding.partition.line[j] <- paste("DNA, ", rownames(noncoding.ranges)[j]," = ",noncoding.ranges[j,1],"-",noncoding.ranges[j,2],sep="")
			}
			write(noncoding.partition.line,file=paste0(dir8,locus.name.temp,"_parts.txt"))
		}
		### Writing the all-data DNA alignment (which contains upstream non-coding, CDS, and downstream non-coding regions)
		if(all(partition.groups)){
			if(partition.groups[1]){
				upstream.dataTable.all    <- data.table(as.matrix(upstream.alignment3),keep.rownames = TRUE)
			} else {
				upstream.dataTable.all    <- NULL
			}
			if(partition.groups[2]){
				cds.datatable.all         <- data.table(as.matrix(cds.temp7),keep.rownames = TRUE)
			} else {
				cds.datatable.all         <- NULL
			}
			if(partition.groups[3]){
				downstream.dataTable.all  <- data.table(as.matrix(downstream.alignment3),keep.rownames = TRUE)
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
			all.alignment2        <- na.replace(dat,"-")
			rn.temp               <- all.alignment2$rn
			all.alignment2        <- apply(all.alignment2[ ,!"rn"], 1, paste, collapse="")
			all.alignment2        <- DNAStringSet(all.alignment2)
			names(all.alignment2) <- rn.temp
			
			if(!all(is.na(ranges2[2,c(1,2)]))){
				ranges2[3,1] <- ranges2[2,1]+1
				ranges2[4,1] <- ranges2[2,1]+2
				ranges2[c(3:4),2] <- ranges2[2,2]
			}
			
			extra.text2       <- c("","\\3","\\3","\\3","")
			ranges3           <- cbind(ranges2,extra.text2)
			rownames(ranges3) <- c("upstream.noncoding","CDS.1","CDS.2","CDS.3","downstream.noncoding")
			ranges3      <- ranges3[c(1,2,2,2,3) %in% which(partition.groups),]
			ape::write.dna(x= all.alignment2,file=paste0(dir1,locus.name.temp,".phy"),format="sequential",append=F,nbcol=1,colw=100000000)
			### updates the alignments.made matrix
			alignments.made[i,1]  <- "yes"
			### preparing to write the partition file
			partition.line       <- NULL
			for(j in 1:nrow(ranges3)){
				partition.line[j] <- paste("DNA, ", rownames(ranges3)[j]," = ",ranges3[j,1],"-",ranges3[j,2],ranges3[j,4],sep="")
			}
			write(partition.line,file=paste0(dir2,locus.name.temp,"_parts.txt"))
		}
		#### Writes table summarizing which alignments were written
		write.table(x=alignments.made,file=paste(output.dir,"partitioned_alignments_made.txt",sep=""))
	} # End i for loop
}
			
			
##############
## plotAlignment function
##############
### Dependencies:
## libraries: stringr
############
## make miniplot of alignment
## alignment must be a DNAStringSet object
##############
plotAlignment <- function(alignment){
	nsamples  <- length(alignment)
	width.al  <- width(alignment[1])
	xvalsA    <- seq(from=1,to=width.al,by=100)
	xvalsB    <- rep(xvalsA,nsamples)
	yvals     <- rep(c(1:nsamples),length(xvalsA))
	plot(xvalsB,yvals,col="white",xlab="position",ylab="sample")
	for(i in 1:length(alignment)){
		gapLocation   <- unique(unlist(str_locate_all(alignment[i],pattern="-")))
		noGapLocation <- c(1:width.al)[-gapLocation]
		points(noGapLocation,rep(i,length(noGapLocation)),pch=15)
	}
}

###########
##### get_ncbi_sequences function
###########
##### outfile = where to save sequences
##### acccessionList <- character vector of accession numbers
##### startList = numeric vector of starting positions; default = 1
##### endList = numeric vector of end positions
##### strandList = character string or vector indicating strand of the sequence; default = "1"
##### 			strandList = "1" | download forward sequence
##### 			strandList = "2" | download reverse complement sequence
##### db = character string of the name of the database to search; changing from "nuccore" will probably cause problems
##### rettype = character string indicating the format to download sequences; default = "fasta", but options include:
#####			rettype = "fasta" | sequences download as fasta files
#####			rettype = "gb" | sequences download as genbank flat files
#####			rettype = "native" | download full sequence record in XML format
#####			rettype = "acc" | download accession numbers
#####			rettype = "seqid" | download SeqID string
#####			rettype = "ft" | download feature table         # doesnt include as much feature information as rettype=gb
##### retmode = string indicating the mode to download sequences; default = "text"; don't change from default
###########
get_ncbi_sequences <- function(outfile="",accessionList,startList=1,endList,strandList="1",db="nuccore",rettype="fasta",retmode="text"){
	URLs <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
		"db=nuccore",
		"&id=",accessionList,
		"&strand=",strandList,
		"&seq_start=",startList,
		"&seq_stop=",endList,
		"&rettype=",rettype,
		"&retmode=",retmode,
		sep=""
	)
	for(i in 1:length(URLs)){
		temp.URL <- getURL(URLs[i])
		write(temp.URL,file=outfile,append=T)
	}
}

####################
#### include.accession.region.gb function
####################
#### this function creates a temporary file (.tmp) containing a modified version of a
#### GenBank flatfile, in which the region information (ie, start and end info)
#### is included as part of the accession name. This modification allows the
#### the region info to be acccessed when using the AA.from.gb function
include.accession.region.gb <- function(gb.filename){
	original.string         <- readChar(gb.filename, file.info(gb.filename)$size)
	new.string <- gsub(pattern=" REGION: ",replacement=":",original.string)
	outfile.tmp <- gsub(".gb",".tmp",gb.filename)
	write(new.string,file=outfile.tmp)
}

####################
#### read.gb function
####################
#### this function imports a GenBank flatflie into R
####
read.gb <- function(gb.filename,progress=T){	
	include.accession.region.gb <- function(gb.filename, outfile){
		original.string  <- readChar(gb.filename, file.info(gb.filename)$size)
		new.string       <- gsub(pattern=" REGION: ",replacement=":",original.string)
		write(new.string,file=outfile)
	}
	if(progress==T){
		progress.val <- T
	} else {
		progress.val <- F
	}
	outfile.tmp    <- paste(gsub("\\..+","",gb.filename),".tmp",sep="")
	include.accession.region.gb(gb.filename, outfile=outfile.tmp)
	gbData         <- biofiles::gbRecord(outfile.tmp,progress=progress.val)   ### read in the GenBank flatfile (which can contain multiple records)
	if(gb.filename!=outfile.tmp){
		file.remove(outfile.tmp)
	}
	gbData
}

####################
### CDS.from.gb function
###################
### obtains CDS portion of a sequence from a genbank flatfile
### also obtains AA sequences (target and full sequence)
### for some loci the start_codon is poorly annotated in ncbi, but this function determines the start codon by comparing the translated CDS to the NCBI AA sequence
### makes lists of the target loci that with unusual or absent annotations
### library(biofiles)
### targetTablePath    <- "~/TargetTable1_28Sep2019.txt"  ### Four column table, with columns:
###                                                                                                           ### [,1]: locus names
###                                                                                                           ### [,2]: genbank ID of contig that contains target locus
###                                                                                                           ### [,3]: start position of target within contig
###                                                                                                           ### [,4]: end position of target locus within contig.
###                                                                                                           ### If start position > end position, then coding sense of target sequence is negative (else sense is positive)
### targetTable        <- fread(targetTablePath)  ### not updated for vision genes and maybe some others
### gb.object          <- gbData
### targetTable.object <- targetTable
### i=1
CDS.from.gb            <- function(gb.object,targetTable.object,output.filename=NA,report.noCDS=T,debug=F){

	#is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
	#	abs(x - round(x)) < tol
	#}
	translate <- Biostrings::translate ### needed to because multiple libraries have a translate function
	
	gbData                  <- gb.object
	loci.gb                 <- gsub("\\.\\.","-",names(gbData))  ### GenBank accession ID plus start and end region
	targetTable             <- targetTable.object
	accession.list          <- gsub("\\..","",targetTable$TargetName_NCBI)
	#start.list              <- apply(X=rbind(targetTable$TargetStart,targetTable$TargetEnd),MARGIN=2,FUN=min) #
	#end.list                <- apply(X=rbind(targetTable$TargetStart,targetTable$TargetEnd),MARGIN=2,FUN=max) #
	#negative.sense          <- which((targetTable$TargetEnd-targetTable$TargetStart)<0)
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
			cds.f2            <- subseq(x=CDS.temp,start=2)
			cds.f3            <- subseq(x=CDS.temp,start=3)
			cds.RC.f1         <- reverseComplement(CDS.temp)
			cds.RC.f2         <- subseq(x=reverseComplement(CDS.temp),start=2)
			cds.RC.f3         <- subseq(x=reverseComplement(CDS.temp),start=3)
			cds.frames        <- DNAStringSet(c(cds.f1,cds.f2,cds.f3,cds.RC.f1,cds.RC.f2,cds.RC.f3))
			
			aa.frames         <- suppressWarnings(translate(cds.frames,no.init.codon=T,if.fuzzy.codon="solve"))
			aa.frames         <- AAStringSet(gsub(pattern="\\*","X",aa.frames))
			aa.frames         <- subseq(aa.frames,start=1,end=(width(aa.frames)-1))  ### removes potential stop codon from search
			
			if(any(width(aa.frames)==0) | (max(width(aa.frames))<3)){
				CDS.test[identifier.temp,2] <- "short"   ### changes this from "yes" to "short"
				next
			}
			
			frame.search      <- str_locate(AA.temp,as.character(aa.frames))[,1]
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

############
### check.performance function
############
## measure how long it takes R to print "checking performance"
####
check.performance <- function(){
	start.time <- Sys.time()
	print("checking performance")
	end.time <- Sys.time()
	#execution.time <- end.time-start.time
	execution.time  <- as.numeric(difftime(time1=end.time,time2=start.time))
	result          <- paste(round((execution.time*1000),digits=3),"milliseconds")
	result
}

####
## collapse.low.support.nodes function
####
## This turns a (usually) bifurcating tree with support values into a multichotomous tree with low-support nodes collapsed
####
## Dependencies
##     packages: ape (required), phangorn (recommended)
####
## rooted.tree = a tree of class phylo, with $node.label indicating the support values
## support.threshold = minimum value required to keep the node (and branch leading to the node)
collapse.low.support.nodes <- function(rooted.tree,support.threshold=80){
	tree1.rooted <- rooted.tree
	tolerance                 <- min(tree1.rooted$edge.length)/2 ### this value is half the length of the shortest edge, and is used as the cuttoff for removing branches that will get set to zero-length
	ntips                     <- length(tree1.rooted$tip.label)
	tree1.nodes.to.drop       <- (ntips + which(as.numeric(tree1.rooted$node.label) <= support.threshold))  ### Also may want to drop nodes without a support value other than the root?
	tree1.crown.edges.to.drop <- which(tree1.rooted$edge[,2] %in% tree1.nodes.to.drop) ### the edges having the descendent node being a weakly supported node that we want to drop
	tree1.rooted2             <- tree1.rooted
	tree1.rooted2$edge.length[tree1.crown.edges.to.drop] <- 0           ### sets edge length to zero for the edges that you want to drop
	tree1.rooted3             <- di2multi(tree1.rooted2,tol=tolerance)  ### deletes the zero-length edges deleted (this is the tree with low-support nodes collapsed!!!)
	tree1.rooted3
}

############
## congruency.test function
####
## Dependencies
# packages = ape, phangorn
# functions = collapse.low.support.nodes
####
## Compares two trees (each of class phylo) and, given a minimum threshold for considering clades as highly supported, returns whether if two trees are congruent (TRUE) or not incongruent (FALSE).
## Use the function multiple.congruency.test if comparing more than two input trees (also works for comparing just two trees)
## Note: If two trees do not share ANY tips, they are treated as not congruent. I may change this in the future.
####
congruency.test <-function(tree1,tree2,min.support=80){
	if(class(tree1)=="multiPhylo"){
		tree1 <- tree1[[1]]
	}
	if(class(tree2)=="multiPhylo"){
		tree2 <- tree2[[1]]
	}
	
	shared.tips    <- intersect(tree1$tip.label,tree2$tip.label)
	
	if(length(shared.tips)<4){
		isCongruent <- NA
	} else {
		arbitrary.root <- intersect(tree1$tip.label,tree2$tip.label)[1]  ### just picks the first shared terminal branch as the root
		tree1.rooted <- root(tree1,arbitrary.root)
		tree2.rooted <- root(tree2,arbitrary.root)

		### collapsing low-support nodes
		tree1.collapsed <- collapse.low.support.nodes(tree1.rooted,support.threshold=min.support)
		tree2.collapsed <- collapse.low.support.nodes(tree2.rooted,support.threshold=min.support)
	
		### dropping tips that are not shared between trees
		shared.tips      <- intersect(tree1.collapsed$tip.label,tree2.collapsed$tip.label)
		tree1.drop.names <- setdiff(tree1.collapsed$tip.label,shared.tips)
		tree2.drop.names <- setdiff(tree2.collapsed$tip.label,shared.tips)
	
		tree1.trimmed <- drop.tip(tree1.collapsed,tree1.drop.names)
		tree2.trimmed <- drop.tip(tree2.collapsed,tree2.drop.names)
	
		### internal node numbers for the collapsed and trimmed trees
		tree1.NodesInternal <- (length(tree1.trimmed$tip.label)+1):max(tree1.trimmed$edge)
		tree2.NodesInternal <- (length(tree2.trimmed$tip.label)+1):max(tree2.trimmed$edge)
	
		### finding the set of clade-sets for each collapsed and trimmed trees
		tree1.cladeSets <- Descendants(tree1.trimmed,node=tree1.NodesInternal,type="tips")
		tree2.cladeSets <- Descendants(tree2.trimmed,node=tree2.NodesInternal,type="tips")
	
		tree1.cladeSets.names <- lapply(X=tree1.cladeSets,FUN=function(input){tree1.trimmed$tip.label[input]})
		tree2.cladeSets.names <- lapply(X=tree2.cladeSets,FUN=function(input){tree2.trimmed$tip.label[input]})
	
		### pairwise comparison of cladeSets to check if the trees are congruent
		pairwise.cladeSets <- list(); length(pairwise.cladeSets) <- length(tree1.cladeSets.names)
		for(i in 1:length(tree1.cladeSets)){
			pairwise.cladeSets[i]  <- all(unlist(lapply(X=tree2.cladeSets.names,FUN=function(input){all(input %in% tree1.cladeSets.names[[i]]) | all(tree1.cladeSets.names[[i]] %in% input) | all(input %notin% tree1.cladeSets.names[[i]])}))) ### checks if every cladeSet of tree1 exists as either (1) a matching set, subset, or does not intersect with each cladeSet of tree2.
		}
		pairwise.cladeSets <- unlist(pairwise.cladeSets)
		isCongruent        <- all(pairwise.cladeSets)
	}
	isCongruent
}

###########
## multiple.congruency.test function
#####
### performs the congruency.test function for every pairwise combination of input trees (of class multiphylo)
### output is a three-column matrix:
###		columns 1 and 2 contain the names of the two trees being compared, respectively
###		column 3 = a logical indicating whether the two trees (of columns 1 and 2) are congruent
#####
multiple.congruency.test <- function(...,min.support=80){
	list.of.trees <- list(...)
	if(length(list.of.trees)==1){
		list.of.trees <- list.of.trees[[1]]
	} else {
		for(i in seq(length(list.of.trees))){
			names(list.of.trees)[i]  <-  as.character(sys.call()[[i+1]])
		}
	}
	if(any(is.null(names(list.of.trees)))){
		names(list.of.trees) <- paste0("tree.name",seq(1:length(list.of.trees)))
	}
	
	list.of.treeNames       <- names(list.of.trees)
	for(i in 1:length(list.of.treeNames)){
		assign(list.of.treeNames[i],list.of.trees[i])
		#class(get(list.of.treeNames[i]))
	}
	pairwise.treesNames.mat <- xprod.combn.mat(list.of.treeNames,list.of.treeNames)
	pairwise.treesNames.mat <- pairwise.treesNames.mat[!duplicated(t(apply(pairwise.treesNames.mat,1,sort))),]
	congruence.result       <- apply(X=pairwise.treesNames.mat,MARGIN=1,FUN=function(input,support.thresh=min.support){tree1=get(input[1]);tree2=get(input[2]);congruency.test(tree1,tree2,min.support=support.thresh)})
	result                  <- cbind(pairwise.treesNames.mat,congruence.result)
	#result                  <- result[!duplicated(t(apply(result,1,sort))),]
	colnames(result) <- c("tree1","tree2","congruent")
	result
}

#########
## congruency.network function
#########
## Dependencies
##  packages:  tidyverse, igraph
##  functions: 
#####
## Description: Plot the result of multiple.congruency.test as a network, with each node representing a tree and each edge connecting two trees that are congruent
#####
## input: the three column matrix produced by the function multiple.congruency.test
#####
congruency.network <- function(congruency.matrix,plot.network=T){
	mat.congruent    <- congruency.matrix[which(congruency.matrix[,"tree1"] != congruency.matrix[,"tree2"] & congruency.matrix[,"congruent"]=="TRUE"),]
	edge_list        <- tibble::tibble(from = mat.congruent[,1], to = mat.congruent[,2])
	node_list        <- tibble::tibble(id = unique(congruency.matrix[,1]))
	result.igraph    <- igraph::graph_from_data_frame(d = edge_list, vertices = node_list, directed = FALSE)
	if(plot.network==T){
		plot(result.igraph)
	}
	result.igraph
}

##########
### unique.tiplabels function
##########
### returns the set of unique tip labels across a multiPhylo object
##########
unique.tiplabels <- function(trees.multiPhylo){
	label.search   <- str_locate(names(unlist(trees.multiPhylo)),pattern="tip.label")
	tip.labels.all <- unlist(trees.multiPhylo)[which(!is.na(label.search[,1]))]
	result         <- unique(tip.labels.all)
	result
}

##########
#### get.all.quad.trees function
##########
#### returns a multiPhylo object containing the set of all quad-trees (4-taxon subclades of a tree) for a set of tip.labels (by default, the tips labels of the input tree).
#### It is often useful to set alt.tips = unique.tiplabels(<set of gene trees>), where the set of gene trees are held in a multiPhylo object
##########
get.all.quad.trees <- function(ref.tree,alt.tips=NULL,support.scaler=1){
	best.tree   <- ref.tree
	if(is.null(alt.tips)){
		taxa = best.tree$tip.label
	} else {
		taxa = alt.tips
	}
	group.index   <- sample(x=c(1:4),size=length(taxa),replace=T)                                                                                 ### assigns an index number to each taxon (i.e., tip label) in taxa vector
	all.quads     <- xprod.combn(list(taxa[group.index==1],taxa[group.index==2],taxa[group.index==3],taxa[group.index==4]))                       ### returns a list of character vectors (each length four), and each is a unique combination of the tip labels
	all.quads.mat <- mat.strsplit(all.quads)                                                                                                      ### each row is a unique 4-taxon set (same as all.quads list, but held in a four column matrix)
	best.tree$edge.length <- rep(1,length(best.tree$edge.length))                                                                                 ### sets all edge lengths equal to 1, because this is a test of topology
	best.tree$node.label[which(best.tree$node.label !="")] <- as.numeric(best.tree$node.label[which(best.tree$node.label !="")])*(support.scaler) ### useful to scale the support values if planning to compare the quad.trees to a set of gene trees with different range of support values compared to ref.tree
	all.quads.counter = 0                                                                                                                         ### empty vector that will contain the set of quad trees within the best tree
	for(i in 1:nrow(all.quads.mat)){
		if(all(all.quads.mat[i,] %in% best.tree$tip.label)){
			drop.temp       <- best.tree$tip.label[!(best.tree$tip.label %in% all.quads.mat[i,])] #### all tips except for those in all.quads.mat[i,]
		} else {
			next
		}
		quad.temp       <- drop.tip(best.tree,drop.temp)
		all.quads.counter = (all.quads.counter +1)        ### increase the counter by 1 if didnt skip to next cycle of loop
		if(all.quads.counter==1){
			all.quads.trees <- c(quad.temp)
		} else {
			all.quads.trees <- c(all.quads.trees,quad.temp)
		}
	}
	all.quads.trees
}

			
##################################################
### Figuring out Functions used to pick UCEs after extracting best-matches... ##
##### Examining the starting set of UCEs. Since SnakeCap loci only selected from the Streicher Micrurus set, then no need to the read them into R
# 
# source.UCEs.dir          <- "~/StreicherWiens2017_UCEs_fasta/"   ### Directory containing fasta file of set of UCEs to be queried in genomes
# source.species           <- c("micrurus")                                                                            ### Micrurus fulvius (even though Streicher also found additional UCEs in python, I only queried the Micrurus set)
# source.species.filenames <- paste(source.UCEs.dir,source.species,"_UCEs.fa",sep="")                                  ### filenames for starting set of UCEs to be queried across species with genomes
# 
# for(i in 1:length(source.species)){                                                           #|reads in UCE alignments   #| Reads in one or more sets of UCEs that will
# 	assign(x=paste(source.species[i],"_UCEs",sep=""),value=readDNAStringSet(filepath=source.species.filenames[i]))          #| be queried within genomes.
# }                                                                                                                         #|

### Comparing the best-match UCEs found in each species genome:
#
#species                 <- c("Thamnophis_sirtalis","Crotalus_horridus","Protobothrops_mucrosquamatus","Ophiophagus_hannah","Vipera_berus","Crotalus_mitchellii","Pantherophis_guttatus","Python_bivittatus")
#UCEs.from.genomes.dir   <- "~/UCEs.In.Snake.Genomes/"
#UCEs.from.genomes.paths <- paste(UCEs.from.genomes.dir,species,"_UCEs.fasta",sep="")
#species.UCEs            <- paste(species,"_UCEs",sep="")
#
#for(i in 1:length(species.UCEs)){                                                              #|reads in UCE sets
#	assign(x=species.UCEs[i],value=readDNAStringSet(filepath=UCEs.from.genomes.paths[i]))      #|for each species
#}
#
#species.UCE.names       <- paste(paste(species,".UCE.names",sep=""))
#for(i in 1:length(species.UCE.names)){                                                            #|assigns a character vector of UCE names to an object for each species
#	assign(x=species.UCE.names[i],value=gsub(pattern=".*_uce","UCE",names(get(species.UCEs[i])))) #|the object names are held in the vector species.UCE.names
#}
#
#shared.UCEs            <- intersect.all(lapply(species.UCE.names,get))                        ### list of UCEs that were found in all of the snake genomes (n=2,968)
#sometimes.present.UCEs <- setdiff(unique(unlist(lapply(species.UCE.names,get))),shared.UCEs)  ### list of UCEs found in one or more but not all snake genomes (n=292)

####################################################################
### These are the alignments produced for the shared UCEs after using MAFFT 
#
# Crotalus.horridus_UCE-containing-contigs.filename <- "~/Crotalus-horridus_UCE-containing-contigs.fasta"
# Crotalus.horridus_UCEs                            <- "~/Crotalus_horridus_UCEs.fasta"
# UCE.alignments.dir                                <- "~/MAFFT-aligned-UCEs"
# UCE.alignment.filenames                           <- list.files(path=UCE.alignments.dir,full.names=T)
# UCE.shortnames                                    <- gsub(".fasta","",list.files(path=UCE.alignments.dir,full.names=F))
# UCE.shortnames                                    <- gsub("uce","UCE",UCE.shortnames)
# 
# for(i in 1:length(UCE.shortnames)){                                                           #|reads in UCE alignments
# 	assign(x=UCE.shortnames[i],value=readDNAStringSet(filepath=UCE.alignment.filenames[i]))   #|and assigns them
# }
# 
# test <- lapply(UCE.shortnames,get)
# test2 <- NULL
# for(i in 1:length(test)){
# 	test2[i] <- length(test[[i]])
# }
# all(test2==8))
#
#

###############			
			
			
##########
# library(ape)
# geneTrees.dir      <- "~/targets_withFlanking_no-trimAL/treefiles_withFlanking_no-trimAL/"
# trees              <- read.tree.multipleFiles(list.files(geneTrees.dir,full.names=T))
# trees              <- trees[1:100] ### just looking at the first 100 for now because there are a lot of trees
# taxa               <- unique.tiplabels(trees)  ### trees = the set of gene trees
# ref.tree           <- read.tree("~/gene-trees_5May2019/targets_withFlanking_no-trimAL/ASTRAL-tree_AllLoci_withFlanking_no-trimAL_Topology2.tre") ### unrooted
# all.quads.trees    <- get.all.quad.trees(ref.tree=ref.tree,alt.tips=unique.tiplabels(c(trees,ref.tree)),support.scaler=100)

#############
# Tests of congruence in R
# * Read IQ-TREE phylogeny into R (with support values)
# * Collapse clades into polytomies if support value is below threshold
# - repeat previous two steps for another tree
# * Find set of shared tips between the two "collapsed" trees
# * Trim each of the two collapsed trees to only include individuals shared between them
# * Test if congruent:
# 	The two trees are congruent if:
# 	• Every pairwise comparison of a clade-sets (i.e., the set of containing individuals of a clade) in Tree 1 and Tree 2 is one of the following:
# 		• the two clade-sets are equal
# 		• one of the two clade-sets is a subset of the other clade-set
# 		• the two clade-sets do not share any individuals
# * Test if equal:
# 	The two trees are equal if:
# 	• The trees are congruent, and,
# 	• Every clade-set in Tree 1 has an equal clade-set in Tree 2 (and vice versa)
# 

#tree1 <- read.tree(file="~/Tropidonophis/IQTREE_Tropidonophis/GeneTrees_PartitionedByCodon/rag1/rag1.tre")
#tree2 <- read.tree(file="~/Tropidonophis/IQTREE_Tropidonophis/GeneTrees_PartitionedByCodon/nt3/nt3.tre")
#
#### may be better to just pick any tip shared between trees
#tree1.rooted <- root(tree1,"Scaphiodontophis-annulatus_KU289943")
#tree2.rooted <- root(tree2,"Scaphiodontophis-annulatus_KU289943")
#
#### collapsing low-support nodes
#tree1.collapsed <- collapse.low.support.nodes(tree1.rooted,support.threshold=80)
#tree2.collapsed <- collapse.low.support.nodes(tree2.rooted,support.threshold=80)
#
#### dropping tips that are not shared between trees
#shared.tips      <- intersect(tree1.collapsed$tip.label,tree2.collapsed$tip.label)
#tree1.drop.names <- setdiff(tree1.collapsed$tip.label,shared.tips)
#tree2.drop.names <- setdiff(tree2.collapsed$tip.label,shared.tips)
#
#tree1.trimmed <- drop.tip(tree1.collapsed,tree1.drop.names)
#tree2.trimmed <- drop.tip(tree2.collapsed,tree2.drop.names)
#
#### internal node numbers for the collapsed and trimmed trees
#tree1.NodesInternal <- (length(tree1.trimmed$tip.label)+1):max(tree1.trimmed$edge)
#tree2.NodesInternal <- (length(tree2.trimmed$tip.label)+1):max(tree2.trimmed$edge)
#
#### finding the set of clade-sets for each collapsed and trimmed trees
#library(phangorn)
#tree1.cladeSets <- Descendants(tree1.trimmed,node=tree1.NodesInternal,type="tips")
#tree2.cladeSets <- Descendants(tree2.trimmed,node=tree2.NodesInternal,type="tips")
#
#tree1.cladeSets.names <- lapply(X=tree1.cladeSets,FUN=function(input){tree1.trimmed$tip.label[input]})
#tree2.cladeSets.names <- lapply(X=tree2.cladeSets,FUN=function(input){tree2.trimmed$tip.label[input]})
#
#### pairwise comparison of cladeSets to check if the trees are congruent
#tree1.cladeSets.in.tree2.cladeSets <- list(); length(tree1.cladeSets.in.tree2.cladeSets) <- length(tree1.cladeSets.names)
#for(i in 1:length(tree1.cladeSets)){
#	tree1.cladeSets.in.tree2.cladeSets[i] <- any(unlist(lapply(X=tree2.cladeSets.names,FUN=function(input){ any(c(all(tree1.cladeSets.names[[i]]  %in% input),(!any(tree1.cladeSets.names[[i]]  %in% input))))}))) ### checks if every cladeSet of tree1 exists as either (1) a set or subset of a cladeSet of tree2, or (2) does not intersect with any cladeSet of tree2.
#}
#tree2.cladeSets.in.tree1.cladeSets <- list(); length(tree2.cladeSets.in.tree1.cladeSets) <- length(tree2.cladeSets.names)
#for(i in 1:length(tree2.cladeSets)){
#	tree2.cladeSets.in.tree1.cladeSets[i] <- any(unlist(lapply(X=tree1.cladeSets.names,FUN=function(input){ any(c(all(tree2.cladeSets.names[[i]]  %in% input),(!any(tree2.cladeSets.names[[i]]  %in% input))))}))) ### checks if the ith tree2 cladeSet exists as a set or subset of any tree2 cladeSet
#}
#tree1.cladeSets.in.tree2.cladeSets <- unlist(tree1.cladeSets.in.tree2.cladeSets)
#tree2.cladeSets.in.tree1.cladeSets <- unlist(tree2.cladeSets.in.tree1.cladeSets)
#isCongruent <- all(tree1.cladeSets.in.tree2.cladeSets,tree2.cladeSets.in.tree1.cladeSets) ### TRUE means that the two trees are congruent (given the threshold level for strong support)
#isCongruent
