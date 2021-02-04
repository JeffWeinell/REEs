#' Remove Poorly Aligned Regions
#' 
#' A sliding window method to replace poorly aligned subsections of an alignment (these are replaced with a string of "-"). Additionally, individuals with less non-missing data than a threshold are removed.
#' Did I or Carl Hutter write this code?
#' 
#' @param input.align Input DNA alignment of class DNAStringSet or DNAMultipleAlignment
#' @param slice.size.bp Window size for the sliding window. Default is 100.
#' @param threshold Maximum fraction of sites that can be parsimony informative sites (mean of pairwise scores calculated for each individual) for the sequence to be considered aligned "good".
#' @param min.length Minimum length non-missing data to keep an individual in the output alignment.
#' @return DNA alignment (class DNAStringSet) with poorly aligned regions replaced with missing data.
#' @export slice.trim
slice.trim <- function(input.align, slice.size.bp = 100, threshold = 0.45,min.length=20){  
  #makes consensus sequence for comparison
  #input.align<-trimal.align
  input.con        <- make.consensus(input.align, method = "majority")   ### Majority-rule consensus sequence of input.align
  names(input.con) <- "Reference_Locus"                                  ### Name of the sequence
  
  comb.align  <- append(input.align, input.con)                          ### Alignment that includes input.align and the consensus reference locus
  
  #Gets slice information ready
  slice.no    <- ceiling(max(width(input.align))/slice.size.bp)          ### Number of slices
  slice.start <- 1
  slice.end   <- slice.size.bp
  
  #checks to see if its out of bounds
  if (slice.end > max(width(input.align))){ 
    slice.end   <-max(width(input.align))
  }#end if check
  output.align <- Biostrings::DNAStringSet()
  
  for (x in 1:slice.no){
    #Slice alignment into number of slices 
    sliced.align <- Biostrings::subseq(comb.align, start = slice.start, end = slice.end)

    #Checks for badly aligned sequences 
    bad.align    <- pairwise.inf.sites(sliced.align, "Reference_Locus")     ### mean pairwise percent of sites that are parsimony informative for each individual in the slice
    
    #Remove bad sequence chunks
    rem.seqs     <- bad.align[bad.align >= threshold]                       ### badly aligned sequences to remove (that have mean pairwise percent parsimony informative greater than the threshold)
    good.align   <- sliced.align[!names(sliced.align) %in% names(rem.seqs)] ### the sequences not considered badly aligned 
    
    #Makes replacement gap seqs for the bad ones
    blank.align <- Biostrings::DNAStringSet()
    if (length(rem.seqs) != 0){
      for (y in 1:length(rem.seqs)){
        blank.align<-append(blank.align, Biostrings::DNAStringSet(paste0(rep("-", slice.end-slice.start+1), collapse = "")) )
      }
      names(blank.align)<-names(rem.seqs)
    }#end rem seqs if
    
    #Saves the slices and cats
    save.slice          <- append(good.align, blank.align)
    save.slice          <- save.slice[order(names(save.slice))]
    save.names          <- names(save.slice)
    output.align        <- Biostrings::DNAStringSet(paste0(as.character(output.align), as.character(save.slice)))
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
        save.slice          <- Biostrings::subseq(comb.align, start = slice.start, end = slice.end)
        save.slice          <- save.slice[order(names(save.slice))]
        save.names          <- names(save.slice)
        output.align        <- Biostrings::DNAStringSet(paste0(as.character(output.align), as.character(save.slice))) ### concatenates previous slice with current slice
        names(output.align) <- save.names
		break
      }
    }#end if
  }#end x loop
  
  #Removes reference
  output.align <- output.align[names(output.align) != "Reference_Locus"]
  
  #removes taxa with too little non-missing data.
  str.splitted <- strsplit(as.character(output.align), "")
  x.align      <- as.matrix(ape::as.DNAbin(str.splitted))
  len.temp     <- as.character(as.list(x.align))
  len.loci     <- lapply(len.temp, function (x) x[x != "-"])
  spp.len      <- unlist(lapply(len.loci, function (x) length(x)))         ### Number of non-gap characters in each sequence
  spp.rem      <- spp.len[spp.len <= min.length]                           ### Individuals to remove from the alignment because they have to few non-gap characters
  return.align <- output.align[!names(output.align) %in% names(spp.rem)]   ### slice.trim alignment without spp.rem individuals
  return(return.align)
}





