#' Run MAFFT
#' 
#' R wrapper for running MAFFT alignment software.
#' Need to check if there is any reason to use this function instead of the REEs::mafft
#' 
#' @param unaligned.contigs Input DNA sequences in an object of class DNAStringSet.
#' @param add.contigs Default is NULL. This argument not currently implemented in this function.
#' @param algorithm Alignment algorithm to use. Default is "localpair". Alternatives are "genafpair" and "globalpair". See MAFFT documentation for more details: https://mafft.cbrc.jp/alignment/software/manual/manual.html.
#' @param rev.dir Should MAFFT try to align reverse complement of sequences? Default is TRUE.
#' @param save.name Filename for DNA alignment. This file is temporary unless the argument delete.files is changed to FALSE.
#' @param threads How many threads should be used (default is 6)
#' @param delete.files If files should be deleted (default is T)
#' @return DNAStringSet object holding the aligned DNA sequences
run.mafft<-function(unaligned.contigs, add.contigs = NULL, algorithm = "localpair", rev.dir = T, save.name = NULL, threads = 6, delete.files = T){
  save.contigs<-as.list(as.character(unaligned.contigs))
  if (is.null(save.name) == T) {
       save.name <- paste(sample(LETTERS, 5, replace = T), collapse = "")
  }
  if (rev.dir == T){
  		adjust.direction <- "--adjustdirection"
  	} else {
  		adjust.direction <- ""
  	}
  
  #Adds a sequence into the alignment. Saves much computation.
  #if (algorithm == "add"){
  #  #Saves to folder to run with mafft
  #  write.fasta(sequences = save.contigs, names = names(save.contigs), paste0(save.name, ".fa"), nbchar = 1000000, as.string = T)
  #  
  #  #Saves to folder to run with mafft
  #  add.save<-as.list(as.character(add.contigs))
  #  write.fasta(sequences = add.save, names = names(add.save), "add_sequences.fa", nbchar = 1000000, as.string = T)
  #  
  #  #Runs MAFFT to align
  #  system(paste0("mafft --",algorithm, " add_sequences.fa ", adjust.direction, " --maxiterate 1000 ", save.name, ".fa > ", save.name, "_align.fa"), ignore.stderr = T)
  #  
  #  alignment<-scanFa(FaFile(paste(save.name, "_align.fa", sep = "")))   # loads up fasta file
  #  unlink(paste(save.name, ".fa", sep = ""))
  #  unlink("add_sequences.fa")    
  #}#end -add
  
  #Does Regular MAFFT Local Pair
  #if (algorithm == "localpair"){
    #Saves to folder to run with mafft
    seqinr::write.fasta(sequences = save.contigs, names = names(save.contigs), paste0(save.name, ".fa"), nbchar = 1000000, as.string = T)
    
    #Runs MAFFT to align
    system(paste0("mafft --",algorithm, " --maxiterate 1000 ", adjust.direction, " --quiet --op 3 --ep 0.123"," --thread ", threads, " ", save.name, ".fa > ", save.name, "_align.fa"))
    # Loads up fasta file
    alignment <- Rsamtools::scanFa(Rsamtools::FaFile(paste0(save.name, "_align.fa")))
    unlink(paste0(save.name, ".fa"))
  #}#end local pair
  
  if (delete.files == T){
    unlink(paste0(save.name, "_align.fa"))
    return(alignment)
  } else {
      return(alignment)
    }
}#function end
