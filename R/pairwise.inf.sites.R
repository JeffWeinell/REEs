#' Mean Pairwise Informative Sites
#' 
#' This looks like something Carl wrote. The function is used in the function slice.trim, which might also be from Carl; slice.trim is not used in any of the other functions.
#' 
#' @param x input sequence
#' @param y second input sequence
#' @return 
#'
pairwise.inf.sites <- function(x, y) {
  temp.align <- strsplit(as.character(x), "")
  mat.align  <- lapply(temp.align, tolower)
  m.align    <- as.matrix(ape::as.DNAbin(mat.align))
  
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
