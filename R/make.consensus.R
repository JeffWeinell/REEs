#' Get Consensus Sequence
#' 
#' This function returns the consensus sequence of an input DNA alignment; it is a slightly modified version of the consensus() function in the seqinr package. The only difference is that this function can additionally take as input alignments of class DNAStringSet and DNAMultipleAlignments, and the result is an alignment of class DNAStringSet.
#' 
#' @param  input.alignment Input DNA alignment of class DNAStringSet or DNAMultipleAlignment
#' @param  method One of the following methods to determine consensus sequences: "majority","threshold","IUPAC","profile". These are defined in seqinr package.
#' @param  threshold When method = threshold, this is the minimum relative frequency for a character to be returned as the consensus character (default = 0.6).
#' @param  warn.non.IUPAC Warn if sequence data includes non-IUPAC standard characters. Default is FALSE.
#' @param  type Only available option is "DNA"
#' @return DNAStringSet object containing the consensus sequence of the input DNA alignment
#' @export 
make.consensus <- function (input.alignment, method = c("majority","threshold","IUPAC","profile"), threshold = 0.6, warn.non.IUPAC = FALSE, type = "DNA") {
  
  if(class(input.alignment) %in% c("DNAStringSet","DNAMultipleAlignment")){
    new.align <- strsplit(as.character(input.alignment), "")                           ### input.alignment converted to a list of vectors
    align.in  <- matrix(unlist(new.align), ncol = length(new.align[[1]]), byrow = T)   ### input.alignment converted to a character matrix
  }

  if(class(input.alignment)=="alignment"){
    align.in <- input.alignment
  }

  ## Internal function to apply when method is "majority"
  ## Returns the most common character in a character vector.
  majority.fun <- function(x) {
      names(which.max(table(x)))
  }

  #Does based on method
  method <- match.arg(method) ### This only works when excuting the function (i.e., not if debugging)
  
  if (method == "IUPAC") {
    type <- match.arg(type)
    res  <- apply(X = align.in, MARGIN = 2, FUN = REEs::bma, warn.non.IUPAC = warn.non.IUPAC, type = type)
    names(res) <- NULL
  }
  
  if (method == "majority") {
    res <- apply(X=align.in, MARGIN=2, FUN=majority.fun)
    names(res) <- NULL
  }
  
  if (method == "profile") {
    obsvalue <- levels(factor(align.in))
    nrow     <- length(obsvalue)
    row.names(align.in) <- NULL
    res <- apply(X=align.in, MARGIN=2, FUN=function(x){table(factor(x, levels = obsvalue))})
  }
  
  if (method == "threshold") {
    profile    <- consensus(align.in, method = "profile")
    profile.rf <- apply(X=profile, MARGIN=2, FUN=function(x) {x/sum(x)})
    res        <- rownames(profile.rf)[apply(X=profile.rf, MARGIN=2, FUN=which.max)]
    res        <- ifelse(apply(X=profile.rf, MARGIN=2, FUN=max) >= threshold, res, NA)
    names(res) <- NULL
  }
  
  out.consensus<-Biostrings::DNAStringSet(paste0(res, collapse = ""))
  names(out.consensus) <- "Consensus_Sequence"
  return(out.consensus)
}
