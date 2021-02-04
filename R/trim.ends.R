#' Trim ends of DNA alignment
#' 
#' This function is PROBABLY used to create an alignment for input into a function that translates CDS DNA.
#' For other (most) applications, the functions "trim.alignment" or "filter.alignment" should be used instead.
#' Need to look further into when/if this function is used in my REEs pipeline.
#' Not sure if I wrote this or if this is Carl Hutter code.
#' 
#' @param x DNA alignment held in an object of class DNAStringSet.
#' @param min.n.seq Minimum number of individuals with non-missing data at site required to keep the site (default is 4).
#' @param codon.trim Trim alignment so that the first site is codon 1 and last site is codon 3. (default is TRUE).
#' @return Trimmed alignment of class DNAStringSet
#' @export
trim.ends <- function (x, min.n.seq = 4, codon.trim = T){
  new.align  <- strsplit(as.character(x), "")    ### Stores the input alignment as a list of character vectors.
  mat.align  <- lapply(new.align, tolower)       ### Same as new.align but with lowercase DNA characters.
  x          <- as.matrix(as.DNAbin(mat.align))  ### DNAbin format of mat.align.

  if (!inherits(x, "DNAbin")) {
    stop("'x' is not of class 'DNAbin'")
  }
  if (!is.matrix(x)) {
    stop("'x' must be a matrix")
  }
  
  replaceWithN <- function(x) {
    id <- (x == as.raw(4))    ### A matrix of logicals, having the same size as x. TRUE for each gap character, FALSE for each non-gap character.
    ## Tests if the first or last alignment column of has any gaps
    if (length(id) > 0 & any(id[c(1, length(id))])) {
      id <- which(id)   ### a vector containing the element index/numbers of gap characters of the alignment
      
      ### Internal function not used in any other function.
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
  
  ### Internal function not used in any other REEs package function.
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
    for (q in 1:length(x3[[y]])){
      if (x3[[y]][q] == "n"){
        x3[[y]][q]<-"-"
        } else {
          break
          }
        }
    for (q in length(x3[[y]]):1){
     if (x3[[y]][q] == "n"){
        x3[[y]][q]<-"-"
      } else {
          break
      }
    }
  }#end x loop
  #Saves final stuff
  temp.align <- lapply(x3, FUN = function(x) paste(x, collapse = ""))
  align.out  <- Biostrings::DNAStringSet(unlist(temp.align))
  names(align.out) <- save.names
  return(align.out)
}
