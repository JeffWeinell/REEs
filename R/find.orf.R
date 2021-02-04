#' Find Open Reading Frames
#' 
#' Given an input DNA sequence, returns either a table of stop codon positions and frame relative to input sequence, or a table of open indicating the positions, lengths, and frames of open reading frames.
#' 
#' @param input.seq DNAStringSet object containing one sequence.
#' @param codons Whether or not the function should return a table with the start and stop positions and frame (relative to sequence) of all potential stop codons. Default is FALSE.
#' @param min.size Minimum open reading frame size to consider (default = 80).
#' @return A data frame holding one of the following:
#'  (1) Three column data frame indicating position and frame of each possible stop codons.
#'  (2) Four column data frame indicating position, nucleotide length, and frame relative to input sequence, of all open reading frames at least as long as specified by the min.size argument.
#' @export
find.orf <- function(input.seq, codons = F, min.size = 80){
  #Sets up data
  # input.seq<-trimmed[j]
  codon.table <- data.frame(Start = rep(0,6), End = rep(0,6), Frame = c("F1", "F2", "F3", "R1", "R2", "R3"))
  for.seq     <- as.character(input.seq)
  
  #Locations of stop codons in forward direction sequence.
  TAA <- Biostrings::matchPattern("TAA", for.seq)
  TGA <- Biostrings::matchPattern("TGA", for.seq)
  TAG <- Biostrings::matchPattern("TAG", for.seq)
  
  #Forward Frame 1 stop codons
  result1 <- TAA[(TAA@ranges@start+2) %% 3 == 0]
  result2 <- TGA[(TGA@ranges@start+2) %% 3 == 0]
  result3 <- TAG[(TAG@ranges@start+2) %% 3 == 0]
  
  starts <- c(result1@ranges@start, result2@ranges@start, result3@ranges@start)
  ends   <- c(result1@ranges@start+2, result2@ranges@start+2, result3@ranges@start+2)

  # Update codon.table with the frame 1 forward stop codons
  if (length(starts) != 0){
    codon.table <- codon.table[codon.table$Frame != "F1",]
    temp.table  <- data.frame(Start = starts, End = ends, Frame = "F1")
    codon.table <- rbind(codon.table, temp.table)
  } 
  
  #Forward Frame 2
  result1 <- TAA[(TAA@ranges@start+1) %% 3 == 0]
  result2 <- TGA[(TGA@ranges@start+1) %% 3 == 0]
  result3 <- TAG[(TAG@ranges@start+1) %% 3 == 0]
  
  starts <- c(result1@ranges@start-1, result2@ranges@start-1, result3@ranges@start-1)
  ends   <- c(result1@ranges@start+1, result2@ranges@start+1, result3@ranges@start+1)
  
  if (length(starts) != 0){
    codon.table <- codon.table[codon.table$Frame != "F2",]
    temp.table  <- data.frame(Start = starts, End = ends, Frame = "F2")
    codon.table <- rbind(codon.table, temp.table)
  }
  
  #Forward Frame 3
  result1 <- TAA[(TAA@ranges@start) %% 3 == 0]   
  result2 <- TGA[(TGA@ranges@start) %% 3 == 0]    
  result3 <- TAG[(TAG@ranges@start) %% 3 == 0]    
  
  starts <- c(result1@ranges@start-2, result2@ranges@start-2, result3@ranges@start-2)
  ends   <- c(result1@ranges@start, result2@ranges@start, result3@ranges@start)
  
  if (length(starts) != 0){
    codon.table <- codon.table[codon.table$Frame != "F3",]
    temp.table  <- data.frame(Start = starts, End = ends, Frame = "F3")
    codon.table <- rbind(codon.table, temp.table)
  }
  
  #Sets up data
  rev.seq <- as.character(Biostrings::reverseComplement(input.seq))
  
  #Gets codon stuff
  TAA <- matchPattern("TAA", rev.seq)
  TGA <- matchPattern("TGA", rev.seq)
  TAG <- matchPattern("TAG", rev.seq)
  
  #Rev Frame 1
  result1<-TAA[(TAA@ranges@start+2) %% 3 == 0]   
  result2<-TGA[(TGA@ranges@start+2) %% 3 == 0]    
  result3<-TAG[(TAG@ranges@start+2) %% 3 == 0]    
  
  starts<-c(result1@ranges@start, result2@ranges@start, result3@ranges@start)
  ends<-c(result1@ranges@start+2, result2@ranges@start+2, result3@ranges@start+2)
  if (length(starts) != 0){
    codon.table <- codon.table[codon.table$Frame != "R1",]
    temp.table  <- data.frame(Start = starts, End = ends, Frame = "R1")
    codon.table <- rbind(codon.table, temp.table)
  }
  
  #Rev Frame 2
  result1 <- TAA[(TAA@ranges@start+1) %% 3 == 0]   
  result2 <- TGA[(TGA@ranges@start+1) %% 3 == 0]    
  result3 <- TAG[(TAG@ranges@start+1) %% 3 == 0]    
  
  starts <- c(result1@ranges@start-1, result2@ranges@start-1, result3@ranges@start-1)
  ends   <- c(result1@ranges@start+1, result2@ranges@start+1, result3@ranges@start+1)
  if (length(starts) != 0){
    codon.table <- codon.table[codon.table$Frame != "R2",]
    temp.table  <- data.frame(Start = starts, End = ends, Frame = "R2")
    codon.table <- rbind(codon.table, temp.table)
  }
  
  #Rev Frame 3
  result1<-TAA[(TAA@ranges@start) %% 3 == 0]   
  result2<-TGA[(TGA@ranges@start) %% 3 == 0]    
  result3<-TAG[(TAG@ranges@start) %% 3 == 0]    
  
  starts<-c(result1@ranges@start-2, result2@ranges@start-2, result3@ranges@start-2)
  ends<-c(result1@ranges@start, result2@ranges@start, result3@ranges@start)
  if (length(starts) != 0){
    codon.table <- codon.table[codon.table$Frame != "R3",]
    temp.table  <- data.frame(Start = starts, End = ends, Frame = "R3")
    codon.table <- rbind(codon.table, temp.table)
  } #end if 
  
  if (codons == T) {
    return(codon.table)
  }
  
  if (codons == F) {  
    frames<-unique(codon.table$Frame)
    orf.frame<-data.frame()
    for (x in 1:length(frames)){
      temp.codon<-codon.table[codon.table$Frame == frames[x],]
      temp.codon<-temp.codon[order(temp.codon$Start),]
      
      if (temp.codon$Start[1] == 0){
        temp.start<- as.numeric(gsub("F|R", "", temp.codon$Frame))
        add.frame <- data.frame(FrameStart = temp.start, FrameEnd = width(input.seq), Size = (width(input.seq)-temp.start)+1, Frame = frames[x])
        orf.frame <- rbind(orf.frame, add.frame)
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
}#
