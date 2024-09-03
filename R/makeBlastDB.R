#' @title Make Local Blast Database
#' 
#' Wrapper for the NCBI makeblastdb command. This function makes an NCBI blast database (which is needed to run blast against).
#' 
#' @param makeblastdb.path Path to makeblastdb (included in the bin folder of blast)
#' @param subject.path Path to the file containing the sequences to include in the database
#' @return Creates a blast database file.
#' @export makeBlastDB
makeBlastDB  <- function(makeblastdb.path,subject.path){
	system(sprintf("%s -n %s -parse_seqids -dbtype nucl -max_file_sz 4GB",makeblastdb.path,subject.path))
}
