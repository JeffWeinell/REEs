#' Left Justify Alignment
#' 
#' 
#' This function is used for asthetics only. Uses the linux commands "column" and "sed" to ensure that columns of the input alignment are lined up (left justified) in the text file
#' even if sequence names are of different lengths. NOTE: This function doesn't quite work yet because there is a limit on the character length of columns, and DNA seqs are often too long.
#' 
#' @param phylip.alignment phylip.alignment = filename or path to sequence alignment in sequential phylip format
#' @param nspaces Minimum number of spaces separating sequences from sequence names (default is 6), i.e., the number of spaces between the longest sequence name and the first alignment site.
#' @export
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
