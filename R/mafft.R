#' MAFFT Multiple sequence alignment
#' 
#' Wrapper for MAFTT. Perform. multiple sequence alignment on DNA, RNA, or AA sequences.
#' This function was lightly modified from the mafft function written by Michael Hahsler and Anurag Nagar.
#'
#' @param x Object of class XStringSet (e.g., DNAStringSet) that contains the sequences to be aligned.
#' @param param Character string containing the mafft command line parameters.
#' @param mafft.path Either NULL (default) or a character string specifying the path to the mafft executable. If NULL, the mafft executable must be on your system PATH. This is very likely to be the case because by default mafft installs to usr/local/bin.
#' A possible exception is when mafft is used on a community cluster. In this case, specify the path to the mafft executable instead of using "auto".
#' @param return.as Of the following: "XStringSet" (default) or "MultipleAlignment"
#' @return Either an object of class XStringSet (if return.as= "XStringSet") or DNAMultipleAlignment (or RNA MultipleAlignment or AAMultipleAlignment) if return.as = MultipleAlignment (see BioStrings).
#' @export mafft
mafft <- function(x, param="--auto",mafft.path=NULL,return.as="XStringSet") {
    #######################################################################
    # BiostringsTools - Interfaces to several sequence alignment and
    # classification tools
    # Copyright (C) 2012 Michael Hahsler and Anurag Nagar
    #
    # This program is free software; you can redistribute it and/or modify
    # it under the terms of the GNU General Public License as published by
    # the Free Software Foundation; either version 2 of the License, or
    # any later version.
    #
    # This program is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    # GNU General Public License for more details.
    #
    # You should have received a copy of the GNU General Public License along
    # with this program; if not, write to the Free Software Foundation, Inc.,
    # 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
    ########################################################################

    #### Prepare the path to the executables. Jeff Weinell added this to avoid using the .findExecutable function.
    if(is.null(mafft.path)){
        blast.exe.path <- "mafft"
    } else {
        blast.exe.path <- mafft.path
    }
    ### Define temporary directory to work in.
    wd        <- tempdir()
    dir       <- getwd()
    temp_file <- basename(tempfile(tmpdir = wd))
    ### Tells R to remove the temporary files once the function is done running
    on.exit({
		file.remove(Sys.glob(paste(temp_file, ".*", sep="")))
		setwd(dir)
	})
    #Change working directory
    setwd(wd)
    ### Define name to use for the mafft input and output files
    infile   <- paste(temp_file, ".in", sep="")
    outfile  <- paste(temp_file, ".aln", sep="")
    ### Check that input data is the correct format. Stop and print a warning if not.
    if(!class(x) %in% c("RNAStringSet","DNAStringSet","AAStringSet")){
        stop("Unknown sequence type!")
    }
    ### Write the unaligned sequence data to a temporary file.
    Biostrings::writeXStringSet(x, infile, append=FALSE, format="fasta")
    ### Run MAFFT
    system(paste(blast.exe.path,param,"--clustalout",infile,">",outfile))
    ### Read the output of mafft
    if(is(x, "RNAStringSet")){
       result <- Biostrings::readRNAMultipleAlignment(outfile, format="clustal")
    }
    if(is(x, "DNAStringSet")){
       result <- Biostrings::readDNAMultipleAlignment(outfile, format="clustal")
    }
    if(is(x, "AAStringSet")){
       result <- Biostrings::readAAMultipleAlignment(outfile, format="clustal")
    }
    if(!(return.as %in% c("XStringSet", "MultipleAlignment"))){
        stop("return.as must be equal to 'XStringSet' or 'MultipleAlignment'")
    } else {
        if(return.as=="XStringSet"){
            if(is(x, "RNAStringSet")){
               result <- Biostrings::RNAStringSet(result)
            }
            if(is(x, "DNAStringSet")){
                result <- Biostrings::DNAStringSet(result)
            }
            if(is(x, "AAStringSet")){
                result <- Biostrings::AAStringSet(result)
            }
        }
    }
    ### Reset the sequence names, because mafft truncates them.
    names(result) <- names(x)
    result
}


