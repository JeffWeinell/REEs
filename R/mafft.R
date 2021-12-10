#' MAFFT Multiple sequence alignment
#' 
#' Wrapper for MAFTT. Perform. multiple sequence alignment on DNA, RNA, or AA sequences.
#' This function was lightly modified from the mafft function written by Michael Hahsler and Anurag Nagar.
#'
#' @param x Object of class XStringSet (e.g., DNAStringSet) that contains the sequences to be aligned.
#' @param param Character string containing the mafft command line parameters.
#' @param mafft.path Either "auto" (default), NULL, or a character string with full path to the mafft executable.
#' If "auto", then the mafft executable must be located somewhere within the blast-mafft/mafft/ directory of the REEs package library (i.e., the default install location if mafft was installed with mafft.install() function).
#' If NULL, the mafft executable must be on your system PATH.
#' A possible exception is when mafft is used on a community cluster. In this case, specify the path to the mafft executable instead of using "auto".
#' @param return.as Of the following: "XStringSet" (default) or "MultipleAlignment"
#' @return Either an object of class XStringSet (if return.as= "XStringSet") or DNAMultipleAlignment (or RNA MultipleAlignment or AAMultipleAlignment) if return.as = MultipleAlignment (see BioStrings).
#' @export mafft
mafft <- function(x, param="--auto",mafft.path="auto",return.as="XStringSet") {
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
	if(mafft.path=="auto"){
		### Path to "/blast-mafft/mafft" directory in REEs library
		REEs.mafft.dir <- file.path(find.package("REEs"),"blast-mafft","mafft")
		### Path to directory containing the mafft executable
		mafft.dir.path <- list.dirs(REEs.mafft.dir)[grep("bin$",list.dirs(REEs.mafft.dir))]
		### Path to the MAFFT executable
		mafft.exe.path <- paste0(mafft.dir.path,"/mafft")
	} else {
		if(is.null(mafft.path)){
			mafft.exe.path <- "mafft"
		} else {
			mafft.exe.path <- mafft.path
		}
	}
	
	#### Verify that mafft.exe.path is executable. Stop and warn if not.
	test.mafft.exe       <- check.if.executable(exe.path=mafft.exe.path)
	if(test.mafft.exe!=0){
		stop(paste0("'",mafft.exe.path," is not executable. Aborting. Run mafft.install() with argument defaults to install MAFFT to REEs package.'"))
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
	infile   <- paste0(temp_file, ".in")
	outfile  <- paste0(temp_file, ".aln")
	### Check that input data is the correct format. Stop and print a warning if not.
	if(!class(x) %in% c("RNAStringSet","DNAStringSet","AAStringSet")){
		stop("Unknown sequence type!")
	}
	### Write the unaligned sequence data to a temporary file.
	Biostrings::writeXStringSet(x, infile, append=FALSE, format="fasta")
	### Run MAFFT
#	system(paste(mafft.exe.path,param,"--clustalout",infile,">",outfile))
	system(sprintf("'%s' %s '%s' > '%s' ",mafft.exe.path,param,infile,outfile))
	### Read the output of mafft
	
	if(is(x, "RNAStringSet")){
	   #result <- Biostrings::readRNAMultipleAlignment(outfile, format="clustal")
	   result <- Biostrings::readRNAStringSet(file.path(wd,outfile))
	}
	if(is(x, "DNAStringSet")){
	   #result <- Biostrings::readDNAMultipleAlignment(outfile, format="clustal")
	   result <- Biostrings::readDNAStringSet(file.path(wd,outfile))
	}
	if(is(x, "AAStringSet")){
	   #result <- Biostrings::readAAMultipleAlignment(outfile, format="clustal")
	   result <- Biostrings::readAAStringSet(file.path(wd,outfile))
	}
	
	if(!(return.as %in% c("XStringSet", "MultipleAlignment"))){
		stop("return.as must be equal to 'XStringSet' or 'MultipleAlignment'")
	} else {
		#if(return.as=="XStringSet"){
		#	if(is(x, "RNAStringSet")){
		#	   result <- Biostrings::RNAStringSet(result)
		#	}
		#	if(is(x, "DNAStringSet")){
		#		result <- Biostrings::DNAStringSet(result)
		#	}
		#	if(is(x, "AAStringSet")){
		#		result <- Biostrings::AAStringSet(result)
		#	}
		#}
		if(return.as=="MultipleAlignment"){
			if(is(x, "RNAMultipleAlignment")){
			   result <- Biostrings::RNAMultipleAlignment(result)
			}
			if(is(x, "readDNAMultipleAlignment")){
				result <- Biostrings::readDNAMultipleAlignment(result)
			}
			if(is(x, "readAAMultipleAlignment")){
				result <- Biostrings::readAAMultipleAlignment(result)
			}
		}
	}
	### Reset the sequence names, because mafft truncates them.
	names(result) <- names(x)
	result
}


