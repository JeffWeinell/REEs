#' Install MAFFT
#' 
#' This function installs MAFFT from source to the REEs directory or a user-specified location.
#' 
#' @param install.loc Path where MAFFT should be installed, or "auto" (Default), which installs MAFFT to blast-mafft subdirectory of REEs directory, or "PATH", which installs to "/usr/local" if it is writeable.
#' @return NULL. MAFFT is installed to install.loc
#' @export mafft.install
mafft.install <- function(install.loc="auto",source=F){
	### This is MAFFT version 7.475 (without extensions)
	source.url        <- "https://mafft.cbrc.jp/alignment/software/mafft-7.475-without-extensions-src.tgz"
	if(install.loc=="PATH"){
		install.loc       <- "/usr/local"
		if(is.writeable!=0){
				stop("install.loc is not writeable")
		} else {
			### download MAFFT source code
			download.file(url=source.url,destfile=paste0(install.loc,"/",basename(source.url)))
			# unpack tarball
			system(paste("cd",paste0("'",install.loc,"'"),"&& tar -xzvf",paste0("'",basename(source.url),"'")))
			# Define path to core directory
			core.dir <- list.dirs(install.loc)[grep("core$",list.dirs(install.loc))]
			### Compile and install
			system(paste("cd",core.dir,"&& make clean"))
			system(paste("cd",core.dir,"&& make"))
			system(paste("cd",core.dir,"&& make install"))
		}
	} else {
		if(install.loc=="auto"){
			install.loc       <- paste0(find.package("REEs"),"/blast-mafft/mafft")
			#dir.check.create(install.loc)
			dir.create(install.loc,recursive=T)
			bindir <- paste0(install.loc,"/bin")
			# dir.check.create(bindir)
			dir.create(bindir,recursive=T)
		}
		if(!(install.loc %in% c("auto","PATH"))) {
			### Check for write access to install.loc
			is.writeable <- file.access(install.loc,mode=2)
			if(is.writeable!=0){
				stop("install.loc is not writeable")
			}
			bindir <- paste0(install.loc,"/bin")
			# dir.check.create(bindir)
			dir.create(bindir,recursive=T)
		}
		### download MAFFT source code
		download.file(url=source.url,destfile=paste0(install.loc,"/",basename(source.url)))
		# unpack tarball
		system(paste("cd",paste0("'",install.loc,"'"),"&& tar -xzvf",paste0("'",basename(source.url),"'")))
		# Define path to core directory
		core.dir <- list.dirs(install.loc)[grep("core$",list.dirs(install.loc))]
		# read Makefile into R
		lines.makefile        <- readLines(paste0(list.dirs(install.loc)[grep("core$",list.dirs(install.loc))],"/Makefile"))
		lines.makefile.new    <- lines.makefile
		lines.makefile.new[1] <- paste0("PREFIX = ",install.loc) # gsub("/usr/local$",install.loc,lines.makefile[1])
		lines.makefile.new[3] <- paste0("BINDIR = ",bindir)
		### write the new lines to Makefile
		writeLines(lines.makefile.new,paste0(list.dirs(install.loc)[grep("core$",list.dirs(install.loc))],"/Makefile"))
		### Compile and install.
		system(paste("cd",core.dir,"&& make clean"))
		system(paste("cd",core.dir,"&& make"))
		system(paste("cd",core.dir,"&& make install"))
	}
	# delete tarball
	system(paste("rm -R",paste0(install.loc,"/",basename(source.url))))
}



