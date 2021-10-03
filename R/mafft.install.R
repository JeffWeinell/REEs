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
	### Backup URL to use if the primary URL becomes unavailable: "https://osf.io/qtdfr/download"
	sourcename <- "mafft-7.475-without-extensions-src.tgz"
	if(install.loc=="PATH"){
		install.loc       <- "/usr/local"
		if(is.writeable!=0){
				stop("install.loc is not writeable")
		} else {
			### download MAFFT source code
			download.file(url=source.url,destfile=file.path(install.loc,sourcename))
			# unpack tarball
			#system(paste("cd",paste0("'",install.loc,"'"),"&& tar -xzvf",paste0("'",basename(source.url),"'")))
			system(sprintf("cd '%s' && tar -xzvf '%s'", install.loc, sourcename))
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
			dir.check.create(install.loc)
			bindir <- paste0(install.loc,"/bin")
			# dir.check.create(bindir)
			dir.check.create(bindir)
		}
		if(!(install.loc %in% c("auto","PATH"))) {
			### Check for write access to install.loc
			is.writeable <- file.access(install.loc,mode=2)
			if(is.writeable!=0){
				stop("install.loc is not writeable")
			}
			bindir <- paste0(install.loc,"/bin")
			# dir.check.create(bindir)
			dir.check.create(bindir)
		}
		### download MAFFT source code
		download.file(url=source.url,destfile=file.path(install.loc,sourcename))
		# unpack tarball
		# system(paste("cd",paste0("'",install.loc,"'"),"&& tar -xzvf",paste0("'",basename(source.url),"'")))
		system(sprintf("cd '%s' && tar -xzvf '%s'", install.loc, sourcename))
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
	#system(paste("rm -R",paste0(install.loc,"/",basename(source.url))))
	system(sprintf("rm -R '%s/%s'",install.loc,sourcename))
}

#' Initialize REEs
#' 
#' This function installs some dependency programs from source to a location where REEs knows where to find it.
#'
#' @return Currently NULL. Soon this function will return a data frame showing where REEs is looking for the programs it depends on.
#' @export initialize_REEs
initialize_REEs <- function(){
	mafft.install()
	blast.install()
	pblat.install()
	# cap3.install()
	# pblat.install() : wget -O 'icebert-pblat-e26bf6b.zip' 'https://osf.io/pu82t/download' --no-check-certificate
	# unzip 'icebert-pblat-e26bf6b.zip'
	# 
	# mafft: wget "https://mafft.cbrc.jp/alignment/software/mafft-7.475-without-extensions-src.tgz"
	# wget -O "mafft-7.475-without-extensions-src.tgz" "https://osf.io/qtdfr/download"
}


#' Install PBLAT
#' 
#' This function installs PBLAT from source to a directory where REEs knows to look
#' 
#' @return NULL.
#' @export pblat.install
pblat.install <- function(){
	source.url   <- "https://osf.io/pu82t/download"
	sourcename   <- "icebert-pblat-e26bf6b.zip"
	install.loc  <- file.path(find.package("REEs"),"pblat")
	# dir.check.create(install.loc)
	dir.check.create(install.loc)
	#bindir <- file.path(install.loc,"bin")
	# dir.check.create(bindir)
	#dir.check.create(bindir)
	# download pblat source code
	download.file(url=source.url,destfile=file.path(install.loc,sourcename))
	# unzip the source file
	system(sprintf("cd '%s' && unzip '%s'", install.loc, sourcename))
	# path to directory containing the makefile
	core.dir <- file.path(install.loc,gsub(".zip$","",sourcename))
	# install
	system(sprintf("cd '%s' && make",core.dir))
	# delete zipfile
	system(sprintf("rm -R '%s/%s'",install.loc,sourcename))
}





