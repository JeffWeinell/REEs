#' Install NCBI BLAST+
#' 
#' This function installs BLAST to one of several default directories or to user-specified directory.
#' 
#' @param install.loc Path to directory where source files should be downloaded and blast should be installed. Default is to install blast in a directory called called blast-mafft, which is in the REEs package directory.
#' @param source Whether or not blast should be built and installed from source. Default is FALSE, in which case the function installs from either linux, windows, or macOS tarball depending on the result of Sys.info()
#' @return NULL. BLAST is installed to install.loc be installed to 
#' @export blast.install
blast.install <- function(install.loc="auto",source=F){
	if(install.loc=="auto"){
			install.loc       <- paste0(find.package("REEs"),"/blast-mafft")
			if(!(dir.exists(install.loc))){
				dir.create(install.loc)
			}
	} else {
		is.writeable <- file.access(install.loc,mode=2)
		if(is.writeable!=0){
			stop("Directory defined by install.loc is not writeable.")
		}
	}
	download.dir      <- install.loc
	blast.index.url   <- "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+"
	file.prefix       <- "ncbi-blast-"
	latest.version    <- "2.11.0"
	file.extensions   <- c("+-1.src.rpm","+-1.src.rpm.md5","+-1.x86_64.rpm","+-1.x86_64.rpm.md5","+-src.tar.gz","+-src.tar.gz.md5","+-src.zip","+-src.zip.md5","+-win64.exe","+-win64.exe.md5","+-x64-linux.tar.gz","+-x64-linux.tar.gz.md5","+-x64-macosx.tar.gz","+-x64-macosx.tar.gz.md5","+-x64-win64.tar.gz","+-x64-win64.tar.gz.md5","+.dmg","+.dmg.md5")
	blast.urls        <- paste0(blast.index.url,"/",latest.version,"/",file.prefix,latest.version,file.extensions)
	linux.url         <- blast.urls[grep(pattern="+-x64-linux.tar.gz$",file.extensions)]
	windows.url       <- blast.urls[grep(pattern="+-x64-win64.tar.gz$",file.extensions)]
	macOS.url         <- blast.urls[grep(pattern="+-x64-macosx.tar.gz$",file.extensions)]
	source.url        <- blast.urls[grep(pattern="+-src.tar.gz$",file.extensions)]
	if(source){
		# download tarball
		download.file(url=source.url,destfile=paste0(download.dir,"/",basename(source.url)),method="wget")
		# unpack tarball
		system(paste("cd",paste0("'",download.dir,"'"),"&& tar -xzvf",paste0("'",basename(source.url),"'")))
		# move into c++ subdirectory and configure
		system(paste("cd",paste0("'",download.dir,"/",file.prefix,latest.version,"+-src/c++","'"),"&& ./configure"))
		# move into build subdirectory and then install
		system(paste("cd",paste0("'",download.dir,"/",file.prefix,latest.version,"+-src/c++/ReleaseMT/build","'"),"&& make all_r"))
		# delete tarball
		system(paste("rm -R",paste0(download.dir,"/",basename(source.url))))
	} else{
		if(!(Sys.info()["sysname"] %in% c("Linux","Darwin","Windows"))){
			stop("unrecognized operating system detected using Sys.info()")
		}
		if(Sys.info()["sysname"]=="Linux"){
			# download tarball
			download.file(url=linux.url,destfile=paste0(download.dir,"/",basename(linux.url)),method="wget")
			# unpack tarball
			system(paste("cd",paste0("'",download.dir,"'"),"&& tar -xzvf",paste0("'",basename(linux.url),"'")))
			# delete tarball
			system(paste("rm -R",paste0(download.dir,"/",basename(linux.url))))
		}
		if(Sys.info()["sysname"]=="Darwin"){
			# download tarball
			download.file(url=macOS.url,destfile=paste0(download.dir,"/",basename(macOS.url)),method="wget")
			# move to directory where tarball was downloaded and then unpack it
			system(paste("cd",paste0("'",download.dir,"'"),"&& tar -xzvf",paste0("'",basename(macOS.url),"'")))
			# delete tarball
			system(paste("rm -R",paste0(download.dir,"/",basename(macOS.url))))
		}
		if(Sys.info()["sysname"]=="Windows"){
			# download tarball
			download.file(url=windows.url,destfile=paste0(download.dir,"/",basename(windows.url)),method="wget")
			# move to directory where tarball was downloaded and then unpack it
			system(paste("cd",paste0("'",download.dir,"'"),"&& tar -xzvf",paste0("'",basename(windows.url),"'")))
			# delete tarball
			system(paste("rm -R",paste0(download.dir,"/",basename(windows.url))))
		}
		return(paste0(file.prefix,latest.version," installed to ",download.dir))
	}
}





