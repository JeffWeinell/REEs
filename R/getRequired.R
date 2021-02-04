#' State Required Packages
#' 
#' Get a list of all packages required (imports and depends) for a vector of packages.
#' Package versions are ignored.
#'
#' @param packages Character vector of package names
#' @param include.input Should the output vector also include the input packages? Default FALSE.
#' @param include.Rcore Should the output vector include R-core packages? Default TRUE.
#' @return Character vector of names of packages required (directly or indirectly) by the input set of packages.
#' @export required.packages
required.packages <- function(packages,include.input=F,include.Rcore=T){
	# Rcore.packages <- c("foreign","base","compiler","datasets","graphics","grDevices","grid","nlme","methods","parallel","splines","stats","stats4","tcltk","tools","utils","Matrix")
	Rcore.packages <- names(which(!is.na(installed.packages()[,"Priority"])))
	required <- packages
	length.temp <- 0
	while(length.temp  < length(required)){
		length.temp    < length(required)
		length.temp    <- length(required)
		imports.temp   <- gsub(" .+","",unlist(lapply(required,FUN=function(x){strsplit(mgsub(c("\n",", +"),c("",","),packageDescription(x,fields="Imports")),split=",")})))
		depends.temp   <- gsub(" .+","",unlist(lapply(required,FUN=function(x){strsplit(mgsub(c("\n",", +"),c("",","),packageDescription(x,fields="Depends")),split=",")})))
		linkingto.temp <- gsub(" .+","",unlist(lapply(required,FUN=function(x){strsplit(mgsub(c("\n",", +"),c("",","),packageDescription(x,fields="LinkingTo")),split=",")})))
		required.temp  <- gsub("\\(.+","",c(required,imports.temp,depends.temp,linkingto.temp))
		required       <- unique(required.temp[which(required.temp!="R")])
	#	required      <- setdiff(required.temp,c("R"))

	}
	if(include.input==F){
		required <- setdiff(required,packages)
	}
	if(include.Rcore==F){
		required <- setdiff(required,Rcore.packages)
	}
	required <- sort(required)
	required
}

#REEs.required        <- required.packages(packages="REEs",include.input=F,include.Rcore=T)
#### names of R-core packages
#base                 <- intersect(REEs.required,Rcore.packages)
#### names of non R-core packages
#remainder            <- setdiff(REEs.required,Rcore.packages)
#### For each non R-core package, a vector of package names imported, dependent, and required (imported or dependent), respectively
#imports.list         <- lapply(remainder,FUN=function(x){mgsub(c(" .+","\\(.+"),c("",""),unlist(strsplit(mgsub(c("\n",", +"),c("",","),packageDescription(x,fields="Imports")),split=",")))})
#depends.list         <- lapply(remainder,FUN=function(x){mgsub(c(" .+","\\(.+"),c("",""),unlist(strsplit(mgsub(c("\n",", +"),c("",","),packageDescription(x,fields="Depends")),split=",")))})
#required.list        <- mapply(FUN=function(A,B){result.temp=unique(c(A,B));result = result.temp[which(result.temp!="R")];result},A=imports.list,B=depends.list)
#names(required.list) <- remainder
#
#### This is the list of non R-core packages that should be loaded first
#group1               <- names(required.list)[which(sapply(required.list,FUN=function(x){length(setdiff(x,base))})==0)]
#### Remaining packages to be loaded
#remainder2           <- setdiff(REEs.required,c(Rcore.packages,group1))
#### For each remainder2 package, a vector of package names required (imported or dependent)
#required.list2       <- required.list[which(names(required.list) %in% remainder2)]
#### This is the list of non R-core packages that should be loaded second
#group2               <- names(required.list2)[which(sapply(required.list2,FUN=function(x){length(setdiff(x,c(Rcore.packages,group1)))})==0)]
#### Remaining packages to be loaded
#remainder3           <- setdiff(REEs.required,c(Rcore.packages,group1,group2))
#### For each remainder2 package, a vector of package names required (imported or dependent)
#required.list3       <- required.list[which(names(required.list) %in% remainder3)]
#### This is the list of non R-core packages that should be loaded second
#group3               <- names(required.list3)[which(sapply(required.list3,FUN=function(x){length(setdiff(x,c(Rcore.packages,group1,group2)))})==0)]
#
###########

#' State Required Packages
#' 
#' Get a list of all packages required (imports and depends) for a vector of packages.
#' Package versions are ignored.
#'
#' @param pkgs Character vector of package names.
#' @param search.for.more Should the output include dependencies of dependencies, etc.? Default TRUE.
#' @param order.by.repos Not yet implemented.
#' @return A list of character vectors with the package names. Index of list indicates relative order of groups of packages that should be installed or loaded.
#' @export package.load.order
package.load.order   <- function(pkgs,search.for.more=T,order.by.repos=F){
	if(order.by.repos){
		warning("order.by.repos argument not yet implemented; output is not ordereded by package repository")
	}
	# Rcore.packages <- c("foreign","base","compiler","datasets","graphics","grDevices","grid","nlme","methods","parallel","splines","stats","stats4","tcltk","tools","utils","Matrix")
	### Character vector containing the names of the packages that were automatically installed when R was installed.
	Rcore.packages <- names(which(!is.na(installed.packages()[,"Priority"])))
	if(search.for.more){
		pkgs.all         <- required.packages(packages=pkgs,include.input=T,include.Rcore=T)
	} else {
		pkgs.all         <- pkgs
	}
	
	i=1
	group.temp           <- intersect(pkgs.all,Rcore.packages)
	remainder.temp       <- setdiff(pkgs.all,Rcore.packages)
	#result               <- list(NULL)
	result               <- list(group.temp)
	#result[[i]]          <- group.temp
	imports.list         <- lapply(remainder.temp,FUN=function(x){mgsub(c(" .+","\\(.+"),c("",""),unlist(strsplit(mgsub(c("\n",", +"),c("",","),packageDescription(x,fields="Imports")),split=",")))})
	depends.list         <- lapply(remainder.temp,FUN=function(x){mgsub(c(" .+","\\(.+"),c("",""),unlist(strsplit(mgsub(c("\n",", +"),c("",","),packageDescription(x,fields="Depends")),split=",")))})
	link.list            <- lapply(remainder.temp,FUN=function(x){mgsub(c(" .+","\\(.+"),c("",""),unlist(strsplit(mgsub(c("\n",", +"),c("",","),packageDescription(x,fields="LinkingTo")),split=",")))})
	required.list        <- mapply(FUN=function(A,B,C){D=unique(c(A,B,C));D[which(D!="R")]},A=imports.list,B=depends.list,C=link.list)
	names(required.list) <- remainder.temp
	required.list.temp   <- required.list
	
	while(length(remainder.temp)>0 & length(group.temp)>0){
		### This is the list of non R-core packages that should be loaded first
		i=i+1
		result             <- c(result,list(NULL))
		group.temp         <-  names(required.list.temp)[which(sapply(required.list.temp,FUN=function(x){length(setdiff(x,unlist(result)))})==0)]
		result[[i]]        <- group.temp
		remainder.temp     <- setdiff(pkgs.all,unlist(result))
		required.list.temp <- required.list.temp[which(names(required.list.temp) %in% remainder.temp)]
	}
	result
}
#' @example
#' ### The output of this example includes R-core packages that should already be loaded
#' ### The second vector includes packages that only require R-core packages to be loaded
#' ### All packages in the ith vector should be loaded before loaded any packages in vectors i+n
#' package.load.order(pkgs="REEs")
