#' Filter Feature Annotation Table
#' 
#' Filters a GFF-format feature annotation table to include only features of a particular feature type and minimum feature length.
#' The output of this function is a GFF file containing the filtered an annotation table, which can be used as input to the functions get.loci.from.annotationTable or get.exome.from.annotationTable
#' 
#' @param input.gff Either a path to input feature annotation table (gff format) or an object of class data.table. See "ref_Thamnophis_sirtalis-6.0_top_level_JLW.gff3" as an example input.
#' @param output.gff Where to write the output (filtered) annotation table. Default is NULL, in which case the output is not written. See "CDS_ref_Thamnophis_sirtalis-6.0_top_level_JLW_longer120bp.gff3" as an example output.
#' @param feature.type Type of sequences to keep. Options include all unique values in the "feature" column of input.gff. The usually include the following: "region", "gene", "mRNA", "exon", "CDS", "lnc_RNA", "cDNA_match", "transcript", "pseudogene", "tRNA", "V_gene_segment", "C_gene_segment".
#' @param min.length Minimum nucleotide length to include a feature in the output annotation table. Default is 120.
#' @param write.table If TRUE, then the filtered feature table is written to the value of 
#' @return Writes the filtered annotation table to the path indicated by the output.gff argument. Additionally returns a data.table object. 
#' @export filter.gff
filter.gff <- function(input.gff,output.gff=NULL,feature.type,min.length=0) {
#	if(is.character(input.gff)){
#		#unfiltered.gff   <- data.table::fread(input.gff)
#		unfiltered.gff    <- ape::read.gff(input.gff)
#	}
#	if("data.frame" %in% class(input.gff)){
#		unfiltered.gff   <- input.gff
#	}
	# Set x to be input.gff because load.gff has an argument also called input.gff
	x <- input.gff
	# Load the GFF as a data.table object and then coerce to a data.frame
	unfiltered.gff   <- as.data.frame(load.gff(x))
	##### Set column modes
	# Identify which columns are named "start" or "end", because there are the columnd that should be mode "numeric"
	numeric.columns <- which(colnames(unfiltered.gff) %in% c("start","end"))
	# Set mode to numeric for those columns that should be numeric
	unfiltered.gff[, numeric.columns] <- sapply(unfiltered.gff[, numeric.columns], as.numeric)
	# Get column indices for all of the columns other than the columns with names "start" and "end".
	character.columns <- which(!(colnames(unfiltered.gff) %in% c("start","end")))
	# Set mode to "character" for the columns indexed in the character.columns vector
	unfiltered.gff[, character.columns] <- sapply(unfiltered.gff[, character.columns], as.character)
	# Filter rows by feature.type argument. Column three holds the feature type and is usually named "type"
	filtered.gff1A   <- unfiltered.gff[which(unfiltered.gff[,3]==feature.type),]
	# Calculate feature lengths
	widths1A         <- ((abs(filtered.gff1A$start-filtered.gff1A$end))+1)
	# Filter rows by feature lengths and the min.length argument
	filtered.gff1B   <- filtered.gff1A[which(widths1A>=min.length),]
	if(!is.null(output.gff)){
		utils::write.table(x=filtered.gff1B,file=output.gff,sep="\t",row.names=F)
	}
	# Coerce to data table object because it can be printed safely.
	result <- data.table::as.data.table(filtered.gff1B)
	result
}
#' @examples
#' ### Load GFF table from NCBI repository.
#' Thamnophis.sirtalis_GFF.url <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.gff.gz"
#' Thamnophis.sirtalis_GFF     <- load.gff(input=Thamnophis.sirtalis_GFF.url,local=F)
#' 
#' # Filter Thamnophis.sirtalis_GFF to only include CDS features with length at least 120bp
#' Thamnophis.sirtalis_GFF_CDS_longer120bp <- filter.gff(input.gff=Thamnophis.sirtalis_GFF,feature.type="CDS",min.length=120) ## This object is now automatically included with REEs package
#' 