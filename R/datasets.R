#' REEs Package Datasets
#' 
#' This function is used to return data associated with the REEs package
#' 
#' @param x A character string, character vector, or NULL (default). If NULL, a character vector of possible values to use for this argument is returned.
#' @return Character vector containing the URLs of the fasta-formatted genomes for the species used to design REEs in the SnakeCap study. The names attribute holds the corresponding species names.
#' @export datasets
datasets <- function(x=NULL){
	
	#### URLs to some genomes (mostly snakes) available from NCBI genome database.
	Thamnophis.sirtalis_genome.url                  <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.fna.gz"
	Anolis.carolinensis_genome.url                  <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/090/745/GCF_000090745.1_AnoCar2.0/GCF_000090745.1_AnoCar2.0_genomic.fna.gz"
	Gekko.japonicus_genome.url                      <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/447/785/GCF_001447785.1_Gekko_japonicus_V1.1/GCF_001447785.1_Gekko_japonicus_V1.1_genomic.fna.gz"
	Pogona.vitticeps_genome.url                     <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/067/755/GCF_900067755.1_pvi1.1/GCF_900067755.1_pvi1.1_genomic.fna.gz"
	Crotalus.horridus_genome.url                    <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/625/485/GCA_001625485.1_ASM162548v1/GCA_001625485.1_ASM162548v1_genomic.fna.gz"
	Crotalus.mitchellii_genome.url                  <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/737/285/GCA_000737285.1_CrotMitch1.0/GCA_000737285.1_CrotMitch1.0_genomic.fna.gz"
	Ophiophagus.hannah_genome.url                   <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/516/915/GCA_000516915.1_OphHan1.0/GCA_000516915.1_OphHan1.0_genomic.fna.gz"
	Pantherophis.guttatus_genome.url_JTLQ00000000.2 <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/185/365/GCF_001185365.1_UNIGE_PanGut_3.0/GCF_001185365.1_UNIGE_PanGut_3.0_genomic.fna.gz"
	Protobothrops.mucrosquamatus_genome.url         <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/527/695/GCF_001527695.2_P.Mucros_1.0/GCF_001527695.2_P.Mucros_1.0_genomic.fna.gz"
	Python.bivittatus_genome.url                    <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/186/305/GCF_000186305.1_Python_molurus_bivittatus-5.0.2/GCF_000186305.1_Python_molurus_bivittatus-5.0.2_genomic.fna.gz"
	Vipera.berus_genome.url                         <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/800/605/GCA_000800605.1_Vber.be_1.0/GCA_000800605.1_Vber.be_1.0_genomic.fna.gz"
	Thermophis.baileyi_genome.url                   <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/457/575/GCA_003457575.1_DSBC_Tbai_1.0/GCA_003457575.1_DSBC_Tbai_1.0_genomic.fna.gz"
	Crotalus.viridis_genome.url                     <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/400/415/GCA_003400415.2_UTA_CroVir_3.0/GCA_003400415.2_UTA_CroVir_3.0_genomic.fna.gz"
	Pantherophis.obsoletus_genome.url               <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/654/085/GCA_012654085.1_UNIGE_PanObs_1.0/GCA_012654085.1_UNIGE_PanObs_1.0_genomic.fna.gz"
	Hydrophis.melanocephalus_genome.url             <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/320/005/GCA_004320005.1_hydMel_1.0/GCA_004320005.1_hydMel_1.0_genomic.fna.gz"
	Emydocephalus.ijimae_genome.url                 <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/319/985/GCA_004319985.1_emyIji_1.0/GCA_004319985.1_emyIji_1.0_genomic.fna.gz"
	Hydrophis.hardwickii_genome.url                 <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/023/765/GCA_004023765.1_ASM402376v1/GCA_004023765.1_ASM402376v1_genomic.fna.gz"
	Hydrophis.cyanocinctus_genome.url               <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/023/725/GCA_004023725.1_ASM402372v1/GCA_004023725.1_ASM402372v1_genomic.fna.gz"
	Pseudonaja.textilis_genome.url                  <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/518/735/GCF_900518735.1_EBS10Xv2-PRI/GCF_900518735.1_EBS10Xv2-PRI_genomic.fna.gz"
	Laticauda.laticaudata_genome.url                <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/320/025/GCA_004320025.1_latLat_1.0/GCA_004320025.1_latLat_1.0_genomic.fna.gz"
	Ptyas.mucosa_genome.url                         <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/654/045/GCA_012654045.1_UNIGE_Pmuc_v1.0/GCA_012654045.1_UNIGE_Pmuc_v1.0_genomic.fna.gz"
	Notechis.scutatus_genome.url                    <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/518/725/GCF_900518725.1_TS10Xv2-PRI/GCF_900518725.1_TS10Xv2-PRI_genomic.fna.gz"
	Thamnophis.elegans_genome.url                   <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/769/535/GCF_009769535.1_rThaEle1.pri/GCF_009769535.1_rThaEle1.pri_genomic.fna.gz"
	Naja.naja_genome.url                            <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/733/165/GCA_009733165.1_Nana_v5/GCA_009733165.1_Nana_v5_genomic.fna.gz"
	Podarcis_muralis_genome.url                     <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/329/235/GCF_004329235.1_PodMur_1.0/GCF_004329235.1_PodMur_1.0_genomic.fna.gz"
	Aspidoscelis_marmoratus_genome.url              <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/337/955/GCA_014337955.1_AspMar1.0/GCA_014337955.1_AspMar1.0_genomic.fna.gz"
	Varanus_komodoensis_genome.url                  <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/798/865/GCA_004798865.1_ASM479886v1/GCA_004798865.1_ASM479886v1_genomic.fna.gz"
	Salvator_merianae_genome.url                    <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/586/115/GCA_003586115.2_HLtupMer6/GCA_003586115.2_HLtupMer6_genomic.fna.gz"
	Lacerta_viridis_genome.url                      <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/245/905/GCA_900245905.1_ASM90024590v1/GCA_900245905.1_ASM90024590v1_genomic.fna.gz"
	Lacerta_bilineata_genome.url                    <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/245/895/GCA_900245895.1_L._bilineata_genome_assembly/GCA_900245895.1_L._bilineata_genome_assembly_genomic.fna.gz"
	Paroedura_picta_genome.url                      <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/118/565/GCA_003118565.1_Ppicta_assembly_v1/GCA_003118565.1_Ppicta_assembly_v1_genomic.fna.gz"
	Zootoca_vivipara_genome.url                     <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/800/845/GCF_011800845.1_UG_Zviv_1/GCF_011800845.1_UG_Zviv_1_genomic.fna.gz"
	Lacerta_agilis_genome.url                       <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/819/535/GCF_009819535.1_rLacAgi1.pri/GCF_009819535.1_rLacAgi1.pri_genomic.fna.gz"
	Pantherophis.guttatus_genome.url_JTLQ00000000.1 <- "https://osf.io/r2xa3/download"

	#### GFF feature tables associated with annotated genomes.
	Thamnophis.sirtalis_genome.GFF.url                   <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/077/635/GCF_001077635.1_Thamnophis_sirtalis-6.0/GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.gff.gz"
	Anolis.carolinensis_genome.GFF.url                   <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/090/745/GCF_000090745.1_AnoCar2.0/GCF_000090745.1_AnoCar2.0_genomic.gff.gz"
	Gekko.japonicus_genome.GFF.url                       <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/447/785/GCF_001447785.1_Gekko_japonicus_V1.1/GCF_001447785.1_Gekko_japonicus_V1.1_genomic.gff.gz"
	Pogona.vitticeps_genome.GFF.url                      <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/067/755/GCF_900067755.1_pvi1.1/GCF_900067755.1_pvi1.1_genomic.gff.gz"
	Crotalus.horridus_genome.GFF.url                     <- ""
	Crotalus.mitchellii_genome.GFF.url                   <- ""
	Ophiophagus.hannah_genome.GFF.url                    <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/516/915/GCA_000516915.1_OphHan1.0/GCA_000516915.1_OphHan1.0_genomic.gff.gz"
	Pantherophis.guttatus_genome.GFF.url_JTLQ00000000.2  <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/185/365/GCF_001185365.1_UNIGE_PanGut_3.0/GCF_001185365.1_UNIGE_PanGut_3.0_genomic.gff.gz"
	Protobothrops.mucrosquamatus_genome.GFF.url          <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/527/695/GCF_001527695.2_P.Mucros_1.0/GCF_001527695.2_P.Mucros_1.0_genomic.gff.gz"
	Python.bivittatus_genome.GFF.url                     <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/186/305/GCF_000186305.1_Python_molurus_bivittatus-5.0.2/GCF_000186305.1_Python_molurus_bivittatus-5.0.2_genomic.gff.gz"
	Vipera.berus_genome.GFF.url                          <- ""
	Thermophis.baileyi_genome.GFF.url                    <- ""
	Crotalus.viridis_genome.GFF.url                      <- ""
	Pantherophis.obsoletus_genome.GFF.url                <- ""
	Hydrophis.melanocephalus_genome.GFF.url              <- ""
	Emydocephalus.ijimae_genome.GFF.url                  <- ""
	Hydrophis.hardwickii_genome.GFF.url                  <- ""
	Hydrophis.cyanocinctus_genome.GFF.url                <- ""
	Pseudonaja.textilis_genome.GFF.url                   <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/518/735/GCF_900518735.1_EBS10Xv2-PRI/GCF_900518735.1_EBS10Xv2-PRI_genomic.gff.gz"
	Laticauda.laticaudata_genome.GFF.url                 <- ""
	Ptyas.mucosa_genome.GFF.url                          <- ""
	Notechis.scutatus_genome.GFF.url                     <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/518/725/GCF_900518725.1_TS10Xv2-PRI/GCF_900518725.1_TS10Xv2-PRI_genomic.gff.gz"
	Thamnophis.elegans_genome.GFF.url                    <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/769/535/GCF_009769535.1_rThaEle1.pri/GCF_009769535.1_rThaEle1.pri_genomic.gff.gz"
	Naja.naja_genome.GFF.url                             <- ""
	Podarcis_muralis_genome.GFF.url                      <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/329/235/GCF_004329235.1_PodMur_1.0/GCF_004329235.1_PodMur_1.0_genomic.gff.gz"
	Aspidoscelis_marmoratus_genome.GFF.url               <- ""
	Varanus_komodoensis_genome.GFF.url                   <- ""
	Salvator_merianae_genome.GFF.url                     <- ""
	Lacerta_viridis_genome.GFF.url                       <- ""
	Lacerta_bilineata_genome.GFF.url                     <- ""
	Paroedura_picta_genome.GFF.url                       <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/118/565/GCA_003118565.1_Ppicta_assembly_v1/GCA_003118565.1_Ppicta_assembly_v1_genomic.gff.gz"
	Zootoca_vivipara_genome.GFF.url                      <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/800/845/GCF_011800845.1_UG_Zviv_1/GCF_011800845.1_UG_Zviv_1_genomic.gff.gz"
	Lacerta_agilis_genome.GFF.url                        <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/819/535/GCF_009819535.1_rLacAgi1.pri/GCF_009819535.1_rLacAgi1.pri_genomic.gff.gz"
	Pantherophis.guttatus_genome.GFF.url_JTLQ00000000.1  <- ""

	genome.url      <- c(Thamnophis.sirtalis_genome.url,Anolis.carolinensis_genome.url,Gekko.japonicus_genome.url,Pogona.vitticeps_genome.url,Crotalus.horridus_genome.url,Crotalus.mitchellii_genome.url,Ophiophagus.hannah_genome.url,Pantherophis.guttatus_genome.url_JTLQ00000000.2,Protobothrops.mucrosquamatus_genome.url,Python.bivittatus_genome.url,Vipera.berus_genome.url,Thermophis.baileyi_genome.url,Crotalus.viridis_genome.url,Pantherophis.obsoletus_genome.url,Hydrophis.melanocephalus_genome.url,Emydocephalus.ijimae_genome.url,Hydrophis.hardwickii_genome.url,Hydrophis.cyanocinctus_genome.url,Pseudonaja.textilis_genome.url,Laticauda.laticaudata_genome.url,Ptyas.mucosa_genome.url,Notechis.scutatus_genome.url,Thamnophis.elegans_genome.url,Naja.naja_genome.url,Podarcis_muralis_genome.url,Aspidoscelis_marmoratus_genome.url,Varanus_komodoensis_genome.url,Salvator_merianae_genome.url,Lacerta_viridis_genome.url,Lacerta_bilineata_genome.url,Paroedura_picta_genome.url,Zootoca_vivipara_genome.url,Lacerta_agilis_genome.url,Pantherophis.guttatus_genome.url_JTLQ00000000.1)
	species         <- c("Thamnophis sirtalis","Anolis carolinensis","Gekko japonicus","Pogona vitticeps","Crotalus horridus","Crotalus mitchellii","Ophiophagus hannah","Pantherophis guttatus","Protobothrops mucrosquamatus","Python bivittatus","Vipera berus","Thermophis baileyi","Crotalus viridis","Pantherophis obsoletus","Hydrophis melanocephalus","Emydocephalus ijimae","Hydrophis hardwickii","Hydrophis cyanocinctus","Pseudonaja textilis","Laticauda laticaudata","Ptyas mucosa","Notechis scutatus","Thamnophis elegans","Naja naja","Podarcis muralis","Aspidoscelis marmoratus","Varanus komodoensis","Salvator merianae","Lacerta viridis","Lacerta bilineata","Paroedura picta","Zootoca vivipara","Lacerta agilis","Pantherophis guttatus")
	notes           <- c("","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Genome version JTLQ00000000.1, which is a different individial than JTLQ00000000.2. The URL provided here contains the sequences in the three files: https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/JT/LQ/JTLQ01/JTLQ01.1.fsa_nt.gz, https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/JT/LQ/JTLQ01/JTLQ01.2.fsa_nt.gz, https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/JT/LQ/JTLQ01/JTLQ01.3.fsa_nt.gz")
	GFF.url         <- c(Thamnophis.sirtalis_genome.GFF.url,Anolis.carolinensis_genome.GFF.url,Gekko.japonicus_genome.GFF.url,Pogona.vitticeps_genome.GFF.url,Crotalus.horridus_genome.GFF.url,Crotalus.mitchellii_genome.GFF.url,Ophiophagus.hannah_genome.GFF.url,Pantherophis.guttatus_genome.GFF.url_JTLQ00000000.2,Protobothrops.mucrosquamatus_genome.GFF.url,Python.bivittatus_genome.GFF.url,Vipera.berus_genome.GFF.url,Thermophis.baileyi_genome.GFF.url,Crotalus.viridis_genome.GFF.url,Pantherophis.obsoletus_genome.GFF.url,Hydrophis.melanocephalus_genome.GFF.url,Emydocephalus.ijimae_genome.GFF.url,Hydrophis.hardwickii_genome.GFF.url,Hydrophis.cyanocinctus_genome.GFF.url,Pseudonaja.textilis_genome.GFF.url,Laticauda.laticaudata_genome.GFF.url,Ptyas.mucosa_genome.GFF.url,Notechis.scutatus_genome.GFF.url,Thamnophis.elegans_genome.GFF.url,Naja.naja_genome.GFF.url,Podarcis_muralis_genome.GFF.url,Aspidoscelis_marmoratus_genome.GFF.url,Varanus_komodoensis_genome.GFF.url,Salvator_merianae_genome.GFF.url,Lacerta_viridis_genome.GFF.url,Lacerta_bilineata_genome.GFF.url,Paroedura_picta_genome.GFF.url,Zootoca_vivipara_genome.GFF.url,Lacerta_agilis_genome.GFF.url,Pantherophis.guttatus_genome.GFF.url_JTLQ00000000.1)
	genome.urls.mat <- cbind(species,genome.url,GFF.url,notes)

	#### Details about blast output format 6
	column.names  <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
	column.descriptions <- c("query (e.g., unknown gene) sequence id","subject (e.g., reference genome) sequence id","percentage of identical matches","alignment length (sequence overlap)","number of mismatches","number of gap openings","start of alignment in query","end of alignment in query","start of alignment in subject","end of alignment in subject","expect value","bit score")
	blast.format6 <- cbind(column.names,column.descriptions)

	RE.name               <- c("PstI","HpaII","SbfI","EcoRI","MluCI","NlaIII")
	recognition.sequence  <- c("CTGCAG","CCGG","CCTGCAGG","GAATTC","AATT","CATG")
	rec.seq.length        <- c(6,4,8,6,4,4)
	cut.after.length      <- c(5,1,6,1,0,0)
	restriction.enzymes   <- cbind(RE.name,recognition.sequence,rec.seq.length,cut.after.length)

	# list of datasets
	result        <- list("genome.urls.mat","blast.format6","alignment.DNAStringSet1","alignment.DNAStringSet2","alignment.DNAStringSet3","fastestExonPerGene","fastestExonPerGene.best","Thamnophis.sirtalis_GFF","restriction.enzymes")
	
	# names of datasets
	names(result) <- c("SnakeCap.Genome.URLs","BLAST.Output.Format.6","Sample.DNA.Alignment1","Sample.DNA.Alignment2","Sample.DNA.Alignment3","SnakeCap.Fastest.Exon.Per.Gene","SnakeCap.Fastest.Exon.Per.Gene.Best","Thamnophis.sirtalis.genome.GFF","restriction.enzymes")

	if(is.null(x)){
		return(names(result))
	}
	if(is.character(x)){
		return(get(result[[which(names(result) %in% x)]]))
	}
	if(is.numeric(x)){
		return(get(result[[x]]))
	}
}


