#' Characterization of protein groups identified from AP-MS data using MaxQuant 
#'
#' A dataset containing the characterization of protein groups identified using MaxQuant from AP-MS samples.
#' Samples correspond to pull-down experiments performed on Cd4+ T cells from
#' CBL-OST mice (expressing a tagged version of the bait protein Cbl) and WT mice (expressing the endogeneous Cbl protein).
#' Cd4+ T cells from CBL-OST and WT backgrounds were stimulated by cross-linking of Cd3e and Cd4 for different times (0, 30, 120, 300 and 600 seconds).
#' Samples were then lysed and subjected to MS-MS analysis. 
#' Experiments were conducted in three biological replicates (S1, S2 and S3).
#' Sample were analysed by MS-MS in three technical replicates (R1, R2 and R3).
#' 
#' @format A data frame with 1540 rows (protein groups) and 358 variables:
#' 
#' See the \href{ftp://ftp.lrz.de/transfer/proteomics/Phospho_Library/archive/MaxQuant/LibraryResults/EtdNoFilter/combined/txt/tables.pdf}{Maxquant documentation}
#' for the description of the different variables
#' 
"proteinGroups_Cbl"


#' An example Interactome computed from the 'proteinGroups_Cbl' file. 
#' This interactome was created using the following code :
#' \code{InteRact(proteinGroups_Cbl, bait_gene_name = "Cbl")}
#' @format An InteRactome
"Interactome_Cbl"