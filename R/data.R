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

#' Proteome of antigen-experienced conventional CD4+ T cells.
#'
#' The cellular abundance of proteins extracted from lysates 
#' of CD4+ T cells lysates are shown as numbers of copies per cell
#' 
#' @format A data frame with 6355 rows (one per protein ID) and 23 variables.
#' 
"proteome_CD4"

#' Proteome of CD4+ T cells from Cas9-GFP mice after expansion and transfection with a CRISPR gRNA.
#'
#' The cellular abundance of proteins extracted from lysates 
#' of CD4+ T cells lysates are shown as numbers of copies per cell
#' 
#' @format A data frame with 19690 rows (one per protein ID) and 20 variables.
#' 
"proteome_CD4_expanded"

#' Proteome of Jurkat T cells
#'
#' The cellular abundance of proteins extracted from lysates 
#' of Jurkat T cells lysates are shown as numbers of copies per cell
#' 
#' @format A data frame with 4399 rows (one per protein ID) and 20 variables.
#' 
"proteome_Jurkat"