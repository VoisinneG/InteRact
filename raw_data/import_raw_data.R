
# import uniprot data

uniprot_data<-read.table( paste("~/Google_Drive/++Work/++Research/Resources-Databases/",
                        "Uniprot/uniprot-mus_musculus-SwissProt+TrEMBL.txt", sep=""),
                  sep="\t", header=TRUE, fill=TRUE, quote=c("\""),comment.char="")

# import proteome data

proteome_data <- read.table(paste("~/Google_Drive/++Work/++Research/++Projects/",
                          "Proteomes/Proteome_Comparison_OST_vs_kinetics/",
                          "Copy_number_hist_all_aligned_short.txt",
                          sep="" ),
                    sep="\t",header=TRUE)

devtools::use_data( uniprot_data, proteome_data, pkg=".", internal = TRUE, overwrite = TRUE)
