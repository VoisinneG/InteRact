
# import uniprot data

uniprot_data_mouse <- read.table( paste("~/Google_Drive/++Work/++Research/Resources-Databases/",
                        "Uniprot/uniprot-mus+musculus.txt", sep=""),
                  sep="\t", header=TRUE, fill=TRUE, quote=c("\""),comment.char="")

# Note that only reviewd human protein are imported
uniprot_data_human <- read.table( paste("~/Google_Drive/++Work/++Research/Resources-Databases/",
                                      "Uniprot/uniprot-homo+sapiens+AND+reviewed.txt", sep=""),
                                sep="\t", header=TRUE, fill=TRUE, quote=c("\""),comment.char="")

# import proteome data

proteome_data <- read.table(paste("~/Google_Drive/++Work/++Research/++Projects/",
                          "Proteomes/Proteome_Comparison_OST_vs_kinetics/",
                          "Copy_number_hist_all_aligned_short.txt",
                          sep="" ),
                    sep="\t",header=TRUE)


devtools::use_data( uniprot_data_mouse, uniprot_data_human, proteome_data, pkg=".", internal = TRUE, overwrite = TRUE)
