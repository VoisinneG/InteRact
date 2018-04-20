import_KEGG_pathway_info <- function( dest_dir="~/Google_Drive/++Work/++Research/Resources-Databases/KEGG/", organism="mouse"){
  
  url_adres <- switch(organism,
                      "mouse" = "http://rest.kegg.jp/link/mmu/pathway",
                      "human" = "http://rest.kegg.jp/link/hsa/pathway")
  
  file_name <- switch(organism,
                      "mouse" = "KEGG_mmu.txt",
                      "human" = "KEGG_hsa.txt")
  
  dest_file <- paste(dest_dir, file_name, sep="")
  
  if(!file.exists(dest_file)){
    download.file(url_adres, destfile = dest_file)
  }
  
  KEGG <- read.table(dest_file, header=FALSE)
  names(KEGG) <- c("pathway", "id")
  
  u_pathway <- unique(KEGG$pathway)
  name_pathway <- rep("", length(u_pathway))
  gene_pathway <- rep("", length(u_pathway))
  
  for ( i in 1:length(u_pathway) ){
    
    gene_pathway[i] <- paste(as.character(KEGG$id[which(KEGG$pathway == u_pathway[i])]), collapse=";")
    
    #for ( i in 1:100 ){
    url_adres <- paste("http://rest.kegg.jp/get/", u_pathway[i],sep="")
    dest_file <- paste("~/Google_Drive/++Work/++Research/Resources-Databases/KEGG/",u_pathway[i],".txt")
    if(!file.exists(dest_file)){
      download.file(url_adres, destfile = dest_file)
    }
    KEGG_pthw <- readLines(dest_file)
    s <- strsplit(KEGG_pthw[2], split = " ")[[1]]
    s <- s[s!=""]
    name_pathway[i] <- paste(s[2:(which(s=="-")-1)], collapse=" ")
  }
  
  df <- data.frame(pathway = u_pathway, name = name_pathway, IDs = gene_pathway)
  
  return(df)
}

KEGG_mouse <- import_KEGG_pathway_info(organism="mouse")

KEGG_human <- import_KEGG_pathway_info(organism="human")

# import uniprot data ------------------------------------------------------------------------------------------

uniprot_data_mouse <- read.table( paste("~/Google_Drive/++Work/++Research/Resources-Databases/",
                        "Uniprot/uniprot-mus+musculus.txt", sep=""),
                  sep="\t", header=TRUE, fill=TRUE, quote=c("\""),comment.char="")

# Note that only reviewd human protein are imported
uniprot_data_human <- read.table( paste("~/Google_Drive/++Work/++Research/Resources-Databases/",
                                      "Uniprot/uniprot-homo+sapiens+AND+reviewed.txt", sep=""),
                                sep="\t", header=TRUE, fill=TRUE, quote=c("\""),comment.char="")

# import Reactome data ------------------------------------------------------------------------------------------

reactome <- read.table( paste("~/Google_Drive/++Work/++Research/Resources-Databases/Reactome/",
                              "UniProt2Reactome_All_Levels.txt", sep=""),
                        sep="\t", header=FALSE, fill=TRUE, quote=c("\""),comment.char="", as.is=TRUE)

names(reactome) <- c("Protein.ID", "ID", "link", "Name", "Type", "Organism")

reactome_mouse <- reactome[ which(reactome$Organism == "Mus musculus"), ]

reactome_human <- reactome[ which(reactome$Organism == "Homo sapiens"), ]

# import PFAM data ----------------------------------------------------------------------------------------------

pfam_mouse <- read.table( paste("~/Google_Drive/++Work/++Research/Resources-Databases/PFAM/",
                                "10090.tsv", sep=""),
                          sep="\t", header=FALSE, skip=3, fill=TRUE, quote=c("\""),comment.char="", as.is=TRUE)
names(pfam_mouse) <- c("seq id", 
                       "alignment.start",
                       "alignment.end",
                       "envelope.start",
                       "envelope.end",
                       "hmm.acc",
                       "hmm.name",
                       "type",
                       "hmm.start",
                       "hmm.end",
                       "hmm.length",
                       "bit.score",
                       "E.value",
                       "clan")

pfam_human <- read.table( paste("~/Google_Drive/++Work/++Research/Resources-Databases/PFAM/",
                                "9606.tsv", sep=""),
                          sep="\t", header=FALSE, skip=3, fill=TRUE, quote=c("\""),comment.char="", as.is=TRUE)

names(pfam_human) <- names(pfam_mouse)

# import proteome data ------------------------------------------------------------------------------------------

proteome_data <- read.table(paste("~/Google_Drive/++Work/++Research/++Projects/",
                          "Proteomes/Proteome_Comparison_OST_vs_kinetics/",
                          "Copy_number_hist_all_aligned_short.txt",
                          sep="" ),
                    sep="\t",header=TRUE)

# save in ./R/sysdata.rda  ---------------------------------------------------------------------------------------

devtools::use_data( uniprot_data_mouse, 
                    uniprot_data_human, 
                    reactome_mouse, 
                    reactome_human, 
                    pfam_mouse,
                    pfam_human,
                    KEGG_mouse,
                    KEGG_human,
                    proteome_data, 
                    pkg=".", 
                    internal = TRUE, 
                    overwrite = TRUE)

