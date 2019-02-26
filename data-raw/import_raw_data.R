import_KEGG_pathway_info <- function( dest_dir="~/ownCloud//++Work/++Research/Resources-Databases/KEGG/", organism="mouse"){
  
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
    dest_file <- paste("~/ownCloud/++Work/++Research/Resources-Databases/KEGG/",u_pathway[i],".txt")
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

# # import uniprot data ------------------------------------------------------------------------------------------

uniprot_data_mouse <- read.table( paste("~/ownCloud/++Work/++Research/Resources-Databases/",
                                        "Uniprot/uniprot-mus+musculus.txt", sep=""),
                                  sep="\t", header=TRUE, fill=TRUE, quote=c("\""),comment.char="")

# Note that only reviewd human protein are imported

uniprot_data_human <- read.table( paste("~/ownCloud/++Work/++Research/Resources-Databases/",
                                        "Uniprot/uniprot-homo+sapiens+AND+reviewed.txt", sep=""),
                                  sep="\t", header=TRUE, fill=TRUE, quote=c("\""),comment.char="")

# import Reactome data ------------------------------------------------------------------------------------------

reactome <- read.table( paste("~/ownCloud/++Work/++Research/Resources-Databases/Reactome/",
                              "UniProt2Reactome_All_Levels.txt", sep=""),
                        sep="\t", header=FALSE, fill=TRUE, quote=c("\""),comment.char="", as.is=TRUE)

names(reactome) <- c("Protein.ID", "ID", "link", "Name", "Type", "Organism")

reactome_mouse <- reactome[ which(reactome$Organism == "Mus musculus"), ]

reactome_human <- reactome[ which(reactome$Organism == "Homo sapiens"), ]

# import PFAM data ----------------------------------------------------------------------------------------------

pfam_mouse <- read.table( paste("~/ownCloud/++Work/++Research/Resources-Databases/PFAM/",
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

pfam_human <- read.table( paste("~/ownCloud/++Work/++Research/Resources-Databases/PFAM/",
                                "9606.tsv", sep=""),
                          sep="\t", header=FALSE, skip=3, fill=TRUE, quote=c("\""),comment.char="", as.is=TRUE)

names(pfam_human) <- names(pfam_mouse)

# import Hallmark data ------------------------------------------------------------------------------------------

Hallmark_input <- readLines("~/ownCloud/++Work/++Research/Resources-Databases/Hallmark_database/h.all.v6.1.symbols.gmt.txt")
Hallmark_split <- strsplit(Hallmark_input, split="\t")

Hallmark_name <- rep("", length(Hallmark_split))
Hallmark_genes <- rep("", length(Hallmark_split))

for ( i in 1:length(Hallmark_split) ){
  Hallmark_name[i] <- Hallmark_split[[i]][1]
  Hallmark_genes[i] <- paste(Hallmark_split[[i]][3:length( Hallmark_split[[i]] )], collapse=";")
}

Hallmark <- data.frame(name = Hallmark_name,  gene = Hallmark_genes)


# import GO data -----------------------------------------------------------------------------------------------

library("ontologyIndex")

onto <- get_ontology("~/ownCloud/++Work/++Research/Resources-Databases/GO/go.obo", propagate_relationships = "is_a",
                     extract_tags = "everything")

names_gaf <- c(
  "DB",
  "DB_Object_ID",
  "DB_Object_Symbol",
  "Qualifier",
  "GO_ID",
  "DB:Reference",
  "Evidence Code",
  "With (or) From",
  "Aspect",
  "DB_Object_Name",
  "DB_Object_Synonym",
  "DB_Object_Type",
  "Taxon and Interacting taxon",
  "Date",
  "Assigned_By",
  "Annotation_Extension",
  "Gene_Product_Form_ID"
)

#Import GOA annotations for mouse uniprot proteome
GOA_mouse <- read.table("~/ownCloud/++Work/++Research/Resources-Databases/GO/goa_mouse.gaf", sep="\t", skip=12, quote="\"")

names(GOA_mouse) <- names_gaf
idx_match <- match(GOA_mouse$GO_ID, onto$id)
GOA_mouse$GO_type <- onto$namespace[idx_match]
GOA_mouse$GO_name <- onto$name[idx_match]

#Import GOA_slim annotations for mouse uniprot proteome
GOA_mouse_slim <- read.table("~/ownCloud/++Work/++Research/Resources-Databases/GO/goa_mouse_mapped_to_goslim_generic.gaf", sep="\t", skip=12, quote="\"")

names(GOA_mouse_slim) <- names_gaf
idx_match <- match(GOA_mouse_slim$GO_ID, onto$id)
GOA_mouse_slim$GO_type <- onto$namespace[idx_match]
GOA_mouse_slim$GO_name <- onto$name[idx_match]

#Import GOA annotations for mouse uniprot proteome
GOA_human <- read.table("~/ownCloud/++Work/++Research/Resources-Databases/GO/goa_human.gaf", sep="\t", skip=12, quote="\"")

names(GOA_human) <- names_gaf
idx_match <- match(GOA_human$GO_ID, onto$id)
GOA_human$GO_type <- onto$namespace[idx_match]
GOA_human$GO_name <- onto$name[idx_match]

#Import GOA_slim annotations for mouse uniprot proteome
GOA_human_slim <- read.table("~/ownCloud/++Work/++Research/Resources-Databases/GO/goa_human_mapped_to_goslim_generic.gaf", sep="\t", skip=12, quote="\"")

names(GOA_human_slim) <- names_gaf
idx_match <- match(GOA_human_slim$GO_ID, onto$id)
GOA_human_slim$GO_type <- onto$namespace[idx_match]
GOA_human_slim$GO_name <- onto$name[idx_match]



# import HPRD database

THPRD <- read.table(paste("~/ownCloud/++Work/++Research/Resources-Databases/HPRD/",
                          "HPRD_Release9_062910/BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt",
                          sep=""),
                    sep="\t",header = TRUE,quote="\"")

# import proteome data ------------------------------------------------------------------------------------------

load("~/ownCloud/++Work/++Research/++Projects/Proteomes/Merge_Proteome_OST_and_CD4_kinetics/proteome_CD4.rda")
load("~/ownCloud/++Work/++Research/++Projects/Proteomes/Proteome_CD4+_transfected_expanded/proteome_CD4_expanded.rda")
load("~/ownCloud/++Work/++Research/++Projects/Proteomes/Proteome_Jurkat_Itsn2_KO/proteome_Jurkat.rda")


# proteome_data <- read.table(paste("~/ownCloud/++Work/++Research/++Projects/",
#                                   "Proteomes/Proteome_Comparison_OST_vs_kinetics/",
#                                   "Copy_number_hist_all_aligned_short.txt",
#                                   sep="" ),
#                             sep="\t",header=TRUE)




# save in ./R/sysdata.rda  ---------------------------------------------------------------------------------------

devtools::use_data(
  #uniprot_data_mouse,
  #uniprot_data_human,
  #reactome_mouse,
  #reactome_human,
  #pfam_mouse,
  #pfam_human,
  #KEGG_mouse,
  #KEGG_human,
  #Hallmark,
  #GOA_mouse,
  #GOA_mouse_slim,
  #GOA_human,
  #GOA_human_slim,
  proteome_CD4,
  proteome_CD4_expanded,
  proteome_Jurkat,
  THPRD,
  pkg=".",
  internal = TRUE,
  overwrite = TRUE)

# devtools::use_data( uniprot_data_mouse,
#                     #uniprot_data_human,
#                     proteome_data,
#                     pkg=".",
#                     internal = TRUE,
#                     overwrite = TRUE)


proteinGroups_Cbl <- read.csv("~/Documents/Data/BMMlab/CBL-CBLB_paper/C-CBL/data/ProteinGroups_Cbl.txt", sep="\t", nrows=-1, fill=TRUE, na.strings="", dec=",")

devtools::use_data(proteinGroups_Cbl,
                   proteome_CD4,
                   proteome_CD4_expanded,
                   proteome_Jurkat,
                   pkg=".",
                   internal = FALSE,
                   overwrite = TRUE)
