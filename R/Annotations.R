#' Perform enrichment analysis
#' @description Perform enrichment analysis for protein annotations stored in a formatted data.frame
#' @param df a data.frame with annotations corresponding to each row. Types of annotations are organized by columns. 
#' For a given type of annotations, annotations are separated by \code{sep}.
#' @param sep Character string separating different annotations of a given type
#' @param idx_detect indexes of the foreground set.
#' @param annotation_selected set of annotations on which to perform the analysis
#' @param names row names. Used in the output data.frame
#' @param organism organism for which the analysis is to be performed ("mouse" or "human")
#' @param updateProgress logical, function to show progress in shiny app
#' @param showProgress logical, show progress in console
#' @param orderOutput logical, order annotations by enrichment p-values in the output data.frame
#' @return a data.frame
#' @import utils
#' @importFrom stats phyper p.adjust
#' @export
annotation_enrichment_analysis <- function( df, 
                                            sep = ";",
                                            idx_detect, 
                                            annotation_selected = c("Keywords", "Protein.families") , 
                                            names = df$Gene.names...primary.., 
                                            organism = "mouse",
                                            updateProgress = NULL, 
                                            showProgress = TRUE,
                                            orderOutput = TRUE){
  # df : data frame with annotation data
  # idx_detect : indices of the subset of proteins in df for which the enrichment analysis is performed
  # (list of indices are supported. The output will then be a list of data-frames)
  # against the background formed by all proteins in df
  # annotation_selected : set of annotation terms to consider. 
  # Annotations supported are stored in varaiable "supported_annotations:
  
  #supported_annotations <- c( 
  # "Protein.families",  
  # "Pfam", 
  # "Keywords", 
  # "Reactome", 
  # "GO",
  # "Hallmark",
  # "KEGG",
  # "GO_molecular_function",
  # "GO_biological_process",
  # "GO_cellular_component",
  # "GOslim_molecular_function",
  # "GOslim_biological_process",
  # "GOslim_cellular_component",
  # "Motif",
  # "Kinase")
  # 
  
  if(  is.null(df) |  (sum(annotation_selected %in% names(df)) != length(annotation_selected)) ){
    stop("Annotations not available. Import annotations first.")
  }else if ( length(annotation_selected) == 0) {
    stop("No annotations selected. Change selected annotations")
  }
  
  if (showProgress) cat("Perform annotation enrichment analysis...\n")
  
  #list annotation terms found in the dataset ------------------------------------------------
  
  df_int <- df
  
  # if("Protein.families" %in% annotation_selected){
  #   df_int[["Protein.families"]] <- as.character(df_int[["Protein.families"]])
  #   for (i in 1:length(df_int[["Protein.families"]])){
  #     s <- strsplit(df_int[["Protein.families"]][i], split=", ")[[1]]
  #     s_develop <- s[1]
  #     if (length(s) > 1) {
  #       for (j in 2:length(s)){
  #         s_develop <- c(s_develop, paste(s_develop[length(s_develop)], s[j], sep=", "))
  #       }
  #     }
  #     df_int[["Protein.families"]][i] <- paste(s_develop, collapse="; ")
  #   }
  # }
  
  annot_terms <- NULL
  annot_type <- NULL
  annot_names <- NULL
  
  for (annot_type_sel in annotation_selected){
    
    collapse_sep <- sep
    
    # if( annot_type_sel %in% c("Protein.families", "Keywords") ) collapse_sep <- "; "
    
    
    u_annot<-paste(unique(df_int[[annot_type_sel]]), collapse = collapse_sep)
    terms <- unique(strsplit(u_annot, split = collapse_sep)[[1]])
    
    annot_names_int <- terms
    # if(annot_type_sel == "Reactome"){
    #   reactome <- switch(organism, "mouse" = reactome_mouse, "human" = reactome_human)
    #   annot_names_int <- paste(
    #                         reactome$Name[match(terms, reactome$ID)],
    #                         ", [",
    #                         reactome$ID[match(terms, reactome$ID)],
    #                         "]",
    #                         sep = "")
    # }
    # if(annot_type_sel == "Pfam"){
    #   pfam <- switch(organism, "mouse" = pfam_mouse, "human" = pfam_human)
    #   annot_names_int <- paste(
    #                           pfam$hmm.name[match(terms, pfam$hmm.acc)],
    #                           ", ",
    #                           pfam$type[match(terms, pfam$hmm.acc)],
    #                           ", [",
    #                           pfam$hmm.acc[match(terms, pfam$hmm.acc)],
    #                           "]",
    #                           sep = "")
    # }
    
    annot_terms <- c(annot_terms,  terms)
    annot_type <- c(annot_type, rep(annot_type_sel, length(terms)))
    annot_names <- c(annot_names, annot_names_int)
    
  }
  
  df.annot <- data.frame(annot_terms = annot_terms, annot_type = annot_type, annot_names = annot_names)
  df.annot <- df.annot[which(df.annot$annot_terms != ""), ]
  
  
  
  # Compute Background -------------------------------------------------------------------------------------
  
  nodes_tot <- as.character(names);
  
  u_annot_nodes_collapse <- rep("", length(nodes_tot));
  idx_tot <- rep(0, length(nodes_tot));
  
  for ( i in 1:length(nodes_tot) ){
    s <- NULL
    for (annot_type in annotation_selected) {
      s <- c(s, as.character(df[[annot_type]][ i ]))
    }
    u_annot_nodes_collapse[i] <- paste(s, collapse = ";")
  }
  
  N_background = length(nodes_tot);
  
  n_annot <- dim(df.annot)[1]
  
  N_annot_background <- rep(0, n_annot);
  freq_annot_background <- rep(0, n_annot);
  nodes_annot_background <- rep("", n_annot);
  
  if (showProgress & typeof(idx_detect)!="list") pb <- txtProgressBar(min = 0, max = 2*n_annot, style = 3)
  count<-0
  
  for ( k in 1:dim(df.annot)[1] ){
    
    idx_annot <- grep(df.annot$annot_terms[k], u_annot_nodes_collapse, fixed=TRUE)
    N_annot_background[k] = length(idx_annot);
    nodes_annot_background[k] = paste(nodes_tot[idx_annot], collapse=";")
    freq_annot_background[k] = N_annot_background[k]/N_background;
    
    count <- count +1
    if (showProgress & typeof(idx_detect)!="list") setTxtProgressBar(pb, count)
    # progress bar
    if (is.function(updateProgress)) {
      text <- paste0( round(count/(2*n_annot)*100, 0), " %")
      updateProgress(value = count/(2*n_annot)*100, detail = text)
    }
  }
  
  N_annotation_test <- length(which(N_annot_background>0))
  
  # Perform enrichment test for each annotation and each subset of indices ---------------------------------------
  
  if (typeof(idx_detect)=="list"){
    n_sets <- length(idx_detect)
  } else {
    n_sets = 1
  }
  
  df.annot.tot <- list()
  
  if (showProgress  & typeof(idx_detect)=="list") pb <- txtProgressBar(min = 0, max = n_sets, style = 3)
  
  for (i in 1:n_sets){
    
    if (showProgress & typeof(idx_detect)=="list") setTxtProgressBar(pb, i)
    
    if (typeof(idx_detect)=="list"){
      idx_d <- idx_detect[[i]]
    } else {
      idx_d <- idx_detect
    }
    
    N_annot <- rep(0, n_annot);
    freq_annot <- rep(0, n_annot);
    nodes_annot <- rep("", n_annot);
    p_value <- rep(0, n_annot);
    fold_change <- rep(0, n_annot);
    p_value_adjust <- rep(0, n_annot);
    
    for(k in 1:n_annot ){
      
      
      idx_annot <- idx_d[ grep(df.annot$annot_terms[k], u_annot_nodes_collapse[idx_d],fixed=TRUE) ]
      N_annot[k]=length(idx_annot);
      N_sample = length(idx_d);
      
      freq_annot[k] = N_annot[k]/N_sample;
      nodes_annot[k]=paste(nodes_tot[idx_annot], collapse=";")
      
      p_value[k] = 1-phyper(N_annot[k]-1, 
                            N_annot_background[k],  
                            N_background-N_annot_background[k],  
                            N_sample);
      
      fold_change[k] = freq_annot[k]/freq_annot_background[k];
      
      count <- count +1
      if (showProgress & typeof(idx_detect)!="list") setTxtProgressBar(pb, count)
      # progress bar
      if (is.function(updateProgress)) {
        text <- paste0( round(count/(2*n_annot)*100, 0), " %")
        updateProgress(value = count/(2*n_annot)*100, detail = text)
      }
    }
    
    if (showProgress & typeof(idx_detect)!="list"){
      close(pb)
      cat("Done.\n")
    }
    
    
    idx_annot_exist <-  which(N_annot>0);
    p_value_adjust_fdr <- rep( 1,length(p_value) );
    p_value_adjust_bonferroni <- rep( 1,length(p_value) );
    p_value_adjust_bonferroni[idx_annot_exist] <- p.adjust(p_value[idx_annot_exist], method = "bonferroni");
    p_value_adjust_fdr[idx_annot_exist] <- p.adjust(p_value[idx_annot_exist], method = "fdr");
    
    df.annot.set <- data.frame(
      N_annot,
      freq_annot,
      fold_change, 
      p_value, 
      p_value_adjust_fdr,
      nodes_annot,
      p_value_adjust_bonferroni,
      N_annot_background, 
      freq_annot_background,
      nodes_annot_background)
    
    df.annot.set <- cbind(df.annot, df.annot.set)
    
    if (orderOutput) df.annot.set <- df.annot.set[ order(df.annot.set$p_value, decreasing = FALSE), ]
    
    df.annot.tot[[i]] <- df.annot.set
    
  }
  
  if (showProgress & typeof(idx_detect)=="list") close(pb)
  
  if (typeof(idx_detect)=="list"){
    return(df.annot.tot)
  } else {
    return(df.annot.tot[[1]])
  }
  
}


#' Append annotations to an \code{InteRactome}
#' @param res an \code{InteRactome}
#' @param annotations type of annotations to append
#' @param name_id column name used to map protein identifiers
#' @param organism organism for which the annotations have to be appended
#' @return an \code{InteRactome}
#' @export
append_annotations <- function( res, annotations=NULL, name_id = "Protein.IDs", organism = "mouse"){
  
  res_int<-res
  if(is.null(annotations)){
    warning("No annotations to append")
    return(res_int)
  }else if(is.data.frame(annotations)){
    df_annot <- annotations
  } else{
    
    df_annot <- get_annotations(res, name_id = name_id, organism = organism)
    
    for (annot in annotations){
      df_annot <- switch(annot,
                         "GO" =  add_GO_data(df_annot, GO_type = "molecular_function")
      )
    }
    
  }
  
  n_annot <- 0
  for (annot_var in names(df_annot)){
    n_annot <- n_annot + length(df_annot[[annot_var]])
  }
  
  if( is.null(df_annot) | n_annot == 0){
    warning("No annotations to append")
  }
  else{
    cat("Append annotation to interactome...\n")
    idx_match<-rep(NA,length(res$names))
    for(i in 1:length(res$names) ){
      idx_match[i] <- which(as.character(df_annot[[name_id]]) == as.character(res[[name_id]][i]) )
    }
    
    for( var_names in setdiff(names(df_annot), name_id) ){
      res_int[[var_names]] <- as.character(df_annot[[var_names]][idx_match])
    }
    cat("Done.\n")
  }
  
  return(res_int)
  
}

#' @export
append_annotations_enrichr <- function(res, annotations, name_id = "names"){
  
  res_int<-res
  if(is.null(annotations)){
    warning("No annotations to append")
    return(res_int)
  }else if(is.data.frame(annotations)){
    df_annot <- annotations
  } else{
    df_annot <- get_annotations_enrichr(res["names"], name_id = name_id, dbs = annotations)
    
  }
  
  n_annot <- 0
  for (annot_var in names(df_annot)){
    n_annot <- n_annot + length(df_annot[[annot_var]])
  }
  
  if( is.null(df_annot) | n_annot == 0){
    warning("No annotations to append")
  }
  else{
    cat("Append annotation to interactome...\n")
    idx_match<-rep(NA,length(res$names))
    for(i in 1:length(res$names) ){
      idx_match[i] <- which(as.character(df_annot[[name_id]]) == as.character(res[[name_id]][i]) )
    }
    
    for( var_names in setdiff(names(df_annot), name_id) ){
      res_int[[var_names]] <- as.character(df_annot[[var_names]][idx_match])
    }
    cat("Done.\n")
  }
  
  return(res_int)
}

#' @import enrichR
#' @export
get_annotations_enrichr <- function(data, name_id = "names", dbs = "GO_Biological_Process_2018"){
  #library(enrichR)
  #enrichR::listEnrichrDbs()
  #dbs<-"GO_Biological_Process_2017"
  
  df <- data
  
  if( length(setdiff(dbs, names(data)))==0 ){
    warning("Annotations already loaded")
    return(df)
  }
  
  dbs_int <- setdiff(dbs, names(data))
  enriched <- enrichR::enrichr(as.character(data[[name_id]]), dbs_int)
  
  annot <- vector("list", length(dbs_int))
  names(annot) <- dbs_int
  
  
  
  for(i in 1:length(dbs_int)){
    
    annot[[i]] <- rep("", length(data[[name_id]]) )
    
    for(j in 1:length(enriched[[dbs_int[i]]]$Term)){
      genes <- strsplit(enriched[[dbs_int[i]]]$Genes[j], split = ";")[[1]]
      idx_match <- match(genes, toupper(data[[name_id]]))
      if(length(idx_match)>0){
        for(k in 1:length(idx_match)){
          if(nchar(annot[[i]][idx_match[k]])>0){
            annot[[i]][idx_match[k]] <- paste( c(annot[[i]][idx_match[k]], enriched[[dbs_int[i]]]$Term[j]), collapse = ";")
          }else{
            annot[[i]][idx_match[k]] <- enriched[[dbs_int[i]]]$Term[j]
          }
          
        }
      }
      
    }
    
    
  }
  
  annot <- as.data.frame(annot)
  annot[[name_id]] <- data[[name_id]]
  
  df <- merge(df, annot, by = name_id)
  
  return(df)
}

#' Get annotations from uniprot for a set of protein identifiers
#' @description Get annotations from uniprot for a set of protein identifiers.
#' From a set of IDs, keep the first that correspond to a "reviewed" protein, 
#' or by default the first ID of the set
#' name_id : column containing the set of protein identifiers separated by "split_param"
#' organism = c("mouse", "human")
#' @param data a data.frame with protein IDs in column \code{name_id}
#' @param name_id column name used to map protein identifiers
#' @param split_param split character used to separate different protein IDs
#' @param organism organism for which the annotations have to be appended
#' @param updateProgress logical, function to show progress in shiny app
#' @return an annotated data.frame
#' @import utils
#' @export
get_annotations <- function( data, name_id = "Protein.IDs", split_param = ";", organism = "mouse", updateProgress = NULL ){
  # Get annotations from uniprot for a set of protein identifiers
  # From a set of IDs, keep the first that correspond to a "reviewed" protein, 
  # or by default the first ID of the set
  # name_id : column containing the set of protein identifiers separated by "split_param"
  # organism = c("mouse", "human")
  
  df_annot <- switch(organism, "mouse" = uniprot_data_mouse, "human" = uniprot_data_human)
  
  df<-NULL
  
  nodes_IDs<- as.character(data[[name_id]]);
  Nnodes<-length(nodes_IDs)
  
  if( length( strsplit( paste(nodes_IDs, collapse=" "), split = split_param, fixed = TRUE)[[1]] ) == 1){
    warning(paste("The split parameter '",split_param, "' was not detected in '", name_id,"' : you might need to change 'split_param'.", sep = "") )
  }
  
  nodes_ID<-rep("", Nnodes);
  imatch<-rep(NA, Nnodes);
  
  
  # create progress bar
  cat("Creating annotation table...\n")
  pb <- txtProgressBar(min = 0, max = Nnodes, style = 3)
  
  
  for (i in 1:Nnodes ){
    
    list_nodes_ID<- unique(strsplit(nodes_IDs[i], split=split_param, fixed = TRUE)[[1]]);
    for(j in 1:length(list_nodes_ID) ){
      nID_clean <- strsplit(list_nodes_ID[j], split="-")[[1]][1];
      idx_entry<-which(df_annot$Entry == nID_clean)
      if(length(idx_entry>0)) {
        
        if(is.na(imatch[i])){
          imatch[i]<-idx_entry;
        }
        
        if(df_annot$Status[idx_entry]=="reviewed"){
          imatch[i]<-idx_entry;
          break;
        }
        
      }
    }
    setTxtProgressBar(pb, i)
    
    if (is.function(updateProgress)) {
      text <- paste0( round(i/Nnodes*100, 0), " %")
      updateProgress(value = i/Nnodes*100, detail = text)
    }
    
  }
  close(pb)
  
  df<-df_annot[imatch, ]
  df <- cbind(data[[name_id]], df)
  names(df)[1] <- name_id
  
  names(df)[which(names(df) == "Cross.reference..Pfam.")]<-"Pfam"
  names(df)[which(names(df) == "Cross.reference..Reactome.")]<-"Reactome"
  names(df)[which(names(df) == "Gene.ontology..GO.")]<-"GO"
  
  cat("Done.\n")
  
  return(df)
  
}

#' Add GO annotations
#' @description Add GO annotations corresponding to a set of protein identifiers
#' @param df a data.frame with protein IDs in column \code{map_id}
#' @param map_id column name used to map protein identifiers
#' @param GO_type type of GO term ("molecular_function", "biological_process" or "cellular_component")
#' @param organism organism for which the annotations have to be mapped
#' @param slim logical, use GO_slim annotations
#' @param updateProgress logical, function to show progress in shiny app
#' @return a data.frame with an additional GO annotation column
#' @import utils
#' @export
add_GO_data <- function(df, map_id = "Entry", GO_type="molecular_function", organism = "mouse", slim = FALSE, updateProgress = NULL){
  # GO_type = "molecular_function", biological_process" or "cellular_component"
  # organism = "mouse" or "human
  # slim = TRUE or FALSE (use GO slim annotations)
  
  df_int <- df
  n <- length(df_int[[map_id]])
  
  GO_terms <- rep("", n)
  
  name_GOA <- paste("GOA_", organism, sep = "");
  if (slim) {
    name_GOA <- paste(name_GOA, "_slim", sep = "")
  }
  GOA <- get(name_GOA)
  
  idx_type <- as.vector(which(GOA$GO_type == GO_type))
  
  # create progress bar
  cat("Add GO annotation data...\n")
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  
  for(i in 1:n){
    # progress bar
    if (is.function(updateProgress)) {
      text <- paste0( round(i/n*100, 0), " %")
      updateProgress(value = i/n*100, detail = text)
    }
    setTxtProgressBar(pb, i)
    
    idx_GO <- idx_type[ grep(df_int[[map_id]][i], GOA$DB_Object_ID[idx_type], fixed = TRUE) ]
    if(length(idx_GO) > 0){
      term <- paste(GOA$GO_name[idx_GO], " [", GOA$GO_ID[idx_GO], "]", sep = "")
      GO_terms[i] <- paste(as.character(unique(term)), collapse = ";")
    }
    
  }
  close(pb)
  
  name_annot <- "GO_"
  if( slim ) name_annot <- paste(name_annot, "slim_", sep ="")
  name_annot <- paste(name_annot, GO_type, sep= "")
  df_int[[name_annot]] <- GO_terms
  
  return(df_int)
}

#' Add KEGG pathway annotations
#' @description Add KEGG pathway annotations corresponding to a set of KEGG IDs
#' @param df a data.frame with KEGG IDs in column \code{map_id}
#' @param map_id column name used to map KEGG IDs
#' @param organism organism for which the annotations have to be mapped
#' @param updateProgress logical, function to show progress in shiny app
#' @return a data.frame with an additional KEGG pathway annotation column
#' @import utils
#' @export
add_KEGG_data <- function(df, map_id = "Cross.reference..KEGG.", organism="mouse", updateProgress = NULL){
  
  df_int <- df
  n <- length(df_int[[map_id]])
  
  # Add KEGG data
  KEGG_pathways <- rep("", n)
  KEGG <- switch(organism, "mouse" = KEGG_mouse, "human" = KEGG_human)
  for(i in 1:n){
    # progress bar
    if (is.function(updateProgress)) {
      text <- paste0( round(i/n*100, 0), " %")
      updateProgress(value = i/n*100, detail = text)
    }
    
    if(!is.na(df_int[[map_id]][i])){
      if( nchar(as.character(df_int[[map_id]][i])) >0 ){
        idx_KEGG <- grep(df_int[[map_id]][i], as.character(KEGG$IDs), fixed = TRUE)
        pathways <- paste(KEGG$name[idx_KEGG]," [",KEGG$pathway[idx_KEGG],"]",sep="")
        KEGG_pathways[i] <- paste(as.character(pathways), collapse = ";")
      }
    }
    
  }
  df_int$KEGG <- KEGG_pathways
  
  return(df_int)
}

#' Add Hallmark annotations
#' @description Add Hallmark annotations corresponding to a set of protein identifiers
#' @param df a data.frame with protein IDs in column \code{map_id}
#' @param map_id column name used to map protein identifiers
#' @param updateProgress logical, function to show progress in shiny app
#' @return a data.frame with an additional Hallmark annotation column
#' @import utils
#' @export
add_Hallmark_data <- function(df, map_id="Gene.names...primary..", updateProgress = NULL){
  
  df_int <- df
  
  n <- length(df_int[[map_id]])
  hallmark_set <- rep("", n )
  
  
  # create progress bar
  cat("Add Hallmark annotation data...\n")
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  
  for(i in 1:n){
    # progress bar
    if (is.function(updateProgress)) {
      text <- paste0( round(i/n*100, 0), " %")
      updateProgress(value = i/n*100, detail = text)
    }
    setTxtProgressBar(pb, i)
    idx_set <- NULL
    for(j in 1:dim(Hallmark)[1]){
      idx_in_geneset <- match(toupper(df_int[[map_id]][i]), strsplit(as.character(Hallmark$gene[j]), split=";")[[1]])
      if(!is.na(idx_in_geneset)){
        idx_set <- c(idx_set, j)
      }
    }
    hallmark_set[i] <- paste( as.character(Hallmark$name[idx_set]), collapse=";")
  }
  close(pb)
  
  df_int$Hallmark <- hallmark_set
  
  return(df_int)
}

#' Retrieve protein-protein interaction information using PSICQUIC
#' @param gene_name the gene name for which to retrieve PPI
#' @param tax_ID taxon ID for which to retrieve PPI
#' @param provider database from which to retrieve PPI
#' @return a data.frame PPI information
#' @import PSICQUIC
get_PPI_from_psicquic <- function( gene_name, tax_ID = c(9606,10090) , provider = c("IntAct","MINT") ){
  
  psicquic <- PSICQUIC::PSICQUIC()
  
  for (k in 1:length(tax_ID) ){
    
    tbl <- PSICQUIC::interactions(psicquic, 
                                  gene_name, 
                                  species = tax_ID[k] , 
                                  provider = provider )
    
    
    s<-strsplit(tbl$aliasA, split="|", fixed = TRUE);
    gene_name_A <- rep("",length(s))
    if(length(s)>0){
      for (i in 1:length(s) ){
        ign<-grep("(gene name)",s[[i]],fixed=TRUE)
        if(length(ign)>0){
          gene_name_A[i] <- strsplit( strsplit(s[[i]][ign],split=":")[[1]][2], 
                                      split="(" ,fixed=TRUE )[[1]][1]
          
        }
      }
    }
    
    s<-strsplit(tbl$aliasB,split="|",fixed = TRUE);
    gene_name_B <- rep("",length(s))
    if(length(s)>0){
      for (i in 1:length(s) ){
        ign<-grep("(gene name)",s[[i]],fixed=TRUE)
        if(length(ign)>0){
          gene_name_B[i] <- strsplit( strsplit(s[[i]][ign],split=":")[[1]][2], 
                                      split="(" ,fixed=TRUE )[[1]][1]
        }
      }
    }
    
    s<-strsplit(tbl$A,split = ":")
    uniprot_A <- rep("",length(s))
    if(length(s)>0){
      for (i in 1:length(s) ){
        uniprot_A[i] <- s[[i]][2]
      }
    }
    
    s<-strsplit(tbl$B,split = ":")
    uniprot_B <- rep("",length(s))
    if(length(s)>0){
      for (i in 1:length(s) ){
        uniprot_B[i] <- s[[i]][2]
      }
    }
    
    s<-strsplit(tbl$publicationID,split = "|",fixed=TRUE)
    Pubmed_ID <- rep("",length(s))
    if(length(s)>0){
      for (i in 1:length(s) ){
        ip<-grep("pubmed",s[[i]],fixed=TRUE);
        Pubmed_ID[i] <- strsplit(s[[i]][ip],split=":")[[1]][2];
      }
    }
    
    s<-strsplit(tbl$type,split = "(",fixed=TRUE)
    Int_type <- rep("",length(s))
    if(length(s)>0){
      for (i in 1:length(s) ){
        Int_type[i] <- strsplit(s[[i]][2],split=")",fixed=TRUE)[[1]][1];
      }
    }
    
    s<-strsplit(tbl$detectionMethod,split = "(",fixed=TRUE)
    Detection_method <- rep("",length(s))
    if(length(s)>0){
      for (i in 1:length(s) ){
        Detection_method[i] <- strsplit(s[[i]][2],split=")",fixed=TRUE)[[1]][1];
      }
    }
    
    s<-strsplit(tbl$firstAuthor,split = " ",fixed=TRUE)
    Author <- rep("",length(s))
    if(length(s)>0){
      for (i in 1:length(s) ){
        Author[i] <- paste(s[[i]][1], s[[i]][length(s[[i]])],sep=" ");
      }
    }
    
    taxon<-rep(tax_ID[k],length(s) );
    
    if(k>1){
      #df2<-data.frame(gene_name_A, gene_name_B, uniprot_A, uniprot_B, taxon, Int_type, Detection_method, Author=Author, pubmed_ID=pubmed_ID, Database=tbl$provider)
      df2<-data.frame(gene_name_A, gene_name_B, taxon, Int_type, 
                      Detection_method, Author=Author, Pubmed_ID=Pubmed_ID, Database=tbl$provider)
      df1<-rbind(df1,df2);
      
    }
    else{
      #df1<-data.frame(gene_name_A, gene_name_B, uniprot_A, uniprot_B, taxon, Int_type, Detection_method, Author=Author, pubmed_ID=pubmed_ID, Database=tbl$provider)
      df1<-data.frame(gene_name_A, gene_name_B, taxon, Int_type, 
                      Detection_method, Author=Author, Pubmed_ID=Pubmed_ID, Database=tbl$provider)
      
    }
    
    df1 <- df1[nchar(as.character(df1$gene_name_A))>0 & nchar(as.character(df1$gene_name_B))>0 , ]
    
  }
  
  return(df1)
}

#' Retrieve protein-protein interaction information from BioGRID
#' @param gene_name the gene name for which to retrieve PPI
#' @param tax_ID taxon ID for which to retrieve PPI
#' @return a data.frame PPI information
get_PPI_from_BioGRID <- function( gene_name, tax_ID = c(9606,10090) ){
  
  access_key <- "7ad36061b7644111aa9f5b3948429fb2"
  
  for (k in 1:length(tax_ID) ){
    
    url_adress <- paste("http://webservice.thebiogrid.org/interactions?searchNames=true&geneList=",
                        gene_name,"&includeInteractors=true&format=tab2&includeHeader=true&taxId=",
                        tax_ID[k],"&accesskey=",
                        access_key,sep="");
    
    Tbiogrid <- read.table(url_adress,header=TRUE,fill=TRUE,sep="\t",comment.char="", quote="\"")
    
    Tbiogrid <- Tbiogrid[Tbiogrid$Organism.Interactor.A == Tbiogrid$Organism.Interactor.B,]
    
    taxon_biogrid <- rep(tax_ID[k],dim(Tbiogrid)[1] );  
    
    s<-strsplit(as.character(Tbiogrid$Author),split = " ",fixed=TRUE)
    Author_Biogrid <- rep("",length(s))
    
    if(length(s)>0){
      for (i in 1:length(s) ){
        Author_Biogrid[i] <- paste(s[[i]][1], s[[i]][length(s[[i]])],sep=" ");
      }
    }
    
    
    
    if(k>1){
      df_biogrid_2 <- data.frame(gene_name_A=Tbiogrid$Official.Symbol.Interactor.A, 
                                 gene_name_B = Tbiogrid$Official.Symbol.Interactor.B, 
                                 taxon=taxon_biogrid, 
                                 Int_type=Tbiogrid$Experimental.System.Type, 
                                 Detection_method=Tbiogrid$Experimental.System, 
                                 Author=Author_Biogrid, 
                                 Pubmed_ID=Tbiogrid$Pubmed.ID,  
                                 Database=Tbiogrid$Source.Database )
      df_biogrid_1 <- rbind(df_biogrid_1,df_biogrid_2);
    }
    else{
      df_biogrid_1 <- data.frame(gene_name_A=Tbiogrid$Official.Symbol.Interactor.A, 
                                 gene_name_B = Tbiogrid$Official.Symbol.Interactor.B, 
                                 taxon=taxon_biogrid, 
                                 Int_type=Tbiogrid$Experimental.System.Type, 
                                 Detection_method=Tbiogrid$Experimental.System, 
                                 Author=Author_Biogrid, 
                                 Pubmed_ID=Tbiogrid$Pubmed.ID,  
                                 Database=Tbiogrid$Source.Database )
      
    }
    
  }
  
  return(df_biogrid_1)
  
}

#' Retrieve protein-protein interaction information from HPRD
#' @param gene_name the gene name for which to retrieve PPI
get_PPI_from_HPRD <- function( gene_name ){
  
  THPRD <- THPRD[which(THPRD$Gene_symbol_1 == toupper(gene_name) | THPRD$Gene_symbol_2 == toupper(gene_name)), ]
  
  df_HPRD <- data.frame(gene_name_A = THPRD$Gene_symbol_1, 
                        gene_name_B = THPRD$Gene_symbol_2, 
                        taxon = rep(9606, dim(THPRD)[1] ), 
                        Int_type = rep("NA", dim(THPRD)[1] ), 
                        Detection_method = THPRD$Experiment_type, 
                        Author = rep("NA", dim(THPRD)[1] ), 
                        Pubmed_ID = THPRD$Pubmed_id, 
                        Database = rep("HPRD", dim(THPRD)[1] ));
  return(df_HPRD)
  
}

#' Retrieve protein-protein interaction information from databses 
#' IntAct, MINT, BioGRID and HPRD
#' @param gene_name the gene name for which to retrieve PPI
#' @return a data.frame PPI information
#' @import utils
#' @export
create_summary_table_PPI <- function(gene_name){
  
  cat("Fetching PPi from databases...\n")
  pb <- txtProgressBar(min = 0, max = 3, style = 3)
  
  df_psicquic <- try(get_PPI_from_psicquic(gene_name = gene_name), silent = FALSE)
  setTxtProgressBar(pb, 1)
  df_biogrid <- get_PPI_from_BioGRID(gene_name = gene_name)
  setTxtProgressBar(pb, 2)
  df_HPRD <- get_PPI_from_HPRD(gene_name = gene_name)
  setTxtProgressBar(pb, 3)
  close(pb)
  
  df_tot <- rbind(df_biogrid, df_psicquic, df_HPRD)
  
  uInteractors <- sort(unique(toupper(c(as.character(df_tot$gene_name_A), as.character(df_tot$gene_name_B) ) )))
  uInteractors <- uInteractors[ uInteractors != toupper(gene_name) ];
  
  N_pub_0 <- rep(0,length(uInteractors));
  Authors_0 <- rep("",length(uInteractors));
  Pubmed_ID_0 <- rep("",length(uInteractors));
  Detection_method_0 <- rep("",length(uInteractors));
  Int_type_0 <- rep("",length(uInteractors));
  Database_0 <- rep("",length(uInteractors));
  
  for (i in 1:length(uInteractors) ){
    
    i_int <- which( toupper(as.character(df_tot$gene_name_A)) == uInteractors[i] | toupper(as.character(df_tot$gene_name_B)) == uInteractors[i]  )
    
    Authors_0[i] <- paste(as.character(unique(df_tot$Author[i_int])), collapse="|")
    
    Pubmed_ID_0[i] <- paste(as.character(unique(df_tot$Pubmed_ID[i_int])), collapse=",")
    spl <- strsplit(Pubmed_ID_0[i], split=",");
    Pubmed_ID_0[i] <- paste(spl[[1]], collapse="|");
    N_pub_0[i] <- length(unique(spl[[1]]));
    
    #N_pub[idx_int] <- length(unique(df_tot$pubmed_ID[i_int]));
    
    Detection_method_0[i] <- paste(as.character(unique(df_tot$Detection_method[i_int])), collapse="|")
    Int_type_0[i] <- paste(as.character(unique(df_tot$Int_type[i_int])), collapse="|")
    Database_0[i] <- paste(as.character(unique(df_tot$Database[i_int])), collapse="|")
    
  }
  
  df_summary <-data.frame(gene_name_A=rep(toupper(gene_name),length(uInteractors)), 
                          gene_name_B=uInteractors, 
                          N_pub = N_pub_0, 
                          Authors = Authors_0, 
                          Pubmed_ID = Pubmed_ID_0,
                          Detection_method = Detection_method_0,
                          Int_type = Int_type_0,
                          Database= Database_0)
  
  return(df_summary)
}

#' Append protein-protein interaction 
#' @description Append protein-protein interaction information to an \code{InteRactome}.
#' PPI are retrieved from databases IntAct, MINT, BioGRID and HPRD
#' @param res an \code{InteRactome}
#' @param mapping name of the \code{InteRactome}'s variable containing gene names
#' @param df_summary data.frame with PPI information obatined from a call to \code{create_summary_table_PPI()}
#' @return an \code{InteRactome}
#' @export
append_PPI <- function( res, mapping = "names", df_summary = NULL){
  
  res_int <- res
  
  if(is.null(df_summary)){
    df_ppi <- create_summary_table_PPI( res$bait )
  }else{
    if(length(unique(df_summary$gene_name_A))>1 | toupper(df_summary$gene_name_A[1]) != toupper(res$bait)){
      warning("incorrect PPI data")
      return(res)
    }else{
      df_ppi <- df_summary
    }
    
  }
  
  uInteractors <- df_ppi$gene_name_B
  
  sh_preys<-intersect(toupper(res[[mapping]]), uInteractors)
  sh_preys<-setdiff(sh_preys, toupper(res$bait))
  
  N_pub <- rep(0,length(res$names));
  Authors <- rep("",length(res$names));
  Pubmed_ID <- rep("",length(res$names));
  Detection_method <- rep("",length(res$names));
  Int_type <- rep("",length(res$names));
  Database <- rep("",length(res$names));
  
  if(length(sh_preys)>0){
    
    for (i in 1:length(sh_preys) ){
      
      idx_int <- which(toupper(res[[mapping]]) == sh_preys[i] );
      i_int <- which( toupper(as.character(df_ppi$gene_name_A)) == sh_preys[i] | toupper(as.character(df_ppi$gene_name_B)) == sh_preys[i]  )
      
      Authors[idx_int] <- paste(as.character(unique(df_ppi$Authors[i_int])), collapse="|")
      
      Pubmed_ID[idx_int] <- paste(as.character(unique(df_ppi$Pubmed_ID[i_int])), collapse="|")
      spl <- strsplit(Pubmed_ID[idx_int], split="|", fixed=TRUE);
      
      N_pub[idx_int] <- length(unique(spl[[1]]));
      
      #N_pub[idx_int] <- length(unique(df_ppi$pubmed_ID[i_int]));
      
      Detection_method[idx_int] <- paste(as.character(unique(df_ppi$Detection_method[i_int])), collapse="|")
      Int_type[idx_int] <- paste(as.character(unique(df_ppi$Int_type[i_int])), collapse="|")
      Database[idx_int] <- paste(as.character(unique(df_ppi$Database[i_int])), collapse="|")
      
    }
    
  }
  
  res_int$N_pub <- N_pub
  res_int$Authors <- Authors
  res_int$Pubmed_ID <- Pubmed_ID
  res_int$Detection_method <- Detection_method
  res_int$Int_type <- Int_type
  res_int$Database <- Database
  
  return(res_int)
  
}