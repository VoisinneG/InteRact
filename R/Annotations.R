#' Append annotations to an \code{InteRactome}
#' @param res an \code{InteRactome}
#' @param annotations type of annotations to append
#' @param name_id column name from \code{res} used to map protein identifiers
#' @param name_id_annot column name from the annotation data.frame used to map protein identifiers
#' @param sep Character separating different proteoin ids
#' @import pannot
#' @return an \code{InteRactome}
#' @export
append_annotations <- function(res, 
                               annotations=c("keywords", "families", "go"), 
                               name_id = "Protein.IDs",
                               sep = ";",
                               name_id_annot = "Protein.IDs"){
  
  res_int<-res
  if(is.null(annotations)){
    warning("No annotations to append")
    return(res_int)
  }else if(is.data.frame(annotations)){
    df_annot <- annotations
  } else{
    df_annot <- pannot::get_annotations_uniprot(id = res[[name_id]], 
                                                sep = ";",
                                                columns = annotations)
    df_annot[[name_id_annot]] <- res[[name_id]]
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
      idx_match[i] <- which(as.character(df_annot[[name_id_annot]]) == as.character(res[[name_id]][i]) )
    }
    
    for( var_names in setdiff(names(df_annot), name_id_annot) ){
      res_int[[var_names]] <- as.character(df_annot[[var_names]][idx_match])
    }
    cat("Done.\n")
  }
  
  return(res_int)
  
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