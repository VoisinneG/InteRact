#' Append annotations to an \code{InteRactome}
#' @param res an \code{InteRactome}
#' @param annotations type of annotations to append
#' @param name_id column name used to map protein identifiers
#' @param organism organism for which the annotations have to be appended
#' @import pannot
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
    
    df_annot <- pannot::get_annotations(res, name_id = name_id, organism = organism)
    
    for (annot in annotations){
      df_annot <- switch(annot,
                         "GO" =  add_GO_data(df_annot, GO_type = "molecular_function"),
                          df_annot
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