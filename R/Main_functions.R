#' Analysis of AP-MS data
#' @param df A dataframe containing protein intensities. By default, protein intensity column names start by "Intensity." 
#' (use parameter \code{Column_intensity_pattern} to change)
#' @param updateProgress function to show progress bar in shiny app
#' @param N_rep Number of iterations for the replacement of missing values
#' @param method Method to replace missing values. Methods from the "mice" package are supported. By default, 
#' missing values are sampled from a normal distribution centered on the quantile of ctrl intensities defined by parameter \code{quantile_rep}
#' with the standard deviation set to the mean SD of ctrl intensities across all proteins.
#' @param quantile_rep Numeric value between 0 and 1. Quantile of the distribution of mean intensities 
#' in the control background used to replace missing values.
#' @param pool_background option to use all control background conditions as one control group for all conditions
#' @param log_test logical, perform t-test on log transform intensities
#' @param log_stoichio logical, use the geometric mean instead of the arithmetic mean to compute stoichiometries
#' @param log_mean logical, use the geometric mean instead of the arithmetic mean to compute the mean \code{InteRactome}
#' @param substract_ctrl logical, substract ctrl intensities in the calculation of stoichiometries
#' @param use_mean_for_bait logical, average bait intensities across all conditions to compute interaction stoichiometries
#' @param by_conditions option to perform the comparison between bait and control group for each condition
#' @param preprocess_df list obtained by  the function \code{preprocess_data()}
#' @param ... Additional parameters passed to function \code{preprocess_data()} and \code{identify_conditions}
#' @return a list containing the preprocessed data and on object of class \code{InteRactome}, i.e a list including the following elements :
#' @return \code{conditions} : a vector of experimental conditions.
#' @return \code{names} : a vector of names (by default gene names are used).
#' @return \code{p_val} : a list of vectors containing the p values associated to each experimental condition.
#' @return \code{fold_change} : a list of vectors containing the fold change associated to each experimental condition. 
#' @return \code{...} : other variables.
#' @importFrom stats quantile rnorm
#' @import mice
#' @export
#' @author Guillaume Voisinne
#' @examples
#' #load data :
#' data("proteinGroups_Cbl")
#' #Run InteRact with default parameters
#' res <- InteRact(proteinGroups_Cbl, bait_gene_name = "Cbl")
#' 
#' #You now have an \code{InteRactome}. See its elements.
#' class(res)
#' names(res)
#' #Generate volcano plots
#' plot_volcanos(res)
#' 
#' #Identify specific interactors
#' res <- identify_interactors(res, p_val_thresh = 0.05, fold_change_thresh = 2)
#' 
#' #Visualize interaction kinetics
#' plot_per_condition(res)
#' 
#' # Append protein abundance information
#' res <- merge_proteome(res)
#' # Append annotations
#' annot <- get_annotations(res)
#' res <- append_annotations(res,  annot)
#' #Create a summary data frame
#' sum_tbl <- summary_table(res)
InteRact <- function(
  df,
  updateProgress = NULL,
  N_rep=3,
  method = "default",
  quantile_rep  = 0.05,
  pool_background = FALSE, 
  log_test = TRUE,
  log_stoichio = TRUE,
  log_mean = TRUE,
  by_conditions = TRUE,
  substract_ctrl = TRUE,
  use_mean_for_bait = FALSE,
  preprocess_df = NULL,
  ...
){
  if(is.null(preprocess_df)){
    avg <- preprocess_data(df=df, ...)
  } else {
    avg <- preprocess_df
  }
  
  # identify missing values
  
  log10_I_norm_mean <- log10(avg$Intensity);
  q <- stats::quantile(log10_I_norm_mean[ , avg$conditions$bckg == avg$bckg_ctrl], na.rm=TRUE, probs=quantile_rep);
  s <- mean( row_sd(log10_I_norm_mean[ , avg$conditions$bckg == avg$bckg_ctrl]), na.rm=TRUE);
  
  log10_I_norm_mean_rep <- log10_I_norm_mean
  
  # replace missing values N_rep times
  
  if(N_rep>0){
    
    res <- vector("list", N_rep)
    names(res)<-paste( rep('Rep_',N_rep), 1:N_rep, sep="" )
    cat(paste("Replace missing values and perform interactome analysis for",N_rep,"replicates\n",sep=" "))
    
    n_replace <- length(which(is.na(log10_I_norm_mean)));
    
    for(i in 1:N_rep){
      
      if (is.function(updateProgress)) {
        text <- paste0( i, "/", N_rep)
        updateProgress(value = i/N_rep*100, detail = text)
      }
      
      cat(paste("Nrep=",i,"\n",sep=""));
      
      if(method == "default"){
        log10_I_norm_mean_rep[is.na(log10_I_norm_mean)] <- stats::rnorm( n_replace, mean=q, sd=s) 
      }
      else{
        imp <- mice::mice(log10_I_norm_mean, printFlag = FALSE, method = method)
        log10_I_norm_mean_rep <- mice::complete(imp)
      }
      
      res[[i]] <- analyse_interactome(Intensity_na_replaced = 10^log10_I_norm_mean_rep,
                                      Intensity = avg$Intensity,
                                      conditions = avg$conditions,
                                      bait_gene_name = avg$bait_gene_name, 
                                      Npep = avg$Npep,
                                      Protein.IDs = avg$Protein.IDs, 
                                      names = avg$names,
                                      bckg_bait = avg$bckg_bait, 
                                      bckg_ctrl = avg$bckg_ctrl,
                                      pool_background = pool_background, 
                                      log_test = log_test, 
                                      log_stoichio = log_stoichio,
                                      by_conditions = by_conditions,
                                      substract_ctrl = substract_ctrl,
                                      use_mean_for_bait =use_mean_for_bait)
      
    }
    
    res_mean = mean_analysis(res, log = log_mean, na.rm = TRUE);
    
    
  }else{
    
    res_mean <- analyse_interactome(Intensity_na_replaced = 10^log10_I_norm_mean_rep,
                                    Intensity = avg$Intensity,
                                    conditions = avg$conditions,
                                    bait_gene_name = avg$bait_gene_name, 
                                    Npep = avg$Npep,
                                    Protein.IDs = avg$Protein.IDs,
                                    names = avg$names,
                                    bckg_bait = avg$bckg_bait, 
                                    bckg_ctrl = avg$bckg_ctrl,
                                    pool_background = pool_background, 
                                    log_test = log_test, 
                                    log_stoichio = log_stoichio,
                                    by_conditions = by_conditions,
                                    substract_ctrl = substract_ctrl,
                                    use_mean_for_bait = use_mean_for_bait)
    
    
  }
  
  res_mean <- global_analysis(res_mean)
  
  res_mean$params <- c( list(N_rep=N_rep,
                             method = method,
                             quantile_rep  = quantile_rep),
                        res_mean$params)
  
  res_mean$data <- c(res_mean$data, 
                     list(Intensity_filter = avg$Intensity_filter, 
                          conditions_filter = avg$conditions_filter))
  
  return(res_mean)
  
}


#' Identify indirect interactions by comparing stoichiometries between two interactomes.
#' @param resA reference interactome
#' @param resB intermediate interactome.
#' @param conditions set of conditions for which stoichiometries will be compared. If \code{NULL} then all shared conditions will be used.
#' @export

identify_indirect_interactions <- function(resA, resB, conditions = NULL){
  
  shared_cond <- intersect(resA$conditions, resB$conditions)
  common_interactors <- setdiff(intersect(resA$interactor, resB$interactor), c(resA$bait, resB$bait) )
  
  # Verify that the analysis can be performed
  if(!(resB$bait %in% resA$interactor)){
    cat("B is not an interactor of bait A\n")
    return(NULL)
  }
  if(length(common_interactors) == 0){
    cat("Empty set of shared interactors\n")
    return(NULL)
  }
  
  #identify conditions for which interaction stoichiometries will be compared
  
  if(is.null(conditions)){
    if(length(shared_cond)==0){
      cond <- "max"
      warning("No shared condition. Analysis will be performed using the maximum stoichiometry across conditions.")
    }else{
      cond <- shared_cond
    }
  }else{
    if( length(intersect(conditions, shared_cond)) == 0){
      cond <- "max"
      warning("Specified conditions are not shared by both interactomes. Analysis will be performed using the maximum stoichiometry across conditions.")
    }else{
      cond <- conditions
    }
  }
  
  #cat(common_interactors)
  #cat(cond)
  
  idxA <- match(common_interactors, resA$names)
  idxB <- match(common_interactors, resB$names)
  idxB_in_A <- match(resB$bait, resA$names)
  
  
  stA <- data.frame(do.call(cbind, resA$stoichio),  row.names = resA$names)
  stB <- data.frame(do.call(cbind, resB$stoichio),  row.names = resB$names)
  names(stA)<- resA$conditions
  names(stB)<- resB$conditions
  
  stoichio_direct <- sapply(shared_cond, function(x){stA[idxA, x]} )
  if(length(common_interactors)==1){
    stoichio_direct <- t(stoichio_direct)
  }
  rownames(stoichio_direct) <- common_interactors
  stoichio_direct <- data.frame(stoichio_direct, check.names = FALSE)
  
  stoichio_indirect <- sapply(shared_cond, function(x){stB[idxB, x]*stA[idxB_in_A, x]} )
  if(length(common_interactors)==1){
    stoichio_indirect <- t(stoichio_indirect)
  }
  rownames(stoichio_indirect) <- common_interactors
  stoichio_indirect <- data.frame(stoichio_indirect, check.names = FALSE)
  
  stoichio_C_in_B <- sapply(shared_cond, function(x){stB[idxB, x]} )
  if(length(common_interactors)==1){
    stoichio_C_in_B <- t(stoichio_C_in_B)
  }
  rownames(stoichio_C_in_B) <- common_interactors
  stoichio_C_in_B <- data.frame(stoichio_C_in_B, check.names = FALSE)
  
  perc_C_in_B <- sapply(shared_cond, function(x){stB[idxB, x]*resB$Copy_Number[resB$names == resB$bait]/resB$Copy_Number[idxB]} )
  if(length(common_interactors)==1){
    perc_C_in_B <- t(perc_C_in_B)
  }
  rownames(perc_C_in_B) <- common_interactors
  perc_C_in_B <- data.frame(perc_C_in_B, check.names = FALSE)
  
  delta_stoichio_log <- log10(stoichio_direct) - log10(stoichio_indirect)
  max_delta_stoichio_log <- apply(abs(delta_stoichio_log), MARGIN = 1, max)
  sum_delta_stoichio_log <- apply(abs(delta_stoichio_log), MARGIN = 1, sum)
  
  return(list(interactors = common_interactors,
              conditions = shared_cond,
              stoichio_direct = stoichio_direct, 
              stoichio_indirect=stoichio_indirect,
              stoichio_C_in_B=stoichio_C_in_B,
              perc_C_in_B = perc_C_in_B,
              delta_stoichio_log = delta_stoichio_log,
              max_delta_stoichio_log = max_delta_stoichio_log,
              sum_delta_stoichio_log = sum_delta_stoichio_log,
              bait_A = resA$bait,
              bait_B = resB$bait)
  )
  
}

#' Preprocessing of raw data
#' @param df Data.frame with protein intensities
#' @param bait_gene_name The gene name of the bait
#' @param Column_score Column with protein identification score
#' @param Column_ID Column with protein IDs
#' @param Column_Npep Column with number of theoretically observable peptides per protein
#' @param Column_gene_name Column with gene names
#' @param Column_intensity_pattern Pattern (regular exrpression) used to identfy df's columns containing protein intensity values
#' @param condition data.frame with columns "column", bckg", "bio", "time" and "tech" indicating 
#' for each intensity column ("sample") its corresponding background ("bckg"), biologicla replicate ("bio), 
#' experimental condition ("tine) and technical replicate ("tech).
#' @param bckg_bait Name of the bait background as found in \code{condition$bckg} (see below)
#' @param bckg_ctrl Name of the control background as found in \code{condition$bckg} (see below)
#' @param log logical, use geometric mean to average technical replicates
#' @param filter_time vector of experimental conditions to exclude from analysis
#' @param filter_bio vector of biological replicates to exclude from analysis
#' @param filter_tech vector of technical replicates to exclude from analysis
#' @param min_score threshold on identification score
#' @param filter_gene_name logical, filter out proteins withy empty gene name
#' @param ... Additional parameters passed to function \code{identify_conditions}
#' @export
preprocess_data <- function(df,
                            Column_gene_name = "Gene.names",
                            Column_score = "Score",
                            Column_ID = "Protein.IDs",
                            Column_Npep = NULL,
                            Column_intensity_pattern = "^Intensity.",
                            bait_gene_name,
                            condition = NULL,
                            bckg_bait = bait_gene_name,
                            bckg_ctrl = "WT",
                            log = TRUE,
                            filter_time=NULL,
                            filter_bio=NULL,
                            filter_tech=NULL,
                            min_score = 0, 
                            filter_gene_name = TRUE, 
                            ...            
){
  
  #if( sum( sapply( grep(Column_intensity_pattern,names(df)), function(x) is.factor( df[, x] ) ) ) >0 ){
  #  stop("Some intensity columns are factors, try changing the decimal separator (most likely '.' or ',') used for importing the data")
  #}
  
  if(! Column_ID %in% names(df)){
    warning(paste("Column ", Column_ID, " could not be found", sep=""))
  }
  if(! Column_gene_name %in% names(df)){
    stop(paste("Column ", Column_gene_name, " could not be found", sep=""))
  }
  
  # Identify conditions corresponding to intensity columns
  
  if( is.null(condition) ){
    
    cond <- identify_conditions(df, 
                                Column_intensity_pattern = Column_intensity_pattern, 
                                ...)
    
    
  } else if( length(setdiff(c("column", "bckg", "bio", "time", "tech"), names(condition))) > 0 ) {
    stop("incorrect dimensions for data.frame condition")
  } else {
    cond <- dplyr::tibble(column = condition$column,
                          bckg = condition$bckg,
                          bio = condition$bio,
                          time = condition$time,
                          tech = condition$tech)
  }
  
  # filter out some experimental conditions
  
  cond_filter <- cond[!is.na(cond$time), ]
  
  idx_filter <- c( unlist( lapply(filter_time, function(x) l=which(cond_filter$time==x) ) ) , 
                   unlist( lapply(filter_bio, function(x) l=which(cond_filter$bio==x) ) ),
                   unlist( lapply(filter_tech, function(x) l=which(cond_filter$tech==x) ) ) )
  
  if(!is.null(idx_filter) && length(idx_filter)>0 ){
    cat("Filter following intensity columns :\n")
    cat(idx_filter)
    cat("\n")
    cond_filter <- cond[-idx_filter,] 
  }
  
  # match conditions to samples
  idx_match <- match(cond_filter$column, names(df))
  if( sum(is.na(idx_match)) == length(idx_match) ){
    stop("Could not match conditions to samples")
  }
  cond_filter <- cond_filter[!is.na(idx_match), ]
  
  #print(cond_filter$column[is.na(idx_match)])
  
  #discard columns that are factors
  is_factor <- sapply(match(cond_filter$column, names(df)), function(x){is.factor(df[[x]])})
  
  if( sum(is_factor) == length(is_factor) ){
    stop("All intensity columns are factors, try changing the decimal separator (most likely '.' or ',') used for importing the data")
  }else if(sum(is_factor) > 0) {
    warning("Some intensity columns identified as factors were discarded. Check imported data")
    df <- df[ , -match(cond_filter$column, names(df))[is_factor] ]
    cond_filter <- cond_filter[!is_factor, ]
  }
  
  #format gene name
  #df$gene_name <- sapply(df[[Column_gene_name]], function(x) strsplit(as.character(x),split=";")[[1]][1] )
  
  # Filter protein with no gene name and a low score
  df <- filter_Proteins(df, 
                        min_score = min_score, 
                        filter_gene_name = filter_gene_name, 
                        Column_ID = Column_ID, 
                        Column_gene_name = Column_gene_name, 
                        Column_score = Column_score)
  
  # Compute number of theoretically observable peptides
  df$Npep <- estimate_Npep(df, Column_Npep = Column_Npep)
  
  #Merge protein groups with the same gene name
  idx_col = match(cond_filter$column, names(df))
  idx_col <- idx_col[!is.na(idx_col)]
  df <- merge_duplicate_groups(df, idx_col = idx_col, merge_column = "gene_name")
  
  
  
  #identify bait
  ibait <- which(df$gene_name == bait_gene_name);
  if(length(ibait)==0){
    stop(paste("Could not find bait '",bait_gene_name,"' in column '",Column_gene_name,"'", sep="")) 
  }
  
  #Select intensity columns corresponding to selected conditions
  T_int <- df[ , idx_col];
  T_int[T_int==0] <- NA;
  
  #Discard proteins with NA values for all conditions
  idx_all_na <- which( rowSums(!is.na(T_int)) == 0 )
  if(length(idx_all_na)>0){
    T_int <- T_int[-idx_all_na, ]
    df <- df[-idx_all_na, ]
  }
  row.names(T_int) <- df$gene_name
  
  #Normalize on median intensity across conditions
  T_int_norm <- rescale_median(T_int);
  
  cat("Rescale median intensity across conditions\n")
  
  avg <- average_technical_replicates(T_int_norm, cond_filter, log = log)
  
  avg$Intensity_filter <- T_int_norm
  avg$conditions_filter <- cond
  avg$Npep <- df$Npep
  avg$Protein.IDs <- df[[Column_ID]]
  avg$names <- df$gene_name
  row.names(avg$Intensity) <- avg$names
  avg$bckg_bait <- bckg_bait
  avg$bckg_ctrl <- bckg_ctrl
  avg$bait_gene_name <- bait_gene_name
  
  return(avg)
  
}

#' Identify conditions (background, time of stimulation, biological and technical replicates) 
#' from column names
#' @param df A dataframe containing protein intensities. By default, protein intensity column names start by "Intensity." 
#' (use parameter \code{Column_intensity_pattern} to change)
#' @param Column_intensity_pattern Pattern (regular exrpression) used to identfy df's columns containing protein intensity values
#' @param split Character used to split column names into substrings
#' @param bckg_pos Position of the sample background in splitted column names
#' @param bio_pos Position of the sample biological replicate in splitted column names
#' @param time_pos Position of the sample experimental condition in splitted column names
#' @param tech_pos Position of the sample technical replicate in splitted column names
#' @return a data frame describing experimental samples in terms of background, 
#' biological and technical replicates, and experimental conditions
#' @importFrom dplyr tibble
#' @export
#' @examples
#' #load data :
#' data("proteinGroups_Cbl")
#' # You can identify columns and their description separately using \code{identify_conditions()}
#' cond <- identify_conditions(proteinGroups_Cbl)
#' print.data.frame(cond)
#' # and use it as parameters for function InteRact()
#' res <- InteRact(proteinGroups_Cbl, bait_gene_name = "Cbl", condition = cond)
identify_conditions <- function(df,
                                Column_intensity_pattern = "^Intensity.",
                                split = "_",
                                bckg_pos = 1,
                                bio_pos = 3,
                                time_pos = 2,
                                tech_pos = 4
){
  
  idx_col<-grep(Column_intensity_pattern,colnames(df))
  if(length(idx_col)==0){
    stop("Couldn't find pattern in column names")
  }
  
  T_int <- df[ ,idx_col];
  
  col_I <- colnames(T_int)
  
  s0 <- sapply(strsplit(col_I, Column_intensity_pattern), function(x){x[2]})
  s <- strsplit(s0, split=split, fixed=TRUE)
  
  bckg <- rep("", length(col_I))
  bio <- rep("", length(col_I))
  tech <- rep("", length(col_I))
  time <- rep("", length(col_I))
  #s <- strsplit(col_I, split=split, fixed=TRUE)
  
  for (i in 1:length(col_I)){
    n <- length(s[[i]])
    
    if(bckg_pos <= n){ bckg[i] <- s[[i]][bckg_pos]}
    if(bio_pos <= n){ bio[i] <-s[[i]][bio_pos] }
    if(tech_pos <= n){ tech[i] <- s[[i]][tech_pos] }
    if(time_pos <= n){ time[i] <- s[[i]][time_pos] }
    
  }
  
  bckg[nchar(bckg)==0] <- paste("bckg", "0", sep=split)
  bio[nchar(bio)==0] <- paste("bio", "0", sep=split)
  tech[nchar(tech)==0] <- paste("tech", "0", sep=split)
  time[nchar(time)==0] <- paste("time", "0", sep=split)
  
  
  cond <- dplyr::tibble(column=col_I, bckg, time, bio, tech)
  
}


#' Average protein intensities over technical replicates
#' @param df A data frame of protein intensities
#' @param cond A data frame containing the description of df's columns (i.e "bckg", "time", "bio"  and "tech") 
#' as returned by function \code{identify_conditions()}
#' @param log use geometric mean
#' @return A list containing :
#' @return \code{Intensity}, a data frame of protein intensities averaged over technical replicates;
#' @return \code{conditions}, a data frame containing the description of \code{Intensity}'s columns
#' @importFrom dplyr group_by summarise
#' @export
#' @examples
#' #load data :
#' data("proteinGroups_Cbl")
#' df <- proteinGroups_Cbl
#' cond <- identify_conditions(df)
#' Column_intensity_pattern <- "^Intensity."
#' df_int <- df[ , grep(Column_intensity_pattern, colnames(df))]
#  avg <- average_technical_replicates(df_int, cond)
average_technical_replicates<-function(df, cond, log = TRUE){
  
  cond$idx_match <- match(cond$column, names(df))
  cond_group <- dplyr::group_by(cond, bckg, time, bio)
  idx_cond <-  dplyr::summarise(cond_group, idx_tech=list(idx_match))
  
  cond_name <- vector("character", dim(idx_cond)[1])
  df_mean = data.frame( matrix( NA, nrow = dim(df)[1], ncol=dim(idx_cond)[1] ) );
  
  for(j in 1:dim(idx_cond)[1]){
    cond_name[j] = paste( idx_cond$bckg[j], 't', idx_cond$time[j], 'rep', idx_cond$bio[j], sep="_");
    df_mean[[j]] <- row_mean(df[ idx_cond$idx_tech[[j]] ], na.rm=TRUE, log = log);
  }
  colnames(df_mean)=cond_name;
  idx_cond$column <- cond_name
  idx_cond$tech <- rep("tech_0", dim(idx_cond)[1])
  
  output = list(Intensity=df_mean, conditions=idx_cond)
  
}

#' Filtering of a data frame using a threshold on protein identification score and 
#' gene names
#' @param df A data frame
#' @param min_score Threshold for protein identification score
#' @param filter_gene_name logical, filter out proteins withy empty gene name
#' @param Column_gene_name The name of df's column containing gene names
#' @param Column_ID Column with protein IDs
#' @param Column_score The name of df's column containing protein identification score
#' @param split_param Character used to split gene names into substrings. 
#' @return A filtered data frame. Contains an extra column with the first substring of the column \code{Column_gene_name}
#' @export
filter_Proteins <- function( df, 
                             min_score=0, 
                             filter_gene_name =TRUE, 
                             Column_ID = "Protein.IDs", 
                             Column_gene_name = "Gene.names", 
                             Column_score= "Score", 
                             split_param=";"){
  
  idx_row = 1:dim(df)[1]
  if( nchar(Column_score)>0 & Column_score %in% colnames(df)){
    if(class(df[[Column_score]]) == "numeric"){
      idx_row = which( df[[Column_score]] > min_score )
      df<-df[idx_row, ]
      cat("Data Filtered based on portein identification score\n")
    }else{
      warning("Column 'Score' is not numeric : Data NOT Filtered based on portein identification score\n")
    }
  }else{
    warning("Column 'Score' not available : Data NOT Filtered based on portein identification score\n")
  }
  
  #Remove contaminants from dataset
  
  if( Column_gene_name %in% colnames(df) & Column_ID %in% colnames(df)){
    
    idx_cont <- unique( c(grep("^KRT",toupper(df[[Column_gene_name]])),
                          grep("^CON_",toupper(df[[Column_ID]])),
                          grep("^REV_",toupper(df[[Column_ID]]))
    ))
    
    if (length(idx_cont) > 0){
      df<-df[ -idx_cont, ]
      cat("Contaminant proteins discarded\n")
    }
    
    df$gene_name <- as.character(df[[Column_gene_name]])
    
    length_name <- nchar(df$gene_name)    
    
    idx_no_name <- which( length_name == 0  | is.na(length_name) )
    
    if(length(idx_no_name) > 0){
      df$gene_name[idx_no_name] <- as.character(df[[Column_ID]][idx_no_name])
    }
    
    
    if(filter_gene_name){
      if(length(idx_no_name) > 0){
        df <- df[ -idx_no_name, ]
        cat("Proteins with no gene name available discarded\n")
      }
    }
    
    df$gene_name <- sapply(df$gene_name, function(x){ strsplit(as.character(x),split=split_param)[[1]][1]} )
    
  }else{
    warning(paste("Column_gene_name '", Column_gene_name, "' or Column_ID '", Column_ID,"' not available",sep=""))
  }
  
  
  
  return(df)
}

#' Merge protein groups with the same gene name.
#' @param df A data frame
#' @param idx_col idx of columns for which values will be merged across protein groups
#' @param merge_column column to identify rows to be be merged
#' @param sum_intensities Logical. Sum intensities across merged groups.
#' @return A merged data frame 
merge_duplicate_groups <- function(df, idx_col = NULL, merge_column = "gene_name", sum_intensities = TRUE){
  
  cat("Merge protein groups associated to the same gene name (sum of intensities) \n")
  
  df_int <- df
  
  ugene <- unique(df[[merge_column]]);
  
  idx_merge = rep(1, dim(df)[1]); # rows with idx_merge = 0 will be filtered out
  
  for(i in 1:length(ugene) ){
    
    idx_u <- which( df[[merge_column]] == ugene[i] ) 
    
    if (length(idx_u) > 1) {
      max_I <- rep(0, length(idx_u));
      for (j in 1:length(idx_u) ){
        idx_merge[idx_u[j]] = 0;
        max_I[j] = max(df[ idx_u[j], idx_col], na.rm = TRUE)
      }
      jmax = which(max_I == max(max_I) );
      idx_merge[ idx_u[jmax[1]] ] = 1;
      
      for (k in idx_col ){
        s=0;
        for(j in idx_u ){
          s= s + df[j, k];
        }
        if(sum_intensities){
          df_int[idx_u[jmax[1]], k] = s;
        }
        
      }
      
    }
    
  }
  
  return(df_int[idx_merge>0, ])
  
}

#' Get the number of theoretically observable peptides per protein
#' @param df A data frame
#' @param Column_Npep column containing the number of theoretically observable peptides per protein.
#' If NULL try to compute the number of theoretically observable peptides using iBAQ values, 
#' or use molecular weight.
#' @return A data frame with the column 'Npep'
estimate_Npep <- function(df, Column_Npep = NULL){
  
  if ( is.null(Column_Npep) ){
    if( "Intensity" %in% colnames(df) & "iBAQ" %in% colnames(df)){
      Npep <- round( as.numeric( as.character(df$Intensity))/as.numeric( as.character(df$iBAQ)) )
      cat("Number of theoretically observable peptides computed using iBAQ values\n")
    }else if("Mol..weight..kDa." %in% colnames(df) ){ 
      Npep <- df$Mol..weight..kDa.
      cat("Number of theoretically observable peptides unavailable : used MW instead\n")
    }else{
      Npep <- rep(1,dim(df)[1])
      cat("Number of theoretically observable peptides unavailable : set to 1\n")
    }
  } else {
    Npep <- df[[Column_Npep]]
  }
  
  output<-Npep
}

#' Normalize data frame by columns using the median
#' @param df A data frame
#' @importFrom stats median
#' @return A normalized data frame
rescale_median <- function(df){
  df_out<-df;
  for (i in seq_along(df)){
    df_out[[i]] <- df[[i]]/median(df[[i]], na.rm=TRUE);
  }
  df_out
}

#' Perform the geometric mean of a numeric vector
#' @param x A numeric vector
#' @param na.rm remove NA values
#' @return A numeric value
geom_mean = function(x, na.rm = TRUE){
  idx_na <- which(is.na(x))
  x <- x[x>0]
  if(sum(!is.na(x)) == 0){
    return(NA)
  }else if(!na.rm & sum(is.na(x)) >0){
    return(NA)
  }else{
    return( exp(sum(log(x[!is.na(x)])) / sum(!is.na(x)) ))
  }
  
}

#' Compute the standard deviation by row
#' @param df a data frame
#' @importFrom stats sd
#' @return A numeric vector
row_sd <- function(df){
  output<-vector("double", dim(df)[1] )
  for(i in 1:dim(df)[1] ){
    output[i] <- sd(as.numeric(df[i,]),na.rm=TRUE);
  }
  output
}

#' Compute the mean by row
#' @param df a data frame
#' @param log logical, use geometric mean instead of arithmetic mean
#' @param na.rm logical, remove NA values
#' @return A numeric vector
#' @export
row_mean <- function(df, na.rm = TRUE, log = FALSE){
  output <- vector("double", dim(df)[1] )
  for(i in 1:dim(df)[1] ){
    if(log){
      output[i] <- geom_mean(as.numeric(df[i,]), na.rm = na.rm);
    } else{
      output[i] <- mean(as.numeric(df[i,]), na.rm = na.rm);
    }
    
  }
  output
}

#' Perform a t-test comparison between two groups by row
#' 
#' @param df a data frame
#' @param idx_group_1 column indexes corresponding to the first group
#' @param idx_group_2 column indexes corresponding to the second group
#' @param log option to perform the t-test on log transformed data
#' @importFrom stats t.test
#' @return A data frame with columns 'p_val' and 'fold_change' (group_1 vs group_2)
row_ttest <- function(df, idx_group_1, idx_group_2, log = TRUE){
  
  p_val <- rep(NaN,dim(df)[1]);
  fold_change <- rep(NaN,dim(df)[1]);
  output <- data.frame(p_val=p_val, fold_change=fold_change);
  
  if(log){
    df_test<-log10(df)
  }else{
    df_test<-df
  }
  
  for(i in (1:dim(df)[1]) ){
    res<-try(t.test(df_test[i, idx_group_1], df_test[i, idx_group_2]), silent=TRUE)
    if(!inherits(res,"try-error")){
      p_val[i] <- res$p.value;
      if(log){
        fold_change[i] <- 10^(res$estimate[1] - res$estimate[2])
      } else {
        fold_change[i] <- res$estimate[1] / res$estimate[2]
      }
    }
  }
  output$p_val <- p_val;
  output$fold_change <- fold_change;
  output
  
}

#' Compute the stoichiometry of interaction using the method described in ...
#' 
#' @param df a data frame
#' @param idx_group_1 column indexes corresponding to the first group (bait background)
#' @param idx_group_2 column indexes corresponding to the second group (ctrl background)
#' @param idx_bait row index for the bait protein
#' @param Npep numeric vector containing the number of theoretically observable peptides for each protein
#' @param log logical, use the geometric mean instead of the arithmetic mean
#' @param substract_ctrl logical, substract ctrl intensities in the calculation of stoichiometries
#' @param use_mean_for_bait logical, average bait intensities across all conditions to compute interaction stoichiometries
#' @param idx_group_1_mean if \code{use_mean_for_bait} is TRUE, column indexes for the bait protein corresponding to the first group (bait background)
#' @param idx_group_2_mean if \code{use_mean_for_bait} is TRUE, column indexes for the bait protein corresponding to the second group (ctrl background)
#' @return A numeric vector of interaction stoichiometries
row_stoichio <- function(df, 
                         idx_group_1, 
                         idx_group_2, 
                         idx_bait, 
                         Npep, 
                         log = TRUE, 
                         substract_ctrl = TRUE,
                         use_mean_for_bait = TRUE,
                         idx_group_1_mean,
                         idx_group_2_mean){
  # compute stoichiometry of interaction for each row (protein).
  # idx_group_1 : indexes of columns for group #1 (OST bait)
  # idx_group_2 : indexes of columns for group #2 (WT)
  # idx_bait : row index for the bait
  # N_pep : vector with the number of theoretical observable peptides for each row (protein).
  
  stoichio <- rep(NaN,dim(df)[1]);
  
  if(use_mean_for_bait){
    xbait1<-df[idx_bait, idx_group_1_mean]
    xbait2<-df[idx_bait, idx_group_2_mean]
  }else{
    xbait1<-df[idx_bait, idx_group_1]
    xbait2<-df[idx_bait, idx_group_2]
  }
  
  
  for(i in (1:dim(df)[1]) ){
    x1<-df[i, idx_group_1]
    x2<-df[i, idx_group_2]
    if (log) {
      if (substract_ctrl){
        stoichio[i] <- ( geom_mean(x1) - geom_mean(x2) ) / ( geom_mean(xbait1) - geom_mean(xbait2) )*Npep[idx_bait]/Npep[i]
      } else {
        stoichio[i] <- ( geom_mean(x1) ) / ( geom_mean(xbait1) )*Npep[idx_bait]/Npep[i]
      }
      
    } else {
      if (substract_ctrl){
        stoichio[i] <- ( mean(x1) - mean(x2) ) / ( mean(xbait1) - mean(xbait2) )*Npep[idx_bait]/Npep[i]
      } else {
        stoichio[i] <- ( mean(x1) ) / ( mean(xbait1)  )*Npep[idx_bait]/Npep[i]
      }
      
    }
    
  }
  stoichio
}

#' Construct an interactome by comparing bait and control background across experimental conditions
#'
#' @param Intensity a data frame of protein intensities (without replacement of missing values). columns are experimental samples and rows are proteins
#' @param Intensity_na_replaced a data frame of protein intensities (with replacement of missing values). columns are experimental samples and rows are proteins
#' @param conditions : data frame describing the conditions corresponding to the columns of the data frame \code{Intensity}
#' @param bait_gene_name : The gene name of the bait
#' @param Npep : vector containing the number of theoretically observable peptide per protein (same length as \code{dim(df)[1]})
#' @param Protein.IDs : vector containing protein IDs (same length as \code{dim(df)[1]})
#' @param names : vector containing protein names (same length as \code{dim(df)[1]})
#' @param bckg_bait : Name of the bait background as found in \code{conditions$bckg}
#' @param bckg_ctrl : Name of the control background as found in \code{conditions$bckg}
#' @param pool_background option to use all control background conditions as one control group for all conditions
#' @param log_test logical, perform t-test on log transform intensities
#' @param log_stoichio logical, use the geometric mean instead of the arithmetic mean to compute stoichiometries
#' @param substract_ctrl logical, substract ctrl intensities in the calculation of stoichiometries
#' @param use_mean_for_bait logical, average bait intensities across all conditions to compute interaction stoichiometries
#' @param by_conditions option to perform the comparison between bait and control group for each condition
#' @return an object of class \code{InteRactome}, i.e a list including the following elements :
#' @return \code{conditions} : a vector of experimental conditions.
#' @return \code{names} : a vector of names (by default gene names are used).
#' @return \code{p_val} : a list of vectors containing the p values associated to each experimental condition.
#' @return \code{fold_change} : a list of vectors containing the fold change associated to each experimental condition. 
#' @return \code{...} : other variables.
analyse_interactome <- function(Intensity,
                                Intensity_na_replaced = Intensity,
                                conditions,
                                bait_gene_name,
                                Npep,
                                Protein.IDs,
                                names,
                                bckg_bait, 
                                bckg_ctrl,
                                by_conditions = TRUE, pool_background = TRUE, 
                                log_test = TRUE, log_stoichio = TRUE,
                                substract_ctrl = TRUE,
                                use_mean_for_bait = TRUE){
  
  
  if(!by_conditions){
    conds<-rep(-1,length(conditions$time))
  }else{
    conds<-conditions$time
  }
  cond<-unique(conds)
  background <- conditions$bckg
  
  p_val <- vector("list",length(cond))
  fold_change <- vector("list",length(cond))
  stoichio <- vector("list",length(cond))
  
  names(p_val)<-as.character(cond)
  names(fold_change)<-as.character(cond)
  names(stoichio)<-as.character(cond)
  
  replicates <- conditions$bio
  ubio <- unique(replicates)
  stoichio_bio <- vector("list",length(ubio))
  names(stoichio_bio) <- as.character(ubio)
  
  for (i_bio in 1:length(ubio)){
    stoichio_bio[[i_bio]] <- vector("list",length(cond))
    names(stoichio_bio[[i_bio]])<-as.character(cond)
  }
  
  for( i in seq_along(cond) ){
    
    idx_ctrl <- which( background == bckg_ctrl )
    if(!pool_background) {
      idx_ctrl <- which( background == bckg_ctrl & conds == cond[[i]] )
    } 
    
    ttest <- row_ttest(Intensity_na_replaced, 
                       idx_group_1 = which( background == bckg_bait & conds==cond[[i]]), 
                       idx_group_2 = idx_ctrl, 
                       log = log_test)
    p_val[[i]] <- ttest$p_val;
    fold_change[[i]] <- ttest$fold_change;
    
    stoichio[[i]] <- row_stoichio(Intensity_na_replaced, 
                                  idx_group_1 = which( background==bckg_bait & conds==cond[[i]] ), 
                                  idx_group_2 = idx_ctrl, 
                                  idx_bait = which(names == bait_gene_name),
                                  Npep=Npep,
                                  log = log_stoichio,
                                  substract_ctrl = substract_ctrl,
                                  use_mean_for_bait = use_mean_for_bait,
                                  idx_group_1_mean = which( background == bckg_bait ),
                                  idx_group_2_mean = which( background == bckg_ctrl )
    )
    for (i_bio in 1:length(ubio)){
      
      idx_ctrl_bio <- which( background == bckg_ctrl & replicates == ubio[i_bio])
      if(!pool_background) {
        idx_ctrl_bio <- which( background == bckg_ctrl & conds == cond[[i]] & replicates == ubio[i_bio])
      } 
      
      idx_bait_bio <- which( background==bckg_bait & conds==cond[[i]] & replicates==ubio[i_bio])
      stoichio_bio[[i_bio]][[i]] <- row_stoichio(Intensity_na_replaced, 
                                                 idx_group_1 = idx_bait_bio, 
                                                 idx_group_2 = idx_ctrl_bio, 
                                                 idx_bait = which(names == bait_gene_name),
                                                 Npep=Npep,
                                                 log = log_stoichio,
                                                 substract_ctrl = substract_ctrl,
                                                 use_mean_for_bait = use_mean_for_bait,
                                                 idx_group_1_mean = which( background == bckg_bait ),
                                                 idx_group_2_mean = which( background == bckg_ctrl )
      )
      # intensity_bait[[i]][[i_bio]] <- as.numeric(df[ , idx_bait_bio])
      # intensity_ctrl[[i]][[i_bio]]<- as.numeric(df[ , idx_ctrl_bio])
      # intensity_na_bait[[i]][[i_bio]] <- as.numeric(Intensity[ , idx_bait_bio])
      # intensity_na_ctrl[[i]][[i_bio]]<- as.numeric(Intensity[ , idx_ctrl_bio])
    }
  }
  
  res = list(bait = bait_gene_name, 
             bckg_bait = bckg_bait,
             bckg_ctrl = bckg_ctrl,
             conditions = as.character(cond),
             replicates = as.character(ubio),
             names = names,
             Protein.IDs = Protein.IDs,
             Npep = Npep,
             p_val=p_val,
             fold_change=fold_change,
             stoichio=stoichio,
             stoichio_bio = stoichio_bio,
             # intensity_bait = intensity_bait,
             # intensity_ctrl = intensity_ctrl,
             # intensity_na_bait = intensity_na_bait,
             # intensity_na_ctrl = intensity_na_ctrl,
             data = list(Intensity_na_replaced = Intensity_na_replaced,
                         Intensity = Intensity,
                         conditions = conditions
             ),
             params = list(by_conditions = by_conditions,
                           pool_background = pool_background,
                           log_test = log_test, 
                           log_stoichio = log_stoichio,
                           substract_ctrl = substract_ctrl,
                           use_mean_for_bait = use_mean_for_bait
             )
  )
  
  class(res) <- 'InteRactome'
  
  return(res)
}

#' Performs a running average on a numeric vector
#' @param x a numeric vector
#' @param n integer, radius of the moving avergae (number of points extending on each side of the 
#' center point on which the average is computed)
#' @return a smoothed numeric vector
moving_average <- function(x, n){
  
  x_smooth <- x
  
  for (i in 1:length(x)){
    idx <- max(c(i-n, 1)):min(c(i+n, length(x))) 
    x_smooth[i] <- mean(x[idx], na.rm=TRUE)
  }
  
  return(x_smooth)
}

#' Smooth, using a moving average across conditions, selected variables of an \code{InteRactome}
#' @param res an \code{InteRactome}
#' @param n integer, radius of the moving avergae (number of points extending on each side of the 
#' center point on which the average is computed)
#' @param order_conditions a numeric vector ordering conditions in \code{res$conditions}
#' @param var_smooth variables on which the moving average will be computed
#' @return an smoothed \code{InteRactome}
#' @export
smooth_interactome <- function( res,  n = 1, order_conditions = NULL, var_smooth = c("fold_change","p_val") ){
  
  res_smooth <- res
  
  if(!is.null(order_conditions)){
    idx_order <- order_conditions
  } else {
    idx_order <- 1:length(res$conditions)
  }
  
  for (var in var_smooth){
    M <- do.call(cbind, res[[var]])
    M_smooth <- M
    
    for (i in 1:dim(M)[1] ){
      M_smooth[i, idx_order] <- 10^(moving_average(log10(M[i, idx_order]), n))
    }
    
    for (i in 1:length(res$conditions)){
      res_smooth[[var]][[res$conditions[i]]] <- M_smooth[ , i]
    }
  }
  
  
  if( "norm_stoichio" %in% names(res) ){
    res_smooth <- global_analysis(res_smooth)
  }
  
  return(res_smooth)
}


#' Normalize the log fold change by its standard deviation for each condition
#' @param res an \code{InteRactome}
#' @return an \code{InteRactome} with the additional variable \code{norm_log_fold_change}
#' @export
normalize_interactome <- function( res ){
  
  res_norm <- res
  res_norm[["norm_log_fold_change"]] <- vector("list", length = length(res$conditions))
  
  for (i in 1:length(res$conditions)){
    sd_norm <- sd(log(res[["fold_change"]][[res$conditions[i]]]))
    res_norm[["norm_log_fold_change"]][[res$conditions[i]]] <- log(res[["fold_change"]][[res$conditions[i]]])/sd_norm
  }
  
  return(res_norm)
}

#' Merge different conditions from different interactomes into a single data.frame 
#' @param res a list of \code{InteRactomes}
#' @param selected_conditions a character vector containing names of conditions to merge
#' @return a data.frame with columns bait, names, Protein.IDs, conditions, p_val, 
#' @return fold_change and norm_log_fold_change
#' @export
merge_conditions <- function( res,  selected_conditions = NULL){
  
  
  df_merge <- NULL
  
  if (class(res) == "list"){
    
    for (i in 1:length(res)){
      if(class(res[[i]]) == "InteRactome"){
        if(is.null(selected_conditions)){
          conditions <- res[[i]]$conditions
        } else {
          conditions <- selected_conditions
        }
        
        res_int <- normalize_interactome(res[[i]])
        
        for (cond in conditions){
          names <- res_int$names
          df <- data.frame(
            bait = rep(res_int$bait, length(names)),
            names = names,
            Protein.IDs = res_int$Protein.IDs,
            conditions = rep(cond, length(names)), 
            p_val = res_int$p_val[[cond]],
            fold_change = res_int$fold_change[[cond]],
            norm_log_fold_change = res_int$norm_log_fold_change[[cond]]
          )
          
          df_merge <- rbind(df_merge, df)
        }
      } else {
        stop("Input is not of class 'InteRactome'")
      }
    }
  } else{
    if(class(res) == "InteRactome"){
      if(is.null(selected_conditions)){
        conditions <- res$conditions
      } else {
        conditions <- selected_conditions
      }
      
      res_int <- normalize_interactome(res)
      
      for (cond in conditions){
        names <- res_int$names
        df <- data.frame(
          bait = rep(res_int$bait, length(names)),
          names = names,
          Protein.IDs = res_int$Protein.IDs,
          conditions = rep(cond, length(names)), 
          p_val = res_int$p_val[[cond]],
          fold_change = res_int$fold_change[[cond]],
          norm_log_fold_change = res_int$norm_log_fold_change[[cond]]
        )
        
        df_merge <- rbind(df_merge, df)
      }
    } else {
      stop("Input is not of class 'InteRactome'")
    }
    
  }
  return(df_merge)
}

#' Compute the FDR from the asymmetry of the volcano plot
#' @description Compute the FDR (False Discovery Rate) using the asymmetry of the volcano plot.
#' It uses the fonction f(x) = c / (|x|-x0) with x = log10(fold_change), y=-log10(p_value) by default.
#' Otherwise, custom x and y vectors can be provided.
#' Points with x>x0 and y>f(x) are taken as true positive (TP)
#' Points with x<x0 and y>f(x) are taken as false positive (FP)
#' For a given set of parameters (c,x0), the FDR is given by TP/(TP+FP)
#' @param df : a data.frame containing columns \code{p_val} and \code{fold_change}
#' @param x : numeric vector of x values. Only used if \code{df} is NULL.
#' @param y : numeric vector of y values. Only used if \code{df} is NULL.
#' @param c : numeric vector
#' @param x0 : numeric vector
#' @return  a data.frame with a extra column \code{FDR} if \code{df} is not NULL. A vector of FDR values otherwise.  
#' @return If parameters \code{c} and \code{x0} are vectors, \code{FDR} is taken as the minimum FDR value across all sets of parameters
#' @import utils
#' @export
#' @examples
#' #' #load data :
#' data("proteinGroups_Cbl")
#' #Run InteRact with default parameters
#' res <- InteRact(proteinGroups_Cbl, bait_gene_name = "Cbl")
#' df_merge <- merge_conditions(res)
#' df_FDR <- compute_FDR_from_asymmetry(df_merge)
#' Interactome <- append_FDR(res, df_FDR)
compute_FDR_from_asymmetry <- function( df = NULL,
                                        x = NULL,
                                        y = NULL,
                                        c = seq(from = 0, to = 10, by = 0.1),
                                        x0 = seq(from = 0, to =10, by = 0.1)){
  if(!is.null(df)){
    df_int<-df
  }else if(is.null(x) & is.null(y)){
    stop("No input data provided.")
  }
  
  
  if(is.null(x)){
    if("fold_change" %in% names(df)){
      x <- log10(df$fold_change)
    }else{
      stop("Variable 'fold_change' not available")
    }
    
  }
  if(is.null(y)){
    if("p_val" %in% names(df)){
      y <- -log10(df$p_val)
    }else{
      stop("Variable 'p_val' not available")
    }
  }
  
  if(length(x) != length(y)){
    stop("Lengths of x and y differ.")
  }
  
  FDR <- rep(1, length(x))
  cat("Compute FDR...\n")
  pb <- txtProgressBar(min = 0, max = length(c)*length(x0), style = 3)
  count <- 0
  
  c_tot <- NULL
  x0_tot <- NULL
  TP_tot <- NULL
  FP_tot <- NULL
  FDR_tot <- NULL
  
  for (i in 1:length(c)){
    for (j in 1:length(x0)){
      idx_TP <- which(x>x0[j] & y>c[i]/(x-x0[j]) )
      idx_FP <- which(x<(-x0[j]) & y>c[i]/(-x-x0[j]) )
      TP_int <- length(idx_TP)
      FP_int <- length(idx_FP)
      FDR_int <- FP_int/(FP_int + TP_int)
      
      FDR[idx_TP[ FDR[idx_TP] >= FDR_int ]] <- FDR_int
      
      c_tot <- c(c_tot, c[i])
      x0_tot <- c(x0_tot, x0[j])
      TP_tot <- c(TP_tot, TP_int)
      FP_tot <- c(FP_tot, FP_int)
      FDR_tot <- c(FDR_tot, FDR_int)
      
      count <- count +1
      setTxtProgressBar(pb, count)
    }
  }
  
  df_params = data.frame(c = c_tot, 
                         x0 = x0_tot, 
                         TP = TP_tot, 
                         FP = FP_tot, 
                         FDR = FDR_tot)
  close(pb)
  
  uFDR <- unique(df_params$FDR)
  uFDR <- uFDR[order(uFDR)]
  idx_max <- rep(NA, length(uFDR))
  for(i in 1:length(uFDR)){
    if(!is.na(uFDR[i])){
      idx_FDR <- which(df_params$FDR <= uFDR[i] & !is.na(df_params$FDR) )
      idx_max[i] <- idx_FDR[which.max(df_params$TP[idx_FDR])]
    }
  }
  df_max_TP <- data.frame( FDR = uFDR, 
                           c = df_params$c[idx_max], 
                           x0 = df_params$x0[idx_max],
                           TP = df_params$TP[idx_max],
                           FP = df_params$FP[idx_max]) 
  
  return(list(FDR = FDR, 
              parameters = df_params,
              max_TP_parameters = df_max_TP ))
}


#' Append a FDR column to an \code{InteRactome}
#' @param res an \code{InteRactome}
#' @param df a data.frame containing (at least) columns 'bait', names', 'FDR' and 'conditions'
#' @return an \code{InteRactome}
#' @export
#' @examples
#' #' #load data :
#' data("proteinGroups_Cbl")
#' #Run InteRact with default parameters
#' res <- InteRact(proteinGroups_Cbl, bait_gene_name = "Cbl")
#' df_merge <- merge_conditions(res)
#' df_FDR <- compute_FDR_from_asymmetry(df_merge)
#' Interactome <- append_FDR(res, df_FDR)
append_FDR <- function(res, df){
  
  
  res_int <- res
  df_int <- df
  
  df_int <- df_int[which(df_int$bait == res$bait), ]
  
  FDR <- list()
  
  for (cond in res$conditions){
    
    df_cond <- df_int[df_int$conditions == cond, ]
    idx_match <- match(res$names, df_cond$names)
    FDR[[cond]] <- df_cond$FDR[idx_match]
    
  }
  
  res_int$FDR <-FDR
  res_int <- global_analysis(res_int)
  
  return(res_int)
  
}

#' Filter conditions from an interactome
#' @param res an \code{InteRactome}
#' @param conditions_to_filter_out character vector with names of conditions to filter out
#' @return an \code{InteRactome}
#' @export
filter_conditions <- function( res, conditions_to_filter_out ){
  
  idx_conditions <- which(is.element(res$conditions, conditions_to_filter_out ))
  res_filter<-res
  res_filter$conditions <- res$conditions[-idx_conditions]
  
  for(i in 1:length(res)){
    if( setequal( names(res[[i]]), res$conditions ) ){
      res_filter[[i]] <- res[[i]][-idx_conditions]
    }
  }
  if( "norm_stoichio" %in% names(res) ){
    res_filter <- global_analysis(res_filter)
  }
  return(res_filter)
}

#' Compute the mean \code{InteRactome} (on variables 'p_val', 'fold_cahnge', 'stoichio' and 'stoichio_bio')
#' from a list of \code{InteRactomes}
#' @param res a list of \code{InteRactome}
#' @param log logical, use the geometric mean instead of the arithmetic mean
#' @param na.rm logical, remove NA values
#' @export
mean_analysis <- function(res, log = TRUE, na.rm = TRUE ){
  
  # Performs the average of different interactomes
  
  cat(paste("Averaging", length(res) ,"interactomes\n",sep=" ") )
  res_mean = res[[1]];
  
  for ( cond in seq_along(res_mean$conditions) ){
    p_val=vector("list",length(res));
    fold_change=vector("list",length(res));
    stoichio = vector("list",length(res));
    for ( i in seq_along(res)){
      p_val[[i]] <- res[[i]]$p_val[[cond]]
      fold_change[[i]] <- res[[i]]$fold_change[[cond]]
      stoichio[[i]] <- res[[i]]$stoichio[[cond]]
    }
    if (log){
      res_mean$p_val[[cond]] <- 10^( rowMeans( log10(do.call(cbind, p_val)), na.rm =na.rm) )
      res_mean$fold_change[[cond]] <- 10^( rowMeans( log10(do.call(cbind, fold_change)), na.rm =na.rm) )
      res_mean$stoichio[[cond]] <- 10^ (rowMeans( log10(do.call(cbind, stoichio)), na.rm =na.rm) )
    } else {
      res_mean$p_val[[cond]] <-  rowMeans( do.call(cbind, p_val), na.rm =na.rm) 
      res_mean$fold_change[[cond]] <- rowMeans( do.call(cbind, fold_change), na.rm =na.rm) 
      res_mean$stoichio[[cond]] <- rowMeans( do.call(cbind, stoichio), na.rm =na.rm)
    }
    
    
    for (bio in res_mean$replicates){
      stoichio_bio <- vector("list",length(res));
      for ( i in seq_along(res)){
        stoichio_bio[[i]] <- res[[i]]$stoichio_bio[[bio]][[cond]]
      }
      if (log){
        res_mean$stoichio_bio[[bio]][[cond]] <- 10^( rowMeans( log10(do.call(cbind, stoichio_bio)), na.rm =na.rm) )
      } else {
        res_mean$stoichio_bio[[bio]][[cond]] <- rowMeans( do.call(cbind, stoichio_bio), na.rm =na.rm)
      }
    }
    
  }
  
  res_mean
  
}

#' Adds global variables by analysing values
#' across all conditions of an \code{InteRactome}
#' @param res an \code{InteRactome}
#' @return an \code{InteRactome} with global varaiables
#' @export
global_analysis <- function( res ){
  
  res_int <- res;
  max_stoichio<- apply( do.call(cbind, res$stoichio), 1, function(x)  ifelse( sum(!is.na(x))>0, max(x,na.rm=TRUE),NaN) )
  max_stoichio[!is.finite(max_stoichio)]<-NaN
  res_int$max_stoichio <- max_stoichio
  
  matrix_fold_change<-do.call(cbind, res$fold_change)
  res_int$max_fold_change <- apply(matrix_fold_change , 1, function(x)  ifelse( sum(!is.na(x))>0, max(x,na.rm=TRUE),NaN) )
  
  matrix_p_val<-do.call(cbind, res$p_val)
  matrix_p_val[matrix_fold_change<1]<-NA # keep only p_values corresponding to a fold_change >1 
  res_int$min_p_val <- apply( matrix_p_val, 1, function(x)  ifelse( sum(!is.na(x))>0, min(x,na.rm=TRUE),NaN) )
  
  if( "FDR" %in% names(res) ){
    matrix_FDR<-do.call(cbind, res$FDR)
    res_int$min_FDR <- apply( matrix_FDR, 1, function(x)  ifelse( sum(!is.na(x))>0, min(x,na.rm=TRUE),NaN) )
  }
  
  norm_stoichio = vector("list",length(res$conditions))
  names(norm_stoichio)<-res$conditions
  for ( i in seq_along(res$conditions) ){
    norm_stoichio[[i]] <- res$stoichio[[i]] / max_stoichio
  }
  res_int$norm_stoichio <- norm_stoichio
  
  output=res_int
  
}

#' Add protein abundance to an \code{InteRactome}
#' @description Add protein abundance to an \code{InteRactome}. 
#' Protein abundance are obtained from CD4+ effector T cells.
#' @param res an \code{InteRactome}
#' @export
merge_proteome <- function( res ){
  
  res_int <- res
  
  gene_name_prot <- proteome_data$Gene.names;
  
  ibait <- match(res$bait, res$names)
  
  ######### Retrieve protein abundance and compute related quantities
  
  Copy_Number <- rep(0, length(res$names));
  
  for( i in 1:length(res$names) ){
    idx_prot <- which(gene_name_prot==as.character(res$names[i]));
    idx_ID_x <-  which(proteome_data$Protein.IDs.x==as.character(res$Protein.IDs[i]));
    idx_ID_y <-  which(proteome_data$Protein.IDs.y==as.character(res$Protein.IDs[i]));
    
    #idx_ID_x <-  grep(as.character(res[[Interactome_ID_name]][i]),  as.character(proteome_data$Protein.IDs.x), fixed=TRUE);
    #idx_ID_y <-  grep(as.character(res[[Interactome_ID_name]][i]),  as.character(proteome_data$Protein.IDs.y), fixed=TRUE);
    # 
    abund_x <- NA;
    abund_y <- NA;
    if( length(idx_ID_x)>0 ){ 
      abund_x <- proteome_data$mean.x[idx_ID_x] 
    }
    if( length(idx_ID_y)>0 ){ 
      abund_y <- proteome_data$mean.y[idx_ID_y] 
    }
    
    #Copy_Number[i] <- mean(c(abund_x, abund_y), na.rm = TRUE)
    
    if(length(idx_prot)>0){
      Copy_Number[i] = proteome_data$mean[ idx_prot[1] ];
    }else if( length(idx_ID_x)>0 || length(idx_ID_y)>0 ){
      Copy_Number[i] = mean( c(abund_x, abund_y), na.rm=TRUE);
    }else{
      Copy_Number[i] = NA;
    }
    
  }
  res_int$Copy_Number = Copy_Number
  res_int$stoch_abundance = Copy_Number / Copy_Number[ibait]
  
  #res_int$N_complex= Tsum$max_stoch*Copy_Number[ibait]
  #res_int$percentage_prey_in_complex = Tsum$max_stoch*Copy_Number[ibait]/Copy_Number;
  
  # Tsum$Perc_t_0 = Tsum$Stoch_t_0*Copy_Number[ibait]/Copy_Number;
  # Tsum$Perc_t_30 = Tsum$Stoch_t_30*Copy_Number[ibait]/Copy_Number;
  # Tsum$Perc_t_120 = Tsum$Stoch_t_120*Copy_Number[ibait]/Copy_Number;
  # Tsum$Perc_t_300 = Tsum$Stoch_t_300*Copy_Number[ibait]/Copy_Number;
  # Tsum$Perc_t_600 = Tsum$Stoch_t_600*Copy_Number[ibait]/Copy_Number;
  
  output=res_int
}

#' Identify specific interactors in an \code{InteRactome}
#' @param res an \code{InteRactome}
#' @param var_p_val name of the p-value variable
#' @param p_val_thresh p-value threshold
#' @param fold_change_thresh fold-change threshold
#' @param conditions Select conditions used to identify interactors
#' @param n_success_min minimal number of conditions in which the interactor 
#' must pass the the p-value and the fold-change thresholds
#' @param consecutive_success logical, impose that the interactor must pass selection thresholds 
#' in \code{n_success_min} consecutive conditions.
#' @param ... additionnal paramters passed to function \code{order_interactome()}
#' @return an \code{InteRactome} with extra variables \code{is_interactor}, 
#' \code{n_success} and \code{interactor}
#' @export
identify_interactors <- function(res, 
                                 var_p_val = "p_val", 
                                 p_val_thresh = 0.05, 
                                 fold_change_thresh = 2,
                                 conditions = res$conditions,
                                 n_success_min = 1, 
                                 consecutive_success = FALSE,
                                 p_val_breaks = c(1, min(1,4*p_val_thresh), min(1, 2*p_val_thresh), min(1, p_val_thresh)),
                                 var_min_p_val = paste("min_",var_p_val,sep=""),
                                 ...){
  
  res_int <- res
  
  n_cond <- length(conditions)
  is_interactor <- rep(0, length(res$names))
  n_success <- rep(0, length(res$names))
  
  M <- do.call(cbind, res[[var_p_val]][conditions]) <= p_val_thresh & do.call(cbind, res[["fold_change"]][conditions]) >= fold_change_thresh
  
  for (i in 1:length(res$names)){
    
    M_test <- M[i, ]
    n_success[i] <- sum( M_test, na.rm = TRUE )
    
    if (consecutive_success & n_success_min > 1){
      for (k in 1:(n_success_min-1)){
        idx_mod <- ((1:n_cond) + k) %% n_cond
        idx_mod[ idx_mod == 0] <- n_cond
        M_test <- rbind(M_test, M[i, idx_mod])
      }
      is_interactor[i] <- sum( colMeans( M_test ) == 1 ) > 0
    } else {
      is_interactor[i] <- n_success[i] >= n_success_min
    }
    
  }
  
  res_int$is_interactor <- is_interactor
  res_int$n_success <- n_success
  res_int$interactor <- res$names[is_interactor>0]
  
  res_int <- order_interactome(res_int, 
                               idx = NULL, 
                               p_val_breaks = p_val_breaks,
                               var_min_p_val = var_min_p_val,
                               ...)
  
  res_int$params$var_p_val <- var_p_val
  res_int$params$p_val_thresh <- p_val_thresh
  res_int$params$fold_change_thresh <- fold_change_thresh
  res_int$params$conditions <- conditions
  res_int$params$n_success_min <- n_success_min 
  res_int$params$consecutive_success <- consecutive_success
  
  return(res_int)
  
}

#' Discretize values in a vector based on a finite set of values
#' @param x numeric vector
#' @param breaks numeric vector. Set of discrete values on which \code{x} values will be mapped. 
#' Non-mapped values will be set to NA
#' @param decreasing_order logical. Map \code{beaks} values from the greatest to the smallest
#' @return a numeric vector
#' @export
discretize_values <- function( x, breaks = c(1,0.1,0.05,0.01), decreasing_order = TRUE){
  
  breaks_order <- breaks[order(breaks, decreasing=decreasing_order)]
  
  x_discrete <- rep(NA, length(x));
  
  for( i in 1:length(breaks_order) ){
    x_discrete[ x <= breaks_order[i] ] <- breaks_order[i];
  }
  
  return(x_discrete)
}

#' Order proteins within an \code{InteRactome}
#' @param res an \code{InteRactome}
#' @param idx indices used to order proteins. Overrides ordering using \code{var_p_val}, \code{p_val_breaks} and \code{var_order}
#' @param var_min_p_val name of the p-value variable
#' @param p_val_breaks numeric vector to discretize p-value
#' @param var_order Variable used to order interactors for each p-value
#' @param bait_first logical, puts bait in first position
#' @return an \code{InteRactome}
#' @export
order_interactome <- function(res, idx = NULL, 
                              var_min_p_val = "min_p_val", 
                              p_val_breaks = c(1,0.1,0.05,0.01), 
                              var_order = "max_stoichio",
                              bait_first = TRUE){
  
  min_p_val_discrete <- discretize_values(res[[var_min_p_val]], breaks = p_val_breaks, decreasing_order = TRUE)
  
  if(!is.null(idx)){
    idx_order <- idx
    if(length(idx_order)!=length(res$names)){
      stop("Vector of ordering indexes does not have the proper length")
    }
  } else if( "interactor" %in% names(res)){
    Ndetect<-length(res$interactor)
    idx_order<-order(res$is_interactor, 1/min_p_val_discrete, res[[var_order]], decreasing = TRUE)
  } else{
    stop("idx not provided")
  }
  
  if(bait_first){
    idx_bait <- which(res$names == res$bait)
    idx_order_bait <- which(idx_order == idx_bait)
    idx_order <- c(idx_bait, idx_order[-idx_order_bait])
  }
  
  res_order<-res;
  for( var in setdiff( names(res), c("bait", "bckg_bait", "bckg_ctrl","conditions", "interactor", "replicates", "data", "params") ) ){
    
    
    if(length(res[[var]]) == length(res$names)){
      res_order[[var]] <- res[[var]][idx_order]
    }
    else{
      names_var <- names(res[[var]])
      if( setequal(names_var, res$conditions) ){
        for(i in 1:length(names_var) ){
          res_order[[var]][[i]] <- res[[var]][[i]][idx_order]
        }
      }
      else{
        for(i in 1:length(names_var) ){
          names_var_2 <- names(res[[var]][[i]]) 
          #if( setequal(names_var_2, res$conditions) ){
          for(j in 1:length(names_var_2) ){
            if(length(res[[var]][[i]][[j]]) == length(res$names)){
              res_order[[var]][[i]][[j]] <- res[[var]][[i]][[j]][idx_order]
            } else if (dim(res[[var]][[i]][[j]])[1] == length(res$names)){
              res_order[[var]][[i]][[j]] <- res[[var]][[i]][[j]][idx_order, ]
            } else {
              warning(paste("Couldn't order ",var, sep=""))
            }
          }
          #}
        }
      }
    }
    
  }
  
  
  output = res_order
  
  
}