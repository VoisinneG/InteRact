#' Plot indirect interactions
#' @param score output of the \code{identify_indirect_interactions()} function
#' @param var_threshold Variable of \code{score} on which the threshold will be applied
#' @param threshold maximum difference between observed and predicted stoichiometries (in log10)
#' @param min_range Minimum range (in log10 scale) displayed on the y-axis (centered on the mean)
#' @param save_dir path to the directory where the plot will be saved
#' @param plot_width set plot width
#' @param plot_height set plot height
#' @param show_legend logical, shows plot legend
#' @param theme_name name of the ggplot2 theme function to use ('theme_gray' by default)
#' @export
plot_indirect_interactions <- function(score,
                                       var_threshold = "max_delta_stoichio_log",
                                       threshold = 1.0, 
                                       min_range = NULL,
                                       save_dir = NULL,
                                       plot_width = 2, 
                                       plot_height = 2,
                                       show_legend = FALSE,
                                       theme_name = "theme_gray"
                                       ){
  theme_function <- function(...){
    do.call(theme_name, list(...))
  }
  
  idx_select <- which(!is.na(score[[var_threshold]]))
  if(length(idx_select) > 0){
    if(sum(score[[var_threshold]][idx_select] <= threshold, na.rm = TRUE) == 0){
      #p <- NULL
      plist <- NULL
    }else{
      
      idx_plot <- idx_select[ score[[var_threshold]][idx_select] <= threshold]
      plist <- list()
      
      for(i in 1:length(idx_plot)){
        name_interactor <- score$interactor[idx_plot[i]]
        score_display <- score[[var_threshold]][idx_plot[i]]
        s_direct <- score$stoichio_direct[idx_plot[i], ]
        s_indirect <- score$stoichio_indirect[idx_plot[i], ]
        df_stoichio <- data.frame(type = c( rep("direct", length(s_direct)),
                                            rep("reciprocal", length(s_indirect))),
                                  stoichio = c(as.numeric(s_direct), as.numeric(s_indirect)),
                                  conditions = c(names(s_direct), names(s_indirect)))
        
        
          log10_stoichio <- log10(df_stoichio$stoichio[df_stoichio$stoichio>0 & !is.na(df_stoichio$stoichio)])
          datarange <- range( log10_stoichio )
          
          if(!is.null(min_range)){
            if(abs(datarange[2]-datarange[1]) <= min_range){
              ymin <- mean(log10_stoichio) - min_range/2
              ymax <- mean(log10_stoichio) + min_range/2
            }else{
              ymin <- datarange[1]
              ymax <- datarange[2]
            }
          }else{
            ymin <- datarange[1]
            ymax <- datarange[2]
          }
          
        df_stoichio$log10_stoichio <- log10(df_stoichio$stoichio)
        
        plist[[i]] <- ggplot(df_stoichio, aes_string(x='conditions', y='log10_stoichio', color='type', group='type')) +
          theme_function() +
          theme(axis.text.x = element_text(angle=90, hjust = 1),
                title = element_text(size = 6) ) +
          ggtitle(paste(score$bait_A,"<",score$bait_B,"<", name_interactor," (", signif(score_display, 3), ")", sep="")) +
          geom_point(show.legend = FALSE) +
          geom_line(show.legend = show_legend) +
          coord_cartesian(ylim = c(ymin, ymax))
          
          
        if(!is.null(save_dir)){
          pdf(paste(save_dir,"A_",
                    score$bait_A,"_B_",
                    score$bait_B,"_C_",
                    name_interactor,
                    ".pdf",sep=""), 
              width = plot_width, 
              height = plot_height)
          print(plist[[i]])
          dev.off()
        }
        
      }
    }
  }
  return( plist )
  
}      

#' Plot abundance versus interaction stoichiometries
#' @param res an \code{InteRactome}
#' @param condition condition selected. If "max", the maximum stoichiometry across conditions will be used.
#' @param names names of the proteins to display. If not NULL, supersedes \code{N_display} and
#' \code{only_interactors}
#' @param xlim range of x values
#' @param ylim range of y values
#' @param N_display maximum number of protein to display
#' @param only_interactors display only interactors 
#' (identified using the function \code{identify_interactors()})
#' @param color_values color vector passed to \code{scale_color_manual()}
#' @param fill_values color vector passed to \code{scale_fill_manual()}
#' @param shape Point shape aesthetics passed to \code{geom_point()}
#' @param stroke Point stroke aesthetics passed to \code{geom_point()}
#' @param p_val_thresh Threshold on p-value used to identify regulated interactions
#' @param fold_change_thresh Threshold on fold-change used to identify regulated interactions
#' @param ref_condition Reference condition used to identify regulated interactions
#' @param label_size_min minimum label size (between 0 and \code{label_size_max})
#' @param label_size_max maximum label size (a threshold on log10(fold-change))
#' @param label_size_scale_factor scale label size according to plot range (the higher the bigger the label)
#' @param label_range if NULL, scales labels according to plot range (in log10 scale).
#' @param theme_name name of the ggplot2 theme function to use ('theme_gray' by default)
#' @param show_core logical. Show core interaction area?
#' @param range_factor numeric factor to expand plot range
#' @param ratio_strong ratio of available proteins bound to bait usd to define gray-shaded area 
#' @param n_character_max max number of label characters (ignored if NULL)
#' @param ... parameters passed to \code{geom_text_repel()}
#' @return a plot
#' @import ggplot2
#' @import ggrepel
#' @export
plot_2D_stoichio <- function( res, 
                              condition = "max", 
                              names = NULL,
                              xlim = NULL, 
                              ylim = NULL,
                              N_display=30,
                              only_interactors = FALSE,
                              fill_values = c("not_regulated" = "black",
                                              "induced" = "red",
                                              "repressed" = "blue",
                                              "bait" = "yellow"),
                              color_values = c("not_regulated" = "black",
                                               "induced" = "red",
                                               "repressed" = "blue",
                                               "bait" = "black"),
                              shape = 21,
                              stroke = 1,
                              p_val_thresh = 0.05,
                              fold_change_thresh = 1,
                              ref_condition = res$conditions[1],
                              label_size_min = 1,
                              label_size_max = 3,
                              label_size_scale_factor = 25,
                              label_range = NULL,
                              theme_name = "theme_gray",
                              show_core = TRUE,
                              range_factor = 1.1,
                              ratio_strong = 0.3,
                              n_character_max = 8,
                              ...
                              
){
  
  theme_function <- function(...){
    do.call(theme_name, list(...))
  }
  
  if(!"Copy_Number" %in% names(res)){
    warning("Protein abundances not available. Please import and merge a proteome first.")
    return(NULL)
  }
  
  plist <- list()
  
  for(icond in 1:length(condition)){
    
    cond <- condition[icond]
    
    df<- data.frame( Y=log10(res$stoch_abundance), 
                     names=res$names)
    
    if(cond=="max"){
      df$X <- log10(res$max_stoichio)
      df$size <- res$max_fold_change
    }else if(cond %in% res$conditions){
      df$X <- log10(res$stoichio[[cond]])
      df$size <- res$fold_change[[cond]]
    }else{
      stop("Condition is not defined")
    }
    
    if(!is.null(names)){
      df <- df[df$names %in% names, ]
    }else{
      if(only_interactors){
        df <- df[!is.na(match(df$names, res$interactor)), ]
      }
      
      df<-df[1:min(N_display, dim(df)[1]), ]
    }
    
    
    xc <- -0.5
    yc <- 0
    rc<-1
    
    
    ylow <- NULL
    
    if(sum(is.na(df$Y))>0){
      ylow <- min(df$Y, na.rm = TRUE) - 0.25
      df$Y[is.na(df$Y)] <- ylow - 0.25
    }
   
    
    if(is.null(xlim) & is.null(ylim)){
      max_range <- max( max(df$X,na.rm=TRUE)-min(df$X,na.rm=TRUE),  max(df$Y,na.rm=TRUE)-min(df$Y,na.rm=TRUE))
      center_x <- ( max(df$X,na.rm=TRUE)+min(df$X,na.rm=TRUE) )/2
      center_y <- (max(df$Y,na.rm=TRUE) + min(df$Y,na.rm=TRUE))/2
    }else{
      max_range <- max( xlim[2] - xlim[1],  ylim[2] - ylim[1] )
      center_x <- ( xlim[2] + xlim[1] )/2
      center_y <- ( ylim[2] + ylim[1] )/2
    }

    if(is.null(label_range)){
      label_range <- max_range
    }
    
    xmin<-center_x - max_range/2*range_factor
    xmax<-center_x + max_range/2*range_factor
    ymin<-center_y - max_range/2*range_factor
    ymax<-center_y + max_range/2*range_factor
    
    
    if(is.null(ylow)){
      ylow <- ymin - 1
    }
    
    #df$size_prey <- log10(df$size)/max_range*20
    df$size_prey <- log10(df$size)*5
    df$size_label <- unlist(lapply(log10(df$size), function(x) { 
      ifelse(x>label_size_min, min(c(x,label_size_max)), label_size_min) 
      })
      )/label_range*label_size_scale_factor/label_size_max
    
    df$sat_max_fold_t0 <- rep(1,dim(df)[1])
    
    idx_plot <- which(df$X<=xmax & df$X>=xmin & df$Y<=ymax & df$Y>=ymin)
    
    
    if(cond=="max"){
      test <- compare_stoichio(res, 
                               names = df$names, 
                               ref_condition = ref_condition, 
                               test_conditions = setdiff(res$conditions, ref_condition))
      
      M_p_val <- do.call(cbind, test$p_val)
      M_fold_change <- do.call(cbind, test$fold_change)
      
      df$color <- rep("not_regulated", dim(df)[1])
      for(i in 1:length(df$names)){
        idx_pval <- which( M_p_val[i, ] <= p_val_thresh )
        if(length(idx_pval)>0){
          idx_max <- idx_pval[ which.max( abs(log10(M_fold_change[i, idx_pval]))) ]
          if(M_fold_change[i, idx_max] > fold_change_thresh){
            df$color[i] <- "induced"
          }
          if(M_fold_change[i, idx_max] <= 1/fold_change_thresh){
            df$color[i] <- "repressed"
          }
        }
      }
      df$color[df$names == res$bait] <- "bait"
      
    }else if(cond %in% res$conditions){
      test <- compare_stoichio(res, names = df$names, ref_condition = ref_condition, test_conditions = cond)
      df$p_val <- test$p_val[[cond]]
      df$fold_change <- test$fold_change[[cond]]
      
      df$color <- rep("not_regulated", dim(df)[1])
      df$color[df$p_val <= p_val_thresh & df$fold_change >= fold_change_thresh] <- "induced"
      df$color[df$p_val <= p_val_thresh & df$fold_change <= 1/fold_change_thresh] <- "repressed"
      df$color[df$names == res$bait] <- "bait"
    }
    
    df <- df[idx_plot, ]
    
    p<-ggplot(df,aes_string(x='X', y='Y', label='names')) +
      theme_function() +
      theme(aspect.ratio=1) +
      ggtitle(cond) + 
      geom_polygon(data=data.frame(x=c(ylow, xmax, xmax), 
                                   y=c(ylow, ylow, xmax)), 
                   mapping=aes_string(x='x', y='y'),
                   alpha=0.1,
                   inherit.aes=FALSE)
    
    p <- p + geom_polygon(data=data.frame(y=c(ylow, 0, ymax, ymax, ylow),
                                          x=c(ylow + log10(ratio_strong), 
                                              log10(ratio_strong), 
                                              log10(ratio_strong),
                                              xmax, 
                                              xmax)), 
                          mapping=aes_string(x='x', y='y'),
                          alpha=0.1,
                          inherit.aes=FALSE)
    
    if(show_core){
      p <- p + annotate("path",
                        x=xc+rc*cos(seq(0,2*pi,length.out=100)),
                        y=yc+rc*sin(seq(0,2*pi,length.out=100)), color=rgb(0,0,0,0.5) )
    }
    
    df$label <- df$names
      
    if(!is.null(n_character_max)){
      df$label <- unlist(lapply(as.character(df$names), function(x){
        l <- nchar(x)
        if(l > n_character_max){
          return(paste(substr(x,1,n_character_max), "...", sep = ""))
        }else{
          return(x)
        }
      }))
    }
    
    p <- p + 
      annotate("segment", x = ylow, xend = xmax, y = ylow, yend = xmax, colour = rgb(0,0,0,0.5), linetype = "dashed" ) +
      annotate("segment", x = xmin, xend = xmax, y = ylow, yend = ylow, colour = rgb(0,0,0,0.5) ) +
      xlab(expression(paste('Interaction Stoichiometry (log'[10], ')'))) +
      ylab(expression(paste('Abundance Stoichiometry (log'[10], ')'))) +
      #ylab(bquote('Abundance Stoichiometry ('~log[10]~')')) +
      geom_point(data = df,
                 mapping=aes_string(x='X', y='Y', color='color', fill = 'color'), 
                 size=df$size_prey, 
                 alpha=0.2,
                 shape = shape,
                 stroke = stroke, 
                 inherit.aes = FALSE, 
                 show.legend = FALSE) +
      coord_cartesian(xlim = c(xmin,xmax), ylim = c(ymin,ymax), expand = FALSE)+
      geom_text_repel(data = df,
                      mapping=aes_string(x='X', y='Y', label='label'),
                      size=df$size_label,
                      ...
                      #inherit.aes = FALSE, 
                      #show.legend = FALSE,
                      #force=force, 
                      #segment.size = segment.size,
                      #min.segment.length = unit(0.15, "lines"), 
                      #point.padding = NA,
                      #max.iter = 100000
                      )
    
    if(!is.null(color_values)) {
      p <- p + scale_color_manual( values = color_values) +
        scale_fill_manual( values = fill_values)
    }
    
    plist[[icond]] <- p
  }
  
  return(plist)

}


#' @importFrom grDevices dev.off pdf rgb
#' @importFrom graphics hist lines plot
plot_Intensity_histogram <- function( I, I_rep, breaks=20, save_file=NULL){
  # plot histogram of intensities for all columns in two different datasets 
  # (1st dataset is in in black, 2nd dataset is in red)
  if(!is.null(save_file)){
    pdf( save_file, 4, 4 )
  }
  
  
  for( j in seq_along(I) ){
    h<-hist( I[[j]] , breaks=breaks, plot=FALSE);
    if(j>1){
      lines(h$mids,h$density,col=rgb(0,0,0,0.2));
    }else{
      plot(h$mids, h$density, type="l",col=rgb(0,0,0,0.2),ylim=c(0,1));
    }
  }
  
  for( j in seq_along(I_rep) ){
    h<-hist( I_rep[[j]] , breaks=breaks, plot=FALSE);
    lines(h$mids,h$density,col=rgb(1,0,0,0.2));
  }
  
  if(!is.null(save_file)){
    dev.off()
  }
  
}

#' Plot protein enrichement fold-change versus p-value
#' @param res an \code{InteRactome}
#' @param data data.frame with columns 'names', 'p_val' and 'fold_change'
#' @param names Names of proteins highlighted
#' @param N_print maximum of protein labels to display
#' @param labels labels for proteins in plot. Must the same length as \code{res$names}
#' @param conditions conditions to plot
#' @param p_val_thresh threshold on p-value to display
#' @param fold_change_thresh threshold on fold-change to display
#' @param x0 parameters x0 of the line dividing the volcano plot according to \code{f(x) = c / (|x|-x0)}.
#' Ignored unless parameters \code{p_val_thresh} and \code{fold_change_thresh} are set to \code{NULL}
#' @param c parameters c of the line dividing the volcano plot according to \code{f(x) = c / (|x|-x0)}.
#' Ignored unless parameters \code{p_val_thresh} and \code{fold_change_thresh} are set to \code{NULL}
#' @param save_file path of output file (.pdf)
#' @param xlim range of x values
#' @param ylim range of y values
#' @param asinh_transform logical, display asinh(log10(p-value)) on the y-axis
#' @param norm Use normalized fold-changes
#' @param both_sides logical. Shading on right and left upper graphs.
#' @param show_thresholds Show thresholds using red lines? 
#' @param alpha_segment transparency of threshold segments
#' @param theme_name name of the ggplot2 theme function to use ('theme_gray' by default)
#' @param size dot size
#' @param alpha dot transparency
#' @param color dot color
#' @param label_size size of labels (5 by default)
#' @param n_character_max max number of label characters
#' @param ... parameters passed to \code{geom_text_repel()}
#' @return a plot
#' @importFrom grDevices dev.off pdf rgb
#' @import ggplot2
#' @import ggrepel
#' @export
plot_volcanos <- function( res=NULL,
                           data = NULL,
                           names = NULL,
                           labels=NULL, 
                           N_print=15, 
                           conditions=NULL,
                           p_val_thresh=0.05, 
                           fold_change_thresh=2,
                           x0 = NULL,
                           c = NULL,
                           save_file=NULL,
                           xlim=NULL,
                           ylim=NULL,
                           asinh_transform = TRUE,
                           norm = FALSE,
                           both_sides = FALSE,
                           show_thresholds = TRUE,
                           alpha_segment = 0.2,
                           theme_name = "theme_gray",
                           size = 0.3,
                           alpha = 0.1,
                           color = rgb(0.5, 0.5, 0.5, 0.5),
                           label_size  =5,
                           n_character_max = 8,
                           ...
                           ){
  
  theme_function <- function(...){
    do.call(theme_name, list(...))
  }
  
  if(!is.null(res)){
    res_int <- res
    if (is.null(labels)) labels=res_int$names
    if (is.null(conditions)) conditions=res_int$conditions
  }
  
  
  
  
  
  
  plist <- vector("list",length(conditions));
  
  if(is.null(data)){
    if(norm){
      if("norm_log_fold_change" %in% names(res_int)){
        xval <- do.call(cbind, res_int$norm_log_fold_change)
      }else{
        res_int <- normalize_interactome(res_int)
        xval <- do.call(cbind, res_int$norm_log_fold_change)
      }
    }else{
      xval <- log10(do.call(cbind, res_int$fold_change))
    }
    
    yval <- -log10(do.call(cbind,res_int$p_val))
  }else{
      conditions <- ""
      xval <- log10(data$fold_change)
      yval <- -log10(data$p_val)
  }
  
  ymax <- max(yval[is.finite(yval)])
  if (asinh_transform) ymax <- asinh(ymax)
  xmax <- max(abs(xval[is.finite(xval)]))
  
  if(!is.null(fold_change_thresh) & !is.null(p_val_thresh)){
    if(norm){
      x1 <- fold_change_thresh
    }else{
      x1 <- log10(fold_change_thresh)
    }
    
    y1 <- -log10(p_val_thresh)
    x2 <- xmax
    
    if (asinh_transform) y1 <- asinh(y1)
    y2 <- ymax
  }
  
  
  
  
  if(!is.null(xlim) ){
    xrange <- xlim
  }else{
    xrange <- c(-xmax*1.05, xmax*1.05)
  }
  
  if(!is.null(ylim)){
    yrange <- ylim
  }else{
    yrange <- c(-0.1, ymax*1.05)
  }
  
  for( i in seq_along(conditions) ){
    
    if(is.null(data)){
      if(norm){
        df <- data.frame(p_val=res_int$p_val[[conditions[i]]], 
                         fold_change= res_int$norm_log_fold_change[[conditions[i]]], 
                         names=labels)
        df$X <- df$fold_change
        df$Y <- -log10(df$p_val)
        
      }else{
        df <- data.frame(p_val=res_int$p_val[[conditions[i]]], 
                         fold_change= res_int$fold_change[[conditions[i]]], 
                         names=labels)
        df$X <- log10(df$fold_change)
        df$Y <- -log10(df$p_val)
      }
    }else{
      df <- data
      df$names <- as.character(data$names)
      df$Y <- -log10(data$p_val)
      df$X <- log10(df$fold_change)
    }
    
    if (asinh_transform) df$Y <- asinh(df$Y)
    score_print <- rep(0, dim(df)[1])
    
    if(!is.null(p_val_thresh) & !is.null(fold_change_thresh)){
      is_above_thresh <- rep(0, dim(df)[1])
      is_in_frame <- rep(0, dim(df)[1])
      is_above_thresh[ which(df$p_val <= p_val_thresh & df$fold_change >= fold_change_thresh ) ] <- 1
      is_in_frame[ which(df$X >= xrange[1] & df$X <= xrange[2] & df$Y >= yrange[1] & df$Y <= yrange[2]) ] <- 1
      score_print<- is_above_thresh + is_in_frame
      score_print[is_in_frame==0]<-0
      N_show <- min(N_print, sum(score_print>0))
      if( N_show>0 ){
        idx_print <- order(score_print, df$fold_change, decreasing = TRUE)[ 1 : N_show ]
      }else{
        idx_print <- NULL
      }
    } else if(!is.null(x0) & !is.null(c)){
      is_above_thresh <- rep(0, dim(df)[1])
      is_in_frame <- rep(0, dim(df)[1])
      is_above_thresh[ which(-log10(df$p_val) >= c/(log10(df$fold_change) - x0) & 
                               log10(df$fold_change) >= x0 ) ] <- 1
      is_in_frame[ which(df$X >= xrange[1] & df$X <= xrange[2] & df$Y >= yrange[1] & df$Y <= yrange[2]) ] <- 1
      score_print<- is_above_thresh + is_in_frame
      score_print[is_in_frame==0]<-0
      N_show <- min(N_print, sum(score_print>0))
      if( N_show>0 ){
        idx_print <- order(score_print, df$fold_change, decreasing = TRUE)[ 1 : N_show ]
      }else{
        idx_print <- NULL
      }
    }else{
      if( N_print>0 ){
        idx_print <- order(df$fold_change, decreasing = TRUE)[ 1:N_print ]
      }else{
        idx_print <- NULL
      }
    }
    
    df$label_color <- as.factor(score_print)
    if(norm){
      label_x <- "norm. log(fold_change) [sd units]"
    }else{
      label_x <- "log10(fold_change)"
    }
    
    label_y <- "-log10(p_value)"
    if (asinh_transform) label_y <- "asinh(-log10(p_value))"
    

    df$label <- unlist(lapply(as.character(df$names), function(x){
      l <- nchar(x)
      if(!is.null(n_character_max)){
        if(l > n_character_max){
          return(paste(substr(x,1,n_character_max), "...", sep = ""))
        }else{return(x)}
      }else{
        return(x)
      }
    }))
    
    
    plist[[i]] <- ggplot( df , aes_string( label='label' ) ) +
      coord_cartesian(xlim = xrange, ylim = yrange, expand = FALSE) +
      xlab(label_x ) + 
      ylab(label_y) +
      ggtitle(conditions[i])
    
    
    if(!is.null(p_val_thresh) & !is.null(fold_change_thresh)){
      plist[[i]] <- plist[[i]] +
        geom_polygon(data=data.frame(x=c(x1,xrange[2],xrange[2],x1),
                                     y=c(y1,y1,yrange[2],yrange[2])), 
                     mapping=aes_string(x='x', y='y'),
                     alpha=0.1,
                     inherit.aes=FALSE)
      if(both_sides){
        plist[[i]] <- plist[[i]] +
          geom_polygon(data=data.frame(x=c(-x1,xrange[1],xrange[1],-x1),
                                       y=c(y1,y1,yrange[2],yrange[2])), 
                       mapping=aes_string(x='x', y='y'),
                       alpha=0.1,
                       inherit.aes=FALSE) 
      }
      
      plist[[i]] <- plist[[i]] + scale_x_continuous(limits = xrange)
      
      if(show_thresholds){
        plist[[i]] <- plist[[i]] +
          annotate("segment", x = xrange[1], xend = xrange[2], y = y1, yend = y1, colour = rgb(1,0,0, alpha_segment) ) +
          annotate("segment", x = -x1, xend = -x1, y = 0, yend = yrange[2], colour = rgb(1,0,0, alpha_segment) ) +
          annotate("segment", x = x1, xend = x1, y = 0, yend = yrange[2], colour = rgb(1,0,0, alpha_segment) )
      }
     
    }
    
    if(!is.null(x0) & !is.null(c)){
      xpath_right = seq(x0 + 0.01, xrange[2], by = 0.1)
      ypath_right = c/( xpath_right - x0)
      
      xpath_left = seq(xrange[1], -x0 - 0.01, by = 0.1)
      ypath_left = -c/( xpath_left + x0)
      
      if(asinh_transform){
        ypath_right <- asinh(ypath_right)
        ypath_left <- asinh(ypath_left)
      }
      
      plist[[i]] <- plist[[i]] +
        geom_polygon(data=data.frame(x=c(xpath_right, 
                                         xpath_right[length(xpath_right)] ),
                                     y=c(ypath_right,
                                         ypath_right[1])),
                     mapping=aes_string(x='x', y='y'),
                     alpha=0.1,
                     inherit.aes=FALSE) +
        geom_path(data=data.frame(x=xpath_right, y=ypath_right), mapping=aes_string(x='xpath_right', y='ypath_right'), 
                  colour = rgb(1,0,0,0.5), inherit.aes=FALSE)+
        geom_path(data=data.frame(x=xpath_left, y=ypath_left), mapping=aes_string(x='xpath_left', y='ypath_left'), 
                  colour = rgb(1,0,0,0.5), inherit.aes=FALSE)
    }
    
    
    if(!is.null(names)){
      idx_print <- which(df$names %in% names)
      df$label_color[idx_print] <- "2"
    }
    
    plist[[i]] <- plist[[i]] + 
      geom_point(data=df, mapping=aes_string(x='X', y='Y'), size = size, alpha = alpha, color = color) +
      geom_point(data=df[idx_print, ],
                 aes_string(x='X', y='Y'), colour = "red", size = size, alpha=0.8) +
      geom_text_repel(data=df[idx_print, ],
                      aes_string(x='X', y='Y', label = 'label', colour='label_color'), size = label_size, ...) +
      scale_color_manual(values = c("0" = color, "1" = color, "2" = rgb(1,0,0) ), guide=FALSE) +
      theme_function() +
      theme(legend.position="none")
      
    
  }
  
  if( length(save_file)>0 ){
    pdf( save_file, 6, 6)
    print(plist)
    dev.off()
  }
  
  return(plist)
}

#' Dot plot representation of interaction as a function of experimental conditions
#' @param res an \code{InteRactome}
#' @param names vector of names to be displayed
#' @param idx_cols numeric vector to select and order conditions to be displayed
#' @param idx_rows numeric vector to select proteins to display
#' @param size_var name of the variable corresponding to dot size
#' @param size_range range of dot sizes to display
#' @param size_limits limits used for the dot size scale
#' @param color_var name of the variable corresponding to dot color
#' @param color_breaks vector used to discretize colors
#' @param color_values values parameter passed to \code{scale_color_manual()}
#' @param color_default value corresponding to the default color
#' @param save_file path of output file (.pdf)
#' @param plot_width width of the output .pdf file
#' @param plot_height height of the output .pdf file
#' @param clustering logical or numeric vector. If logical, use hierarchical 
#' @param theme_name name of the ggplot2 theme function to use ('theme_gray' by default)
#' @param n_character_max max number of label characters (ignored if NULL)
#' @param ... additionnal arguments passed to \code{dot_plot()}
#' clustering to order proteins. If numeric, ordering indexes for displayed proteins 
#' (must be the same length as \code{idx_rows})
#' @return a list conataining :
#' @return a plot "plot"
#' @return a numeric vector "idx_order" containing the position of the protein 
#' displayed within the \code{InteRactome}
#' @importFrom grDevices dev.off pdf rgb
#' @importFrom stats dist hclust
#' @export
plot_per_condition <- function( res,
                                names = NULL,
                                idx_cols = 1:length(res$conditions),
                                idx_rows=1:20,
                                size_var="norm_stoichio",
                                size_range=c(0, 5.5),
                                size_limits=c(0, 1),
                                color_var="p_val", 
                                color_breaks=c(1,0.1,0.05,0.01),
                                color_values = c( "black", "blue", "purple", "red"),
                                color_default = 1,
                                save_file=NULL,
                                plot_width=2.5 + length(res$conditions)/5,
                                plot_height=2 + length(idx_rows)/5,
                                clustering = FALSE,
                                theme_name = "theme_gray",
                                n_character_max = 8,
                                ...){
  
  if(length(idx_rows)==1){
    idx_rows<-1:idx_rows
  }
  
  M<-do.call(cbind, res[[size_var]])
  M1<-do.call(cbind, res[[color_var]])
  
  row.names(M) <- unlist(lapply(res$names, function(x){
    l <- nchar(x)
    
    if(!is.null(n_character_max)){
      if(l > n_character_max){
        return(paste(substr(x,1,n_character_max), "...", sep = ""))
      }else{
        return(x)
      }
    }else{return(x)}
    
  }))
  
  if(!is.null(names)){
    row.names(M) <- names
  }  
  
  Mcol<-M
  Mcol[!is.null(M)]<-color_default
  
  if(!is.null(color_var)){
    names(color_values) <- as.character(color_breaks)
    idx_order_col <- order(color_breaks, decreasing = TRUE);
    for(i in seq_along(color_breaks)){
      Mcol[M1<color_breaks[idx_order_col[i]]]<-color_breaks[idx_order_col[i]]
    }
  }
  
  if("interactor" %in% names(res)){
    title_text <- paste(res$bckg_bait," vs ", 
                        res$bckg_ctrl,
                        " (n=",length(res$interactor[res$interactor != res$bait]),")",
                        sep="")
  } else {
    title_text <- paste(res$bckg_bait," vs ", res$bckg_ctrl, sep="")
  }
  
  M <- M[idx_rows, idx_cols]
  Mcol <- Mcol[idx_rows, idx_cols]
  M[is.na(M)]<-0
  
  idx_order <- 1:length(idx_rows)
  if(is.logical(clustering)){
    if(clustering){
      d<-dist(M)
      h<-hclust(d)
      idx_order <- h$order
    }
  }else if(is.numeric(clustering)){
    if(length(clustering) == length(idx_rows)){
      idx_order <- clustering
    }
  } 
  
  M <- M[idx_order, ]
  Mcol <- Mcol[idx_order, ]
  
  p<-dot_plot( as.matrix(M), 
               as.matrix(Mcol), 
               title = title_text,
               size_var = size_var, 
               size_range=size_range,
               size_limits=size_limits,
               color_var=color_var,
               color_values = color_values, 
               theme_name = theme_name,
               ...)
    
  if(!is.null(save_file)){
    pdf(save_file, plot_width, plot_height)
    print(p)
    dev.off()
  }
  
  return(list(plot = p, idx_order = idx_order))
  
}

#' Dot plot representation of matrices
#' @param Dot_Size a matrix of dot sizes
#' @param Dot_Color a matrix of dot colors (optionnal)
#' @param title plot title
#' @param size_range range of dot sizes to display
#' @param size_limits limits for dot size (as used in \code{ggplot2::scale_radius()})
#' @param size_breaks breaks used for dot size scale
#' @param size_var name of the variable corresponding to dot size
#' @param color_var name of the variable corresponding to dot color
#' @param color_values values parameter passed to \code{scale_color_manual()}
#' @param size_label_y size of y-axis labels
#' @param size_label_x size of x-axis labels
#' @param theme_name name of the ggplot2 theme function to use ('theme_gray' by default)
#' @return a plot
#' @import ggplot2
#' @export
dot_plot <- function(Dot_Size, 
                     Dot_Color=NULL, 
                     title="Dot Plot", 
                     size_range = c(0, 5.5) ,
                     size_limits = range(Dot_Size),
                     size_breaks = NULL,
                     size_var ="size", 
                     color_var="color",
                     color_values = c( "red", "purple",  "blue", "black" ),
                     size_label_y = NULL,
                     size_label_x = NULL,
                     theme_name = "theme_gray"){
  
  theme_function <- function(...){
    do.call(theme_name, list(...))
  }
  
  # Dot_Size: matrix of dot size
  
  M<-Dot_Size
  Mcol <- Dot_Color
  
  ylabels <- row.names(M)
  if(length(ylabels)==0){
    ylabels <- 1:dim(M)[1]
  }
  
  xlabels <- colnames(M)
  if(length(xlabels)==0){
    xlabels <- 1:dim(M)[2]
  }
  
  xpos <- vector("list", dim(M)[2] );
  ypos <- vector("list", dim(M)[2] );
  size <- vector("list", dim(M)[2] );
  if(length(Dot_Color)>0){
    if( sum(dim(M) != dim(M))==0 ){
      color <- vector("list", dim(M)[2] );
    }else{
      stop("Dimensions of size and color matrix do not match")
    }
  }
  
  pos <-  - ( 1:dim(M)[1] )
  
  for( k in 1:dim(M)[2] ){
    xpos[[k]] <- rep(k, dim(M)[1] )
    ypos[[k]] <- pos
    size[[k]] <- as.numeric(M[,k]);
    if(length(Dot_Color)>0){
      color[[k]] <- Mcol[,k];
    }
  }
  
  df<-data.frame( xpos=unlist(xpos), ypos=unlist(ypos), size=unlist(size) );
  if(length(Dot_Color)>0){
    df$color = unlist(color)
  }else{
    df$color = rep( 1, dim(df)[1] )
  }
  
  df$color<-as.factor(df$color)
  
  unique_col <- unique(df$color);
  if(is.null(size_label_y)){
    size_label_y <- max(6, 16 - (dim(M)[1] %/% 10)*1.5 )
  }
  if(is.null(size_label_x)){
    size_label_x <- max(6, 16 - (dim(M)[2] %/% 5)*1.5 )
  }
  
  p <- ggplot(df, aes_string(x='xpos', y='ypos', size='size', col='color' ) ) +
    theme_function()+
    theme(
      plot.title = element_text(size=12),
      axis.text.y= element_text(size=size_label_y), 
      axis.text.x = element_text(size=size_label_x, angle = 90, hjust = 1,vjust=0.5) ) +
    ggtitle(title)+
    scale_color_manual( values=color_values , name = color_var) +
    xlab("") +
    ylab("") +
    scale_x_continuous(breaks=1:dim(M)[2],
                       limits=c(0.5, dim(M)[2]+0.5),
                       labels=xlabels) +
    scale_y_continuous(breaks=pos,
                       limits= -c(dim(M)[1]+0.75, 0.25 ),
                       labels=ylabels) +
    geom_point(alpha=0.5, show.legend = TRUE)
  
  if(!is.null(size_breaks)){
    p <- p + scale_radius(limits = size_limits, range = size_range, name=size_var, breaks = size_breaks)
  }else{
    p <- p + scale_radius(limits = size_limits, range = size_range, name=size_var)
  }

  return(p)
  
}

#' Plot interaction stoichiometries per biological replicate
#' @param res an \code{InteRactome}
#' @param name name of the protein to display
#' @param conditions set of conditions to display
#' @param ref_condition name of the reference condition for all \code{test}
#' @param test name of the test function to compare sttoichiometries between conditions
#' @param test.args arguments passed to function \code{test}
#' @param map_signif_level named vector with labels and corresponding significance levels
#' @param save_file path of output file (.pdf)
#' @param theme_name name of the ggplot2 theme function to use ('theme_gray' by default)
#' @return a plot
#' @import ggplot2
#' @import ggsignif
#' @importFrom grDevices dev.off pdf rgb
#' @export
plot_stoichio <- function(res, 
                          name,
                          conditions = res$conditions,
                          ref_condition = conditions[1], 
                          test="t.test", 
                          test.args = list("paired"=TRUE),
                          map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
                          save_file = NULL,
                          theme_name = "theme_gray"){
  
  theme_function <- function(...){
    do.call(theme_name, list(...))
  }
  
  plot_title <- paste(name, " / ",test, sep = "") 
  
  if("paired" %in% names(test.args)){
    if(test.args$paired){
      plot_title <- paste(plot_title, "(paired)", sep="")
    }
  }
  
  
  idx_match <- which(res$names == name)
  df_tot <- NULL
  for ( bio in res$replicates){
    stoichio <- do.call(cbind, res$stoichio_bio[[bio]])[idx_match, conditions]
    df <- data.frame(stoichio = stoichio, cond = factor(names(stoichio), levels=conditions), bio = rep(bio, length(stoichio)))
    df_tot <- rbind(df_tot, df)
  }
  
  comparisons <- list()
  cond_test <- setdiff(conditions, ref_condition)
  for (i in 1:length(cond_test)){
    comparisons[[i]] <- c(ref_condition, cond_test[i])
  }
  
  df_tot$log10_stoichio <- log10(df_tot$stoichio)
  
  p <- ggplot(df_tot, aes_string(x='cond', y='log10_stoichio')) + 
    theme_function() +
    theme(axis.text.x = element_text(size=10, angle = 90, hjust = 1,vjust=0.5),
          axis.text.y = element_text(size=10)
    ) +
    geom_point(size=0, alpha = 0) +  
    ggtitle(plot_title) + 
    xlab("conditions") +
    geom_line( data = df_tot, mapping = aes_string(x='cond', y='log10_stoichio', group='bio', color='bio'), alpha = 0.2) +
    geom_point( data = df_tot, mapping = aes_string(x='cond', y='log10_stoichio', color='bio'), size=3, alpha = 0.8) + 
    geom_signif(comparisons = comparisons, 
                step_increase = 0.1,
                test = test,
                textsize = 2.5, 
                test.args = test.args,
                map_signif_level = map_signif_level
    )
  if(!is.null(save_file)){
    pdf(save_file, 4, 4)
    print(p)
    dev.off()
  }
  return(p)
  
}

#' Plot protein intensities per biological replicate and background
#' @param res an \code{InteRactome}
#' @param names name of the protein to display
#' @param conditions set of conditions to display
#' @param textsize size of labels corresponding to significance levels
#' @param ylims plot limits on the y axis
#' @param mapping name of the plot elemnet on which the aestethics color and alpha are mapped. Can be either "point" or "bar".
#' @param var_x x variable 
#' @param levels_x defines factor levels for x variable
#' @param labels_x defines the labels for the levels of the x variable
#' @param var_facet_x variable used for faceting plots horizontally
#' @param var_facet_y variable used for faceting plots vertically
#' @param var_color variable used for the color aestethic
#' @param var_alpha variable used for the alpha aestethic
#' @param color_values named vector of colors.
#' @param alpha_values named vector of alpha values.
#' @param offset_bar logical, use an offset to ensure all values are positive
#' @param show_bar logical, show bars using \code{geom_bar}
#' @param show_error_bar logical, show error bars using \code{geom_errorbar}
#' @param show_signif logical, show significance of comparison tests using \code{geom_signif}
#' @param show_violin logical, show point distribution using \code{geom_violin}
#' @param comparisons list of comparison pairs (as indices or x variable names)
#' @param test name of the test function to compare intensities between background
#' @param test.args arguments passed to function \code{test()}
#' @param map_signif_level named vector with labels and corresponding significance levels
#' @param position name of the function used to position data points 
#' @param position.args arguments passed to function \code{position()}
#' @param theme_name name of the ggplot2 theme function to use ('theme_gray' by default)
#' @return a plot
#' @import ggplot2
#' @import ggsignif
#' @export
plot_comparison <- function(res,
                            names,
                            conditions = res$conditions, 
                            textsize = 4,
                            ylims= NULL,
                            mapping = "point",
                            var_x = "bckg",
                            var_facet_x = "cond",
                            var_facet_y = "name",
                            var_color = "bio",
                            var_alpha = "missing",
                            levels_x = c(res$bckg_ctrl, res$bckg_bait),
                            labels_x = c(res$bckg_ctrl, res$bckg_bait),
                            color_values = NULL,
                            alpha_values = c("TRUE" = 0.5, "FALSE" = 1),
                            show_bar = FALSE,
                            show_error_bar = FALSE,
                            show_signif = TRUE,
                            show_violin = TRUE,
                            offset_bar = FALSE,
                            comparisons = list(c(1,2)),
                            test="t.test",
                            test.args = list("paired"=FALSE),
                            map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
                            position = "position_jitter",
                            position.args = list(width=0.3, height=0),
                            theme_name = "theme_gray"){
  
  theme_function <- function(...){
    do.call(theme_name, list(...))
  }
  
  plot_title <- paste(test, sep = " / ") 
  
  if("paired" %in% names(test.args)){
    if(test.args$paired){
      plot_title <- paste(plot_title, "(paired)", sep="")
    }
  }
  
  df_tot <- NULL
  for(name in names){
    idx_match <- which(rownames(res$data$Intensity) == name)
    for( cond in conditions){
      for( bio in unique(res$data$conditions$bio) ){
        
        
        idx_cond <- which( res$data$conditions$time == cond &
                             res$data$conditions$bio == bio &
                             res$data$conditions$bckg == res$bckg_bait )
        
        
        
        if(length(idx_cond)>0){
          intensity = as.numeric( res$data$Intensity_na_replaced[idx_match, idx_cond] )
          intensity_na = as.numeric( res$data$Intensity[idx_match, idx_cond] )
          df <- data.frame(intensity = intensity,
                           intensity_na = intensity_na,
                           bckg = rep(res$bckg_bait, length(intensity)), 
                           bio = rep(bio, length(intensity)), 
                           cond=  rep(cond, length(intensity)),
                           name = rep(name, length(intensity)))
          
          df_tot <- rbind(df_tot, df)
        }
        
        if(res$params$pool_background){
          idx_cond <- which(res$data$conditions$bio == bio &
                              res$data$conditions$bckg == res$bckg_ctrl )
          
        }else{
          idx_cond <- which( res$data$conditions$time == cond &
                               res$data$conditions$bio == bio &
                               res$data$conditions$bckg == res$bckg_ctrl )
        }
        
        if(length(idx_cond)>0){
          intensity_na = as.numeric( res$data$Intensity[idx_match, idx_cond] )
          intensity = as.numeric( res$data$Intensity_na_replaced[idx_match, idx_cond] )
          df <- data.frame(intensity_na = intensity_na,
                           intensity = intensity,
                           bckg = rep(res$bckg_ctrl, length(intensity)), 
                           bio = rep(bio, length(intensity)), 
                           cond=  rep(cond, length(intensity)),
                           name = rep(name, length(intensity)))
          df_tot <- rbind(df_tot, df)
        }
        
      }
    }
  }
  
  
  df_tot$x <- df_tot[[var_x]]
  if(!is.null(levels_x)){
    df_tot$x <- factor(df_tot$x, levels = levels_x, labels = labels_x) 
  }
  
  df_tot$facet_x <- df_tot[[var_facet_x]]
  df_tot$facet_y <- df_tot[[var_facet_y]]
  df_tot$color <- df_tot[[var_color]]
  df_tot$missing <- is.na(df_tot[["intensity_na"]])
  df_tot$alpha <- df_tot[[var_alpha]]
  
  label_y <- "log10(Intensity)"
  
  if(offset_bar){
    df_tot$intensity <- df_tot$intensity / min( c(df_tot$intensity, df_tot$intensity_na), na.rm = TRUE ) 
    df_tot$intensity_na <- df_tot$intensity_na / min( c(df_tot$intensity, df_tot$intensity_na), na.rm = TRUE ) 
    label_y <- "log10(Intensity Norm.)"
  }
  
  df_tot$log10_intensity <- log10(df_tot$intensity)
  
  p <- ggplot(df_tot, aes_string(x='x', y='log10_intensity'  )) + 
    theme_function() +
    theme(axis.text = element_text(size=12),
          axis.text.x = element_text(angle=90, hjust = 1)) +
    geom_point(size=0, alpha = 0) + 
    ggtitle(plot_title) +
    xlab(var_x) +
    ylab(label_y)
  
  if(show_bar){
    if(mapping == "bar")  {
      p <- p + geom_bar(data = df_tot,
                        mapping = aes_string(x='x', y='log10_intensity', alpha = 'alpha', fill = 'color',
                        stat = "summary", fun.y = "mean") )
      if(!is.null(color_values)){
        p <- p + scale_fill_manual(values = color_values)
      }
    } else {
      p <- p + geom_bar(data = df_tot, 
                        mapping = aes_string(x='x', y='log10_intensity', stat = "summary", fun.y = "mean"))
    }
  }
  
  if(show_violin){
    p <- p + geom_violin()
  }
  
  if(mapping == "point")  {
    p <- p + geom_point( data = df_tot, 
                         mapping = aes_string(x='x', y='log10_intensity', color='color', alpha = 'alpha'),
                         size=1.5,
                         position = do.call(position, position.args) )
    if(!is.null(color_values)){
      p <- p + scale_color_manual(values = color_values)
    }
    if(!is.null(var_color)){
      p <- p + guides(color=guide_legend(title=var_color))
    }
    
  }else{
    p <- p + geom_point( data = df_tot,
                         mapping = aes_string(x='x', y='log10_intensity'), 
                         size=1.5,
                         alpha = 0.8,
                         position = do.call(position, position.args) )
  }
  
  
  if(!is.null(alpha_values)){
    p <- p + scale_alpha_manual(values = alpha_values)
  }
  if(!is.null(var_alpha)){
    p <- p + guides(alpha=guide_legend(title=var_alpha))
  }   
  
  
  
  
  if(show_error_bar){
    p <- p + geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), width = 0.1)
  }
  
  if(show_signif){
    p <- p + geom_signif(comparisons = comparisons, 
                         step_increase = 0.1,
                         test = test,
                         textsize = textsize, 
                         test.args = test.args,
                         map_signif_level = map_signif_level)
  }
  
  p <- p + scale_y_continuous(limits= c(min(log10(df_tot$intensity)), 
                                        max(log10(df_tot$intensity)) + 
                                          0.2*(max(log10(df_tot$intensity)) - 
                                                 min(log10(df_tot$intensity)) )) )
  
  if(!is.null(ylims)){
    if(ylims[1] <= min(log10(df_tot$intensity)) & ylims[2] >= max(log10(df_tot$intensity)) ){
      p <- p + scale_y_continuous(limits = c(ylims[1], ylims[2]))
    }
  }
  
  p <- p + facet_grid(facet_y ~ facet_x)
  
  return(p)
  
}

#' Quality check plots for preprocessed AP-MS data
#' @param data an \code{Interactome} or preprocessed data as obtained using function \code{preprocess_data()}
#' @param theme_name name of the ggplot2 theme function to use ('theme_gray' by default)
#' @return Several QC plots
#' @import ggplot2
#' @importFrom stats quantile IQR
#' @importFrom Hmisc rcorr
#' @export
plot_QC <- function(data, theme_name = "theme_gray"){
  
  theme_function <- function(...){
    do.call(theme_name, list(...))
  }
  
  p_list <- list()
  
  df <- data
  
  if(class(data) == "InteRactome"){
    df$Intensity <- df$data$Intensity
    df$conditions <- df$data$conditions
  }
  
  ibait <- which(rownames(df$Intensity) == df$bait)
  
  M <- as.matrix(df$Intensity)
  R <- Hmisc::rcorr(M)
  
  
  idx_bait <- df$conditions$bckg == df$bckg_bait
  idx_ctrl <- df$conditions$bckg == df$bckg_ctrl
  
  df$conditions$bait <- rep("", dim(df$conditions)[1])
  df$conditions$bait[idx_bait] <- "Bait"
  df$conditions$bait[idx_ctrl] <- "Ctrl"
  
  Ravg_bait <- cbind(df$conditions[idx_bait, ], 
                     Ravg = row_mean(R$r[idx_bait, idx_bait]))
  Ravg_ctrl <- cbind(df$conditions[idx_ctrl, ], 
                     Ravg = row_mean(R$r[idx_ctrl, idx_ctrl]))
  Ravg <- rbind(Ravg_bait, Ravg_ctrl)
  
  
  x <- Ravg_bait$Ravg[Ravg_bait$bckg==df$bckg_bait]
  qnt <- stats::quantile(x, probs=c(.25, .75), na.rm = TRUE)
  H <- 1.5 * stats::IQR(x, na.rm = T)
  
  outlier_bio_1 <- unique(Ravg_bait$bio[x < (qnt[1] - H)])
  if(length(outlier_bio_1)>0){
    message_outlier_1 <- paste("outliers in", paste(outlier_bio_1, sep=", "), sep=" ")
  }else{
    message_outlier_1 <- "No outliers"
  }
  message_outlier_1 <- paste(message_outlier_1, "(in bait bckg)", sep=" ")
  message_outlier_1 = ""
  
  p1 <- ggplot(Ravg, aes_string(x='bckg', y='Ravg', col='bio', shape = "time")) + 
    theme_function()+
    ggtitle("QC: Intensity Correlation", subtitle = message_outlier_1) +
    ylab("average R")+
    geom_boxplot(data=Ravg, mapping=aes_string(x='bckg', y='Ravg'), inherit.aes = FALSE, outlier.alpha = 0) +
    geom_point( size = 3, position=position_jitter(width=0.25), alpha = 0.8)
  
  p_list[[1]] <- p1
  
  Ibait <- cbind(df$conditions, Ibait = as.numeric(df$Intensity[ibait, ]) )
  Ibait <- Ibait[Ibait$bckg %in% c(df$bckg_bait, df$bckg_ctrl), ]
  x <- Ibait$Ibait[Ibait$bckg==df$bckg_bait]
  qnt <- stats::quantile(x, probs=c(.25, .75), na.rm = TRUE)
  H <- 1.5 * stats::IQR(x, na.rm = TRUE)
  
  outlier_bio_2 <- unique(Ibait$bio[x < (qnt[1] - H) | x > (qnt[2] + H)])
  if(length(outlier_bio_2)>0){
    message_outlier_2 <- paste("outliers in", paste(outlier_bio_2, sep=", "), sep=" ")
  }else{
    message_outlier_2 <- "No outliers"
  }
  message_outlier_2 <- paste(message_outlier_2, "(in bait bckg)", sep=" ")
  message_outlier_2 = ""
  Ibait$log10_Ibait <- log10(Ibait$Ibait)
  
  p2 <- ggplot(Ibait, aes_string(x='bckg', y='log10_Ibait', col='bio', shape = "time")) + 
    theme_function() +
    ggtitle("QC: Bait Purification", subtitle = message_outlier_2) +
    ylab("norm. Intensity (log10)") +
    geom_boxplot(data=Ibait, mapping=aes_string(x='bckg', y='log10_Ibait'), inherit.aes = FALSE, outlier.alpha = 0) +
    geom_point( size = 3, position=position_jitter(width=0.25), alpha = 0.8)
  
  p_list[[2]] <- p2
  
  nNA <- cbind(df$conditions,
               nNA = sapply(1:dim(df$Intensity)[2], FUN=function(x){sum(is.na(df$Intensity[,x]))})
  )
  nNA <- nNA[nNA$bckg %in% c(df$bckg_bait, df$bckg_ctrl), ]
  x <- nNA$nNA[nNA$bckg==df$bckg_bait]
  qnt <- stats::quantile(x, probs=c(.25, .75), na.rm = TRUE)
  H <- 1.5 * stats::IQR(x, na.rm = TRUE)
  
  outlier_bio_3 <- unique(nNA$bio[x < (qnt[1] - H) | x > (qnt[2] + H)])
  if(length(outlier_bio_3)>0){
    message_outlier_3 <- paste("outliers in", paste(outlier_bio_3, sep=", "), sep=" ")
  }else{
    message_outlier_3 <- "No outliers"
  }
  message_outlier_3 <- paste(message_outlier_3, "(in bait bckg)", sep=" ")
  message_outlier_3 = ""
  nNA$log10_nNA <- log10(nNA$nNA)
  
  p3 <- ggplot(nNA, aes_string(x='bckg', y='log10_nNA', col='bio', shape = "time")) +
    theme_function() +
    ggtitle("QC: Missing Values", subtitle = message_outlier_3) +
    ylab("NA counts (log10)") +
    geom_boxplot(data=nNA, mapping=aes_string(x='bckg', y='log10_nNA'), inherit.aes = FALSE, outlier.alpha = 0) +
    geom_point( size = 3, position=position_jitter(width=0.25), alpha = 0.8)
  
  p_list[[3]] <- p3
  
  names(p_list) <- c("Intensity_Correlation", "Bait_Purification", "Missing_values")
  
  return(list(plot = p_list, Ravg = Ravg, Ibait = Ibait, nNA = nNA))
  
}

#' Plot an interactive correlation network with communities highlighted
#' @param res an \code{InteRactome}
#' @param idx indexes of the set of proteins in \code{res} for which correlations will be computed.
#' @param df_corr a data.frame with columns 'r_corr' and 'p_corr'. Has priority over parameters \code{res}.
#' @param source variable of \code{df_corr} with the names of the source protein
#' @param target variable of \code{df_corr} with the names of the target protein
#' @param cluster named vector containing the cluster number for each node
#' @param var_p_val variable for which the threshold \code{p_val_thresh} will be applied. Set to 'p_corr' by default
#' @param var_r_corr variable for which the threshold \code{r_corr_thresh} will be applied. Set to 'r_corr' by default
#' @param r_corr_thresh threshold for variable 'r_corr' (min)
#' @param p_val_thresh threshold for variable \code{var_p_val} (max)
#' @param ... other parameters passed either to function \code{compute_correlations()} if \code{df_corr} is NULL
#' or to function \code{restrict_network_degree()} otherwise.
#' @return an interactive networkD3 plot
#' @import igraph
#' @import networkD3
#' @export
plot_correlation_network <- function(res, 
                                     idx = NULL, 
                                     source = "name_1", 
                                     target = "name_2",
                                     cluster = NULL,
                                     df_corr = NULL, 
                                     var_p_val = "p_corr", 
                                     var_r_corr = "r_corr", 
                                     r_corr_thresh = 0.8, 
                                     p_val_thresh = 0.05,
                                     ...){
  
  if(is.null(df_corr)){
    if(is.null(idx)){
      if("interactor" %in% names(res)){
        idx_filter <- which(res$is_interactor > 0)
      }else{
        idx_filter <- 1:length(res$names)
      }
    }else{
      idx_filter = idx
    }
    df_corr <- compute_correlations(res = res, idx = idx_filter, ...)
    
    if(var_r_corr %in% names(df_corr) & var_p_val %in% names(df_corr)){
      idx_filter_corr <- which(df_corr[[var_r_corr]]>=r_corr_thresh & df_corr[[var_p_val]]<=p_val_thresh)
      df_corr_filtered <- df_corr[idx_filter_corr, ]
    }else{
      df_corr_filtered <- df_corr[ , ]
    }
    
  }else{
    if(is.null(idx)){
      if(var_r_corr %in% names(df_corr) & var_p_val %in% names(df_corr)){
        idx_filter_corr <- which(df_corr[[var_r_corr]]>=r_corr_thresh & df_corr[[var_p_val]]<=p_val_thresh)
        df_corr_filtered <- df_corr[idx_filter_corr, ]
      }else{
        df_corr_filtered <- df_corr
      }
      
    }else{
      idx_filter_corr <- idx
      df_corr_filtered <- df_corr[idx_filter_corr, ]
    }
    df_corr_filtered <- restrict_network_degree(df_corr_filtered, ...)
  }
  
  
  
  net <- igraph::graph.data.frame(df_corr_filtered[ , c(source, target)], directed=FALSE)
  net <- igraph::simplify(net)
  
  
  if(!is.null(cluster)){
    idx_match <- match(vertex.attributes(net)$name, names(cluster))
    group = cluster[idx_match]
  }else{
    cfg <- igraph::cluster_fast_greedy(as.undirected(net))
    group = cfg$membership
    names(group) <- cfg$names
  }
  net_d3 <- networkD3::igraph_to_networkD3(net, group = group)
  
  p <- forceNetwork(Links = net_d3$links, Nodes = net_d3$nodes,
               Source = 'source', Target = 'target',
               fontFamily = "arial",
               NodeID = 'name', Group = 'group',
               colourScale = JS("d3.scaleOrdinal(d3.schemeCategory20);"),
               charge = -10, opacity = 1,
               linkColour = rgb(0.75, 0.75, 0.75),
               fontSize = 12, bounded = TRUE, zoom=TRUE, opacityNoHover = 1
  )
  
  return( list(plot = p, 
               net_d3 = net_d3, 
               df_corr_filtered = df_corr_filtered, 
               cluster = group))
  
}

#' Plot points with density background with correlation coefficient
#' @param df a data.frame
#' @param var_x name of the x variable
#' @param var_y name of the y variable
#' @param theme_name name of the ggplot2 theme function to use ('theme_gray' by default)
#' @return a plot
#' @import ggplot2
#' @import Hmisc
#' @export
plot_density <- function(df, var_x = names(df)[1], var_y = names(df)[2], theme_name = "theme_gray"){
  
  theme_function <- function(...){
    do.call(theme_name, list(...))
  }
  
  df$x <- df[[var_x]]
  df$y <- df[[var_y]]
  Pcorr <- rcorr(x=df[[var_x]], y=df[[var_y]]  )
  
  p <- ggplot(df, aes_string(y='y', x='x' ) ) +
    theme_function() + 
    theme(axis.text=element_text(size=16), plot.title = element_text(size=16))+
    stat_density2d(aes_string(alpha='..level..', fill='..level..'), size=2, geom="polygon", bins=20) + 
    scale_fill_gradient(low = "yellow", high = "red") +
    scale_alpha(range = c(0.00, 0.5), guide = FALSE) +
    geom_point(alpha=0.2, size=1.5)+
    xlab(var_x) +
    ylab(var_y) +
    ggtitle(paste('Pearson R=',signif(Pcorr$r[1,2],5),", n=",dim(df)[1],sep="") )
  
  return(p)
    
}



#' Plot FDR as a function of parameters used to divide the volcano plot
#' @param FDR_res output from the function \code{compute_FDR_from_asymmetry()}
#' @param FDR_bins FDR levels
#' @param xlim x-axis plot limits
#' @param ylim y-axis plot limits
#' @param colors color palette. Should have \code{length(FDR_bins)-1} colors.
#' @importFrom stats reshape
#' @importFrom graphics filled.contour title
#' @importFrom grDevices terrain.colors
#' @export
plot_FDR_map <- function(FDR_res, 
                         FDR_bins = c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2),
                         xlim=c(0,3), 
                         ylim=c(0,3),
                         colors = terrain.colors(length(FDR_bins)-1)
){
  
  T1 <- FDR_res$parameters
  T2 <- FDR_res$max_TP_parameters
  
  T3 <- T1[ , c("x0", "c", "FDR")]
  T4 <- stats::reshape(T3, idvar = "x0", timevar = "c", direction = "wide")
  T5 <- as.matrix(T4[,2:dim(T4)[2]]);
  
  
  filled.contour(x = unique(T3$x0),
                 y = unique(T3$c),
                 z = T5,
                 levels = FDR_bins,
                 xlim = xlim, 
                 ylim = ylim, 
                 col = colors ,
                 plot.title = title(xlab = "x0", 
                                    ylab = "c")
                 
  )
}
