# Load packages ----

library(shiny)
library(ggplot2)
library(ggrepel)
library(grid)
library(Hmisc)
library(igraph)
library(networkD3)
library(ggsignif)
library("InteRact")

#source("./R/InteRact.R")

options(shiny.maxRequestSize = 100*1024^2) #maximum file size is set to 100MB


# User interface ----
ui <- fluidPage(
  titlePanel("InteRact : Analysis of AP-MS data"),
  
  fluidRow(
    column(3,
           br(),
           wellPanel(
             h3("Parameters"),
             textInput("bait_gene_name", "Bait (gene name)", value = "Bait"),
             checkboxInput("pool_background", "pool_background", value = TRUE),
             checkboxInput("substract_ctrl", "substract_ctrl", value = TRUE),
             numericInput("Nrep", "# iterations (missing values replacement)", value = 3),
             numericInput("p_val_thresh", "p-value (maximum)", value = 0.01),
             numericInput("fold_change_thresh", "fold-change (minimum)", value = 2),
             numericInput("n_success_min", "n_success_min", value = 1),
             checkboxInput("consecutive_success", "consecutive_success", value = TRUE),
             verbatimTextOutput("interactors"),
             downloadButton("download_all", "Save analysis")
             #verbatimTextOutput("plotly_print")
           )
    ),
    column(9,
           br(),
           tabsetPanel(id = "inTabset",
                       tabPanel("Import",
                                column(4,
                                       br(),
                                       wellPanel(
                                         h3("Import"),
                                         fileInput("file", h4("Protein intensity file :"), placeholder = "Enter file here"),
                                         checkboxInput("dec", "Use comma as decimal separator", value = FALSE),
                                         textInput("pattern", "Pattern for intensity columns", value = "^Intensity."),
                                         selectInput("column_gene_name",
                                                     "column for gene name",
                                                     choices = list(),
                                                     selected = NULL),
                                         selectInput("column_ID",
                                                     "column for protein ID",
                                                     choices = list(),
                                                     selected = NULL)
                                       )
                                ),
                                column(8,
                                       br(),
                                       dataTableOutput("data_summary")
                                )
                       ),
                       tabPanel("Group",

                                column(4,
                                       br(),
                                       wellPanel(
                                         h3("General"),
                                         textInput("bckg_bait", "Name of Bait background", value = "Bait"),
                                         textInput("bckg_ctrl", "Name of Control background", value = "WT")
                                         
                                       ),
                                       wellPanel(
                                         h3("Map samples"),
                                         checkboxInput("manual_mapping", "manual mapping", value = FALSE),
                                         conditionalPanel(
                                           condition = "input.manual_mapping == false",
                                           textInput("split", "split character", value = "_"),
                                           h4("Enter position of:"),
                                           numericInput("bckg_pos", "background", value = 1),
                                           numericInput("bio_pos", "biological replicates", value = 2),
                                           numericInput("time_pos", "experimental conditions", value = 3),
                                           numericInput("tech_pos", "technical replicates", value = 4)
                                           # textInput("preffix_bio", "biological replicates", value = "S"),
                                           # textInput("preffix_tech", "technical replicates", value = "R"),
                                           # textInput("preffix_time", "experimental conditions", value = "")
                                           
                                         ),
                                         conditionalPanel(
                                           condition = "input.manual_mapping == true",
                                           fileInput("file_cond", h4("Import file :"), placeholder = "Enter file here"),
                                           checkboxInput("dec_cond", "Use comma as decimal separator", value = FALSE),
                                           checkboxInput("sample_by_rows", "Samples by rows", value = FALSE),
                                           selectInput("manual_bckg",
                                                              "background",
                                                              choices = list(),
                                                              selected = NULL),
                                           selectInput("manual_bio",
                                                              "biological replicates",
                                                              choices = list(),
                                                              selected = NULL),
                                           selectInput("manual_tech",
                                                              "technical replicates",
                                                              choices = list(),
                                                              selected = NULL),
                                           selectInput("manual_time",
                                                              "experimental conditions",
                                                              choices = list(),
                                                              selected = NULL)
                                           
                                         )
                                         
                                       )
                                ),
                                column(8,
                                       br(),
                                       dataTableOutput("condTable")
                                )
                       ),
                       tabPanel("QC / Filter",
                                column(4,
                                    br(),
                                    wellPanel(
                                        helpText("check boxes to filter out samples"),
                                        checkboxGroupInput("filter_bio",
                                                           "Biological replicates (bio)",
                                                           choices = list(),
                                                           selected = NULL),
                                        checkboxGroupInput("filter_tech",
                                                           "Technical replicates (tech)",
                                                           choices = list(),
                                                           selected = NULL),
                                        checkboxGroupInput("filter_time",
                                                           "Experimental conditions (time)",
                                                           choices = list(),
                                                           selected = NULL)

                                    )
                                ),
                                column(8,
                                  br(),
                                  plotOutput("QCPlot1", width="250",height="200"),
                                  plotOutput("QCPlot2", width="250",height="200"),
                                  plotOutput("QCPlot3", width="250",height="200")
                                  # plotlyOutput("QCPlot1ly", width="300",height="250"),
                                  # plotlyOutput("QCPlot2ly", width="300",height="250"),
                                  # plotlyOutput("QCPlot3ly", width="300",height="250")
                                  #dataTableOutput("condTable_bis")
                                )
                       ),
                       tabPanel("Volcano",
                                column(4,
                                       br(),
                                       wellPanel(
                                         selectInput("volcano_cond", "Select condition",
                                                     choices = list(), selected = NULL),
                                         numericInput("N_print", "# labels displayed (maximum) ", value = 15),
                                         checkboxInput("asinh_transform", "asinh_transform", value = TRUE)

                                       ),
                                       wellPanel(
                                         helpText("Hover mouse over point to display extra info"),
                                         helpText("Click to select protein"),
                                         helpText("Brush and double-click to zoom")
                                       ),
                                       
                                       br(),
                                       plotOutput("compPlot_volcano",width="200",height="200")
                                       
                                ),
                                column(8,
                                       br(),
                                       fluidRow(
                                        downloadButton("download_volcano", "Download plot"),
                                        downloadButton("download_all_volcanos", "Download all volcanos")
                                       ),
                                       br(),
                                       plotOutput("volcano", width="400",height="400",
                                                  hover = hoverOpts(id ="volcano_hover"),
                                                  click = "volcano_click",
                                                  dblclick = "volcano_dblclick",
                                                  brush = brushOpts(
                                                  id = "volcano_brush",
                                                  resetOnNew = TRUE) ),
                                       verbatimTextOutput("info_volcano_hover")
                                )
                       ),
                       tabPanel("Dot Plot",
                                column(4,
                                       br(),
                                       wellPanel(
                                         numericInput("Nmax", "N display ", value = 30),
                                         checkboxInput("clustering", "Hierarchical clustering", value = FALSE)
                                       ),
                                       wellPanel(
                                         helpText("Hover mouse over point to display extra info"),
                                         helpText("Click to select protein"),
                                         helpText("Brush and double-click to zoom")
                                       ),
                                       br(),
                                       plotOutput("stoichioPlot",width="200",height="200"),
                                       br(),
                                       plotOutput("compPlot",width="200",height="200")
                                ),
                                column(8,
                                       br(),
                                       downloadButton("download_dotPlot", "Download Plot", value = FALSE),
                                       br(),
                                       plotOutput("dotPlot",width="250",height="500",
                                                  hover = hoverOpts(id ="dotPlot_hover"),
                                                  dblclick = "dotPlot_dblclick",
                                                  click = "dotPlot_click",
                                                  brush = brushOpts(
                                                    id = "dotPlot_brush",
                                                    resetOnNew = TRUE) ),
                                       br(),
                                       verbatimTextOutput("info_dotPlot_hover")
                                )
                       ),
                       tabPanel("2D Stoichio",
                                br(),
                                fluidRow(
                                  column(width=4,
                                    wellPanel(
                                      selectInput("Stoichio2D_cond", "Select condition",
                                                  choices = list(), selected = NULL)
                                    )
                                  ),
                                  column(width=4,
                                    wellPanel(
                                      numericInput("Nmax2D", "# displayed (max) ", value = 30)
                                    )
                                  )
                                ),
                                fluidRow(
                                  column(width=5,
                                         helpText("Brush to select zoom area"),
                                         downloadButton("download_Stoichio2D", "Download Plot", value = FALSE),
                                         plotOutput("Stoichio2D", width="300",height="300",
                                                    hover = hoverOpts(id ="Stoichio2D_hover"),
                                                    brush = brushOpts(
                                                      id = "Stoichio2D_brush",
                                                      resetOnNew = TRUE) ),
                                         verbatimTextOutput("info_Stoichio2D_hover")

                                  ),
                                  column(width=5,
                                         helpText("zoom on selected area"),
                                         downloadButton("download_Stoichio2D_zoom", "Download Plot", value = FALSE),
                                         plotOutput("Stoichio2D_zoom", width="300",height="300",
                                                    hover = hoverOpts(id ="Stoichio2D_zoom_hover")),
                                         br(),
                                         verbatimTextOutput("info_Stoichio2D_zoom_hover")
                                 )
                                )


                       ),
                       tabPanel("Correlations",
                                column(4,
                                       br(),
                                       wellPanel(
                                         numericInput("r_corr_thresh", "Correlation Pearson R (min) ", value = 0.8),
                                         numericInput("p_val_corr_thresh", "Associated p-value (max)", value = 0.05)
                                       ),
                                       wellPanel(
                                         helpText("zoom in : dbl click"),
                                         helpText("zoom out : shift + dbl click "),
                                         helpText("move : click + drag ")
                                       )
                                ),
                                column(8,
                                       br(),
                                       forceNetworkOutput("force_net")
                                       
 
                                )
                       ),
                       tabPanel("Annotations",
                                br(),
                                column(4,
                                       wellPanel(
                                         checkboxGroupInput("annotation_selected",
                                                            "Select annoations",
                                                            choices = c( 
                                                                        "Protein.families",
                                                                        "Keywords",
                                                                        "GO",
                                                                        "KEGG",
                                                                        "Reactome", 
                                                                        "Pfam",
                                                                        "Hallmark",
                                                                        "GO_molecular_function",
                                                                        "GO_biological_process",
                                                                        "GO_cellular_component"),
                                                            selected = c("Keywords", "Protein.families")),
                                         actionButton("launch_annot","Launch analysis")
                                       ),
                                       wellPanel(
                                         #checkboxInput("slim","use GO slim", value=FALSE),
                                         selectInput("method_adjust_p_val", "Method to adjust p-values",
                                                     choices = c("none", "fdr", "bonferroni"), selected = "fdr"),
                                         numericInput("p_val_max", "p-value (maximum)", value = 0.05),
                                         numericInput("fold_change_min", "fold-change (minimum)", value = 2),
                                         numericInput("N_annot_min", "Number of annotated proteins (minimum)", value = 2)
                                       )
                                ),
                                column(8,
                                       br(),
                                       downloadButton("download_annotPlot", "Download Plot", value = FALSE),
                                       plotOutput("annotPlot"),
                                       br(),
                                       downloadButton("download_annotTable", "Download Table", value = FALSE),
                                       br(),
                                       br(),
                                       dataTableOutput("annotTable")
                                )


                       ),
                       tabPanel("Summary",
                                column(4,
                                       br(),
                                       wellPanel(
                                         checkboxGroupInput("columns_displayed",
                                                            "Columns displayed",
                                                            choices = c("names", "max_stoichio", "max_fold_change", "min_p_val"),
                                                            selected = c("names", "max_stoichio", "max_fold_change", "min_p_val")
                                                            )
                                       )
                                ),
                                column(8,
                                       br(),
                                       downloadButton("download_summaryTable", "Download summary table", value = FALSE),
                                       br(),
                                       br(),
                                       dataTableOutput("summaryTable")
                                )

                       )
           )
    )
  )
              
)

# Server logic
server <- function(input, output, session) {
  
  #Reactive values ---------------------------------------------------------------------------------
  
  Ninteractors <- reactiveValues(x=0)
  ranges <- reactiveValues(x = c(-1.5,0.5), y = c(-1,1))
  ranges_volcano <- reactiveValues(x = NULL, y = NULL)
  ranges_dotPlot <- reactiveValues(x = NULL, y = NULL)

  saved_df <- reactiveValues(annot = NULL, enrichment = NULL)
 
  annotation<- reactiveValues(selected = NULL, 
                              loaded = NULL, 
                              to_load = NULL, 
                              enrichment_performed = NULL, 
                              enrichment_to_perform = NULL)
  
  var_to_load <- reactiveValues(names = NULL)
  df_corr_plot <- reactiveValues(names = NULL, x = NULL, y = NULL, cluster=NULL)
  select_dotPlot <- reactiveValues(i_prot = 1, i_cond=1)
  select_volcanoPlot <- reactiveValues(i_min = 1, min_dist1 = 0)
  select_stoichioPlot <- reactiveValues(i_min = 1, min_dist1 = 0)
  idx_order <- reactiveValues(cluster = NULL)
  
  #Main reactive functions -------------------------------------------------------------------------
  
  data <- reactive({
    
    df <- read.csv(input$file$datapath, 
             sep="\t", fill=TRUE, 
             na.strings="", 
             dec=ifelse(input$dec,",",".") )
    
    gn_selected <- names(df)[grep("GENE", toupper(names(df)))[1]]
    id_selected <- names(df)[grep("ID", toupper(names(df)))[1]]
    
    updateSelectInput(session, "column_gene_name",
                      choices = as.list(names(df)),
                      selected = gn_selected) 
    
    updateSelectInput(session, "column_ID",
                      choices = as.list(names(df)),
                      selected = id_selected) 
    
    df
    
  })
  
  data_cond <- reactive({
    if(input$manual_mapping){
      df_cond <- read.table(input$file_cond$datapath, 
                 sep="\t", 
                 fill=TRUE, 
                 na.strings="",
                 header=FALSE)
      
      if (input$sample_by_rows) {
        df_cond <- transpose(df_cond)
      }
      
      updateSelectInput(session, "manual_bckg",
                               choices = as.list(df_cond[ , 1]),
                               selected = NULL) 
      
      updateSelectInput(session, "manual_bio",
                               choices = as.list(df_cond[ , 1]),
                               selected = NULL) 
      
      updateSelectInput(session, "manual_tech",
                               choices = as.list(df_cond[ , 1]),
                               selected = NULL) 
      
      updateSelectInput(session, "manual_time",
                               choices = as.list(df_cond[ , 1]),
                               selected = NULL)
      df_cond
    } else {
      cond()
    }
    
    
  })
  
  cond <- reactive({
    
    if(input$manual_mapping){
      
      df_cond <- data_cond()
      
      idx_cond <- grep(input$pattern, names(data()))
      col_I <- names(data())[idx_cond]
      
      idx_match <- match(col_I, df_cond[1, ])
      
      bckg <- rep("", length(col_I))
      bio <- rep("", length(col_I))
      tech <- rep("", length(col_I))
      time <- rep("", length(col_I))
      
      idx_bckg <- which(df_cond[, 1] == input$manual_bckg)
      cat(idx_bckg)
      if(length(idx_bckg)>0){
        bckg <- factor(unlist(df_cond[idx_bckg, idx_match]))
      }
      idx_bio <- which(df_cond[, 1] == input$manual_bio)
      if(length(idx_bio)>0){
        bio <- factor(unlist(df_cond[idx_bio, idx_match]))
      }
      idx_tech <- which(df_cond[,1] == input$manual_tech)
      if(length(idx_tech)>0){
        tech <- factor(unlist(df_cond[idx_tech, idx_match]))
      }
      idx_time <- which(df_cond[,1] == input$manual_time)
      if(length(idx_time)>0){
        time <- factor(unlist(df_cond[idx_time, idx_match]))
      } 
      
      cond_int <- dplyr::tibble(idx=seq_along(col_I), column=col_I, bckg, time, bio, tech)
                      
    } else {
      # cond_int <- identify_conditions(data(),
      #                               Column_intensity_pattern = input$pattern,
      #                               bckg_bait = input$bckg_bait,
      #                               bckg_ctrl = input$bckg_ctrl,
      #                               preffix_time = input$preffix_time,
      #                               preffix_bio = input$preffix_bio, 
      #                               preffix_tech = input$preffix_tech,
      #                               split = input$split)
      
      cond_int <- identify_conditions(data(),
                                      Column_intensity_pattern = input$pattern,
                                      bckg_pos = input$bckg_pos,
                                      bio_pos = input$bio_pos,
                                      time_pos = input$time_pos, 
                                      tech_pos = input$tech_pos,
                                      split = input$split)
    }
    
    
    
    
    updateCheckboxGroupInput(session, "filter_bio",
                             choices = as.list(unique(cond_int$bio)),
                             selected = NULL) 
    
    updateCheckboxGroupInput(session, "filter_tech",
                             choices = as.list(unique(cond_int$tech)),
                             selected = NULL) 
    
    updateCheckboxGroupInput(session, "filter_time",
                             choices = as.list(unique(cond_int$time)),
                             selected = NULL)
    cond_int
  })
  
  prep_data <- reactive({
    
    updateSelectInput(session, "volcano_cond",
                      choices = as.list(setdiff(unique(cond()$time), input$filter_time) ),
                      selected = NULL)
    updateSelectInput(session, "Stoichio2D_cond",
                      choices = as.list(setdiff( c("max", unique(cond()$time)), input$filter_time)),
                      selected = NULL)

    Column_gene_name <- input$column_gene_name #names(data())[grep("GENE", toupper(names(data())))[1]]
    Column_ID <- input$column_ID #names(data())[grep("ID", toupper(names(data())))[1]]
    
    preprocess_data(  df = data(),
                      Column_intensity_pattern = input$pattern,
                      bait_gene_name = input$bait_gene_name,
                      Column_gene_name = Column_gene_name,
                      Column_ID = Column_ID,
                      bckg_bait = input$bckg_bait,
                      bckg_ctrl = input$bckg_ctrl,
                      # preffix_bio = input$preffix_bio,
                      # preffix_tech = input$preffix_tech,
                      # preffix_time = input$preffix_time,
                      filter_bio = input$filter_bio,
                      filter_tech = input$filter_tech,
                      filter_time = input$filter_time,
                      condition = cond()
                      )
  })
  
  res<- reactive({
    
    # Create a Progress object
    progress <- shiny::Progress$new(min = 0, max = 100)
    progress$set(message = "Compute interactome...", value = 0)
    on.exit(progress$close())
    updateProgress <- function(value = NULL, detail = NULL) {
      progress$set(value = value, detail = detail)
    }

    res_int <- InteRact(preprocess_df = prep_data(),
                      N_rep=input$Nrep,
                      pool_background = input$pool_background,
                      updateProgress = updateProgress,
                      substract_ctrl = input$substract_ctrl)
           
    res_int$Interactome <- merge_proteome(res_int$Interactome)
    df_merge <- merge_conditions(res_int$Interactome)
    df_FDR <- compute_FDR_from_asymmetry(df_merge)
    res_int$Interactome <- append_FDR(res_int$Interactome, df_FDR)
    res_int
    
  })

  order_list <- reactive({
    res_int <- identify_interactors (res()$Interactome,
                                     var_p_val = "p_val", 
                                     p_val_thresh = input$p_val_thresh, 
                                     fold_change_thresh = input$fold_change_thresh, 
                                     n_success_min = input$n_success_min, 
                                     consecutive_success = input$consecutive_success)
    order_list_int <- get_order_discrete(res_int, var_p_val = "min_p_val", p_val_breaks=c(1,0.1,0.05,0.01) )
    Ninteractors$x <- order_list_int$Ndetect
    order_list_int
  })

  ordered_Interactome <- reactive({
      res_int <- identify_interactors (res()$Interactome,
                                     var_p_val = "p_val", 
                                     p_val_thresh = input$p_val_thresh, 
                                     fold_change_thresh = input$fold_change_thresh, 
                                     n_success_min = input$n_success_min, 
                                     consecutive_success = input$consecutive_success)
      res_int <- order_interactome(res_int, order_list()$idx_order)
      res_int
  })
  
  annotated_Interactome <- reactive({
      results<-append_annotations(ordered_Interactome(), saved_df$annot )
      
      names_excluded <- c("names","bait", "groups","conditions")
      updateCheckboxGroupInput(session, "columns_displayed",
                             choices = as.list( setdiff(names(results), names_excluded)),
                             selected=c("max_stoichio", "max_fold_change", "min_p_val") )
      results
  })

  #Observe functions -------------------------------------------------------------------
  
  # observe({
  #   
  #     b_name <- input$bait_gene_name
  #     updateTextInput(session, "bckg_bait", value =  b_name)
  #     
  # })
  
  observeEvent(input$launch_annot, {
    
    annotation$selected <- input$annotation_selected
    annotation$loaded <- names(saved_df$annot)
    annotation$to_load <- setdiff(annotation$selected, annotation$loaded)
    
    if(length(annotation$to_load )>0){
      
      # Create a Progress object
      progress <- shiny::Progress$new(min = 0, max = 100)
      on.exit(progress$close())
      updateProgress <- function(value = NULL, detail = NULL) {
        progress$set(value = value, detail = detail)
      }
      
      if( !("Entry" %in% annotation$loaded) ){
        Sys.sleep(1)
        progress$set(message = "Append annotations...", value = 0)
        saved_df$annot <- get_annotations(data(), updateProgress = updateProgress)
      }
      if( "KEGG" %in% annotation$to_load ){
        progress$set(message = "Add KEGG annotations...", value = 0)
        saved_df$annot <- add_KEGG_data(saved_df$annot, updateProgress = updateProgress)
      }
      if( "Hallmark" %in% annotation$to_load ){
        progress$set(message = "Add Hallmark annotations...", value = 0)
        saved_df$annot <- add_Hallmark_data(saved_df$annot, updateProgress = updateProgress)
      }
      if( "GO_molecular_function" %in% annotation$to_load ){
        progress$set(message = "Add GO molecular function annotations...", value = 0)
        saved_df$annot<- add_GO_data(saved_df$annot, GO_type = "molecular_function", slim = FALSE, updateProgress = updateProgress)
      }
      if( "GO_biological_process" %in% annotation$to_load ){
        progress$set(message = "Add GO biological_process annotations...", value = 0)
        saved_df$annot <- add_GO_data(saved_df$annot, GO_type = "biological_process", slim = FALSE, updateProgress = updateProgress)
      }
      if( "GO_cellular_component" %in% annotation$to_load ){
        progress$set(message = "Add GO cellular_component annotations...", value = 0)
        saved_df$annot <- add_GO_data(saved_df$annot, GO_type = "cellular_component", slim = FALSE, updateProgress = updateProgress)
      }
    }
    
    annotation$loaded <- names(saved_df$annot)
    
  })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$Stoichio2D_brush, {
    brush <- input$Stoichio2D_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)

    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })

  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$volcano_dblclick, {
    brush_volcano <- input$volcano_brush
    if (!is.null(brush_volcano)) {
      ranges_volcano$x <- c(brush_volcano$xmin, brush_volcano$xmax)
      ranges_volcano$y <- c(brush_volcano$ymin, brush_volcano$ymax)

    } else {
      ranges_volcano$x <- NULL
      ranges_volcano$y <- NULL
    }
  })

  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$dotPlot_dblclick, {
    brush_dotPlot <- input$dotPlot_brush
    if (!is.null(brush_dotPlot)) {
      ranges_dotPlot$x <- c(brush_dotPlot$xmin, brush_dotPlot$xmax)
      ranges_dotPlot$y <- c(brush_dotPlot$ymin, brush_dotPlot$ymax)

    } else {
      ranges_dotPlot$x <- NULL
      ranges_dotPlot$y <- NULL
    }
  })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$plot_corr_brush, {

    brush <- input$plot_corr_brush
    click <- input$plot_corr_click
    df_brush <- data.frame( x=c(brush$xmin, brush$xmax), y=c(brush$ymin, brush$ymax))
    dist_brush_x <- (click$x - df_brush$x)^2
    idx_brush_x <- which.max(dist_brush_x)
    dist_brush_y <- (click$y - df_brush$y)^2
    idx_brush_y <- which.max(dist_brush_y)

    dist<- (click$x - df_corr_plot$x)^2 + (click$y - df_corr_plot$y)^2
    idx_selected <- which.min(dist)

    df_corr_plot$x[idx_selected] <- df_brush$x[idx_brush_x]
    df_corr_plot$y[idx_selected] <- df_brush$y[idx_brush_y]

  })
  
  #Reactive functions for output ---------------------------------------------------------------
  
  df_corr <- reactive({
    compute_correlations(ordered_Interactome(), 
                         idx=which(ordered_Interactome()$is_interactor > 0))
    
  })
  
  df_corr_filtered <- reactive({
    
    df1 <- df_corr()
    df1 <- df1[df1$r_corr>=input$r_corr_thresh & df1$p_corr<=input$p_val_corr_thresh, ]
    
    net <- graph.data.frame(df1, directed=FALSE);
    net <- igraph::simplify(net)
    layout <- layout_nicely(net)
    cfg <- cluster_fast_greedy(as.undirected(net))
    net_d3 <- igraph_to_networkD3(net, group = cfg$membership)
    
    # vatt <- vertex.attributes(net)
    # vertex_names <- as.character(vatt$name)
    # 
    # #layout <- data.frame(x=rnorm(length(vertex_names)), y=rnorm(length(vertex_names)))
    # #idx_vertex <- as.numeric(vertex_names)
    # #df_corr_plot$names <- names[idx_vertex]
    # 
    # df_corr_plot$names <- vertex_names
    # df_corr_plot$x <- layout[ , 1]
    # df_corr_plot$y <- layout[ , 2]
    # df_corr_plot$cluster <- as.factor(cfg$membership)
    # #df_corr_plot$cluster <- sample(10, length(vertex_names), replace = TRUE)
    # 
    # df1
  })
  
  output$force_net <- renderForceNetwork({
    
    forceNetwork(Links = df_corr_filtered()$links, Nodes = df_corr_filtered()$nodes,
                 Source = 'source', Target = 'target',
                 fontFamily = "arial",
                 NodeID = 'name', Group = 'group',
                 colourScale = JS("d3.scaleOrdinal(d3.schemeCategory20);"),
                 charge = -10, opacity = 1,
                 linkColour = rgb(0.75, 0.75, 0.75),
                 fontSize = 12, bounded = TRUE, zoom=TRUE, opacityNoHover = 1
    )
  })
  
  # corrPlot <- reactive({
  #   
  #   df2 = data.frame(x=df_corr_plot$x, y=df_corr_plot$y, label=df_corr_plot$names, cluster=df_corr_plot$cluster)
  #   
  #   #print(head(df_corr_filtered() ))
  #   
  #   p<-ggplot(df2, aes(x, y, label=label, color=cluster)) +
  #     theme_void() +
  #     geom_point(alpha=0.3, size=10) +
  #     geom_text()
  #   
  #   for (i in 1:length(df_corr_filtered()$name_1)){
  #     idx1 <- which(df_corr_plot$names == df_corr_filtered()$name_1[i])
  #     idx2 <- which(df_corr_plot$names == df_corr_filtered()$name_2[i])
  #     p <- p + annotate("segment",
  #                       x = df_corr_plot$x[idx1],
  #                       xend = df_corr_plot$x[idx2],
  #                       y = df_corr_plot$y[idx1],
  #                       yend = df_corr_plot$y[idx2],
  #                       colour = "gray50",
  #                       alpha=0.25)
  #   }
  #   
  #   p
  # }) 
  
  data_summary <- reactive({
    df <- data()
    columns <- names(df)
    data_class <- sapply(1:dim(df)[2], FUN=function(x){class(df[,x])})
    data_median <- sapply(1:dim(df)[2], FUN=function(x){ 
                                              if(is.numeric(df[,x])){
                                                median(df[,x], na.rm=TRUE)
                                              } else {
                                                NA
                                              }
                                            })
    
    data.frame(columns = columns, class = data_class, median = data_median)
  })
  
  condTable <- reactive({
    data_cond()
  })
    
  summaryTable <- reactive({
    summary_table(annotated_Interactome(),  add_columns = input$columns_displayed)
  })
  
  
  observeEvent(input$launch_annot,{
    
    if(!is.null(annotation$loaded)){
      
      if("annot_type" %in% names(saved_df$enrichment)){
        annotation$enrichment_performed <- unique(saved_df$enrichment$annot_type)
      }
      
      annotation$enrichment_to_perform <- setdiff(annotation$selected, annotation$enrichment_performed)
      
      if( length(annotation$enrichment_to_perform) > 0 ){
        
        # Create a Progress object
        progress2 <- shiny::Progress$new(min = 0, max = 100)
        on.exit(progress2$close())
        updateProgress2 <- function(value = NULL, detail = NULL) {
          progress2$set(value = value, detail = detail)
        }
        progress2$set(message = "Perform enrichment analysis...", value = 0)
        
        df_annot <- annotation_enrichment_analysis( annotated_Interactome(), 
                                                    1:order_list()$Ndetect, 
                                                    annotation_selected = annotation$enrichment_to_perform, 
                                                    names = annotated_Interactome()$names,
                                                    updateProgress = updateProgress2)
        saved_df$enrichment <- rbind(saved_df$enrichment, df_annot)
        annotation$enrichment_performed <- unique(saved_df$enrichment$annot_type)
      }
      
    }
    
  })
  
  annotTable <- reactive({
    
    if(!is.null(annotation$enrichment_performed) & length(setdiff(annotation$selected, annotation$enrichment_performed))==0 ){
      
      idx_annot <- which(saved_df$enrichment$annot_type %in% annotation$selected)
      df<-saved_df$enrichment[idx_annot, ]
      
      idx_annot_exist <-  which(df$N_annot>0);
      p_value_adjust_fdr <- rep( 1,length(df$p_value) );
      p_value_adjust_bonferroni <- rep( 1,length(df$p_value) );
      p_value_adjust_bonferroni[idx_annot_exist] <- p.adjust(df$p_value[idx_annot_exist], method = "bonferroni");
      p_value_adjust_fdr[idx_annot_exist] <- p.adjust(df$p_value[idx_annot_exist], method = "fdr");
      
      df$p_value_adjust_fdr <- p_value_adjust_fdr
      df$p_value_adjust_bonferroni <- p_value_adjust_bonferroni
      df<- df[ order(df$p_value, decreasing = FALSE), ]
      
      output = df
    }else{
      output = NULL
    }
    
  })
  
  Stoichio2D_zoom <- reactive({
    plot_2D_stoichio(ordered_Interactome(),
                     condition = input$Stoichio2D_cond,
                     xlim = ranges$x,
                     ylim = ranges$y,
                     N_display=min(order_list()$Ndetect, input$Nmax2D) )
  })
  
  Stoichio2D <- reactive({
    plot_2D_stoichio(ordered_Interactome(),
                     condition = input$Stoichio2D_cond,
                     N_display = min(order_list()$Ndetect, input$Nmax2D) )
  })
  
  dotPlot <- reactive({
    p <- plot_per_conditions(ordered_Interactome(),
                        idx_rows = min(input$Nmax, order_list()$Ndetect),
                        clustering = input$clustering)
    idx_order$cluster <- p$idx_order
    p$plot + coord_cartesian(xlim = ranges_dotPlot$x, ylim = ranges_dotPlot$y, expand = FALSE)
  })
  
  volcano <- reactive({
      plot_volcanos( ordered_Interactome(),
                     conditions = input$volcano_cond,
                     p_val_thresh = input$p_val_thresh,
                     fold_change_thresh = input$fold_change_thresh,
                     xlim = ranges_volcano$x,
                     ylim = ranges_volcano$y,
                     N_print=input$N_print,
                     asinh_transform = input$asinh_transform)[[1]]
  })

  all_volcanos <- reactive({
    plot_volcanos( ordered_Interactome(),
                   p_val_thresh = input$p_val_thresh,
                   fold_change_thresh = input$fold_change_thresh,
                   N_print=input$N_print )
  })
  
  annotPlot <- reactive({
    plot_annotation_results(annotTable(),
                            method_adjust_p_val = input$method_adjust_p_val,
                            p_val_max = input$p_val_max,
                            fold_change_min = input$fold_change_min,
                            N_annot_min = input$N_annot_min)
  })
  
  observe({
    
    if(!is.null(input$dotPlot_hover)){
      select_dotPlot$i_prot <- round(-input$dotPlot_hover$y)
      select_dotPlot$i_cond <- round(input$dotPlot_hover$x)
    }
    
  })
  
  stoichioPlot <- reactive({
    
    plot_stoichio(ordered_Interactome(), 
                  name = ordered_Interactome()$names[idx_order$cluster[select_dotPlot$i_prot]],
                  test = "t.test",
                  test.args = list("paired" = FALSE))
   
  })
  
  compPlot <- reactive({
    
    plot_comparison(ordered_Interactome(), 
                    name = ordered_Interactome()$names[ idx_order$cluster[select_dotPlot$i_prot]],
                    condition = ordered_Interactome()$conditions[select_dotPlot$i_cond])
  })
  
  observe({
    if(!is.null(input$volcano_hover)){
      hover=input$volcano_hover
      
      dist1=sqrt((hover$x-log10(ordered_Interactome()$fold_change[[input$volcano_cond]]) )^2 +
                   (hover$y+log10(ordered_Interactome()$p_val[[input$volcano_cond]]) )^2)
      
      if(input$asinh_transform) {
        dist1=sqrt((hover$x-log10(ordered_Interactome()$fold_change[[input$volcano_cond]]) )^2 +
                     (hover$y+asinh(log10(ordered_Interactome()$p_val[[input$volcano_cond]])) )^2)
      }
      
      select_volcanoPlot$min_dist1 <- min(dist1, na.rm=TRUE)
      select_volcanoPlot$i_min <- which.min(dist1)
    }
  })
  
  observe({
    
    if(!is.null(input$Stoichio2D_hover)){
      hover=input$Stoichio2D_hover
      
      if(input$Stoichio2D_cond == "max"){
        dist1=sqrt( (hover$x-log10(ordered_Interactome()$max_stoichio[1:min(order_list()$Ndetect, input$Nmax2D)]) )^2 +
                      (hover$y-log10(ordered_Interactome()$stoch_abundance[1:min(order_list()$Ndetect, input$Nmax2D)]) )^2)
      }else{
        dist1=sqrt( (hover$x-log10(ordered_Interactome()$stoichio[[input$Stoichio2D_cond]][1:min(order_list()$Ndetect, input$Nmax2D)]) )^2 +
                      (hover$y-log10(ordered_Interactome()$stoch_abundance[1:min(order_list()$Ndetect, input$Nmax2D)]) )^2)
      }
      
      select_stoichioPlot$min_dist1 <- min(dist1, na.rm=TRUE)
      select_stoichioPlot$i_min <- which.min(dist1)
    }
    
  })
  
  compPlot_volcano <- reactive({
    
    plot_comparison(ordered_Interactome(), 
                    name = ordered_Interactome()$names[select_volcanoPlot$i_min],
                    condition = input$volcano_cond)
  })
  
  QCPlot <- reactive({
      plot_QC(prep_data()) 
  })
    
  #Output Table functions -------------------------------------------------------------------------
  
  output$condTable <- renderDataTable(condTable())
  output$condTable_bis <- renderDataTable(condTable())
  output$data_summary <- renderDataTable(data_summary())
  output$summaryTable <- renderDataTable({summaryTable()[1:order_list()$Ndetect, ] })
  output$annotTable <- renderDataTable(annotTable())
  
  #Output Plot functions -------------------------------------------------------------------------
  
  output$Stoichio2D <- renderPlot( Stoichio2D() )
  output$Stoichio2D_zoom <- renderPlot( Stoichio2D_zoom() )
  output$dotPlot <- renderPlot( dotPlot() )
  output$volcano <- renderPlot( volcano() )
  output$annotPlot <- renderPlot(annotPlot())
  #output$plot_corr <- renderPlot(corrPlot())
  output$stoichioPlot <- renderPlot(stoichioPlot())
  output$compPlot <- renderPlot(compPlot())
  output$compPlot_volcano <- renderPlot(compPlot_volcano())
  output$QCPlot1 <- renderPlot(QCPlot()[[1]])
  output$QCPlot2 <- renderPlot(QCPlot()[[2]])
  output$QCPlot3 <- renderPlot(QCPlot()[[3]])
  #output$QCPlot1ly <- renderPlotly( ggplotly(QCPlot()[[1]], source="QCPlot1ly") )
  #output$QCPlot2ly <- renderPlotly( ggplotly(QCPlot()[[2]]) )
  #output$QCPlot3ly <- renderPlotly( ggplotly(QCPlot()[[3]]) )
  
  
  # output$plotly_print <- renderPrint({
  #   d <- event_data("plotly_hover", source = "QCPlot1ly")
  #   if (is.null(d)) "Hover on a point!" else names(d)
  # })
  
  #Output Download functions ---------------------------------------------------------------------
  
  output$download_summaryTable <- downloadHandler(
    filename = "summary_table.txt",
    content = function(file) {
      write.table(summaryTable(), file, sep = "\t", dec = ".", row.names = FALSE)
    }
  )
  
  output$download_Stoichio2D <- downloadHandler(
    filename = "Stoichio2D_plot.pdf",
    content = function(file) {
      pdf(file, 5, 5)
      print(Stoichio2D())
      dev.off()
    }
  )
  
  output$download_Stoichio2D_zoom <- downloadHandler(
    filename = "Stoichio2D_zoom_plot.pdf",
    content = function(file) {
      pdf(file, 5, 5)
      print(Stoichio2D_zoom())
      dev.off()
    }
  )
  
  output$download_dotPlot <- downloadHandler(
    filename = "dot_plot.pdf",
    content = function(file) {
      plot_width = 2.5 + length(ordered_Interactome()$conditions)/5
      plot_height = 1.5 + input$Nmax/5
      pdf(file,plot_width,plot_height)
      print(dotPlot())
      dev.off()
    }
  )
  
  output$download_volcano <- downloadHandler(
    filename = "volcano_plot.pdf",
    content = function(file) {
      pdf(file, 5, 5)
      print(volcano())
      dev.off()
    }
  )

  output$download_all_volcanos <- downloadHandler(
    filename = "all_volcanos_plot.pdf",
    content = function(file) {
      pdf(file, 5, 5)
      print(all_volcanos())
      dev.off()
    }
  )
  
  output$download_annotPlot <- downloadHandler(
    filename = "annotations.pdf",
    content = function(file) {
      pdf(file, 5, 5)
      print(annotPlot())
      dev.off()
    }
  )
  
  output$download_annotTable <- downloadHandler(
    filename = "annotations.txt",
    content = function(file) {
      write.table(annotTable(), file,  sep = "\t", dec = ".", row.names = FALSE)
    }
  )
  
  output$download_all <- downloadHandler(
    filename = paste("report.tar", sep=""),
    content = function(file) {
      
      dir.create(paste("./report/",sep=""))
      setwd("./report/")
      on.exit(setwd(".."))
      
      pdf(paste("./volcano.pdf", sep=""), 5, 5)
      print(all_volcanos())
      dev.off()
      
      write.table(annotTable(), paste("./summary_table.txt", sep=""),  sep = "\t", dec = ".", row.names = FALSE)
      
      tar(file)
      
    }
  )
  
  #Output Info functions -------------------------------------------------------------------------
  
  output$interactors<- renderPrint({
    s1 <- paste("# interactors : ", Ninteractors$x, sep="")
    cat(s1)
  })
  
  output$info_volcano_hover <- renderPrint({
    if(!is.null(input$volcano_hover)){
      hover=input$volcano_hover
      
      dist1=sqrt((hover$x-log10(ordered_Interactome()$fold_change[[input$volcano_cond]]) )^2 +
                  (hover$y+log10(ordered_Interactome()$p_val[[input$volcano_cond]]) )^2)
      
      if(input$asinh_transform) {
        dist1=sqrt((hover$x-log10(ordered_Interactome()$fold_change[[input$volcano_cond]]) )^2 +
                     (hover$y+asinh(log10(ordered_Interactome()$p_val[[input$volcano_cond]])) )^2)
      }
      
      min_dist1 <- min(dist1, na.rm=TRUE)
      i_min <- which.min(dist1)
      if( min_dist1 < 0.25){
        s1<-paste("name: ", ordered_Interactome()$names[ i_min ],sep="")
        s2<-paste("p_val: ", ordered_Interactome()$p_val[[input$volcano_cond]][ i_min ],sep="")
        s3<-paste("fold_change: ", ordered_Interactome()$fold_change[[input$volcano_cond]][ i_min ],sep="")
        s4<-paste("stoichio: ", ordered_Interactome()$stoichio[[input$volcano_cond]][ i_min ],sep="")
        cat(s1,s2,s3,s4,sep="\n")
      }
    }
  })

  output$info_Stoichio2D_zoom_hover <- renderPrint({
    if(!is.null(input$Stoichio2D_zoom_hover)){
      hover=input$Stoichio2D_zoom_hover
      if(input$Stoichio2D_cond == "max"){
        dist1=sqrt( (hover$x-log10(ordered_Interactome()$max_stoichio[1:min(order_list()$Ndetect, input$Nmax2D)]) )^2 +
                      (hover$y-log10(ordered_Interactome()$stoch_abundance[1:min(order_list()$Ndetect, input$Nmax2D)]) )^2)
      }else{
        dist1=sqrt( (hover$x-log10(ordered_Interactome()$stoichio[[input$Stoichio2D_cond]][1:min(order_list()$Ndetect, input$Nmax2D)]) )^2 +
                      (hover$y-log10(ordered_Interactome()$stoch_abundance[1:min(order_list()$Ndetect, input$Nmax2D)]) )^2)
      }
      min_dist1 <- min(dist1, na.rm=TRUE)
      i_min <- which.min(dist1)

      if(min_dist1 < 0.25){
        s1<-paste("name: ", ordered_Interactome()$names[ i_min ],sep="")
        s2<-paste("min_p_val: ", ordered_Interactome()$min_p_val[ i_min ],sep="")
        s3<-paste("max_fold_change: ", ordered_Interactome()$max_fold_change[ i_min ],sep="")
        cat(s1,s2,s3,sep="\n")
      }
    }
  })

  output$info_Stoichio2D_hover <- renderPrint({
    
        
      if(select_stoichioPlot$min_dist1 < 0.25){
        s1<-paste("name: ", ordered_Interactome()$names[ select_stoichioPlot$i_min ],sep="")
        s2<-paste("min_p_val: ", ordered_Interactome()$min_p_val[ select_stoichioPlot$i_min ],sep="")
        s3<-paste("max_fold_change: ", ordered_Interactome()$max_fold_change[ select_stoichioPlot$i_min ],sep="")
        cat(s1,s2,s3,sep="\n")
      }
    
  })

  output$info_dotPlot_hover <- renderPrint({
    #if(!is.null(input$dotPlot_hover)){
    #  i_prot = round(-input$dotPlot_hover$y)
    #  i_cond = round(input$dotPlot_hover$x)
      prot_idx <- idx_order$cluster[select_dotPlot$i_prot]
      cond_idx <- select_dotPlot$i_cond
      
      s1 <- paste("Name: ", ordered_Interactome()$names[ prot_idx ], sep="")
      s2 <- paste("Condition: ", ordered_Interactome()$conditions[ cond_idx ], sep="")
      s3 <- paste("p-value: ", ordered_Interactome()$p_val[[ cond_idx ]][ prot_idx ], sep="")
      s4 <- paste("fold-change: ", ordered_Interactome()$fold_change[[ cond_idx ]][ prot_idx ], sep="")
      s5 <- paste("stoichio: ", ordered_Interactome()$stoichio[[ cond_idx ]][ prot_idx ], sep="")
      s6 <- paste("norm_stoichio: ", ordered_Interactome()$norm_stoichio[[ cond_idx ]][ prot_idx ], sep="")
      cat(s1, s2, s3, s4, s5, s6,sep="\n")
    #}
  })

  
}

# Run the app
shinyApp(ui, server)
