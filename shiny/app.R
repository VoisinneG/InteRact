#Load packages ----

library(shiny)
library(shinyBS)

library("InteRact")

# source("../R/InteRact.R")
# library(ggrepel)
# library(ggsignif)
# library(grid)
# library(mice)
# library(Hmisc)
# library(igraph)
# library(networkD3)
# library(dplyr)

library(Hmisc)
library(ggplot2)
library(networkD3)
library(data.table)
library(BiocInstaller)

`%then%` <- shiny:::`%OR%`

options(repos = BiocInstaller::biocinstallRepos())
#getOption("repos")
options(shiny.maxRequestSize = 100*1024^2) #maximum file size is set to 100MB

method_choices <- c("default",
                    "pmm", 
                    "midastouch",
                    "sample",
                    "cart",
                    "rf",
                    "mean",
                    "norm",
                    "norm.nob",
                    "norm.boot",
                    "norm.predict",
                    "quadratic",
                    "ri",
                    "logreg",
                    "logreg.boot",
                    "polr",
                    "polyreg",
                    "lda",
                    "2l.norm",
                    "2l.lmer",
                    "2l.pan",
                    "2lonly.mean",
                    "2lonly.norm")

# User interface ----
ui <- fluidPage(
  titlePanel("InteRact : Analysis of AP-MS data"),
  
  fluidRow(
    column(3,
           br(),
           
           #conditionalPanel(
             #condition = "input.mode == 'compute'",
             wellPanel(
               h4("General parameters"),
               textInput("bait_gene_name", "Bait (gene name)", value = "Bait"),
               bsTooltip("bait_gene_name", 
                         "The gene name of the bait protein"),
               checkboxInput("pool_background", "pool ctrl intensities", value = FALSE),
               bsTooltip("pool_background", 
                         "Perform protein enrichment tests using control intensities from all conditions"),
               checkboxInput("merge_conditions", "merge conditions", value = FALSE),
               bsTooltip("merge_conditions", 
                         "Perform protein enrichment tests using control intensities from all conditions against bait intensities from all conditions"),
               checkboxInput("substract_ctrl", "substract ctrl (stoichio)", value = FALSE),
               bsTooltip("substract_ctrl", 
                         "Substract protein intensity from ctrl background to compute interaction stoichiometry"),
               checkboxInput("use_mean_for_bait", "Use mean bait intensity (stoichio)", value = FALSE),
               bsTooltip("use_mean_for_bait", 
                         "Use mean bait intensity across all conditions and biological repliactes to compute interaction stoichiometry"),
               h4("Missing values :"),
               numericInput("N_rep", "# replacements", value = 1),
               bsTooltip("N_rep", 
                         "Number of times missing values will be replaced. Use 0 if you do not want to replace missing values"),
               selectInput("method", "Replacement method", choices = as.list(method_choices), selected = "default"),
               bsTooltip("method", 
                         "Method used for the replacement of missing values.", 
                         placement = "right"),
               actionButton("start", label = "Compute Interactome")
             ),

          #),
          wellPanel(
             h4("Interaction parameters"),
             selectInput("var_p_val", "p-value variable", choices=c("p_val","FDR"), selected = "p_val"),
             bsTooltip("var_p_val", "p_val : p-value corresponding to the enrichment t-test; FDR : false discovery rate calculated using the asymmetry of the volcano plot.", 
                       placement = "top"),
             numericInput("p_val_thresh", "p-value (maximum)", value = 0.01),
             bsTooltip("p_val_thresh", "Threshold on interaction p-value"),
             numericInput("fold_change_thresh", "fold-change (minimum)", value = 2),
             bsTooltip("fold_change_thresh", "Threshold on interaction fold-change"),
             numericInput("n_success_min", "n_success_min", value = 1),
             bsTooltip("n_success_min", 
                       "Minimum number of conditions for which both interaction p-value and fold-change must pass the defined thresholds"),
             checkboxInput("consecutive_success", "consecutive_success", value = FALSE),
             bsTooltip("consecutive_success", 
                       "Should the successful passing of thresholds happen for consecutive conditions?"),
             verbatimTextOutput("interactors"),
             bsTooltip("interactors", 
                       "Number of proteins that pass the detection criteria defined above")
          )
    ),
    column(9,
           br(),
           tabsetPanel(id = "inTabset",
                       tabPanel("Import",
                                column(4,
                                       br(),
                                       wellPanel(
                                         selectInput("mode", "Select data type", choices=c("Raw data","Interactome"), selected = "raw data"),
                                         bsTooltip("mode", "Raw data : file with protein intensity values (.txt or .csv), Interactome : a pre-computed interactome", placement = "top")
                                       ),
                                       conditionalPanel(
                                         condition = "input.mode == 'Raw data'",
                                         wellPanel(
                                           fileInput("file", h4("Import file :"), placeholder = "Enter file here"),
                                           checkboxInput("delim", "Use comma as delimiter", value = FALSE),
                                           checkboxInput("dec", "Use comma as decimal separator", value = FALSE)
                                         ),
                                         wellPanel(
                                           h4("Select columns"),
                                           selectInput("column_gene_name",
                                                       "column for gene name",
                                                       choices = list(),
                                                       selected = NULL),
                                           bsTooltip("column_gene_name", 
                                                     "Choose the column containing gene names. This is where the Bait (gene name) defined in the general parameters panel should be.", 
                                                     placement = "top"
                                           ),
                                           selectInput("column_ID",
                                                       "column for protein ID",
                                                       choices = list(),
                                                       selected = NULL),
                                           bsTooltip("column_ID", 
                                                     "Choose the column containing protein IDs (from uniprot). This information is used to retrieve additional information such as GO annotations.", 
                                                     placement = "top")
                                         )
                                      ),
                                      conditionalPanel(
                                        condition = "input.mode == 'Interactome'",
                                        wellPanel(
                                          fileInput("load", h4("Load interactome :"), placeholder = "Enter file here")
                                        )  
                                      )   
                                  ),
                                  column(8,
                                         br(),
                                         conditionalPanel(
                                           condition = "input.mode == 'Raw data'",
                                           dataTableOutput("data_summary")
                                         )
                                  )  
                       ),
                       tabPanel("Group",
                                
                                  column(4,
                                         br(),
                                         wellPanel(
                                           h4("Select background"),
                                           selectInput("bckg_bait", "Bait background", choices = list(), selected = NULL ),
                                           bsTooltip("bckg_bait", 
                                                     "Enter the name of the bait background (as displayed in the table on the right).", 
                                                     placement = "top"),
                                           selectInput("bckg_ctrl", "Control background", choices = list(), selected = NULL ),
                                           #textInput("bckg_ctrl", "Control background", value = "WT"),
                                           bsTooltip("bckg_ctrl", 
                                                     "Enter the name of the control background (as displayed in the table on the right).", 
                                                     placement = "top")
                                           
                                         ),
                                         conditionalPanel(
                                           condition = "input.mode == 'Raw data'",
                                           wellPanel(
                                             h4("Map samples"),
                                             checkboxInput("manual_mapping", "manual mapping", value = FALSE),
                                             bsTooltip("manual_mapping", "Option to import custom definition of samples"),
                                             conditionalPanel(
                                               condition = "input.manual_mapping == false",
                                               textInput("pattern", "Pattern for intensity columns", value = "^Intensity."),
                                               bsTooltip("pattern", "Columns whose name contains this pattern will be identified as protein intensity columns. Regular expressions are supported."),
                                               textInput("split", "split character", value = "_"),
                                               bsTooltip("split", "split character used to divide column names in multiple substrings"),
                                               h4("Enter position of:"),
                                               numericInput("bckg_pos", "background", value = 1),
                                               bsTooltip("bckg_pos", "Position, within column names, of the substrings containing the background name (id)"),
                                               numericInput("bio_pos", "biological replicates", value = 2),
                                               bsTooltip("bio_pos", "Position, within column names, of the substrings containing the name (id) of the biological replicate"),
                                               numericInput("time_pos", "experimental conditions", value = 3),
                                               bsTooltip("time_pos", "Position, within column names, of the substrings containing the name (id) of the experimental condition"),
                                               numericInput("tech_pos", "technical replicates", value = 4),
                                               bsTooltip("tech_pos", "Position, within column names, of the substrings containing the name (id) of the technical replicate")
                                               
                                             ),
                                             conditionalPanel(
                                               condition = "input.manual_mapping == true",
                                               fileInput("file_cond", h4("Import file :"), placeholder = "Enter file here"),
                                               checkboxInput("sep_cond", "Use comma as separator", value = FALSE),
                                               checkboxInput("transpose", "transpose", value = FALSE),
                                               bsTooltip("transpose", "Invert rows/columns. Rows should correspond to protein intensity column names (corresponding to those shown in the import tab)"),
                                               h4("Choose column for:"),
                                               selectInput("column_name",
                                                           "column name",
                                                           choices = list(),
                                                           selected = NULL),
                                               bsTooltip("column_name", "Choose column containing protein intensity column names (corresponding to those shown in the import tab)", 
                                                         placement = "right"),
                                               selectInput("manual_bckg",
                                                           "background",
                                                           choices = list(),
                                                           selected = NULL),
                                               bsTooltip("manual_bckg", "Choose column containing the background name (id) for each sample", 
                                                         placement = "right"),
                                               selectInput("manual_bio",
                                                           "biological replicates",
                                                           choices = list(),
                                                           selected = NULL),
                                               bsTooltip("manual_bio", "Choose column containing the name (id) of the biological replicate for each sample", 
                                                         placement = "right"),
                                               selectInput("manual_tech",
                                                           "technical replicates",
                                                           choices = list(),
                                                           selected = NULL),
                                               bsTooltip("manual_tech", "Choose column containing the name (id) of the technical replicate for each sample", 
                                                         placement = "right"),
                                               selectInput("manual_time",
                                                           "experimental conditions",
                                                           choices = list(),
                                                           selected = NULL),
                                               bsTooltip("manual_time", "Choose column containing the name (id) of the experimental condition for each sample", 
                                                         placement = "right")
                                             )
                                           ),
                                           conditionalPanel(
                                             condition = "input.manual_mapping == false",
                                             wellPanel(
                                               h4("Format column names"),
                                               selectInput("format_function",
                                                           "Function",
                                                           choices = c("gsub", "sub"),
                                                           selected = "gsub"),
                                               bsTooltip("format_function", "Name of the function used to modify column names. The function gsub replaces all occurences of pattern by replacement while the function sub replaces only the first occurence.", 
                                                         placement = "top"),
                                               textInput("format_pattern", "pattern", value = "."),
                                               bsTooltip("format_pattern", "pattern to be replaced in column names"),
                                               textInput("format_replacement", "replacement", value = "_"),
                                               bsTooltip("format_replacement", "character string replacing pattern"),
                                               actionButton("format_names","Format column names")
                                             )
                                             
                                           )
                                         )
                                         
                                  ),
                                  column(8,
                                         br(),
                                         uiOutput("warning"),
                                         dataTableOutput("condTable")
                                  )

                                
                                # conditionalPanel(
                                #   condition = "input.mode == 'Interactome'",
                                #   br(),
                                #   helpText("Not available for imported interactomes")
                                # )
                       ),
                       tabPanel("QC / Select",
                                #conditionalPanel(
                                #  condition = "input.mode == 'Raw data'",
                                  column(4,
                                         br(),
                                         wellPanel(
                                           h4("Select samples"),
                                           selectizeInput("bio_selected", "Biological replicates", choices = list(), multiple = TRUE),
                                           selectizeInput("tech_selected", "Technical replicates", choices = list(), multiple = TRUE),
                                           selectizeInput("time_selected", "Experimental conditions (ordered)", choices = list(), multiple = TRUE),
                                           actionButton("apply_filter", label = "Apply")
                                         ),
                                         wellPanel(
                                           selectInput("QC_selected", "Select QC plot",
                                                       choices = c("Intensity correlation", "Bait purification", "Missing values"), 
                                                       selected = "Bait purification") 
                                         )
                                  ),
                                  column(8,
                                         br(),
                                         plotOutput("QCPlot", width="400",height="350")
                                  )
                                
                                #)
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
                                         helpText("Brush and double-click to zoom")
                                       ),
                                       verbatimTextOutput("info_volcano_hover")
                                       
                                       
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
                                       br(),
                                       plotOutput("compPlot_volcano",width="450",height="225")
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
                                         helpText("Brush and double-click to zoom"),
                                         helpText("Hover mouse over point to display extra info or select protein below"),
                                         selectizeInput("name_focus", "Select protein", choices = list(), multiple = FALSE)
                                       ),
                                       verbatimTextOutput("info_dotPlot_hover"),
                                       plotOutput("stoichioPlot",width="200",height="200")
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
                                       plotOutput("compPlot",width="450",height="225")

                                       
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

                       ),
                       tabPanel("Save",
                                column(4,
                                       br(),
                                       wellPanel(
                                         checkboxGroupInput("saved_items",
                                                            "Items to save",
                                                            choices = c("Interactome", "preprocessed data", "summary table", "correlation network", "volcano plot", "dot plot", "2D stoichio plot", "enrichment plots", "stoichio plots" ),
                                                            selected = c("Interactome", "preprocessed data", "summary table", "correlation network", "volcano plot", "dot plot", "2D stoichio plot" )
                                         )
                                       ),
                                       wellPanel(
                                         downloadButton("download_all", "Save report"),
                                         bsTooltip("download_all", "Download a report of the analysis")
                                       )
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
  saved_df <- reactiveValues(res = NULL, 
                             cond = NULL,
                             cond_data = NULL,
                             cond_select = NULL,
                             data = NULL, 
                             annot = NULL, 
                             enrichment = NULL, 
                             summary = NULL, 
                             params = NULL)
 
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
  check_var <- reactiveValues(is_numeric = NULL, ncol = NULL)
  
  saved_df$params <- list(N_rep = 3, 
                          method = "default", 
                          #quantile_rep = 0.05, 
                          #by_conditions = TRUE,
                          pool_background = FALSE,
                          #log_test = TRUE,
                          #log_stoichio = TRUE,
                          substract_ctrl = FALSE,
                          use_mean_for_bait = FALSE,
                          var_p_val = "p_val",
                          p_val_thresh = 0.01,
                          fold_change_thresh = 5,
                          #conditions = c(),
                          n_success_min = 1,
                          consecutive_success = FALSE)

  #Main reactive functions -------------------------------------------------------------------------
  
  observe({
    updateNumericInput(session, "N_rep", value = saved_df$params$N_rep)
    updateTextInput(session, "method", value = saved_df$params$method)
    #updateCheckboxInput(session, "by_conditions", value = saved_df$params$by_conditions)
    updateCheckboxInput(session, "pool_background", value = saved_df$params$pool_background)
    updateCheckboxInput(session, "substract_ctrl", value = saved_df$params$substract_ctrl)
    updateCheckboxInput(session, "use_mean_for_bait", value = saved_df$params$use_mean_for_bait)
    updateTextInput(session, "var_p_val", value = saved_df$params$var_p_val)
    updateNumericInput(session, "p_val_thresh", value = saved_df$params$p_val_thresh)
    updateNumericInput(session, "fold_change_thresh", value = saved_df$params$fold_change_thresh)
    updateNumericInput(session, "n_success_min", value = saved_df$params$n_success_min)
    updateCheckboxInput(session, "consecutive_success", value = saved_df$params$consecutive_success)
  })
  
  observeEvent(input$load, {
    Interactome_name <- load(input$load$datapath)
    saved_df$res <- get(Interactome_name)
    saved_df$params <- saved_df$res$params
    saved_df$cond <- saved_df$res$data$conditions
    saved_df$cond_select <- saved_df$cond
    saved_df$cond_data <- saved_df$cond
    
    updateTextInput(session, "bait_gene_name", value = saved_df$res$bait)
    
  })
  
  data_raw <- reactive({
    
    validate(
      need(input$file$datapath, "Please select a file to import")
    )
    
    df <- read.csv(input$file$datapath, 
             sep=ifelse(input$delim,";","\t"), 
             fill=TRUE, 
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
    
    saved_df$data <- df
    
    df
    
  })
  
  data <- reactive({
    saved_df$data
  })
  
  observe({
    
    if(input$manual_mapping){
      
      validate(
        need(input$file_cond$datapath, "Please select a file to import")
      )
      
      df_cond <- read.table(input$file_cond$datapath, 
                 sep=ifelse(input$sep_cond, ",", "\t"), 
                 fill=TRUE, 
                 colClasses = "character",
                 na.strings="",
                 header=TRUE)
      
      if (input$transpose) {
       
        df_cond_int <- data.table::transpose(df_cond)
        names(df_cond_int) <- df_cond[ , 1]
        df_cond_int <- cbind(names(df_cond)[-1], df_cond_int[-1, ])
        names(df_cond_int)[1] <- names(df_cond)[1]
        df_cond <- df_cond_int
        
      }
      
      updateSelectInput(session, "column_name",
                        choices = as.list(names(df_cond)),
                        selected = NULL) 
      
      updateSelectInput(session, "manual_bckg",
                               choices = as.list(names(df_cond)),
                               selected = NULL) 
      
      updateSelectInput(session, "manual_bio",
                               choices = as.list(names(df_cond)),
                               selected = NULL) 
      
      updateSelectInput(session, "manual_tech",
                               choices = as.list(names(df_cond)),
                               selected = NULL) 
      
      updateSelectInput(session, "manual_time",
                               choices = as.list(names(df_cond)),
                               selected = NULL)
      
      saved_df$cond_data <- df_cond
      
    } 

    
    
  })
  
  observe({
    
    if(input$manual_mapping){
      
      
      validate(
        need(dim(saved_df$cond_data)[1]>0, "No data imported for conditions")
      )
      
      df_cond <- saved_df$cond_data
      
      col_I <- df_cond[[input$column_name]]
      bckg <- df_cond[[input$manual_bckg]]
      time <- df_cond[[input$manual_time]]
      bio <- df_cond[[input$manual_bio]]
      tech <- df_cond[[input$manual_tech]]
      
      cond_int <- data.frame(idx=seq_along(col_I), column=col_I, bckg, time = time, bio, tech, stringsAsFactors = FALSE)
      saved_df$cond <- cond_int
      
    } else {
      
      match_pattern <- grep(input$pattern, names(data()))
      
      validate(
        need(match_pattern > 0, "Pattern could not be found in column names. Please enter another pattern for intensity columns") %then%
        need(input$bckg_pos > 0, "Enter position of background") %then%
        need(input$bio_pos > 0, "Enter position of biological replicates") %then%
        need(input$time_pos > 0, "Enter position of technical replicates") %then%
        need(input$tech_pos > 0, "Enter position of experimental conditions")
          
      )
      
      cond_int <- identify_conditions(data(),
                                      Column_intensity_pattern = input$pattern,
                                      bckg_pos = input$bckg_pos,
                                      bio_pos = input$bio_pos,
                                      time_pos = input$time_pos, 
                                      tech_pos = input$tech_pos,
                                      split = input$split)
      saved_df$cond_data <- cond_int
      saved_df$cond <- cond_int
      
    }
    
    
    # ibait <- which(as.character(data()[[input$column_gene_name]]) == input$bait_gene_name);
    # 
    # validate(
    #   #need(!is.factor(data()[[input$column_gene_name]]), "Column for 'Gene.Names' contains factors. Please try changing the decimal separator (most likely '.' or ',') used for importing the data") %then%
    #   need(length(ibait)>0,
    #        paste("Could not find bait '", input$bait_gene_name,
    #              "' in column '",input$column_gene_name,
    #              "'. Please either :
    #               modify bait gene name in General Parameters,
    #               select another column for gene names
    #               or change the decimal separator (most likely '.' or ',') used for importing the data.", 
    #              sep=""))
    # )
    
    #saved_df$cond <- cond_int
    
    #cond_int
    
  })
  
  
  observe({
    check_var$is_numeric <- sapply(match(saved_df$cond$column, names(data())), function(x){is.numeric(data()[[x]])})
    check_var$ncol <- length(saved_df$cond$column)
  })
  
  
  
  # is_factor <- sapply(match(saved_df$cond$column, names(data())), function(x){is.factor(data()[[x]])})
  # 
  # validate(
  #   need( sum( is_factor ) < length(is_factor) , "All intensity columns are factors, try changing:
  #         the decimal separator (most likely '.' or ',') used for importing the data 
  #         the pattern used to identify intensity columns")
  #   )
  
  output$warning <- renderUI({
    if( input$mode == "Raw data"){
      if(!is.null(check_var$is_numeric)){
        if(sum( !check_var$is_numeric ) >0 ){
          wellPanel(
            helpText( paste( "Warning : ",sum( !check_var$is_numeric ), " intensity columns out of ", check_var$ncol, " were not numeric and will be discarded.
                             Try changing the decimal separator used for importing the data", sep="") )
            )
        }
        }
      else{
        NULL
      }
    }else{
      NULL
    }
    
    
  })
  
  observe({
    if(input$mode == "Raw data"){
      idx_ctrl_guess <- grep("WT", toupper(unique(saved_df$cond$bckg)))
      if(length(idx_ctrl_guess) > 0){
        ctrl_guess <- unique(saved_df$cond$bckg)[idx_ctrl_guess[1]]     
      }else{
        ctrl_guess <- unique(saved_df$cond$bckg)[1]
      }
      bait_guess <- setdiff(unique(saved_df$cond$bckg), ctrl_guess)[1]
      
      updateSelectInput(session, "bckg_bait",
                        choices = as.list(unique(saved_df$cond$bckg)),
                        selected = bait_guess)
      
      updateSelectInput(session, "bckg_ctrl",
                        choices = as.list(unique(saved_df$cond$bckg)),
                        selected = ctrl_guess)
      
      if(nchar(input$bait_gene_name) == 0 ){
        updateTextInput(session, "bait_gene_name", value = bait_guess)
      }
    }else if(input$mode == "Interactome"){
      updateSelectInput(session, "bckg_bait",
                        choices = as.list(unique(saved_df$cond$bckg)),
                        selected = saved_df$res$bckg_bait)
      
      updateSelectInput(session, "bckg_ctrl",
                        choices = as.list(unique(saved_df$cond$bckg)),
                        selected = saved_df$res$bckg_ctrl)
    }
    
    
  })
  
  observe({
     updateSelectInput(session, "time_selected",
                       choices = as.list(unique(as.character(saved_df$cond$time))),
                       selected = as.list(unique(as.character(saved_df$cond$time))))
  })
  
  observe({
    updateSelectInput(session, "tech_selected",
                             choices = as.list(unique(as.character(saved_df$cond$tech))),
                             selected = as.list(unique(as.character(saved_df$cond$tech)))) 
  })
  
  observe({
    updateSelectInput(session, "bio_selected",
                             choices = as.list(unique(as.character(saved_df$cond$bio))),
                             selected = as.list(unique(as.character(saved_df$cond$bio))))
  })
  
  observe({
    updateSelectInput(session, "name_focus", choices = ordered_Interactome()$names, selected = NULL)
  })
  
  observeEvent(input$apply_filter,{

    cond_int <- saved_df$cond
    cond_int$time <- factor(cond_int$time, levels=input$time_selected)
    cond_int$bio <- factor(cond_int$bio, levels=input$bio_selected)
    if(length(input$tech_selected)>0){
      cond_int$tech <- factor(cond_int$tech, levels=input$tech_selected)
    }else{
      cond_int$tech <- rep("tech_0", dim(cond_int)[1])
    }
    
    idx_cond_selected <- which(rowSums(is.na(cond_int)) == 0)
    cond_int <- cond_int[idx_cond_selected, ]
    saved_df$cond_select <- cond_int

    if(input$mode == "Interactome"){
      saved_df$cond_data <- saved_df$cond_select
    }
    
  })
  
  observeEvent(saved_df$data,{
    saved_df$cond <- NULL
  })
  
  prep_data <- reactive({
    
    if(input$mode == "Raw data"){
      validate(
        need(dim(saved_df$cond_select)[1]>0, "No sample selected. Please select samples and/or validate your sample selection in the QC / Select tab") %then%
          need( length(setdiff(c("column", "bckg", "time", "bio", "tech"), names(saved_df$cond))) == 0 , "Please identify conditions")
      )
      
      ibait <- which(data()[[input$column_gene_name]] == input$bait_gene_name);
      check_two_bckg <- length(unique(saved_df$cond_select$bckg))>1
      
      bio_sel <- setdiff(unique(saved_df$cond_select$bio), input$filter_bio)
      check_two_bckg_per_cond <- sum(sapply(unique(saved_df$cond_select$time[saved_df$cond_select$bio %in% bio_sel]),
                                            function(x) { 
                                              input$bckg_bait %in% saved_df$cond_select$bckg[saved_df$cond_select$time==x] & 
                                                input$bckg_ctrl %in% saved_df$cond_select$bckg[saved_df$cond_select$time==x] 
                                            }
      )
      ) == length(unique(saved_df$cond_select$time[saved_df$cond_select$bio %in% bio_sel]))
      
      
      check_two_bio_rep <- length(unique(saved_df$cond_select$bio))>1
      found_bait_bckg <- input$bckg_bait %in% saved_df$cond_select$bckg
      found_ctrl_bckg <- input$bckg_ctrl %in% saved_df$cond_select$bckg
      
      validate(
        need( sum(is.na(match(saved_df$cond_select$column, names(data())))) == 0, "Some column names could not be mapped to original sample names. Please check the file used to map samples" )
      )
      
      #is_factor <- sapply(match(saved_df$cond$column, names(data())), function(x){is.factor(data()[[x]])})
      
      validate(
        #need( sum( is_factor ) < length(is_factor) , "All intensity columns are factors, try changing the decimal separator (most likely '.' or ',') used for importing the data")%then%
        need(check_two_bckg, "Could not identify distinct backgrounds. Please verify the mapping of samples") %then%
          need(check_two_bio_rep, "Could not identify distinct biological replicates. Please verify the mapping of samples") %then%
          need(found_bait_bckg, paste("Could not find", input$bckg_bait ," in possible backgrounds. Please change the background name the Group tab")) %then%
          need(found_ctrl_bckg, paste("Could not find", input$bckg_ctrl ," in possible backgrounds. Please change the background name the Group tab")) %then%
          need(check_two_bckg_per_cond, "Could not identify bait and ctrl backgrounds for each experimental condition. Please verify the mapping of samples") %then%
          need(input$column_ID, "Please select the column containing protein IDs in the Import tab") %then%
          need(input$column_gene_name, "Please select the column containing gene names in the Import tab") %then%
          need(!input$bait_gene_name %in% c("", "Bait"), "Please enter the gene name of the bait (in General Parameters)") %then%
          need(length(ibait)>0,
               paste("Could not find bait '", input$bait_gene_name,"' in column '",input$column_gene_name,"'. Please modify bait gene name in General Parameters.", sep=""))
        
      )
      
      Column_gene_name <- input$column_gene_name #names(data())[grep("GENE", toupper(names(data())))[1]]
      Column_ID <- input$column_ID #names(data())[grep("ID", toupper(names(data())))[1]]
      
      
      preprocess_data(  df = data(),
                        Column_intensity_pattern = input$pattern,
                        bait_gene_name = input$bait_gene_name,
                        Column_gene_name = Column_gene_name,
                        Column_ID = Column_ID,
                        bckg_bait = input$bckg_bait,
                        bckg_ctrl = input$bckg_ctrl,
                        condition = saved_df$cond_select
      )
      
    } else if( input$mode == "Interactome"){
      validate(
        need(!is.null(saved_df$res), "Please load an Interactome")
      )
      
      idx_match <- match(rownames(saved_df$res$data$Intensity), saved_df$res$names)
      idx_match_col <- match(saved_df$cond_select$column, colnames(saved_df$res$data$Intensity) )
      list(
        Intensity = saved_df$res$data$Intensity[ , idx_match_col],
        conditions = saved_df$cond_select,
        Npep = saved_df$res$Npep[idx_match],
        Protein.IDs = saved_df$res$Protein.IDs[idx_match],
        names = saved_df$res$names[idx_match],
        bckg_bait = saved_df$res$bckg_bait,
        bckg_ctrl = saved_df$res$bckg_ctrl,
        bait_gene_name = saved_df$res$bait
      )
        
      
    }
    
  })
  
  
  observeEvent(input$start, {
    
    if(input$mode == "Raw data"){
      validate(
        need(!is.null(saved_df$data), "Please import data (see the Import tab)") %then%
          need(!is.null(saved_df$cond), "Please define protein intensity columns and the conditions associated to each column (see the Group and QC / Select tabs).") %then%
          need(input$start, "Please click on the Compute Interactome button in General parameters")
      )
    }
      
      
      # Create a Progress object
      progress <- shiny::Progress$new(min = 0, max = 100)
      progress$set(message = "Compute interactome...", value = 0)
      on.exit(progress$close())
      updateProgress <- function(value = NULL, detail = NULL) {
        progress$set(value = value, detail = detail)
      }
      
      
      preprocess_df <- prep_data()
      
      if(input$merge_conditions){
        preprocess_df$conditions$time <- rep("merge", length(preprocess_df$conditions$time))
      }
      
      res_int <- InteRact(preprocess_df = preprocess_df,
                          N_rep=input$N_rep,
                          method = input$method,
                          pool_background = input$pool_background,
                          updateProgress = updateProgress,
                          substract_ctrl = input$substract_ctrl,
                          use_mean_for_bait = input$use_mean_for_bait)
      
      #res_int <- merge_proteome(res_int)
      df_merge <- merge_conditions(res_int)
      df_FDR <- compute_FDR_from_asymmetry(df_merge)
      res_int <- append_FDR(res_int, df_FDR)
      #res_int <- append_PPI(res_int)
      
      saved_df$res <- res_int
      
  })
  
  
  
  res<- reactive({
    updateSelectInput(session, "volcano_cond",
                      choices = as.list(saved_df$res$conditions),
                      selected = NULL)
    updateSelectInput(session, "Stoichio2D_cond",
                      choices = as.list(c("max", saved_df$res$conditions)),
                      selected = "max")
    saved_df$res
  })
  
  

  ordered_Interactome <- reactive({
    
    if(input$mode == "Interactome"){
      validate(
        need(!is.null(saved_df$res), "No interactome loaded. Please import an interactome (see Import tab)")
      )
    }
    
    if(input$mode == "Raw data"){
      validate(
          need(!is.null(saved_df$cond), "Please define protein intensity columns and the conditions associated to each column (see the Group and QC / Select tabs).") %then%
          need(input$start, "Please click on the Compute Interactome button in General parameters")
      )
    }
    
    res_int <- identify_interactors( res(),
                                     var_p_val = input$var_p_val, 
                                     p_val_thresh = input$p_val_thresh, 
                                     fold_change_thresh = input$fold_change_thresh, 
                                     n_success_min = input$n_success_min, 
                                     consecutive_success = input$consecutive_success)
    
    Ninteractors$x <- length(res_int$interactor)
    #res_int <- order_interactome(res_int)
    res_int
      
  })
  
  annotated_Interactome <- reactive({
    
      results <- append_annotations(ordered_Interactome(), saved_df$annot )
      results
      
  })

  #Observe functions -------------------------------------------------------------------
  
  observeEvent(input$format_names, {
    df<-saved_df$data
    names(saved_df$data) <- do.call(input$format_function, list(x=names(df), pattern = input$format_pattern, replacement = input$format_replacement, fixed=TRUE))
  })
  
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
        saved_df$annot <- get_annotations(res(), updateProgress = updateProgress)
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
  
  # df_corr_filtered <- reactive({
  #   
  #   df1 <- df_corr()
  #   df1 <- df1[df1$r_corr>=input$r_corr_thresh & df1$p_corr<=input$p_val_corr_thresh, ]
  #   
  #   net <- igraph::graph.data.frame(df1, directed=FALSE);
  #   net <- igraph::simplify(net)
  #   layout <- igraph::layout_nicely(net)
  #   cfg <- igraph::cluster_fast_greedy(as.undirected(net))
  #   net_d3 <- networkD3::igraph_to_networkD3(net, group = cfg$membership)
  #   
  #   # vatt <- vertex.attributes(net)
  #   # vertex_names <- as.character(vatt$name)
  #   # 
  #   # #layout <- data.frame(x=rnorm(length(vertex_names)), y=rnorm(length(vertex_names)))
  #   # #idx_vertex <- as.numeric(vertex_names)
  #   # #df_corr_plot$names <- names[idx_vertex]
  #   # 
  #   # df_corr_plot$names <- vertex_names
  #   # df_corr_plot$x <- layout[ , 1]
  #   # df_corr_plot$y <- layout[ , 2]
  #   # df_corr_plot$cluster <- as.factor(cfg$membership)
  #   # #df_corr_plot$cluster <- sample(10, length(vertex_names), replace = TRUE)
  #   # 
  #   # df1
  # })
  
  output$force_net <- renderForceNetwork({
    
    plot_correlation_network(df_corr = df_corr(), 
                             r_corr_thresh = input$r_corr_thresh,
                             p_val_thresh =  input$p_val_corr_thresh)

  })
  
  
  data_summary <- reactive({
    df <- data_raw()
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
    #data_cond()
    #cond_filter()
    saved_df$cond_data
  })
  
  
  observe({
    saved_df$summary <- summary_table(annotated_Interactome())
    updateCheckboxGroupInput(session, "columns_displayed",
                             choices = names(saved_df$summary),
                             selected=c("names", "max_stoichio", "max_fold_change", "min_p_val") )
    
  })  
  summaryTable <- reactive({
    
     saved_df$summary[ , input$columns_displayed]
    
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
                                                    1:Ninteractors$x, 
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
    validate(
      need("Copy_Number" %in% names(ordered_Interactome()), "Protein abundance not available")
    )
    plot_2D_stoichio(ordered_Interactome(),
                     xlim = ranges$x,
                     ylim=ranges$y,
                     condition = input$Stoichio2D_cond,
                     N_display=min(Ninteractors$x, input$Nmax2D) )[[1]]
  })
  
  Stoichio2D <- reactive({
    validate(
      need("Copy_Number" %in% names(ordered_Interactome()), "Protein abundance not available")
    )
    plot_2D_stoichio(ordered_Interactome(),
                     condition = input$Stoichio2D_cond,
                     N_display = min(Ninteractors$x, input$Nmax2D) )[[1]]
  })
  
  dotPlot <- reactive({
    p <- plot_per_condition(ordered_Interactome(),
                            color_var = input$var_p_val,   
                            idx_rows = min(input$Nmax, Ninteractors$x),
                            #idx_cols = match(input$time_selected, ordered_Interactome()$conditions),
                            idx_cols = ordered_Interactome()$conditions,
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
                   N_print=input$N_print,
                   asinh_transform = input$asinh_transform )
  })
  
  annotPlot <- reactive({
    plot_annotation_results(annotTable(),
                            method_adjust_p_val = input$method_adjust_p_val,
                            p_val_max = input$p_val_max,
                            fold_change_min = input$fold_change_min,
                            N_annot_min = input$N_annot_min)
  })
  
  observeEvent(input$dotPlot_hover, {
    
    if(!is.null(input$dotPlot_hover)){
      select_dotPlot$i_prot <- round(-input$dotPlot_hover$y)
      select_dotPlot$i_cond <- round(input$dotPlot_hover$x)
      select_dotPlot$name <- ordered_Interactome()$names[idx_order$cluster[select_dotPlot$i_prot]]
    }
    
  })
  
  observeEvent(input$name_focus, {
    select_dotPlot$name <- input$name_focus
  })
  
  stoichioPlot <- reactive({
    
    plot_stoichio(ordered_Interactome(), 
                  name = select_dotPlot$name,
                  test = "t.test",
                  test.args = list("paired" = FALSE),
                  conditions = ordered_Interactome()$conditions)
   
  })
  
  compPlot <- reactive({
    
    validate(
      need(select_dotPlot$name %in% ordered_Interactome()$names, "Protein not found. Please select another one.")
    )
    plot_comparison(ordered_Interactome(), 
                    names = select_dotPlot$name,
                    condition = ordered_Interactome()$conditions,
                    var_facet_x = "cond",
                    var_facet_y = "name")
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
        dist1=sqrt( (hover$x-log10(ordered_Interactome()$max_stoichio[1:min(Ninteractors$x, input$Nmax2D)]) )^2 +
                      (hover$y-log10(ordered_Interactome()$stoch_abundance[1:min(Ninteractors$x, input$Nmax2D)]) )^2)
      }else{
        dist1=sqrt( (hover$x-log10(ordered_Interactome()$stoichio[[input$Stoichio2D_cond]][1:min(Ninteractors$x, input$Nmax2D)]) )^2 +
                      (hover$y-log10(ordered_Interactome()$stoch_abundance[1:min(Ninteractors$x, input$Nmax2D)]) )^2)
      }
      
      select_stoichioPlot$min_dist1 <- min(dist1, na.rm=TRUE)
      select_stoichioPlot$i_min <- which.min(dist1)
    }
    
  })
  
  compPlot_volcano <- reactive({
    
    plot_comparison(ordered_Interactome(), 
                    name = ordered_Interactome()$names[select_volcanoPlot$i_min],
                    condition = ordered_Interactome()$conditions,
                    var_facet_x = "cond",
                    var_facet_y = "name")
  })
  
  QCPlot <- reactive({
      plot_QC(prep_data()) 
  })
    
  #Output Table functions -------------------------------------------------------------------------
  
  output$condTable <- renderDataTable(condTable())
  output$condTable_bis <- renderDataTable(condTable())
  output$data_summary <- renderDataTable(data_summary())
  output$summaryTable <- renderDataTable({summaryTable()})
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
  output$QCPlot <- renderPlot(
    switch(input$QC_selected,
           "Intensity correlation" = QCPlot()[[1]],
           "Bait purification" = QCPlot()[[2]],
           "Missing values" = QCPlot()[[3]]
           )
  )

  
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
      
      unlink("./*", recursive = TRUE, force = TRUE)
      on.exit(setwd(".."))
      
      if( "volcano plot" %in% input$saved_items){
        # print volcano plots
        pdf(paste("./volcano.pdf", sep=""), 6, 6)
        print(all_volcanos())
        dev.off()
      }
      
      if( "2D stoichio plot" %in% input$saved_items){
        # print stoichio 2D plots
        pdf(paste("./stoichio_2D.pdf", sep=""), 4, 4)
        print( plot_2D_stoichio(ordered_Interactome(),
                                condition = c("max", ordered_Interactome()$conditions),
                                N_display = min(Ninteractors$x, input$Nmax2D) ))
        dev.off()
      
      }
      
      if( "dot plot" %in% input$saved_items){
        # print dot plots
        plot_width = 2.5 + length(ordered_Interactome()$conditions)/5
        plot_height = 1.5 + input$Nmax/5
        pdf(paste("./dot_plot.pdf", sep=""), plot_width, plot_height)
        print(dotPlot())
        dev.off()
        
        # print dot plots (all)
        plot_width = 2.5 + length(ordered_Interactome()$conditions)/5
        plot_height = 1.5 + Ninteractors$x/5
        pdf(paste("./dot_plot_all.pdf", sep=""), plot_width, plot_height)
        print(plot_per_condition(ordered_Interactome(),
                                 color_var = input$var_p_val,
                                 idx_rows = Ninteractors$x,
                                 idx_cols = match(input$time_selected, ordered_Interactome()$conditions),
                                 clustering = input$clustering))
        dev.off()
      }
      
      if( "stoichio plots" %in% input$saved_items){
        # print stoichio plots
        pdf(paste("./stoichio_per_interactor.pdf", sep=""), 3, 3)
        print(lapply(ordered_Interactome()$interactor, 
                     function(x){
                       plot_stoichio(res = ordered_Interactome(), 
                                     name = x,
                                     test="t.test", 
                                     test.args = list("paired"=FALSE))
                     }))
        dev.off()
      }
      
      if( "enrichment plots" %in% input$saved_items){
        # print stoichio 2D plots
        pdf(paste("./enrichment_per_interactor.pdf", sep=""), 2 + length(ordered_Interactome()$conditions), 3)
        print(lapply(ordered_Interactome()$interactor, 
                     function(x){
                       plot_comparison(res = ordered_Interactome(), 
                                     name = x,
                                     conditions =  input$time_selected,
                                     var_facet_x = "cond",
                                     var_facet_y = "name")
                     }))
        dev.off()
      }
      
      # write summary tables
      if( "summary table" %in% input$saved_items){
        write.table(summary_table(annotated_Interactome(),  add_columns = names(annotated_Interactome())),
                  paste("./summary_table.txt", sep=""),  sep = "\t", dec = ".", row.names = FALSE)
      }
      
      # save Interactome
      if("Interactome" %in% input$saved_items){
        Interactome <- annotated_Interactome()
        save(Interactome, file = "./Interactome.Rda")
      }
      
      # if("preprocessed data" %in% input$saved_items){
      #   prep_data <- res()$data
      #   save(prep_data, file = "./prep_data.Rda")
      # }
      
      if("correlation network" %in% input$saved_items){
        write.table(df_corr(), file = "./correlation_network.txt", sep = "\t", dec = ".", row.names = FALSE)
      }
      
      input_name <- rep("", length(names(input)))
      input_value <- rep("", length(names(input)))
      for(i in 1:length(names(input))){
        input_name[i] <- names(input)[i]
        input_value[i] <- paste(input[[names(input)[i]]], collapse = ";")
      }
      df_input <- data.frame(input_parameter = input_name, value = input_value)
      write.table(df_input, file = "./parameters.txt", sep = "\t", dec = ".", row.names = FALSE)
      
      tar(file)
      
    }
  )
  
  #Output Info functions -------------------------------------------------------------------------
  
  output$interactors<- renderPrint({
    s1 <- paste("# interactors : ", Ninteractors$x, sep="")
    cat(s1)
  })
  
  output$info_volcano_hover <- renderPrint({
    s1<-paste("name: ", ordered_Interactome()$names[ select_volcanoPlot$i_min ],sep="")
    s2<-paste("p_val: ", ordered_Interactome()$p_val[[input$volcano_cond]][ select_volcanoPlot$i_min ],sep="")
    s3<-paste("fold_change: ", ordered_Interactome()$fold_change[[input$volcano_cond]][ select_volcanoPlot$i_min ],sep="")
    s4<-paste("stoichio: ", ordered_Interactome()$stoichio[[input$volcano_cond]][ select_volcanoPlot$i_min ],sep="")
    cat(s1,s2,s3,s4,sep="\n")
  })

  output$info_Stoichio2D_zoom_hover <- renderPrint({
    if(!is.null(input$Stoichio2D_zoom_hover)){
      hover=input$Stoichio2D_zoom_hover
      if(input$Stoichio2D_cond == "max"){
        dist1=sqrt( (hover$x-log10(ordered_Interactome()$max_stoichio[1:min(Ninteractors$x, input$Nmax2D)]) )^2 +
                      (hover$y-log10(ordered_Interactome()$stoch_abundance[1:min(Ninteractors$x, input$Nmax2D)]) )^2)
      }else{
        dist1=sqrt( (hover$x-log10(ordered_Interactome()$stoichio[[input$Stoichio2D_cond]][1:min(Ninteractors$x, input$Nmax2D)]) )^2 +
                      (hover$y-log10(ordered_Interactome()$stoch_abundance[1:min(Ninteractors$x, input$Nmax2D)]) )^2)
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

      prot_idx <- idx_order$cluster[select_dotPlot$i_prot]
      cond_idx <- select_dotPlot$i_cond
      
      s1 <- paste("Name: ", ordered_Interactome()$names[ prot_idx ], sep="")
      s2 <- paste("Condition: ", ordered_Interactome()$conditions[ cond_idx ], sep="")
      s3 <- paste("p-value: ", ordered_Interactome()$p_val[[ cond_idx ]][ prot_idx ], sep="")
      s4 <- paste("fold-change: ", ordered_Interactome()$fold_change[[ cond_idx ]][ prot_idx ], sep="")
      s5 <- paste("stoichio: ", ordered_Interactome()$stoichio[[ cond_idx ]][ prot_idx ], sep="")
      s6 <- paste("norm_stoichio: ", ordered_Interactome()$norm_stoichio[[ cond_idx ]][ prot_idx ], sep="")
      cat(s1, s2, s3, s4, s5, s6,sep="\n")

  })

  
}

# Run the app
shinyApp(ui, server)
