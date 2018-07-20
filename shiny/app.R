# Load packages ----

library(shiny)
library(shinyBS)

library("InteRact")

#source("../R/InteRact.R")
# library(ggplot2)
# library(ggrepel)
# library(ggsignif)
# library(grid)
# library(mice)
# library(Hmisc)
# library(igraph)
# library(networkD3)

library(data.table)
library(BiocInstaller)

`%then%` <- shiny:::`%OR%`

options(repos = BiocInstaller::biocinstallRepos())
#getOption("repos")
options(shiny.maxRequestSize = 100*1024^2) #maximum file size is set to 100MB

# User interface ----
ui <- fluidPage(
  titlePanel("InteRact : Analysis of AP-MS data"),
  
  fluidRow(
    column(3,
           br(),
           wellPanel(
             h4("General parameters"),
             textInput("bait_gene_name", "Bait (gene name)", value = "Bait"),
             bsTooltip("bait_gene_name", 
                       "The gene name of the bait protein"),
             checkboxInput("pool_background", "pool ctrl intensities", value = FALSE),
             bsTooltip("pool_background", 
                       "Perform protein enrichment tests using control intensities from all conditions"),
             checkboxInput("substract_ctrl", "substract ctrl (stoichio)", value = FALSE),
             bsTooltip("substract_ctrl", 
                       "Substract protein intensity from ctrl background to compute interaction stoichiometry"),
             numericInput("Nrep", "Missing values : # replacements : ", value = 1),
             bsTooltip("Nrep", 
                       "Number of times missing values will be replaced. Use 0 if you do not want to replace missing values")
          ),
          wellPanel(
             h4("Interaction parameters"),
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
                       "Number of proteins that pass the detection criteria defined above"),
             downloadButton("download_all", "Save report"),
             bsTooltip("download_all", "Download a report of the analysis")
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
                                         fileInput("file", h4("Import protein intensity file :"), placeholder = "Enter file here"),
                                         checkboxInput("dec", "Use comma as decimal separator", value = FALSE)
                                       ),
                                       wellPanel(
                                         h4("Identify columns"),
                                         selectInput("column_gene_name",
                                                     "column for gene name",
                                                     choices = list(),
                                                     selected = NULL),
                                         bsTooltip("column_gene_name", 
                                                   "Choose the column containing gene names. This is where the Bait (gene name) defined in the general parameters panel should be."),
                                         selectInput("column_ID",
                                                     "column for protein ID",
                                                     choices = list(),
                                                     selected = NULL),
                                         bsTooltip("column_ID", 
                                                   "Choose the column containing protein IDs (from uniprot). This information is used to retrieve additional information such as GO annotations.")
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
                                         h4("Enter background"),
                                         selectInput("bckg_bait", "Bait background", choices = list(), selected = NULL ),
                                         bsTooltip("bckg_bait", 
                                                   "Enter the name of the bait background (as displayed in the table on the right)."),
                                         selectInput("bckg_ctrl", "Control background", choices = list(), selected = NULL ),
                                         #textInput("bckg_ctrl", "Control background", value = "WT"),
                                         bsTooltip("bckg_ctrl", 
                                                   "Enter the name of the control background (as displayed in the table on the right).")
                                         
                                       ),
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
                                           # textInput("preffix_bio", "biological replicates", value = "S"),
                                           # textInput("preffix_tech", "technical replicates", value = "R"),
                                           # textInput("preffix_time", "experimental conditions", value = "")
                                           
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
                                           bsTooltip("column_name", "Choose column containing protein intensity column names (corresponding to those shown in the import tab)"),
                                           selectInput("manual_bckg",
                                                              "background",
                                                              choices = list(),
                                                              selected = NULL),
                                           bsTooltip("manual_bckg", "Choose column containing the background name (id) for each sample"),
                                           selectInput("manual_bio",
                                                              "biological replicates",
                                                              choices = list(),
                                                              selected = NULL),
                                           bsTooltip("manual_bio", "Choose column containing the name (id) of the biological replicate for each sample"),
                                           selectInput("manual_tech",
                                                              "technical replicates",
                                                              choices = list(),
                                                              selected = NULL),
                                           bsTooltip("manual_tech", "Choose column containing the name (id) of the technical replicate for each sample"),
                                           selectInput("manual_time",
                                                              "experimental conditions",
                                                              choices = list(),
                                                              selected = NULL),
                                           bsTooltip("manual_time", "Choose column containing the name (id) of the experimental condition for each sample")
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
                                           bsTooltip("format_function", "Name of the function used to modify column names. The function gsub replaces all occurences of pattern by replacement while the function sub replaces only the first occurence."),
                                           textInput("format_pattern", "pattern", value = "."),
                                           bsTooltip("format_pattern", "pattern to be replaced in column names"),
                                           textInput("format_replacement", "replacement", value = "_"),
                                           bsTooltip("format_replacement", "character string replacing pattern"),
                                           actionButton("format_names","Format column names")
                                         )
                                         
                                       )
                                ),
                                column(8,
                                       br(),
                                       dataTableOutput("condTable")
                                )
                       ),
                       tabPanel("QC / Select",
                                column(4,
                                    br(),
                                    wellPanel(
                                        helpText("Select samples"),
                                        # checkboxGroupInput("filter_bio",
                                        #                    "Biological replicates (bio)",
                                        #                    choices = list(),
                                        #                    selected = NULL),
                                        # checkboxGroupInput("filter_tech",
                                        #                    "Technical replicates (tech)",
                                        #                    choices = list(),
                                        #                    selected = NULL),
                                        # checkboxGroupInput("filter_time",
                                        #                    "Experimental conditions (time)",
                                        #                    choices = list(),
                                        #                    selected = NULL),
                                        selectizeInput("bio_selected", "Select biological replicates", choices = list(), multiple = TRUE),
                                        selectizeInput("tech_selected", "Select technical replicates", choices = list(), multiple = TRUE),
                                        selectizeInput("time_selected", "Select experimental conditions (in order)", choices = list(), multiple = TRUE),
                                        actionButton("apply_filter", label = "Apply")
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

  saved_df <- reactiveValues(cond = NULL, data = NULL, annot = NULL, enrichment = NULL)
 
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
  
  data_raw <- reactive({
    
    validate(
      need(input$file$datapath, "Please select a file to import")
    )
    
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
    
    saved_df$data <- df
    
    df
    
  })
  
  data <- reactive({
    saved_df$data
  })
  
  data_cond <- reactive({
    if(input$manual_mapping){
      
      validate(
        need(input$file_cond$datapath, "Please select a file to import")
      )
      
      df_cond <- read.table(input$file_cond$datapath, 
                 sep=ifelse(input$sep_cond, ",", "\t"), 
                 fill=TRUE, 
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
      
      df_cond
    } else {

      cond()
    }
    
    
  })
  
  cond <- reactive({
    
    if(input$manual_mapping){
      
      df_cond <- data_cond()
      col_I <- df_cond[[input$column_name]]
      bckg <- df_cond[[input$manual_bckg]]
      time <- df_cond[[input$manual_time]]
      bio <- df_cond[[input$manual_bio]]
      tech <- df_cond[[input$manual_tech]]
      
      cond_int <- dplyr::tibble(idx=seq_along(col_I), column=col_I, bckg, time, bio, tech)

                      
    } else {
      
      match_pattern <- grep(input$pattern, names(data()))
      n_factor_col <- 0
      if(length(match_pattern) > 0){
        n_factor_col <- sum( sapply( match_pattern, function(x) is.factor(data()[[x]]) ) )
      }
      
      validate(
        need(match_pattern>0, "Pattern could not be found in column names. Please enter another pattern for intensity columns") %then%
          need(n_factor_col == 0,
               "Some intensity columns are factors, try changing the decimal separator (most likely '.' or ',') used for importing the data"
          )
      )
      
      cond_int <- identify_conditions(data(),
                                      Column_intensity_pattern = input$pattern,
                                      bckg_pos = input$bckg_pos,
                                      bio_pos = input$bio_pos,
                                      time_pos = input$time_pos, 
                                      tech_pos = input$tech_pos,
                                      split = input$split)
      
    }
    
    
    # updateCheckboxGroupInput(session, "bio_selected",
    #                          choices = as.list(unique(cond_int$bio)),
    #                          selected = as.list(unique(cond_int$bio))) 
    
    # updateCheckboxGroupInput(session, "filter_bio",
    #                          choices = as.list(unique(cond_int$bio)),
    #                          selected = NULL) 
    
    
                             
    # updateCheckboxGroupInput(session, "filter_bio",
    #                          choices = as.list(unique(cond_int$bio)),
    #                          selected = NULL) 
    # 
    # updateCheckboxGroupInput(session, "filter_tech",
    #                          choices = as.list(unique(cond_int$tech)),
    #                          selected = NULL) 
    
    # updateCheckboxGroupInput(session, "filter_time",
    #                          choices = as.list(unique(cond_int$time)),
    #                          selected = NULL)
    
    
    #saved_df$cond <- cond_int
    
    cond_int
    
  })
  
  observe({
    idx_ctrl_guess <- grep("WT", toupper(unique(cond()$bckg)))
    if(length(idx_ctrl_guess) > 0){
      ctrl_guess <- unique(cond()$bckg)[idx_ctrl_guess[1]]     
    }else{
      ctrl_guess <- unique(cond()$bckg)[1]
    }
    bait_guess <- setdiff(unique(cond()$bckg), ctrl_guess)[1]
    
    updateSelectInput(session, "bckg_bait",
                      choices = as.list(unique(cond()$bckg)),
                      selected = bait_guess)
    
    updateSelectInput(session, "bckg_ctrl",
                      choices = as.list(unique(cond()$bckg)),
                      selected = ctrl_guess)
    
    updateTextInput(session, "bait_gene_name", value = bait_guess)
    
  })
  observe({
    # updateSelectInput(session, "time_selected",
    #                   choices = as.list(unique(cond()$time[cond()$bio %in% input$bio_selected])),
    #                   selected = as.list(unique(cond()$time[cond()$bio %in% input$bio_selected])))
    
     updateSelectInput(session, "time_selected",
                       choices = as.list(unique(cond()$time)),
                       selected = as.list(unique(cond()$time)))
    
  })
  
  observe({
    updateSelectInput(session, "tech_selected",
                             choices = as.list(unique(cond()$tech)),
                             selected = as.list(unique(cond()$tech))) 
  })
  
  observe({
    updateSelectInput(session, "bio_selected",
                             choices = as.list(unique(cond()$bio)),
                             selected = as.list(unique(cond()$bio))) 
  })
  
  observe({
    updateSelectInput(session, "name_focus", choices = ordered_Interactome()$names, selected = NULL)
  })
  
  observeEvent(input$apply_filter,{
    cond_int <- cond()
    cond_int$time <- factor(cond_int$time, levels=input$time_selected)
    cond_int$bio <- factor(cond_int$bio, levels=input$bio_selected)
    cond_int$tech <- factor(cond_int$tech, levels=input$tech_selected)
    idx_cond_selected <- which(rowSums(is.na(cond_int)) == 0)
    cat(rowSums(is.na(cond_int)))
    cat(idx_cond_selected)
    cond_int <- cond_int[idx_cond_selected, ]
    saved_df$cond <- cond_int
  })
  
  # cond_filter <- reactive({
  #   if(input$apply_filter){
  #     cond_int <- cond()
  #     cond_int$time <- factor(cond_int$time, levels=input$time_selected)
  #     cond_int$bio <- factor(cond_int$bio, levels=input$bio_selected)
  #     cond_int$tech <- factor(cond_int$tech, levels=input$tech_selected)
  #     idx_cond_selected <- which(rowSums(is.na(cond_int)) == 0)
  #     cat(rowSums(is.na(cond_int)))
  #     cat(idx_cond_selected)
  #     cond_int <- cond_int[idx_cond_selected, ]
  #     saved_df$cond <- cond_int
  #     cond_int
  #   } else{
  #     cond()
  #   }
  #   
  #   
  # })
  
  prep_data <- reactive({
    
    ibait <- which(data()[[input$column_gene_name]] == input$bait_gene_name);
    check_two_bckg <- length(unique(cond()$bckg))>1
    
    bio_sel <- setdiff(unique(cond()$bio), input$filter_bio)
    check_two_bckg_per_cond <- sum(sapply(unique(cond()$time[cond()$bio %in% bio_sel]),
                                         function(x) { 
                                           input$bckg_bait %in% cond()$bckg[cond()$time==x] & 
                                             input$bckg_ctrl %in% cond()$bckg[cond()$time==x] 
                                         }
                                  )
                               ) == length(unique(cond()$time[cond()$bio %in% bio_sel]))
            
    
    check_two_bio_rep <- length(unique(cond()$bio))>1
    found_bait_bckg <- input$bckg_bait %in% cond()$bckg
    found_ctrl_bckg <- input$bckg_ctrl %in% cond()$bckg
    
    validate(
      need(dim(saved_df$cond)[1]>0, "Please validate your sample selection in the QC / Select tab") %then%
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
                      # preffix_bio = input$preffix_bio,
                      # preffix_tech = input$preffix_tech,
                      # preffix_time = input$preffix_time,
                      #filter_bio = input$filter_bio,
                      #filter_tech = input$filter_tech,
                      #filter_time = input$filter_time,
                      condition = saved_df$cond
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
    
    updateSelectInput(session, "volcano_cond",
                      choices = as.list(res_int$Interactome$conditions),
                      selected = NULL)
    updateSelectInput(session, "Stoichio2D_cond",
                      choices = as.list(c("max", res_int$Interactome$conditions)),
                      selected = "max")
    
    res_int
    
  })

  ordered_Interactome <- reactive({
      res_int <- identify_interactors (res()$Interactome,
                                     var_p_val = "p_val", 
                                     p_val_thresh = input$p_val_thresh, 
                                     fold_change_thresh = input$fold_change_thresh, 
                                     n_success_min = input$n_success_min, 
                                     consecutive_success = input$consecutive_success)
      Ninteractors$x <- length(res_int$interactor)
      res_int <- order_interactome(res_int)
      res_int
  })
  
  annotated_Interactome <- reactive({
      results <- append_annotations(ordered_Interactome(), saved_df$annot )
      names_excluded <- c("names","bait", "groups", "conditions")
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
    data_cond()
    #cond_filter()
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
    plot_2D_stoichio(ordered_Interactome(),
                     condition = input$Stoichio2D_cond,
                     xlim = ranges$x,
                     ylim = ranges$y,
                     N_display=min(Ninteractors$x, input$Nmax2D) )
  })
  
  Stoichio2D <- reactive({
    plot_2D_stoichio(ordered_Interactome(),
                     condition = input$Stoichio2D_cond,
                     N_display = min(Ninteractors$x, input$Nmax2D) )
  })
  
  dotPlot <- reactive({
    p <- plot_per_condition(ordered_Interactome(),
                        idx_rows = min(input$Nmax, Ninteractors$x),
                        idx_cols = match(input$time_selected, ordered_Interactome()$conditions),
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
                  conditions = input$time_selected)
   
  })
  
  compPlot <- reactive({
    
    validate(
      need(select_dotPlot$name %in% ordered_Interactome()$names, "Protein not found. Please select another one.")
    )
    plot_comparison(ordered_Interactome(), 
                    name = select_dotPlot$name,
                    condition = input$time_selected)
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
                    condition = ordered_Interactome()$conditions)
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
