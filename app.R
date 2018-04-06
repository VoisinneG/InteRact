# Load packages ----
library(shiny)
library(ggplot2)
library(ggrepel)
library(grid)
library(DT)

source("./R/InteRact.R")


options(shiny.maxRequestSize = 20*1024^2) #maximum file size is set to 30MB

# Source helpers ----
#source("./helpers.R")

# User interface ----
ui <- fluidPage(
  titlePanel("InteRact : Analysis of AP-MS data"),
  
  fluidRow(
    # column(4, 
    #        wellPanel(
    #           h2("Import"),
    #           fileInput("file", h4("ProteinGroups file :"), placeholder = "Enter file here"),
    #           checkboxInput("dec", "Use comma as decimal separator", value = FALSE)
    #        ),
    #        wellPanel(
    #          h2("Options"),
    #          textInput("bait_gene_name", "Bait (gene name)", value = "Bait"),
    #          textInput("bckg_bait", "Enter name of Bait background", value = "Bait"),
    #          textInput("bckg_ctrl", "Enter name of Control background", value = "WT"),
    #          numericInput("Nrep", "Number of replicates", value = 1)
    #        )
    # ),
    column(12, 
           
           tabsetPanel(id = "inTabset",
                       tabPanel("Options",
                                column(4,
                                  br(),
                                  wellPanel(
                                    h2("Import"),
                                    fileInput("file", h4("ProteinGroups file :"), placeholder = "Enter file here"),
                                    checkboxInput("dec", "Use comma as decimal separator", value = FALSE)
                                  )
                                ),
                                column(4,
                                       br(),
                                       wellPanel(
                                         h2("General"),
                                         textInput("bait_gene_name", "Bait (gene name)", value = "Bait"),
                                         textInput("bckg_bait", "Name of Bait background", value = "Bait"),
                                         textInput("bckg_ctrl", "Name of Control background", value = "WT"),
                                         numericInput("Nrep", "Number of iterations (replacement of missing values)", value = 1)
                                       )
                                ),
                                column(4,
                                       br(),
                                       wellPanel(
                                         h2("Preffix:"),
                                         textInput("preffix_bio", "For biological replicates", value = "S"),
                                         textInput("preffix_tech", "For technical replicates", value = "R"),
                                         textInput("preffix_time", "For experimental conditions", value = "")
                                       )
                                )
                       ),
                       tabPanel("Conditions", 
                                column(4,
                                    br(),
                                    wellPanel(
                                        h2("Filter"),
                                        checkboxGroupInput("filter_bio", 
                                                           "filter_bio", 
                                                           choices = list(),
                                                           selected = NULL),
                                        checkboxGroupInput("filter_tech", 
                                                           "filter_tech", 
                                                           choices = list(),
                                                           selected = NULL),
                                        checkboxGroupInput("filter_time", 
                                                           "filter_time", 
                                                           choices = list(),
                                                           selected = NULL)
                                        
                                    )
                                ),
                                column(8,
                                  br(),
                                  DT::dataTableOutput("contents")
                                )
                       ),
                       tabPanel("Volcano plots", 
                                column(4,
                                       br(),
                                       wellPanel(
                                         selectInput("volcano_cond", h3("Select conditions"), 
                                                     choices = list(), selected = NULL),
                                         #textInput("volcano_cond", "Select conditions", value = ""),
                                         uiOutput("my_output_UI_1"),
                                         numericInput("N_print", "# proteins displayed (maximum) ", value = 15),
                                         actionButton("save_volcano", "Save results", value = FALSE)
                                       )
                                ),
                                column(6,
                                       br(),
                                       plotOutput("VolcanoPlot",width="400",height="400") 
                                )
                       ),
                       tabPanel("Dot Plot", 
                                column(4,
                                       br(),
                                       wellPanel(
                                         uiOutput("my_output_UI_2"),
                                         #numericInput("p_val_thresh", "p-value (maximum)", value = 0.01),
                                         #numericInput("fold_change_thresh", "fold-change (minimum)", value = 1),
                                         numericInput("Nmax", "# proteins displayed (maximum) ", value = 30),
                                         actionButton("save_dot", "Save results", value = FALSE)
                                       )
                                ),
                                column(4,
                                       br(),
                                       plotOutput("plot",width="300",height="700") 
                                )

                       ),
                       tabPanel("2D stoichio", 
                                column(6,
                                       br(),
                                       plotOutput("stoichio_2D",width="400",height="400") 
                                )
                                
                       ),
                       tabPanel("Summary Table", 
                                column(4,
                                       br(),
                                       wellPanel(
                                         checkboxInput("append_annotations", "Append annotations (This can take a few minutes)", value = FALSE),
                                         checkboxGroupInput("columns_displayed", 
                                                            "Columns displayed", 
                                                            choices = c("names", "max_stoichio", "max_fold_change", "min_p_val"), 
                                                            selected = c("names", "max_stoichio", "max_fold_change", "min_p_val")
                                                            ),
                                         actionButton("save_table", "Save summary table", value = FALSE)
                                       )
                                ),
                                column(8,
                                       br(),
                                       DT::dataTableOutput("summary")
                                )
                                
                       )
           )
    )
  )
              
  # sidebarLayout(
  #   sidebarPanel(
  #     
  #     h2("Import"),
  #     
  #     fileInput("file", h4("ProteinGroups file :"), placeholder = "Enter file here"),
  #   
  #     checkboxInput("dec", "Use comma as decimal separator", value = FALSE),
  #     
  #     h2("Options"),
  #     
  #     textInput("bait_gene_name", "Bait (gene name)", value = "Bait"),
  #     textInput("bckg_bait", "Enter name of Bait background", value = "Bait"),
  #     textInput("bckg_ctrl", "Enter name of Control background", value = "WT"),
  #     numericInput("Nrep", "Number of replicates", value = 1) 
  #     
  #     ),
  #   
  #   mainPanel( 
  #     h4("Description of conditions used :"),
  #     DT::dataTableOutput("contents") 
  #   )
  # )
)

# Server logic
server <- function(input, output, session) {
  
  # return a list of UI elements
  output$my_output_UI_1 <- renderUI({
    list(
      numericInput("p_val_thresh_1", "p-value (maximum)", value = p_val_thresh),
      numericInput("fold_change_thresh_1", "fold-change (minimum)", value = fold_change_thresh)
    )
  })
  
  output$my_output_UI_2 <- renderUI({
    list(
      numericInput("p_val_thresh_2", "p-value (maximum)", value = p_val_thresh),
      numericInput("fold_change_thresh_2", "fold-change (minimum)", value = fold_change_thresh)
    )
  })
  
  # initial values
  p_val_thresh <- 0.01
  fold_change_thresh <- 1
  
  observeEvent(input$p_val_thresh_1,{
    p_val_thresh <<- input$p_val_thresh_1
    updateNumericInput(session, "p_val_thresh_2", value=p_val_thresh)})
  observeEvent(input$p_val_thresh_2,{
    p_val_thresh <<- input$p_val_thresh_2
    updateNumericInput(session,"p_val_thresh_1", value=p_val_thresh)})
  observeEvent(input$fold_change_thresh_1,{
    fold_change_thresh <<- input$fold_change_thresh_1
    updateNumericInput(session,"fold_change_thresh_2", value=fold_change_thresh)})
  observeEvent(input$fold_change_thresh_2,{
    fold_change_thresh <<- input$fold_change_thresh_2
    updateNumericInput(session,"fold_change_thresh_1", value=fold_change_thresh)})
  
  observe({
    c_num <- input$bait_gene_name
    updateTextInput(session, "bckg_bait", value =  c_num)
  })
  
  
  data <- reactive({
    read.csv(input$file$datapath, sep="\t", fill=TRUE, na.strings="", dec=ifelse(input$dec,",",".") )
  })
  
  cond <- reactive({
    identify_conditions(data(),
                        bckg_bait = input$bckg_bait,
                        bckg_ctrl = input$bckg_ctrl,
                        preffix_time = input$preffix_time,
                        preffix_bio = input$preffix_bio, 
                        preffix_tech = input$preffix_tech )
  })
  
  res<- reactive({
    
    updateSelectInput(session, "volcano_cond", 
                      choices = as.list(setdiff(unique(cond()$time), input$filter_time) ), 
                      selected = NULL)
    
    res_int<-InteRact(data(), 
           bait_gene_name = input$bait_gene_name, 
           N_rep=input$Nrep, 
           bckg_bait = input$bckg_bait , 
           bckg_ctrl = input$bckg_ctrl,
           preffix_bio = input$preffix_bio,
           preffix_tech = input$preffix_tech,
           preffix_time = input$preffix_time,
           filter_bio = input$filter_bio,
           filter_tech = input$filter_tech,
           filter_time = input$filter_time,
           bckg=cond()$bckg, 
           time=cond()$time, 
           bio=cond()$bio, 
           tech=cond()$tech
           )
    res_int$Interactome <- merge_proteome(res_int$Interactome)
    res_int
    
  })
  
  
  output$contents <- DT::renderDataTable({
    updateCheckboxGroupInput(session, "filter_bio",
                             choices = as.list(unique(cond()$bio)),
                             selected = NULL)
    updateCheckboxGroupInput(session, "filter_tech",
                             choices = as.list(unique(cond()$tech)),
                             selected = NULL)
    updateCheckboxGroupInput(session, "filter_time",
                             choices = as.list(unique(cond()$time)),
                             selected = NULL)
    
    cond()
    #res()$conditions
  })
  
  annotated_Interactome <- reactive({
    append_annotations(res()$Interactome)
  })
  
  ordered_Interactome <- reactive({
    if(input$append_annotations){
      
      order_list <- get_order_discrete(annotated_Interactome(), p_val_thresh, fold_change_thresh )
      results <- order_interactome(annotated_Interactome(), order_list$idx_order)
      names_excluded <- c("names","bait", "groups","conditions")
      updateCheckboxGroupInput(session, "columns_displayed",
                               choices = as.list( setdiff(names(results), names_excluded),selected=c("names", "max_stoichio", "max_fold_change", "min_p_val")) )
      results
      
    }else{
      order_list <- get_order_discrete(res()$Interactome, p_val_thresh, fold_change_thresh )
      results <- order_interactome(res()$Interactome, order_list$idx_order)
      
      names_excluded <- c("names","bait", "groups","conditions")
      updateCheckboxGroupInput(session, "columns_displayed",
                               choices = as.list( setdiff(names(results), names_excluded),selected=c("names", "max_stoichio", "max_fold_change", "min_p_val")))
      results
      
    }
    
  })
  
  
  observe({
    if(input$save_table){
      save_dir = paste("~/desktop/results_",input$bait_gene_name,"/",sep="")
      dir.create(save_dir)
      write.table( summary_table(ordered_Interactome(),  add_columns = input$columns_displayed ), 
                   file=paste(save_dir,"results.txt",sep=""), sep="\t", dec=".", row.names = FALSE)
    } 
  })
  
  
  
  output$summary <- DT::renderDataTable({
        summary_table(ordered_Interactome(),  add_columns = input$columns_displayed )
  })
  
  output$stoichio_2D <- renderPlot(
    plot_2D_stoichio(res()$Interactome)
  )
  
  output$plot <- renderPlot(
    
    if(input$save_dot){
      trigger<-c(input$p_val_thresh_1, input$p_val_thresh_2, input$fold_change_thresh_1, input$fold_change_thresh_2)
      save_dir = paste("~/desktop/results_",input$bait_gene_name,"/",sep="")
      dir.create(save_dir)
      results=res$Interactome()
      
      save( results, file=paste(save_dir,"res.Rdata",sep="") )
      
      plot_per_conditions(Interactome_order, 
                          idx_rows = min(Nmax, order_list$Ndetect), 
                          size_var=size_var, 
                          size_range=size_range, 
                          color_var="p_val", 
                          color_breaks=p_val_breaks, 
                          save_file=save_file )
      
      plot(res()$Interactome, p_val_thresh = p_val_thresh, 
           fold_change_thresh=fold_change_thresh, Nmax = input$Nmax,
           save_file=paste(save_dir,"dot_plot.pdf",sep=""))
    }else{
      trigger<-c(input$p_val_thresh_1, input$p_val_thresh_2, input$fold_change_thresh_1, input$fold_change_thresh_2)
      plot(res()$Interactome, p_val_thresh = p_val_thresh, 
           fold_change_thresh=fold_change_thresh, Nmax = input$Nmax)
    }
    
  )
  
  save_file<-NULL
  
  observeEvent(input$save_volcano, {
    save_dir = paste("~/desktop/results_",input$bait_gene_name,"/",sep="")
    dir.create(save_dir)
    save_file <<- paste(save_dir,"volcano_plot.pdf",sep="")
    plot_volcanos( res()$Interactome, conditions = input$volcano_cond, save_file=save_file, p_val_thresh = p_val_thresh, 
                   fold_change_thresh=fold_change_thresh, N_print=input$N_print )
  })
  
  output$VolcanoPlot <- renderPlot(
    
    if(length(res())>0){
      trigger<-c(input$p_val_thresh_1, input$p_val_thresh_2, input$fold_change_thresh_1, input$fold_change_thresh_2)
      plot_volcanos( res()$Interactome, conditions = input$volcano_cond, p_val_thresh = p_val_thresh, 
                     fold_change_thresh=fold_change_thresh, N_print=input$N_print )
    }
    
     # if(input$save_volcano){
     #     trigger<-c(input$p_val_thresh_1, input$p_val_thresh_2, input$fold_change_thresh_1, input$fold_change_thresh_2)
     #     save_dir = paste("~/desktop/results_",input$bait_gene_name,"/",sep="")
     #     dir.create(save_dir)
     #     results=res()
     #     save( results, file=paste(save_dir,"res.Rdata",sep="") )
     #     save_file=paste(save_dir,"volcano_plot.pdf",sep="")
     #     if(input$volcano_cond %in% res()$Interactome$conditions){
     #       plot_volcanos( res()$Interactome, conditions = input$volcano_cond, save_file=save_file, p_val_thresh = p_val_thresh, 
     #                      fold_change_thresh=fold_change_thresh, N_print=input$N_print )
     #     }
     #     else{
     #       plot_volcanos(res()$Interactome, save_file=save_file, p_val_thresh = p_val_thresh, 
     #                     fold_change_thresh=fold_change_thresh, N_print=input$N_print)
     #     }
     # }else{ 
     #     trigger<-c(input$p_val_thresh_1, input$p_val_thresh_2, input$fold_change_thresh_1, input$fold_change_thresh_2)
     #     save_file=NULL
     #     if(input$volcano_cond %in% res()$Interactome$conditions){
     #       
     #       plot_volcanos( res()$Interactome, conditions = input$volcano_cond, save_file=save_file, p_val_thresh = p_val_thresh, 
     #                      fold_change_thresh=fold_change_thresh, N_print=input$N_print )
     #     }
     #     else{
     #       plot_volcanos(res()$Interactome, save_file=save_file, p_val_thresh = p_val_thresh, 
     #                     fold_change_thresh=fold_change_thresh, N_print=input$N_print)
     #     }
     # }
    
    
  )
  
  
   
  
}

# Run the app
shinyApp(ui, server)
