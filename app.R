# Load packages ----
library(shiny)
library(ggplot2)
library(ggrepel)
library(grid)

source("./R/InteRact.R")


options(shiny.maxRequestSize = 30*1024^2) #maximum file size is set to 30MB


# User interface ----
ui <- fluidPage(
  titlePanel("InteRact : Analysis of AP-MS data"),
  
  fluidRow(
    
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
                                  dataTableOutput("contents")
                                )
                       ),
                       tabPanel("Volcano plots", 
                                column(3,
                                       br(),
                                       wellPanel(
                                         selectInput("volcano_cond", "Select condition", 
                                                     choices = list(), selected = NULL),
                                         #textInput("volcano_cond", "Select conditions", value = ""),
                                         uiOutput("my_output_UI_1"),
                                         numericInput("N_print", "# labels displayed (maximum) ", value = 15)
                                         
                                       )
                                ),
                                column(6,
                                       br(),
                                       fluidRow(
                                        downloadButton("download_volcano", "Download plot"),
                                        downloadButton("download_all_volcanos", "Download all volcanos")
                                       ),
                                       br(),
                                       #helpText("Hover mouse over point to display extra info"),
                                       #helpText("Brush and double-click to zoom"),
                                       plotOutput("volcano", width="400",height="400", 
                                                  hover = hoverOpts(id ="volcano_hover"),
                                                  dblclick = "volcano_dblclick",
                                                  brush = brushOpts(
                                                  id = "volcano_brush",
                                                  resetOnNew = TRUE) ),
                                       verbatimTextOutput("info_volcano_hover") 
                                ),
                                column(3,
                                       br(),
                                       wellPanel(
                                         helpText("Hover mouse over point to display extra info"),
                                         helpText("Brush and double-click to zoom")
                                       )
                                )
                       ),
                       tabPanel("Dot Plot", 
                                column(4,
                                       br(),
                                       wellPanel(
                                         uiOutput("my_output_UI_2"),
                                         #numericInput("p_val_thresh", "p-value (maximum)", value = NULL),
                                         #numericInput("fold_change_thresh", "fold-change (minimum)", value = NULL),
                                         numericInput("Nmax", "# proteins displayed (maximum) ", value = 30)
                                       )
                                ),
                                column(4,
                                       br(),
                                       downloadButton("download_dotPlot", "Download Plot", value = FALSE),
                                       plotOutput("dotPlot",width="250",height="500",
                                                  hover = hoverOpts(id ="dotPlot_hover"),
                                                  dblclick = "dotPlot_dblclick",
                                                  brush = brushOpts(
                                                    id = "dotPlot_brush",
                                                    resetOnNew = TRUE) )
                                ),
                                column(4,
                                       br(),
                                       wellPanel(
                                         helpText("Hover mouse over point to display extra info"),
                                         helpText("Brush and double-click to zoom")
                                       ),
                                       br(),
                                       verbatimTextOutput("info_dotPlot_hover") 
                                )

                       ),
                       tabPanel("2D stoichio", 
                                column(3,
                                       br(),
                                       wellPanel(
                                         selectInput("Stoichio2D_cond", "Select condition", 
                                                     choices = list(), selected = NULL),
                                         uiOutput("my_output_UI_3"),
                                         numericInput("Nmax2D", "# proteins displayed (maximum) ", value = 30)
                                       )
                                ),
                                column(width=4,
                                       br(),
                                       helpText("Brush to select zoom area"),
                                       downloadButton("download_Stoichio2D", "Download Plot", value = FALSE),
                                       plotOutput("Stoichio2D", width="300",height="300",
                                                  hover = hoverOpts(id ="Stoichio2D_hover"),
                                                  #dblclick = "Stoichio2D_dblclick",
                                                  brush = brushOpts(
                                                    id = "Stoichio2D_brush",
                                                    resetOnNew = TRUE) ),
                                       verbatimTextOutput("info_Stoichio2D_hover")
                                       
                                ),
                                column(width=4,
                                       br(),
                                       helpText("zoom on selected area"),
                                       downloadButton("download_Stoichio2D_zoom", "Download Plot", value = FALSE),
                                       plotOutput("Stoichio2D_zoom", width="300",height="300",
                                                  hover = hoverOpts(id ="Stoichio2D_zoom_hover")),
                                       br(),
                                       verbatimTextOutput("info_Stoichio2D_zoom_hover")
                                       
                                       
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
  
  # initial values
  params<-reactiveValues(p_val_thresh = 0.01, fold_change_thresh = 2)
  
  
  # return a list of UI elements
  output$my_output_UI_1 <- renderUI({
    list(
      numericInput("p_val_thresh_1", "p-value (maximum)", value = params$p_val_thresh),
      numericInput("fold_change_thresh_1", "fold-change (minimum)", value = params$fold_change_thresh)
    )
  })
  
  output$my_output_UI_2 <- renderUI({
    list(
      numericInput("p_val_thresh_2", "p-value (maximum)", value = params$p_val_thresh),
      numericInput("fold_change_thresh_2", "fold-change (minimum)", value = params$fold_change_thresh)
    )
  })
  
  output$my_output_UI_3 <- renderUI({
    list(
      numericInput("p_val_thresh_3", "p-value (maximum)", value = params$p_val_thresh),
      numericInput("fold_change_thresh_3", "fold-change (minimum)", value = params$fold_change_thresh)
    )
  })
  
  observeEvent(input$p_val_thresh_1, {params$p_val_thresh <- input$p_val_thresh_1})
  observeEvent(input$p_val_thresh_2, {params$p_val_thresh <- input$p_val_thresh_2})
  observeEvent(input$p_val_thresh_3, {params$p_val_thresh <- input$p_val_thresh_3})
  observeEvent(input$fold_change_thresh_1, {params$fold_change_thresh <- input$fold_change_thresh_1})
  observeEvent(input$fold_change_thresh_2, {params$fold_change_thresh <- input$fold_change_thresh_2})
  observeEvent(input$fold_change_thresh_3, {params$fold_change_thresh <- input$fold_change_thresh_3})
  
  observe({
    b_name <- input$bait_gene_name
    updateTextInput(session, "bckg_bait", value =  b_name)
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
    
    updateSelectInput(session, "Stoichio2D_cond", 
                      choices = as.list(setdiff( c("max", unique(cond()$time)), input$filter_time)), 
                      selected = NULL)
    
    cond()
    
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
  
  
  output$contents <- renderDataTable({
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
  
  
  order_list <- reactive({
    get_order_discrete(res()$Interactome, 
                       p_val_thresh = params$p_val_thresh, 
                       fold_change_thresh = params$fold_change_thresh )
  })
  
  ordered_Interactome <- reactive({
    if(input$append_annotations){
      
      results <- order_interactome(annotated_Interactome(), order_list()$idx_order)
      cat(order_list()$Ndetect)
      cat(params$p_val_thresh)
      names_excluded <- c("names","bait", "groups","conditions")
      updateCheckboxGroupInput(session, "columns_displayed",
                               choices = as.list( setdiff(names(results), names_excluded),selected=c("names", "max_stoichio", "max_fold_change", "min_p_val")) )
      results
      
    }else{
      results <- order_interactome(res()$Interactome, order_list()$idx_order)
      cat(order_list()$Ndetect)
      cat(params$p_val_thresh)
      names_excluded <- c("names","bait", "groups","conditions")
      updateCheckboxGroupInput(session, "columns_displayed",
                               choices = as.list( setdiff(names(results), names_excluded),selected=c("names", "max_stoichio", "max_fold_change", "min_p_val")))
      results
      
    }
    
  })
  
  summaryTable <- reactive({
    summary_table(ordered_Interactome(),  add_columns = input$columns_displayed)
  })
    
  output$summaryTable <- renderDataTable({summaryTable() })
  
  output$download_summaryTable <- downloadHandler(
    filename = "summary_table.txt",
    content = function(file) {
      write.table(summaryTable(), file, sep = "\t", dec = ".", row.names = FALSE)
    }
  )
  
  ranges <- reactiveValues(x = c(-1.5,0.5), y = c(-1,1))
  
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
  
  ranges_volcano <- reactiveValues(x = NULL, y = NULL)
  
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
  
  ranges_dotPlot <- reactiveValues(x = NULL, y = NULL)
  
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
  
  
  Stoichio2D <- reactive({
    plot_2D_stoichio(ordered_Interactome(),
                     condition = input$Stoichio2D_cond,
                     N_display=min(order_list()$Ndetect, input$Nmax2D) )
  })
  
  Stoichio2D_zoom <- reactive({
    plot_2D_stoichio(ordered_Interactome(),
                     condition = input$Stoichio2D_cond,
                     xlim = ranges$x,
                     ylim = ranges$y,
                     N_display=min(order_list()$Ndetect, input$Nmax2D) )
  })
  
  output$Stoichio2D <- renderPlot( Stoichio2D() )
  output$Stoichio2D_zoom <- renderPlot( Stoichio2D_zoom() )
  
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
  
  dotPlot <- reactive({
    plot_per_conditions(ordered_Interactome(), 
                        idx_rows = min(input$Nmax, order_list()$Ndetect))+
    coord_cartesian(xlim = ranges_dotPlot$x, ylim = ranges_dotPlot$y, expand = FALSE)
  })
  
  output$dotPlot <- renderPlot( dotPlot() )
  
  output$download_dotPlot <- downloadHandler(
    filename = "dot_plot.pdf",
    content = function(file) {
      plot_width = 3.4
      plot_height = input$Nmax/(plot_width+1) + 1 
      pdf(file,plot_width,plot_height)
      print(dotPlot())
      dev.off()
    }
  )
  
  
  volcano <- reactive({
      plot_volcanos( ordered_Interactome(),
                     conditions = input$volcano_cond,
                     p_val_thresh = params$p_val_thresh, 
                     fold_change_thresh = params$fold_change_thresh, 
                     N_print=input$N_print)[[1]]+
      coord_cartesian(xlim = ranges_volcano$x, ylim = ranges_volcano$y, expand = FALSE)
  })
  
  all_volcanos <- reactive({
    plot_volcanos( ordered_Interactome(),
                   p_val_thresh = params$p_val_thresh, 
                   fold_change_thresh = params$fold_change_thresh, 
                   N_print=input$N_print )
  })
  
  output$volcano <- renderPlot( volcano() )
  
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
  
  output$info_volcano_hover <- renderPrint({
    
    if(!is.null(input$volcano_hover)){
      hover=input$volcano_hover
      dist1=sqrt((hover$x-log10(ordered_Interactome()$fold_change[[input$volcano_cond]]) )^2 + 
                  (hover$y+log10(ordered_Interactome()$p_val[[input$volcano_cond]]) )^2)
      min_dist1 <- min(dist1, na.rm=TRUE)
      i_min <- which.min(dist1)
      if( min_dist1 < 0.25){
        #print( res()$Interactome$names[which.min(dist)] )
        s1<-paste("name: ", ordered_Interactome()$names[ i_min ],sep="")
        s2<-paste("p_val: ", ordered_Interactome()$p_val[[input$volcano_cond]][ i_min ],sep="")
        s3<-paste("fold_change: ", ordered_Interactome()$fold_change[[input$volcano_cond]][ i_min ],sep="")
        s4<-paste("stoichio: ", ordered_Interactome()$stoichio[[input$volcano_cond]][ i_min ],sep="")
        cat(s1,s2,s3,s4,sep="\n")
        
        #print(list(a=1, b=2))
        #input$volcano_hover$x
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
    
    if(!is.null(input$Stoichio2D_hover)){
      hover=input$Stoichio2D_hover
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
  
  output$info_dotPlot_hover <- renderPrint({
    
    if(!is.null(input$dotPlot_hover)){
      i_prot = round(-input$dotPlot_hover$y)
      i_cond = round(input$dotPlot_hover$x)
      s1 <- paste("Name: ", ordered_Interactome()$names[ i_prot ], sep="")
      s2 <- paste("Condition: ", ordered_Interactome()$conditions[ i_cond ], sep="")
      s3 <- paste("p-value: ", ordered_Interactome()$p_val[[ i_cond ]][ i_prot ], sep="")
      s4 <- paste("fold-change: ", ordered_Interactome()$fold_change[[ i_cond ]][ i_prot ], sep="")
      s5 <- paste("stoichio: ", ordered_Interactome()$stoichio[[ i_cond ]][ i_prot ], sep="")
      s6 <- paste("norm_stoichio: ", ordered_Interactome()$norm_stoichio[[ i_cond ]][ i_prot ], sep="")
      cat(s1, s2, s3, s4, s5, s6,sep="\n")
    }
    
  })
  
}

# Run the app
shinyApp(ui, server)
