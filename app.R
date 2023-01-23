library(shiny)
library(ggplot2)
library(kableExtra)
library(tidyverse)
library(igraph)


source("./R/PCGII.R")
source("./R/Utility.R")


my_ui <- fluidPage(
  
  titlePanel("PCGII, Interactive Shiny App"),
  
  sidebarLayout(position = "left",
                
                sidebarPanel(
                  fileInput(inputId="exp_file", 
                            label="CHOOSE CSV FILE", 
                            multiple=FALSE,
                            accept=".csv",
                            placeholder= "No file selected"),
                  # fileInput(inputId="file1", 
                  #           label="CHOOSE CSV FILE",
                  #           accept = c(
                  #             "text/csv",
                  #             "text/comma-separated-values,text/plain",
                  #             ".csv"),placeholder= "No file selected"
                  # )
                  checkboxInput(inputId="show_data", 
                                label="Show Input Data Table", FALSE),
                  numericInput(inputId = "nominal_fdr", 
                               label = "Nominal FDR", value = .05, min = 0.001, max = 0.2, step = .001),
                  actionButton(inputId = "run_PCGII", label="Run PCGII now")
                  
                ),
                mainPanel(
                  tableOutput(outputId = "InputDataTable"),
                  plotOutput(outputId = "networkplot")
                )
  )
)

my_server <- function(input, output,session) {
  # contents <- reactive({
  #   
  #   shiny::validate(
  #     need(input$exp_file, "Select a csv file!")
  #   )
  #   
  #   expFile <- input$exp_file
  #   
  #   if (is.null(expFile))
  #     return(NULL)
  #   
  #   read.csv(expFile$datapath, header =TRUE)
  # })
  
  output$InputDataTable <- function() {
    shiny::validate(
      need(input$exp_file, "Waiting for file!")
    )
    
    expFile <- input$exp_file
    
    if (is.null(expFile))
      return(NULL)
    
    df <- read.csv(expFile$datapath, header =TRUE) %>% as.data.frame
    
    showDat <- input$show_data
    if(showDat){
    head(df[,1:10]) %>%
      knitr::kable("html") %>%
      kable_styling("striped", full_width = F) 
    }
  }

  
  current_FDR <- eventReactive(input$run_PCGII, {
    input$nominal_fdr
  })
  
  output$networkplot <- renderPlot({
    expFile <- input$exp_file
    
    if (is.null(expFile))
      return(NULL)
    
    df <- read.csv(expFile$datapath, header =TRUE)
    
    set.seed(1234567)
    n=dim(df)[1] # sample size
    p=dim(df)[2] # number of nodes
    
    prior_set=matrix(data=c(9,15, 3,4, 5,24, 16,20, 25,22, 28,8, 11,4), nrow=7, ncol=2, byrow = TRUE)
    colnames(prior_set)=c("row", "col")
    PCGII_out=PCGII(df=df, prior=double_prior(prior_set), lambda = 2*sqrt(log(p)/n))
    inference_out=inference(list=PCGII_out, alpha = current_FDR())
    
    out=inference_out$sigs
    # create the dataframe of edges (edge set E)
    my_link=out  %>%
      transform(row = pmin(row, col), col = pmax(row, col)) %>% 
      arrange(row, col) %>% 
      unique() 
    colnames(my_link)[1:2]=c("from","to")
    # create node set, Gamma
    my_node=cbind.data.frame(id=1:p, gene=colnames(df)) 
    my_net <- graph_from_data_frame(d=my_link, vertices=my_node, directed=F) 
    # plot the network
    plot(my_net, edge.arrow.size=.2, vertex.frame.color="#ffffff", vertex.label=V(my_net)$gene, vertex.label.color="black",
         layout=layout_in_circle(my_net)) 
    
  })
}

shinyApp(ui=my_ui, server=my_server)
