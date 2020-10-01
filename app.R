#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(limma)
library(DESeq2)
library(dplyr)
library(igraph)
library(data.table)
library(Mfuzz)
library(edgeR)
library(stringr)
library(stringi)
library(ClueR)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("pipeR Dual RNA-seq Vizualizer"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput(inputId = "dataset",
                        label = "Choose a dataset:",
                        choices = c("salmonella", "option B", "option C"))
        ),

        # Show a plot of the generated distribution
        mainPanel(
           
            verbatimTextOutput("summary"),
            
            tableOutput("view")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    # Return the requested dataset
    datasetInput <- reactive({
        switch(input$dataset,
               "salmonella" = Alldat,
               "option B" = Alldat,
               "option C" = Alldat)
    })
    
    # Generate a summary of the dataset ----
    output$summary <- renderPrint({
        dataset <- datasetInput()
        summary(dataset)
    })
    
    # Show the first "n" observations ----
    output$view <- renderTable({
        datasetInput()
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
