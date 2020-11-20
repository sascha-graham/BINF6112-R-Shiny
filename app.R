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
library(tools)
library(shinythemes)
library(shinydashboard)
library(shinycssloaders)
library(shinybusy)

# Define UI for application that draws a histogram
ui <- function(request) {
    titlePanel("Dual RNA-seq Vizualizer")
    sidebar <- dashboardSidebar(
        hr(),
        sidebarMenu(id="tabs",
                    menuItem("Upload", tabName="table", icon=icon("table"), selected=TRUE),
                    menuItem("Pipeline Output", tabName = "plot", icon=icon("line-chart")),
                    menuItem("Example Data Formats", tabName = "about", icon = icon("question"))
        ),
        hr(),
        fluidRow(
            column(1),
            column(10,
                   h4("Uploaded files:"),
                   textOutput("pathFile"),
                   textOutput("hostFile"),
                   textOutput("groupFile"),
                   textOutput("annoFile")
            ),
        ),
        conditionalPanel("input.tabs == 'plot'",
                         hr(),
                         fluidRow(
                             column(1),
                             column(10,
                                    # h4("Other goodies"),
                                    # sliderInput("plot_width", "Plot width", value = 3, min = 1, max= 5, step = .5),
                                    # hr(),
                                    # h5("Save current state"),
                                    # bookmarkButton(),
                                    # hr(),
                                    h5("Download Plot"),
                                    downloadButton('downloadData', 'Download')
                             )
                         )
                         
        )
        
    )
    
    body <- dashboardBody(
        tags$head(
            tags$style(HTML('#clickUpload{background:#d6863a; border:#d6863a}'))
        ),
        tabItems(
            tabItem(tabName = "table",
                    box(width = NULL, status = "primary", solidHeader = TRUE, title="Upload Data",
                        fileInput("file1", "Upload RNA-seq Data",
                                  multiple = TRUE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        
                        fileInput("file3", "Upload Pheno Data",
                                  multiple = TRUE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        
                        fileInput("file4", "Upload Annotations",
                                  multiple = TRUE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                    ),
                    box(width = NULL, title = "Optional Cluster Optimisation", solidHeader = TRUE, status = "primary", 
                        actionButton("clickDMIN", "View DMIN Plot"),
                        conditionalPanel(condition = "input.clickDMIN", withSpinner(plotOutput('dmPlot', width = "50%"))),
                        numericInput("numClusters", 
                                     h5("Number of Clusters:"), 
                                     value = 4, min=2, max=12), 
                        numericInput("minMem", 
                                     h5("Minimum Membership:"), 
                                     value=0.3, min=0.1, max=1.0, step=0.1),
                    )
                    
            ),
            tabItem(tabName = "plot",
                    actionButton("clickUpload", "Run the Pipeline"),
                    hr(),
                    tabBox(width = NULL,
                           tabPanel(h5("Cluster Tabulation"),
                                    textOutput("error"),
                                    conditionalPanel(condition = "input.clickUpload", div(style = 'overflow-x: scroll', withSpinner(tableOutput("view"))))
                                    
                           ),
                           tabPanel(
                               # AMY i put this in - the formatting of the plots is terrible though
                               h5("Cluster Plots"),
                               conditionalPanel(condition = "input.clickUpload", div(style = 'overflow-x: scroll', withSpinner(plotOutput("cluster_plots", width = "75%"))))
                           ),
                           tabPanel(
                               h5("Differential Expression"),
                               conditionalPanel(condition = "input.clickUpload", div(style = 'overflow-x: scroll', withSpinner(tableOutput("diffExp"))))
                               # actionButton("clickCluster", "Clustering"),
                               # actionButton("clickDE", "Differential Expression")  
                           )
                    )
                    
            ),
            tabItem(tabName = "about",
                    box(width = NULL, title = "Example Data Format", solidHeader = TRUE, status = "primary", 
                        selectInput(inputId = "dataset",
                                    label = "Select a sample dataset to see formatting:",
                                    choices = c("Pathogen data", "Grouping data", "Annotations")),
                        
                        actionButton("clickExample", "Submit"),
                        p(),
                        div(style = 'overflow-x: scroll', tableOutput("content"))
                    )
            )
        )
    )
    
    dashboardPage(
        dashboardHeader(title = "Dual RNA-seq Vizualiser"),
        sidebar,
        body
    )
    
    
}

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    options(shiny.maxRequestSize=100*1024^2)
    
    samplePathogenData <- read.csv("Alldat.csv",row.names = 1)
    sampleGroupingData <- read.csv("Salmonella.Pheno.csv", row.names = 1)
    sampleAnnotation <- read.csv("Annotation.GO.Function.csv")
    
    
    # variables for num of clusters and min membership
    n_clusters <- reactive({
        n = input$numClusters
    })
    
    min_membership <- reactive({
        min = input$minMem
    })
    
    output$pathFile <- renderText(
        if (!is.null(input$file1$datapath)) {
            # Extract file name
            x <- sub(".csv$", "", basename(input$file1$name))
            return(x)
        })
    
    output$groupFile <- renderText(
        if (!is.null(input$file3$datapath)) {
            # Extract file name
            x <- c(sub(".csv$", "", basename(input$file3$name)))
            return(x)
        })
    
    output$annoFile <- renderText( 
        if (!is.null(input$file4$datapath)) {
            # Extract file name
            x <- c(sub(".csv$", "", basename(input$file4$name)))
            return(x)
        })
    
    # select datasets from dropdowns
    pathogenData <- reactive({
        read.csv(input$file1$datapath, row.names = 1)
    })
    
    customGrouping <- reactive({
        read.csv(input$file3$datapath, row.names = 1)
    })
    
    annotationData <- reactive({
        read.csv(input$file4$datapath)
    })
    
    # save custom groups
    customGroups <- reactive({
        c <- customGrouping()[colnames(pathogenData()), "groups"]
    })
    
    
    # average replicates across groups from raw data
    averageReplicates <- reactive({
        dat = pathogenData()
        pheno = customGrouping()
    #     c = customGroups()
    #     dat = pathogenData()
    #     d <- cbind(rowMeans(dat[,c(1:3)], na.rm = T),
    #                rowMeans(dat[,c(4:6)], na.rm = T),
    #                rowMeans(dat[,c(7:9)], na.rm = T),
    #                rowMeans(dat[,c(10:12)], na.rm = T),
    #                rowMeans(dat[,c(13:15)], na.rm = T),
    #                rowMeans(dat[,c(16:18)], na.rm = T))
    #     # TODO: can't add column names without overwriting the matrix SG
    #     # colnames(averageReplicates) <- c("Mean.0h","Mean.2h","Mean.4h","Mean.8h","Mean.16h","Mean.24h")
        d = aggregate(x=t(dat), by = list(pheno$groups), FUN=mean)
        ord = c()
        for (i in unique(pheno$groups)) { ord = c(ord, which(d[,1]==i)) }
        d = d[ord, ]
        d = data.matrix(t(d[,2:dim(d)[2]]))
        colnames(d) <- unique(pheno$groups)
        d
    })
    
    
    # normalized and filtered data
    adjustedPathogenData <- reactive({
        c = customGroups()
        y <- DGEList(counts=pathogenData(), group=c, genes=rownames(pathogenData()))
        y <- cpm(calcNormFactors(y, method="TMM"), log = TRUE)
        keep <- filterByExpr(y, group = c, min.count = log2(10))
        y <- y[keep,]
        y <- normlise.count.dat <- data.frame(y) # this is a change from original but I think it's correct SG
    })
    
    nGroups <- reactive({
        pheno <- customGrouping()
        n_groups = length(unique(as.vector(pheno$Time))) - 1
    })
    
    # using limma
    differentialExpressionPathogen <- reactive({
        c = customGroups()
        y = adjustedPathogenData()
        pheno = customGrouping()
        n_groups = nGroups()
        
        group = as.factor(c)
        design <- model.matrix(~ 0 + group)
        colnames(design) <- levels(group)
        rownames(design) <- colnames(y)
        fit <- lmFit(y, design = design)
        # cont.matrix <- makeContrasts(WT.02_h - WT.00_h,
        #                              WT.04_h - WT.00_h,
        #                              WT.08_h - WT.00_h,
        #                              WT.16_h - WT.00_h,
        #                              WT.24_h - WT.00_h,levels=design)
        # fit.cont <- eBayes(contrasts.fit(fit, cont.matrix))
        # summa.fit <- decideTests(fit.cont)
        cont.matrix <- aggregate(x=pheno[,(which(colnames(pheno) == 'groups')+1):dim(pheno)[2]], by = list(pheno$groups), FUN=mean)
        cont.matrix[cont.matrix == 0] <- -2
        cont.matrix[cont.matrix == -1] <- 0
        cont.matrix[cont.matrix == -2] <- -1
        ord = c()
        for (i in unique(pheno$groups)) { ord = c(ord, which(cont.matrix[,1]==i)) }
        cont.matrix = cont.matrix[ord, ]
        cont.matrix <- data.matrix(cont.matrix[,2:dim(cont.matrix)[2]])
        rownames(cont.matrix) <- unique(pheno$groups)
        fit.cont <- eBayes(contrasts.fit(fit, cont.matrix))
        summa.fit <- decideTests(fit.cont)
        
        
        de.ppi <- function(fit.cont, coef=1, lfc = 1, adjP =0.05){
            wtdt <- topTable(fit.cont, n = Inf, coef = coef)
            updt = with(wtdt, logFC > lfc & adj.P.Val < adjP)
            downdt = with(wtdt, logFC < lfc & adj.P.Val < adjP)
            wtdt$col="Not sig"
            wtdt$col[updt] = "Up"
            wtdt$col[downdt] = "Down"
            return(wtdt)}
        
        DE <- list()
        
        for (n in 1:n_groups){
            DE[[n]] <- de.ppi(fit.cont, coef = n)
        }
        DE
        
        # Differential.gene.expression.for.2h<-DE[[1]] # differential expression for time 2h as compared to 0h
    })
    
    differentialExpressionIndividual <- reactive({
        DE = differentialExpressionPathogen()
        number_of_hours = nGroups()
        result_of_all <- list()
        for (i in 1:number_of_hours){
            DAT <- DE[[i]]
            Logfc_DAT <- DAT[1]
            Ave_DAT <- DAT[2]
            PVA_DAT <- DAT[4]
            ADJP_DAT <- DAT[5]
            
            resl <- cbind(Logfc_DAT,Ave_DAT)
            
            resl <- cbind(resl,PVA_DAT)
            resl <- cbind(resl, ADJP_DAT)
            # colnames(out)<- c("Cluster", "GO term", "p-value", "Overlap size", "Overlapping genes")
            colnames(resl)[1] <- paste0("LogFC_Group",i)
            colnames(resl)[2] <- paste0("AveExpr_Group",i)
            colnames(resl)[3] <- paste0("P.Value_Group",i)
            colnames(resl)[4] <- paste0("adj.P.Val_Group",i)

            result_of_all[[i]] <- resl
        }
        result_of_all
    })
    
    differentialExpressionAll <- reactive({
        result_of_all <- differentialExpressionIndividual()
        result_in_one <- result_of_all[[1]]
        for (i in 1:(length(result_of_all)-1)){
            aaa <- result_of_all[[i+1]]
            result_in_one <- cbind(result_in_one,aaa)
        }
        result_in_one
    })
    
    # Clustering with Mfuzz
    cluster <- reactive({
        y = adjustedPathogenData()
        dat = pathogenData()
        d = averageReplicates()
        c = customGroups()
        n_clusters = n_clusters()
        min_membership = min_membership()
        pheno = customGrouping()
        
        # browser()
        y.dat<- as.matrix(d)
        y.dat <- y.dat[which(apply(y.dat, 1, var)>2 & apply(y.dat,1,mean)>2), 1:ncol(y.dat)]
        timepoint <- unique(pheno$groups)
        y.dat <- rbind(timepoint, y.dat)
        rownames(y.dat)[1]<- "time"
        tmp <- tempfile()
        write.table(y.dat,file=tmp, sep='\t',quote=FALSE, col.names=NA)
        z.data <- table2eset(tmp)
        data.z <-standardise(z.data)
        class(data.z)
        m1 <-mestimate(data.z)
        c <- mfuzz(data.z, c=n_clusters, m=m1)
    })
    
    # clueR enrichment
    enrich <- reactive({
        c = cluster()
        
        # Make annotation list
        anno <- annotationData()
        GO <- unique(anno$Gene.Ontology.ID)
        Uniprot.ID <- sapply(1:length(GO), function(i) paste(gsub("[[:space:]]", "", anno[which(anno$Gene.Ontology.ID==GO[i]),]$Uniprot.ID),collapse=" "))
        Uniprot.ID <- as.data.frame(Uniprot.ID)
        GO_Pro_ID <- data.frame(GO.ID=unique(anno$Gene.Ontology.ID),
                                Uniprot.ID=Uniprot.ID)
        
        # List of annotations
        Anno <- list()

        GO_IDs <- as.vector(GO_Pro_ID[,1])
        
        for (i in GO_IDs) {
            myindex <- which(GO_Pro_ID == i)
            Anno[i] <- strsplit(as.character(GO_Pro_ID[myindex, 2]), " ")
        }
        
        ce <- clustEnrichment(c, annotation=Anno, effectiveSize=c(2,100), pvalueCutoff=0.01)
        out <- c()
        i <- 1
        for (clus in ce$enrich.list) {
            clus <- cbind(rep(paste0("Cluster_",i), nrow(clus)), clus)
            i = i+1
            out <- rbind(out,clus)
        }
        
        #TODO: DO WE NEED TO CHANGE THIS?
        colnames(out)<- c("Cluster", "GO term", "p-value", "Overlap size", "Overlapping genes")
        out
    })
    
    # running the pipeline and outputting tables/plots
    observeEvent(input$clickUpload, {
        
        # outputs cluster tables
        output$view <- renderTable(
            rownames = TRUE,
            {
                tryCatch({
                    enrich()
                }, 
                warning=function(w) { 
                    output$error <- renderText("Genes in RNA-seq data are not included in annotation file.")
                    return(NA)
                }, 
                error=function(e) {
                    output$error <- renderText("Genes in RNA-seq data are not included in annotation file.")
                    return(NULL)
                })
            })
        
        # differential expression (NOT WORKING)
        output$diffExp <- renderTable(
            rownames = TRUE,
            {
                differentialExpressionAll()
            })
        
        # cluster plot output
        output$cluster_plots <- renderPlot({
            
            y = adjustedPathogenData()
            dat = pathogenData()
            d = averageReplicates()
            c = customGroups()
            n_clusters = n_clusters()
            min_membership = min_membership()
            colum <- ceiling(sqrt(n_clusters))
            rows <- ceiling(n_clusters/colum)
            pheno <- customGrouping()
            
            # browser()
            y.dat<- as.matrix(d)
            y.dat <- y.dat[which(apply(y.dat, 1, var)>2 & apply(y.dat,1,mean)>2), 1:ncol(y.dat)]
            timepoint <- unique(pheno$groups)
            y.dat <- rbind(timepoint, y.dat)
            rownames(y.dat)[1]<- "time"
            tmp <- tempfile()
            write.table(y.dat,file=tmp, sep='\t',quote=FALSE, col.names=NA)
            z.data <- table2eset(tmp)
            data.z <-standardise(z.data)
            class(data.z)
            m1 <-mestimate(data.z)
            c <- mfuzz(data.z, c=n_clusters, m=m1)
            mfuzz.plot(data.z,cl=c,mfrow=c(rows,colum),min.mem=min_membership,time.labels=timepoint, new.window = FALSE)
        })
        
    })
    
    # download tables
    output$downloadData <- downloadHandler(
        filename = function() {
            paste("data-", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
            write.csv(enrich(), file)
        }
    )
    
    # outputs dmin plot if requested
    # the purpose is to find the right number of clusters
    observeEvent(input$clickDMIN, {
        output$dmPlot <- renderPlot({
            y = adjustedPathogenData()
            d = averageReplicates()
            
            # browser()
            y.dat<- as.matrix(d)
            y.dat <- y.dat[which(apply(y.dat, 1, var)>2 & apply(y.dat,1,mean)>2), 1:6]
            timepoint <- c(0,2,4,8,16,24) # 8 hour timepoint was originally 4 SG
            y.dat <- rbind(timepoint, y.dat)
            rownames(y.dat)[1]<- "time"
            tmp <- tempfile()
            write.table(y.dat,file=tmp, sep='\t',quote=FALSE, col.names=NA)
            z.data <- table2eset(tmp)
            data.z <-standardise(z.data)
            class(data.z)
            m1 <-mestimate(data.z)
            Dmin(data.z, m=m1, crange=seq(2,16,1), repeats = 2, visu = TRUE)
        })
    })
    
    # allows users to view example formatting
    observeEvent(input$clickExample, {
        # Return the requested dataset
        datasetInput <- reactive({
            switch(input$dataset,
                   "Pathogen data" = samplePathogenData,
                   "Pheno data" = sampleGroupingData,
                   "Annotations" = sampleAnnotation)
        })
        
        
        # Show the first "n" observations ----
        output$content <- renderTable(
            rownames = TRUE, {
                datasetInput()
            })
        
    })
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
