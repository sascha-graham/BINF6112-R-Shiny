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

# Define UI for application that draws a histogram
ui <- function(request) {
    sidebar <- dashboardSidebar(
        hr(),
        sidebarMenu(id="tabs",
                    menuItem("Upload", tabName="table", icon=icon("table"), selected=TRUE),
                    menuItem("Pipeline Output", tabName = "plot", icon=icon("line-chart")),
                    menuItem("About", tabName = "about", icon = icon("question"))
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
        # tags$head(
        #   tags$style(HTML('#clickCluster{margin-right: 20px}'))
        # ),
        tabItems(
            tabItem(tabName = "table",
                    box(width = NULL, status = "primary", solidHeader = TRUE, title="Upload Data",
                        fileInput("file1", "Upload Pathogen Data",
                                  multiple = TRUE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        
                        fileInput("file2", "Upload Host Data",
                                  multiple = TRUE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        
                        fileInput("file3", "Upload Grouping Data",
                                  multiple = TRUE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        
                        fileInput("file4", "Upload Annotations",
                                  multiple = TRUE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        
                        actionButton("clickUpload", "Run Pipeline")  
                    ),
                    box(width = NULL, title = "Example Data Format", solidHeader = TRUE, status = "primary", 
                        selectInput(inputId = "dataset",
                                    label = "Select a sample dataset to see formatting:",
                                    choices = c("Pathogen data", "Grouping data", "Annotations")),
                        
                        actionButton("clickExample", "Submit"),
                        p(),
                        div(style = 'overflow-x: scroll', tableOutput("content"))
                    )
            ),
            tabItem(tabName = "plot",
                    tabBox(width = NULL,
                           tabPanel(h5("Clustering"),
                                    div(style = 'overflow-x: scroll', tableOutput("view"))
                                    
                           ),
                           tabPanel(
                               h5("Differential Expression"),
                               # actionButton("clickCluster", "Clustering"),
                               # actionButton("clickDE", "Differential Expression")  
                           )
                    )
                    
            ),
            tabItem(tabName = "about",
                    box(width = NULL, plotOutput("plot3", height="650px", brush = "plot_brush"), 
                        title = "About the Pipeline", solidHeader = TRUE, status = "primary")
                    
            )
            
        )
    )
    
    dashboardPage(
        dashboardHeader(title = "Host and Pathogen Dual-transcriptomic Profiles"),
        sidebar,
        body
    )
    
    
}

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    options(shiny.maxRequestSize=100*1024^2)
    
    samplePathogenData <- read.csv("Alldat.csv",row.names = 1)
    sampleGroupingData <- read.csv("pheno.csv", row.names = 1)
    sampleAnnotation <- read.csv("Annotation.csv")
    
    observeEvent(input$clickUpload, {
        
        output$pathFile <- renderText(
            if (!is.null(input$file1$datapath)) {
                # Extract file name
                x <- sub(".csv$", "", basename(input$file1$name))
                return(x)
            })
        
        output$hostFile <- renderText(
            if (!is.null(input$file2$datapath)) {
                # Extract file name
                x <- c(sub(".csv$", "", basename(input$file2$name)))
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
            c = customGroups()
            dat = pathogenData()
            d <- cbind(rowMeans(dat[,c(1:3)], na.rm = T),
                       rowMeans(dat[,c(4:6)], na.rm = T),
                       rowMeans(dat[,c(7:9)], na.rm = T),
                       rowMeans(dat[,c(10:12)], na.rm = T),
                       rowMeans(dat[,c(13:15)], na.rm = T),
                       rowMeans(dat[,c(16:18)], na.rm = T))
            # TODO: can't add column names without overwriting the matrix SG
            # colnames(averageReplicates) <- c("Mean.0h","Mean.2h","Mean.4h","Mean.8h","Mean.16h","Mean.24h")
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
        
        
        # using limma
        differentialExpressionPathogen <- reactive({
            c = customGroups()
            y = adjustedPathogenData()
            
            group = as.factor(c)
            design <- model.matrix(~ 0 + group)
            colnames(design) <- levels(group)
            rownames(design) <- colnames(y)
            fit <- lmFit(y, design = design)
            cont.matrix <- makeContrasts(WT.02_h - WT.00_h,
                                         WT.04_h - WT.00_h,
                                         WT.08_h - WT.00_h,
                                         WT.16_h - WT.00_h,
                                         WT.24_h - WT.00_h,levels=design)
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
            for (n in 1:5){
                DE[[n]] <- de.ppi(fit.cont, coef = n)
            }
            
            # Differential.gene.expression.for.2h<-DE[[1]] # differential expression for time 2h as compared to 0h
        })
        
        differentialExpression <- reactive({
            Logfc2h <- DE[[1]]
            Logfc2h <- Logfc2h[1]
            
            Logfc4h <- DE[[2]]
            Logfc4h <- Logfc4h[1]
            
            Logfc8h <- DE[[3]]
            Logfc8h <- Logfc8h[1]
            
            Logfc16h <- DE[[4]]
            Logfc16h <- Logfc16h[1]
            
            Logfc24h <- DE[[5]]
            Logfc24h <- Logfc24h[1]
            
            reslfc <- cbind(Logfc2h,Logfc4h)
            reslfc <- cbind(reslfc,Logfc8h)
            reslfc <- cbind(reslfc,Logfc16h)
            reslfc <- cbind(reslfc,Logfc24h)
            #reslfc is the final result like:
            #           logFC     logFC     logFC     logFC     logFC
            # E1WFR6  7.061139  7.326080  7.727751  7.431650  5.570146
            # print(reslfc)
            
            
            up1 <- subset(DE[[1]],col == "Up")
            down1 <- subset(DE[[1]],col == "Down")
            
            up2 <- subset(DE[[2]],col == "Up")
            down2 <- subset(DE[[2]],col == "Down")
            
            up3 <- subset(DE[[1]],col == "Up")
            down3 <- subset(DE[[1]],col == "Down")
            
            up4 <- subset(DE[[1]],col == "Up")
            down4 <- subset(DE[[1]],col == "Down")
            
            up5 <- subset(DE[[1]],col == "Up")
            down5 <- subset(DE[[1]],col == "Down")
            
            totalup <- rbind(up1,up2)
            totalup <- rbind(totalup,up3)
            totalup <- rbind(totalup,up4)
            totalup <- rbind(totalup,up5)
            
            totaldown <- rbind(up1,up2)
            totaldown <- rbind(totaldown,up3)
            totaldown <- rbind(totaldown,up4)
            totaldown <- rbind(totaldown,up5)
            
            # print(totalup)
            
            # From high to low the up&down data
            up2hour <- up1[order(up1$logFC),]
            #print(up2hour)
            up4hour <- up2[order(up2$logFC ),]
            up8hour <- up3[order(up3$logFC ),]
            up16hour <- up4[order(up4$logFC ),]
            up24hour <- up5[order(up5$logFC ),]
            
            
            down2hour <-  down1[order(down1$logFC ),]
            down4hour <-  down2[order(down2$logFC ),]
            down8hour <-  down3[order(down3$logFC ),]
            down16hour <- down4[order(down4$logFC ),]
            down24hour <- down5[order(down5$logFC ),]
            
            reslfc
        })
        
        
        # Clustering with Mfuzz
        cluster <- reactive({
            y = adjustedPathogenData()
            dat = pathogenData()
            d = averageReplicates()
            c = customGroups()
            
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
            # next line commented because it takes ages 
            # the purpose is to find the right number of clusters which we know to be 8 SG
            # Dmin(data.z, m=m1, crange=seq(2,22,1), repeats = 3, visu = TRUE)
            clust=8
            c <- mfuzz(data.z, c=clust, m=m1)
        })
        
        
        # we needed the new c frame above so it was made a separate function SG
        clusterList <- reactive({
            y = adjustedPathogenData()
            dat = pathogenData()
            d = averageReplicates()
            c = customGroups()
            
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
            # next line commented because it takes ages 
            # the purpose is to find the right number of clusters which we know to be 8 SG
            # Dmin(data.z, m=m1, crange=seq(2,22,1), repeats = 3, visu = TRUE)
            clust=8
            c <- mfuzz(data.z, c=clust, m=m1)
            mfuzz.plot(data.z,cl=c,mfrow=c(4,4),min.mem=0.5,time.labels=c(0,2,4,8,16,24),new.window=FALSE)
            membership <- c$membership
            membership <- data.frame(membership)
            fd <- data.frame(cor(t(c[[1]])))
            acore <- acore(data.z,c,min.acore = 0.5)
            acore_list<-do.call(rbind,lapply(seq_along(acore), function(i){data.frame(CLUSTER=i, acore[[i]])}))
            # colnames(acore_list)[2]<-"gene_name" # again don't do this it will break SG
            
            genelist<- acore(data.z,cl=c,min.acore=0.7)
            temp <- do.call("rbind", lapply(genelist, FUN = function(x){
                return(paste0(as.character(x$NAME), collapse = ","))
            }))
            # TODO: this stuffs up because Cluster_list is not formatted correctly??
            Cluster_list<-as.data.frame(temp)
            # colnames(Cluster_list) <-"gene_name"
            Cluster_list<-str_split_fixed(Cluster_list$NAME,",", n=Inf)
            Cluster_list<-t(Cluster_list)
            # colnames(Cluster_list)<- c("Cluster1", "Cluster2","Cluster3","Cluster4","Cluster5","Cluster6","Cluster7","Cluster8")
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
            groupSize <- 422
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
            out
        })
        
        
        # Generate a summary of a dataset
        
        # output$summary <- renderPrint({
        #    dataToShow <- pathogenData()
        #    summary(dataToShow)
        # })
        
        # Show the first "n" observations ----
        output$view <- renderTable(
            rownames = TRUE,
            {
                enrich()
            })
        
        output$downloadData <- downloadHandler(
            filename = function() {
                paste("data-", Sys.Date(), ".csv", sep="")
            },
            content = function(file) {
                write.csv(enrich(), file)
            }
        )
        
    })
    
    observeEvent(input$clickExample, {
        # Return the requested dataset
        datasetInput <- reactive({
            switch(input$dataset,
                   "Pathogen data" = samplePathogenData,
                   "Grouping data" = sampleGroupingData,
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
