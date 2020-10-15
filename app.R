#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# load libraries
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

# Define UI for application
ui <- fluidPage(

    # Application title
    titlePanel("pipeR Dual RNA-seq Visualizer"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput(inputId = "pathogenDataInput",
                        label = "Pathogen dataset:",
                        choices = c("Alldat.csv")),
            
            selectInput(inputId = "hostDataInput",
                        label = "Host dataset:",
                        choices = c("n/a")),
            
            selectInput(inputId = "customGroupingInput",
                        label = "Custom grouping:",
                        choices = c("pheno.csv"))
        ),
        

        # Show the input data
        mainPanel(
           
            # verbatimTextOutput("summary"),
            
            tableOutput("view")
        )
    )
)

# Define server logic 
server <- function(input, output) {
    
    samplePathogenData <- read.csv("Alldat.csv",row.names = 1)
    sampleGroupingData <- read.csv("pheno.csv", row.names = 1)
    

    # select datasets from dropdowns
    pathogenData <- reactive({
        switch(input$pathogenDataInput,
               "Alldat.csv" = samplePathogenData)
    })
    
    customGrouping <- reactive({
        switch(input$customGroupingInput,
               "pheno.csv" = sampleGroupingData)
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
    
    
    # Clustering with Mfuzz
    cluster <- reactive({
        y = adjustedPathogenData()
        dat = pathogenData()
        d = averageReplicates()
        
        y.dat<- as.matrix(d)
        y.dat <- y.dat[which(apply(y.dat, 1, var)>2 & apply(y.dat,1,mean)>2), 1:6]
        timepoint <- c(0,2,4,4,16,24)
        y.dat <- rbind(timepoint, y.dat)
        rownames(y.dat)[1]<- "time"
        tmp <- tempfile()
        write.table(y.dat,file=tmp, sep='\t',quote=FALSE, col.names=NA)
        z.data <- table2eset(tmp)
        data.z <-standardise(z.data)
        class(data.z)
        m1 <-mestimate(data.z)
        Dmin(data.z, m=m1, crange=seq(2,22,1), repeats = 3, visu = TRUE)
        clust=8
        c <- mfuzz(data.z, c=clust, m=m1)
        mfuzz.plot(data.z,cl=c,mfrow=c(4,4),min.mem=0.5,time.labels=c(0,2,4,4,16,24),new.window=FALSE)
        membership<-c$membership
        membership<-data.frame(membership)
        fd<-data.frame(cor(t(c[[1]])))
        acore<-acore(data.z,c,min.acore = 0.5)
        acore_list<-do.call(rbind,lapply(seq_along(acore), function(i){data.frame(CLUSTER=i, acore[[i]])}))
        # colnames(acore_list)[2]<-"gene_name" # again don't do this it will break SG
        browser()
        genelist<- acore(data.z,cl=c,min.acore=0.7)
        # temp <- do.call("rbind", lapply(genelist, FUN = function(x){
            # return(paste0(as.character(x$NAME), collapse = ","))
        # }))
        # TODO: this stuffs up because ClusterList is not formatted correctly
        # Cluster_list<-as.data.frame(temp)
        # colnames(Cluster_list) <-"gene_name"
        # Cluster_list<-str_split_fixed(Cluster_list$gene_name,",", n=Inf)
        # Cluster_list<-t(Cluster_list)
        # colnames(Cluster_list)<- c("Cluster1", "Cluster2","Cluster3","Cluster4","Cluster5","Cluster6","Cluster7","Cluster8")
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
        cluster()
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
