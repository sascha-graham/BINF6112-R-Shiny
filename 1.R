
# Load library ------------------------------------------------------------

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


# Load data and statistics analysis---------------------------------------------------------------

dat   <- read.csv("Alldat.csv",row.names = 1)
pheno <- read.csv("pheno.csv", row.names = 1)

# Normalise ---------------------------------------------------------------

c <- pheno[colnames(dat), "groups"]
y <- DGEList(counts=dat, group=c, genes=rownames(dat))
y <- cpm(calcNormFactors(y, method="TMM"), log = TRUE)


# Filter ------------------------------------------------------------------

keep <- filterByExpr(y, group = c, min.count = log2(10))
y <- y[keep,]
normlise.count.dat<-data.frame(y)


# Differential gene expression analysis using "limma" ---------------------

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



#code my

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
print(reslfc)


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

print(totalup)

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

print(down2hour)

