
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
print(DE)

print(DE[[1]])
print(nrow(DE[[1]]))
print(row.names(DE[[1]])[1])
#code my

dataofgene <- read.csv("Annotation.csv")
#print(dataofgene)
#print(dataofgene[3,1])



d1 <- DE[[1]]
d2 <- DE[[2]]
d3 <- DE[[3]]
d4 <- DE[[4]]
d5 <- DE[[5]]

#print(row.names(d1))
#print(row.names(d1)[2])


#for (i in 1:nrow(d1)){
  #for (j in 1:nrow(dataofgene)) {
    #print(i)
    #print(j)
    #print(row.names(d1)[i])
    #print(dataofgene[j,3])
    
    #if (row.names(d1)[i] == dataofgene[j,3] ){
      #print(row.names(d1)[i])
      #row.names(d1)[i] <- dataofgene[j,2]
      #print(row.names(d1)[i])
    #}
    
  #}
#}
#print(d1)

Logfc2h <- DE[[1]]
Logfc2h <- Logfc2h[1]
ave2h <- DE[[1]]
ave2h <- ave2h[2]
pvalue2h <- DE[[1]]
pvalue2h <- pvalue2h[4]
adjp2h <- DE[[1]]
adjp2h <- adjp2h[5]

Logfc4h <- DE[[2]]
Logfc4h <- Logfc4h[1]
ave4h <- DE[[2]]
ave4h <- ave4h[2]
pvalue4h <- DE[[2]]
pvalue4h <- pvalue4h[4]
adjp4h <- DE[[2]]
adjp4h <- adjp4h[5]

Logfc8h <- DE[[3]]
Logfc8h <- Logfc8h[1]
ave8h <- DE[[3]]
ave8h <- ave8h[1]
pvalue8h <- DE[[3]]
pvalue8h <- pvalue8h[4]
adjp8h <- DE[[3]]
adjp8h <- adjp8h[5]

Logfc16h <- DE[[4]]
Logfc16h <- Logfc16h[1]
ave16h <- DE[[4]]
ave16h <- ave16h[2]
pvalue16h <- DE[[4]]
pvalue16h <- pvalue16h[4]
adjp16h <- DE[[4]]
adjp16h <- adjp16h[5]

Logfc24h <- DE[[5]]
Logfc24h <- Logfc24h[1]
ave24h <- DE[[5]]
ave24h <- ave24h[2]
pvalue24h <- DE[[5]]
pvalue24h <- pvalue24h[4]
adjp24h <- DE[[5]]
adjp24h <- adjp24h[5]

reslfc <- cbind(Logfc2h,ave2h)
reslfc <- cbind(reslfc,pvalue2h)
reslfc <- cbind(reslfc,adjp2h)
reslfc <- cbind(reslfc,Logfc4h)
reslfc <- cbind(reslfc,ave4h)
reslfc <- cbind(reslfc,pvalue4h)
reslfc <- cbind(reslfc,adjp4h)
reslfc <- cbind(reslfc,Logfc8h)
reslfc <- cbind(reslfc,ave8h)
reslfc <- cbind(reslfc,pvalue8h)
reslfc <- cbind(reslfc,adjp8h)
reslfc <- cbind(reslfc,Logfc16h)
reslfc <- cbind(reslfc,ave16h)
reslfc <- cbind(reslfc,pvalue16h)
reslfc <- cbind(reslfc,adjp16h)
reslfc <- cbind(reslfc,Logfc24h)
reslfc <- cbind(reslfc,ave24h)
reslfc <- cbind(reslfc,pvalue24h)
reslfc <- cbind(reslfc,adjp24h)
names(reslfc) <- c("logFc2h","AveExpr2h","P.Value2h","adj.P.Val2h","logFc4h","AveExpr4h","P.Value4h","adj.P.Val4h","logFc8h","AveExpr8h","P.Value8h","adj.P.Val8h","logFc16h","AveExpr16h","P.Value16h","adj.P.Val16h","logFc24h","AveExpr24h","P.Value24h","adj.P.Val24h")
#reslfc is the final result like:
#logFc2h AveExpr2h    P.Value2h  adj.P.Val2h   logFc4h AveExpr4h    P.Value4h
#E1WFR6  7.061139 11.149384 5.315688e-16 2.352724e-12  7.326080 11.149384 2.987711e-16
#E1WFR4  6.831966  9.968890 2.080575e-15 4.604313e-12  7.099757  9.968890 1.141955e-15
#E1WA29  6.579758 11.562188 3.777049e-15 5.572406e-12  6.553658 11.562188 4.018204e-15
#E1WFS5  7.198709  9.978317 1.177589e-14 1.303002e-11  7.272092  9.978317 1.005940e-14
#E1WFR1  6.786238  8.859809 1.866366e-14 1.498380e-11  6.193095  9.861629 3.202606e-14
#print(reslfc)
#**********see:   reslfc is the total result of the 2h to 24h.


#res2h is the result of 2h not the total one
res2h <- cbind(Logfc2h,ave2h)
res2h <- cbind(res2h,pvalue2h)
res2h <- cbind(res2h,adjp2h)
names(res2h) <- c("logFc2h","AveExpr2h","P.Value2h","adj.P.Val2h")

res4h <- cbind(Logfc4h,ave4h)
res4h <- cbind(res4h,pvalue4h)
res4h <- cbind(res4h,adjp4h)
names(res4h) <- c("logFc4h","AveExpr4h","P.Value4h","adj.P.Val4h")

res8h <- cbind(Logfc8h,ave8h)
res8h <- cbind(res8h,pvalue8h)
res8h <- cbind(res8h,adjp8h)
names(res8h) <- c("logFc8h","AveExpr8h","P.Value8h","adj.P.Val8h")

res16h <- cbind(Logfc16h,ave16h)
res16h <- cbind(res16h,pvalue16h)
res16h <- cbind(res16h,adjp16h)
names(res16h) <- c("logFc16h","AveExpr16h","P.Value16h","adj.P.Val16h")

res24h <- cbind(Logfc24h,ave24h)
res24h <- cbind(res24h,pvalue24h)
res24h <- cbind(res24h,adjp24h)
names(res24h) <- c("logFc2h","AveExpr2h","P.Value2h","adj.P.Val2h")



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

