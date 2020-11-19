
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
print(length(DE))

print(DE[[1]])
print(nrow(DE[[1]]))
print(row.names(DE[[1]])[1])
#code my
#while we donnot knwo the number of how many hoursdo we have we first
#assume H is like a list [0,2,4,8,16,24........]
H <- c(0,2,4,8,16,24)
#print(length(H))
number_of_hours = length(H)-1
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
  print(resl)
  result_of_all[[i]] <- resl
  
  
}
result_in_one <-result_of_all[[1]]
for (i in 1:(length(result_of_all)-1)){
  aaa <- result_of_all[[i+1]]
  result_in_one <- cbind(result_in_one,aaa)
  
}

print(result_in_one)






