
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
pheno <- read.csv("Salmonella Pheno.csv", row.names = 1)

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
cont.matrix <- aggregate(x=pheno[,(which(colnames(pheno) == 'groups')+1):dim(pheno)[2]], by = list(pheno$groups), FUN=mean)
cont.matrix[cont.matrix == 0] <- -2
cont.matrix[cont.matrix == -1] <- 0
cont.matrix[cont.matrix == -2] <- -1
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
for (n in 1:5){
  DE[[n]] <- de.ppi(fit.cont, coef = n)
}

Differential.gene.expression.for.2h<-DE[[1]]


# Average replecates across each time -------------------------------------

d = aggregate(x=t(dat), by = list(pheno$groups), FUN=mean)
d = t(d[,2:dim(d)[2]])
colnames(d) <- unique(pheno$groups)


# Clustering using "Mfuzz" ------------------------------------------------

y.dat<- as.matrix(d)
y.dat <- y.dat[which(apply(y.dat, 1, var)>2 & apply(y.dat,1,mean)>2), 1:6]
timepoint <- c(0,2,4,4,16,24)
y.dat <- rbind(timepoint, y.dat)
rownames(y.dat)[1]<- "time"
tmp<- tempfile()
write.table(y.dat,file=tmp, sep='\t',quote=FALSE, col.names=NA)
z.data <- table2eset(tmp)
data.z <-standardise(z.data)
class(data.z)
m1 <-mestimate(data.z)
Dmin(data.z, m=m1, crange=seq(2,22,1), repeats = 3, visu = TRUE)
clust=8
c<- mfuzz(data.z, c=clust, m=m1)
mfuzz.plot(data.z,cl=c,mfrow=c(4,4),min.mem=0.5,time.labels=c(0,2,4,4,16,24),new.window=FALSE)
membership<-c$membership
membership<-data.frame(membership)
fd<-data.frame(cor(t(c[[1]])))
acore<-acore(data.z,c,min.acore = 0.5)
acore_list<-do.call(rbind,lapply(seq_along(acore), function(i){data.frame(CLUSTER=i, acore[[i]])}))
colnames(acore_list)[2]<-"gene_name"
genelist<- acore(data.z,cl=c,min.acore=0.7)
temp <- do.call("rbind", lapply(genelist, FUN = function(x){
  return(paste0(as.character(x$NAME), collapse = ","))
}))
Cluster_list<-as.data.frame(temp)
colnames(Cluster_list) <-"gene_name"
Cluster_list<-str_split_fixed(Cluster_list$gene_name,",", n=Inf)
Cluster_list<-t(Cluster_list)
colnames(Cluster_list)<- c("Cluster1", "Cluster2","Cluster3","Cluster4","Cluster5","Cluster6","Cluster7","Cluster8")


# Make list ---------------------------------------------------------------

anno<-read.csv("Annotation.csv")
GO<-unique(anno$Gene.Ontology.ID)
Uniprot.ID <- sapply(1:length(GO), function(i) paste(gsub("[[:space:]]", "", anno[which(anno$Gene.Ontology.ID==GO[i]),]$Uniprot.ID),collapse=" "))
Uniprot.ID<-as.data.frame(Uniprot.ID)
GO_Pro_ID<-data.frame(GO.ID=unique(anno$Gene.Ontology.ID),
                      Uniprot.ID=Uniprot.ID)


# Make list of list -------------------------------------------------------

Anno <- list()
groupSize <- 422
GO_IDs <- as.vector(GO_Pro_ID[,1])

for (i in GO_IDs) {
  myindex <- which(GO_Pro_ID == i)
  Anno[i] <- strsplit(as.character(GO_Pro_ID[myindex, 2]), " ")
}


# Enrichment using "ClueR" ------------------------------------------------

ce <- clustEnrichment(c, annotation=Anno, effectiveSize=c(2,100), pvalueCutoff=0.01)

out <- c()
i <- 1
for (clus in ce$enrich.list) {
  clus<- cbind(rep(paste0("Cluster_",i), nrow(clus)), clus)
  out <- rbind(out,clus)
  i = i+1
}

write.csv(out, file = "./Enrich.csv", quote = F, row.names = F)




