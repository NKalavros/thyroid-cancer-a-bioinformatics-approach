library("devtools")
library(arrayConnector)
install_github("https://github.com/metaOmics/MetaQC")
library("MetaQC")
library("GEOquery")
library("arrayQualityMetrics")
BiocManager::install("hgu133a.db")
BiocManager::install("hgu133a2.db")
GSE58545 <- getGEO("GSE58545",GSEMatrix=T,AnnotGPL=FALSE)[[1]]
arrayQualityMetrics::arrayQualityMetrics(GSE58545)
GSE33630 <- getGEO("GSE33630",GSEMatrix=T,AnnotGPL=FALSE)[[1]]
GSE29265 <- getGEO("GSE29265",GSEMatrix=T,AnnotGPL=FALSE)[[1]]
GSE27155 <- getGEO("GSE27155",GSEMatrix=T,AnnotGPL=FALSE)[[1]]
GSE3678 <- getGEO("GSE3678",GSEMatrix=T,AnnotGPL=FALSE)[[1]]

mean(apply(exprs(GSE58545),2,max,na.rm = TRUE))
mean(apply(exprs(GSE33630),2,max,na.rm = TRUE))
mean(apply(exprs(GSE29265),2,max,na.rm = TRUE))
mean(apply(exprs(GSE27155),2,max,na.rm = TRUE))
mean(apply(exprs(GSE3678),2,max,na.rm = TRUE))

exprs(GSE27155) <- exprs(GSE27155)/4*14
exprs(GSE3678) <- log2(exprs(GSE3678))

annotation(GSE58545) <- "hgu133a"
annotation(GSE33630) <- "hgu133a2"
annotation(GSE29265) <- "hgu133a2"
annotation(GSE27155) <- "hgu133a"
annotation(GSE3678) <- "hgu133a2"

#Get only the samples we care for
GSE33630 <- GSE33630[,11:ncol(GSE33630)-1]
GSE29265 <- GSE29265[,10:ncol(GSE29265)]
GSE27155 <- GSE27155[,c(14:17,38:ncol(GSE27155) - 17)]

data(pathway)
GList <- pathway[[i]]
filterGenes <- TRUE
cutRatioByMean <- 0.3
cutRatioByVar <- 0.3
all_datasets <- list()
all_datasets$data <- list(GSE58545,GSE33630,GSE27155,GSE29265,GSE3678)
names(all_datasets$data) <- c("GSE58545","GSE33630","GSE27155","GSE29265","GSE3678")
all_datasets$data <- lapply(all_datasets$data,exprs)
samples <- lapply(all_datasets$data,ncol)
reversed_rep <- function(integer){return(rep(0,integer))}
samples <- lapply(samples,reversed_rep)
all_datasets$dataLabel <- samples
QCresult <- MetaQC(all_datasets$data, all_datasets$dataLabel, GList,filterGenes,cutRatioByMean,cutRatioByVar)

my_virtualArrays <- NULL

my_virtualArrays$noBatchEffect <- virtualArrayExpressionSets(covars = "all", supervised = TRUE)

arrayQualityMetrics(my_virtualArrays$noBatchEffect)
group <- as.factor(my_virtualArrays$noBatchEffect@phenoData@data$Covariate.1)


###################################################
### code chunk number 12: virtualArray.Rnw:246-250
###################################################
pData(my_virtualArrays$iPSC_hESC_noBatchEffect)[5] <- 
	c(as.character(pData(GSE23402)[,8]),as.character(pData(GSE26428)[,1]))
pData(my_virtualArrays$iPSC_hESC_noBatchEffect)[6] <- 
	c(rep("red",24),rep("blue1",3))


###################################################
### code chunk number 13: virtualArray.Rnw:268-271
###################################################
dist_iPSC_hESC_noBatchEffect <- 
	dist(t(exprs(my_virtualArrays$iPSC_hESC_noBatchEffect)), 
	method="euclidian")


###################################################
### code chunk number 14: virtualArray.Rnw:283-286
###################################################
hc_iPSC_hESC_noBatchEffect <- 
	hclust(dist_iPSC_hESC_noBatchEffect, method="average")
hc_iPSC_hESC_noBatchEffect$call <- NULL


###################################################
### code chunk number 15: virtualArray.Rnw:309-314
###################################################
virtualArrayHclust(hc_iPSC_hESC_noBatchEffect,
	lab.col=pData(my_virtualArrays$iPSC_hESC_noBatchEffect)[,6],
	lab=pData(my_virtualArrays$iPSC_hESC_noBatchEffect)[,5],
	main="batch effect removed",cex=0.7,
	xlab="sample names")


###################################################
### code chunk number 16: virtualArray.Rnw:371-372 (eval = FALSE)
###################################################
## my_virtualArrays$iPSC_hESC_supervised <- virtualArrayExpressionSets(supervised=TRUE)


###################################################
### code chunk number 17: virtualArray.Rnw:375-376
###################################################
my_virtualArrays$iPSC_hESC_supervised <- virtualArrayExpressionSets(supervised=TRUE,sampleinfo=system.file("extdata","sample_info_mod.txt",package="virtualArray"))


###################################################
### code chunk number 18: virtualArray.Rnw:379-382
###################################################
dist_iPSC_hESC_supervised <- 
	dist(t(exprs(my_virtualArrays$iPSC_hESC_supervised)), 
	method="euclidian")


###################################################
### code chunk number 19: virtualArray.Rnw:384-387
###################################################
hc_iPSC_hESC_supervised <<- 
	hclust(dist_iPSC_hESC_supervised, method="average")
hc_iPSC_hESC_supervised$call <- NULL


###################################################
### code chunk number 20: virtualArray.Rnw:390-395
###################################################
virtualArrayHclust(hc_iPSC_hESC_supervised,
	lab.col=pData(my_virtualArrays$iPSC_hESC_noBatchEffect)[,6],
	lab=pData(my_virtualArrays$iPSC_hESC_noBatchEffect)[,5],
	main="batch effect removed - supervised mode",cex=0.7,
	xlab="sample names")


###################################################
### code chunk number 21: virtualArray.Rnw:400-403
###################################################
pca_supervised <- prcomp(t(exprs(my_virtualArrays$iPSC_hESC_supervised)))
plot(pca_supervised$x, pch=19, cex=2, col=c(rep("red",24),rep("blue",3),pch=17))
legend("topleft",c("GSE23402","GSE26428"),col=c("red","blue"),pch=19,cex=1)


