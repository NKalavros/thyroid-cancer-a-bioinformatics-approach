limmavoom <- function(tcgadataframe,lfcfortreat = 0, moleculetype = "mRNA"){

  print("Loading required packages")
  
  controlmatrix <- t(subset(t(tcgadataframe),
                          substring(colnames(tcgadataframe),14,15) == "11"))
  
  conditionmatrix <- t(subset(t(tcgadataframe),
                            substring(colnames(tcgadataframe),14,15) != "11"))
  #Required packages
  library(limma)
  library(Glimma)
  library(edgeR)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  
  genes <- rownames(tcgadataframe)
  
  ensembl <- rownames(tcgadataframe)
  
  controlnames <- colnames(controlmatrix)
  
  conditionnames <- colnames(conditionmatrix)
  
  group <- as.factor(c(rep("control",length(controlnames)),
                       rep("condition",length(conditionnames))))
  
  print("Creating count matrix and performing independent filtering")
  
  counts <- cbind(controlmatrix,conditionmatrix)

  counts <- DGEList(counts)
  
  counts$samples$group <- group
  
  keep <- filterByExpr(counts, group = group)
  
  counts <- counts[keep,,
                   keep.lib.sizes = FALSE]
  
  genes <- genes[keep]
  
  ensembl <- ensembl[keep]
  
  print("Identifying batch variables through TCGA nomelcature")
  
  tss <- as.factor(substring(colnames(counts),6,7))
  
  batchid <- as.factor(substring(colnames(counts),22,25))
  
  center <- as.factor(substring(colnames(counts),27,28))
  
  portion <- as.factor(substring(colnames(counts),18,19))
  
  counts$samples$tss <- tss
  
  counts$samples$batchid <- batchid
  
  counts$samples$center <- center
  
  counts$samples$portion <- portion
  
  if(moleculetype == "mRNA"){
  
  print("Performing gene annotation, since the molecule is mRNA")
    
  genes <- substring(genes,0,15)
  
  genes <- mapIds(org.Hs.eg.db, keys = genes,
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first")
  }
  
  print("Performing TMM normalization")
  
  counts <- calcNormFactors(counts,
                            method = "TMM")
  #MDS plot
  print("Creating MDS plot, this might take a while, integrating batch variables for examination")
  
  lcpm <- cpm(counts, log=TRUE)
  
  glMDSPlot(lcpm, labels= c(controlnames,conditionnames), top = 500,
            log = FALSE,
            gene.selection = "pairwise",
            groups= counts$samples[,c(1,4,5,6,7)], main = 'Control vs Condition - no adjustment for batch variables',
            path = getwd(),
            folder = paste("Plots",project,sep = "_"),launch = FALSE,
            html = paste("MDS-plot",moleculetype,sep = "_"))
  
  
  print("MDS plot creation complete, proceeding to differential expression analysis")
  
  removetss <- 0
  removebatchid <- 0
  removecenter <- 0
  removeportion <- 0
  
  if(length(levels(center))< 2){
    removecenter <- 1
  }
  
  if(length(levels(batchid))< 2){
    removebatchid <- 1
  }
  
  if(length(levels(tss))< 2){
    removetss <- 1
  }
  
  if(length(levels(portion))< 2){
    removeportion <- 1
  }
  
  if((removeportion | removetss | removebatchid | removecenter)){
    if(removeportion == 1){
      design <- model.matrix(~0+group+tss+batchid+center)
    }
    if(removecenter == 1){
      design <- model.matrix(~0+group+tss+batchid+portion)
    }
    if((removecenter + removeportion) == 2){
      design <- model.matrix(~0+group+tss+batchid)
    }
  }else{
  design <- model.matrix(~0+group+tss+batchid+center+portion)
  }
  
  colnames(design) <- gsub("group", "", colnames(design))
  
  print("Creating constant matrix and performing voom transformation")
  
  contr.matrix <- makeContrasts(
    controlvscondition = condition-control,
    levels = colnames(design))
  
  rownames(counts) <- genes
  
  counts$genes <- ensembl
  
  voomobject <- voom(counts, design, plot = FALSE)
  
  print("Voom transformation complete, fitting to linear model")
  
  vfit <- lmFit(voomobject, design)
  
  print("Fitting to contrasts")
  
  vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
  
  efit <- eBayes(vfit)
  
  print("Differential expression analysis complete, checking for LFC cutoff and creating MD plot")
  
  if(lfcfortreat != 0){
  tfit <- treat(vfit, lfc = lfcfortreat)
  
  dt <- decideTests(tfit)
  
  glMDPlot(tfit, counts = counts$counts, status = dt,
           anno = as.data.frame(cbind(genes,ensembl)),
           transform = TRUE,
           samples = colnames(counts),
           main = colnames(tfit)[1], side.main = "Genes",
           groups = group, launch = FALSE, path = getwd(),
           folder = paste("Plots",project,sep = "_"),
           html = paste("MD-plot",moleculetype,sep = "_"))
  
  print("MD plot creation complete, returning files.")
  
  return(topTreat(tfit, n = Inf))
  }
  
  dt <- decideTests(efit)
  
  glMDPlot(efit, counts = counts$counts, status = dt,
           anno = as.data.frame(cbind(genes,ensembl)),
           transform = TRUE,
           samples = colnames(counts),
           main = colnames(efit)[1], side.main = "Genes",
           groups = group, launch = FALSE, path = getwd(),
           folder = paste("Plots",project,sep = "_"),
           html = paste("MD-plot",moleculetype,sep = "_"))

  return(topTreat(efit,n = Inf))
}
