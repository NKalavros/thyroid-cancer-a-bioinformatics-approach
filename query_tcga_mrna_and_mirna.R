obtain_data_for_mrna_and_mir <- function(project){
  #Load packages
  require(TCGAbiolinks)
  require(data.table)
  #Begin with mRNA
  
  if(!file.exists(paste(project,"mrnacounts",sep = "_"))){
    

    projectmrnaquery <- GDCquery(project = project,
                                 data.category = "Transcriptome Profiling",
                                 data.type = "Gene Expression Quantification",
                                 workflow.type = "HTSeq - Counts")
    
    #Download the data
    GDCdownload(projectmrnaquery,
                method = "client")
    
    #Prepare the data
    mrnacounts <- GDCprepare(projectmrnaquery, save = FALSE, summarizedExperiment = FALSE)
    
    #Edit the mRNA sample matrix as needed
    mrnacounts <- mrnacounts[6:nrow(mrnacounts),]
    
    rownames(mrnacounts) <- mrnacounts$X1
    
    mrnacounts <- mrnacounts[,2:ncol(mrnacounts)]
    
    #Write the new matrices in the working directory
    write.table(mrnacounts,
                file = paste(project,"mrnacounts",sep = "_"),
                quote = TRUE, sep = "\t")
  }else{
    mrnacounts <- fread(paste(project,"mrnacounts",sep = "_"))
    rownames <- mrnacounts$V1
    mrnacounts <- mrnacounts[,2:ncol(mrnacounts)]
    rownames(mrnacounts) <- rownames
  }
  
  assign("mrnacountmatrix", mrnacounts, envir = .GlobalEnv)
  
  #Same for miRs
  
  if(!file.exists(paste(project,"mircounts",sep = "_"))){
    projectmirquery <- GDCquery(project = project,
                                data.category = "Transcriptome Profiling",
                                data.type = "miRNA Expression Quantification",
                                workflow.type = "BCGSC miRNA Profiling")
    
    GDCdownload(projectmirquery,
                method = 'client')
    
    mircounts <- GDCprepare(projectmirquery, save = FALSE, summarizedExperiment = FALSE)
    
    #Edit the miR sample matrix as needed
    rownames(mircounts) <- mircounts$miRNA_ID
    
    mircounts <- mircounts[,2:ncol(mircounts)]
    
    keep <- grep("read_count_",colnames(mircounts))
    
    mircounts <- mircounts[,keep]
    
    colnames(mircounts) <- gsub("read_count_","",colnames(mircounts))
    
    write.table(mircounts,
                file = paste(project,"mircounts",sep = "_"),
                quote = TRUE, sep = "\t")
  }else{
    mircounts <- fread(paste(project,"mircounts",sep = "_"))
    head(colnames(mircounts))
    rownames <- mircounts$V1
    mircounts <- mircounts[,2:ncol(mircounts)]
    rownames(mircounts) <- rownames
  }
  
  
  assign("mircountmatrix", mircounts, envir = .GlobalEnv)
  
  return("Function complete, the countmatrices have been assigned to your global enviroment")
}

    
