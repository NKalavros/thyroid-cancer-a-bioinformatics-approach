obtain_data_for_mrna_and_mir <- function(project){
  #Load packages
  require(TCGAbiolinks)
  require(data.table)
  #Begin with mRNA
  
  if(!file.exists(paste(project,"mrnacounts.tsv",sep = "_"))){
    

    projectmrnaquery <- GDCquery(project = project,
                                 data.category = "Transcriptome Profiling",
                                 data.type = "Gene Expression Quantification",
                                 workflow.type = "HTSeq - Counts")
    
    #Download the data
    GDCdownload(projectmrnaquery,
                method = "client")
    
    #Prepare the data
    mrnacounts <- as.data.frame(GDCprepare(projectmrnaquery, save = FALSE, summarizedExperiment = FALSE))
    
    #Edit the mRNA sample matrix as needed
    mrnacounts <- mrnacounts[6:nrow(mrnacounts),]
    
    rownames <- mrnacounts[,1]
    
    mrnacounts <- mrnacounts[,2:ncol(mrnacounts)]
    
    rownames(mrnacounts) <- rownames
    
    #Write the new matrices in the working directory
    write.table(mrnacounts,
                file = paste(project,"mrnacounts.tsv",sep = "_"),
                quote = TRUE, sep = "\t", row.names = TRUE, col.names = TRUE)
  }else{
    mrnacounts <- fread(paste(project,"mrnacounts.tsv",sep = "_"))
    rownames <- as.character(mrnacounts[,1])
    mrnacounts <- mrnacounts[,2:ncol(mrnacounts)]
    rownames(mrnacounts) <- rownames
  }
  
  assign("mrnacountmatrix", mrnacounts, envir = .GlobalEnv)
  
  #Same for miRs
  
  if(!file.exists(paste(project,"mircounts.tsv",sep = "_"))){
    projectmirquery <- GDCquery(project = project,
                                data.category = "Transcriptome Profiling",
                                data.type = "miRNA Expression Quantification",
                                workflow.type = "BCGSC miRNA Profiling")
    
    GDCdownload(projectmirquery,
                method = 'client')
    
    mircounts <- as.data.frame(GDCprepare(projectmirquery, save = FALSE, summarizedExperiment = FALSE))
    
    #Edit the miR sample matrix as needed
    rownames <- mircounts$miRNA_ID
    
    mircounts <- mircounts[,2:ncol(mircounts)]
    
    keep <- grep("read_count_",colnames(mircounts))
    
    mircounts <- mircounts[,keep]
    
    colnames(mircounts) <- gsub("read_count_","",colnames(mircounts))
    
    rownames(mircounts) <- rownames
    
    write.table(mircounts,
                file = paste(project,"mircounts.tsv",sep = "_"),
                quote = TRUE, sep = "\t", row.names = TRUE, col.names = TRUE)
  }else{
    mircounts <- fread(paste(project,"mircounts.tsv",sep = "_"))
    rownames <- as.character(mircounts[,1])
    mircounts <- mircounts[,2:ncol(mircounts)]
    rownames(mircounts) <- rownames
  }
  
  
  assign("mircountmatrix", mircounts, envir = .GlobalEnv)
  
  return("Function complete, the countmatrices have been assigned to your global enviroment")
}
