mirlab_tcga_analysis <- function(mrna_matrix,microrna_matrix,cancer = 1,mrna_keep = 2000, microrna_keep = 100,top_targets = 20){
  #Packages needed
  require("miRLAB")
  require("edgeR")
  #mrna_matrix <- mrnacountmatrix
  #microrna_matrix <- mircountmatrix
  #cancer <- 1
  #mrna_keep <- 2000
  #microrna_keep <- 100
  #Get common patients
  #Keep cancer or normals only, by default keeps only samples of cancer patients
  #Subset matrices
  if(cancer == 1){
    mrna_matrix <- mrna_matrix[,substring(colnames(mrna_matrix),14,15) != "11"]
    microrna_matrix <- microrna_matrix[,substring(colnames(microrna_matrix),14,15) != "11"]
  }else{
    microrna_matrix <- microrna_matrix[,substring(microrna_matrix,14,15) == "11"]
    mrna_matrix <- mrna_matrix[,substring(mrna_matrix,14,15) == "11"]
  }
  mir_in_mrna <- match(substring(colnames(mrna_matrix),1,15),substring(colnames(microrna_matrix),1,15))
  mir_in_mrna <- mir_in_mrna[!is.na(mir_in_mrna)]
  microrna_matrix <- microrna_matrix[,mir_in_mrna]
  mrna_in_mir <- match(substring(colnames(microrna_matrix),1,15),substring(colnames(mrna_matrix),1,15))
  mrna_in_mir <- mrna_in_mir[!is.na(mrna_in_mir)]
  mrna_matrix <- mrna_matrix[,mrna_in_mir]
  if(!all(substring(colnames(mrna_matrix),1,15) == substring(colnames(microrna_matrix),1,15))){
    print("Could not match all patients")
    return(1)
  }
  #Perform TMM normalization
  #Turn into DGEList
  mrna_matrix <- DGEList(mrna_matrix)
  microrna_matrix <- DGEList(microrna_matrix)
  #Normalize, obtain counts
  mrna_matrix <- calcNormFactors(mrna_matrix, method = "TMM")
  mrna_matrix <- mrna_matrix$counts
  microrna_matrix <- calcNormFactors(microrna_matrix, method = "TMM")
  microrna_matrix <- microrna_matrix$counts
  print("Normalization for both datasets complete")
  #In the interest of time, use only the 100 most highly expressed miRs
  #and the 2000 most highly expressed mRNA
  #Compute row averages
  average_expr_mrna <- rowSums(mrna_matrix)/ncol(mrna_matrix)
  average_expr_microrna <- rowSums(microrna_matrix)/ncol(microrna_matrix)
  #Get their order
  order_mrna <- order(average_expr_mrna,decreasing = TRUE)
  order_microrna <- order(average_expr_microrna, decreasing = TRUE)
  #Get only the top expressed ones, up to K and L (set in function definition)
  if(mrna_keep < nrow(mrna_matrix)){
    mrna_matrix <- mrna_matrix[order_mrna[1:mrna_keep],]
  }
  if(microrna_keep < nrow(microrna_matrix)){
    microrna_matrix <- microrna_matrix[order_microrna[1:microrna_keep],]
  }
  print(paste("Completed subsetting of matrices. Final dimensions for mRNA are:",paste(dim(mrna_matrix), collapse = "x"),"and for miR are:",paste(dim(microrna_matrix),collapse="x")))
  #Set cause and effect
  cause <- 1:nrow(microrna_matrix)
  effect <- (cause[length(cause)]+1):(cause[length(cause)] + nrow(mrna_matrix))
  #Create cause and effect matrix
  dataset <- t(rbind(microrna_matrix,mrna_matrix))
  write.csv(dataset,"dataset.csv",row.names = FALSE,quote = FALSE)
  #predict miRNA targets using Mutual Information
  mi = MI("dataset.csv", cause, effect)
  print("Mutual information analysis complete")
  #predict miRNA targets using causal inference
  ida=IDA("dataset.csv", cause, effect, "stable", 0.01)
  print("IDA analysis complete")
  #predict miRNA targets using linear regression
  lasso=Lasso(dataset, cause, effect)
  print("Lasso analysis completed, merging results")
  #Extract the top targets per miR
  results <- BordaTopk(list(mi,ida,lasso),top_targets)
  
  return(list(mi,ida,lasso))
}
