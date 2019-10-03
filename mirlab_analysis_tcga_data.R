mirlab_tcga_analysis <- function(mrna_matrix,microrna_matrix,cancer = 1,mrna_keep = 2000, microrna_keep = 100){
  #Packages needed
  require("miRLAB")
  require("edgeR")
  #Get common patients
  common_patients <- intersect(substring(colnames(mrna_matrix),1,15),substring(colnames(microrna_matrix),1,15))
  #Keep cancer or normals only, by default keeps only samples of cancer patients
  if(cancer == 1){
    common_patients <- common_patients[substring(common_patients,14,15) != "11"]
  }
  else{
    common_patients <- common_patients[substring(common_patients,14,15) == "11"]
  }
  #Subset matrices
  common_patients_mrna <- mrna_matrix[,common_patients]
  common_patients_microrna <- mrna_matrix[,common_patients]
  #Perform TMM normalization
  #Turn into DGEList
  common_patients_mrna <- DGEList(common_patients_mrna)
  common_patients_microrna <- DGEList(common_patients_microrna)
  #Normalize
  common_patients_mrna <- calcNormFactors(common_patients_mrna, method = "TMM")
  common_patients_microrna <- calcNormFactors(common_patients_microrna, method = "TMM")
  #In the interest of time, use only the 100 most highly expressed miRs
  #and the 2000 most highly expressed mRNA
  #Compute row averages
  average_expr_mrna <- rowSums(common_patients_mrna)/ncol(common_patients_mrna)
  average_expr_microrna <- rowSums(common_patients_microrna)/ncol(common_patients_microrna)
  #Get their order
  order_mrna <- order(average_expr_mrna)
  order_microrna <- order(average_expr_microrna)
  #Get only the top expressed ones, up to K and L (set in function definition)
  if(mrna_keep < nrow(common_patients_mrna)){
    common_patients_mrna <- common_patients_mrna[order_mrna[1:mrna_keep],]
  }
  if(microrna_keep < nrow(common_patients_microrna)){
    common_patients_microrna <- common_patients_microrna[order_microrna[1:microrna_keep],]
  }
  #Set cause and effect
  cause <- 1:nrow(common_patients_microrna)
  effect <- cause+1:(cause+1 + nrow(common_patients_mrna))
  #Create cause and effect matrix
  dataset <- rbind(common_patients_mrna,common_patients_microrna)
  #predict miRNA targets using Mutual Information
  mi=MI(dataset, cause, effect)
  
  #predict miRNA targets using causal inference
  ida=IDA(dataset, cause, effect, "stable", 0.01)
  
  #predict miRNA targets using linear regression
  lasso=Lasso(dataset, cause, effect)
  
  return(list(mi,ida,lasso))
}
  
    
