### compute av. and median corr of the regions

# methyl_mat = a dataframe with row = sample id, col name = cpg ids

GetMedCorr <- function (methyl_mat) {
  
  
  cpg.index <- grep("cg", colnames(methyl_mat))
  
  cpgs_df <- methyl_mat[, cpg.index]
  
  corr <- cor(cpgs_df, method = "spearman")
  
  print (corr)
  
  corr_up <- corr [upper.tri (corr)]
  
  #corr_mean <- mean (corr_up)
  
  corr_median <- median (corr_up)
  
  return (corr_median)
}

GetMinCorr <- function (methyl_mat) {
  
  cpg.index <- grep("cg", colnames(methyl_mat))
  
  cpgs_df <- methyl_mat[, cpg.index]
  
  corr <- cor(cpgs_df, method = "spearman")
  
  corr_up <- corr [upper.tri (corr)]
  
  #corr_mean <- mean (corr_up)
  
  corr_min <- min (corr_up)
  
  return (corr_min)
}

GetMaxCorr <- function (methyl_mat) {
  
  cpg.index <- grep("cg", colnames(methyl_mat))
  
  cpgs_df <- methyl_mat[, cpg.index]
  
  corr <- cor(cpgs_df, method = "spearman")
  
  corr_up <- corr [upper.tri (corr)]
  
  #corr_mean <- mean (corr_up)
  
  corr_max <- max (corr_up)
  
  return (corr_max)
}




