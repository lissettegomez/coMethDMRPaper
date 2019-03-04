
# next 
# change function addrandomprobes - change distance of random probes so they are close




#### add random probes to a region
# cpgsPheno_df = data frame, row = sample ids, column = cpg ids
# beta_all = beta values, genomewide
# nRandomProbes = number of random probes to be added

AddRandomProbes <- function (cpgsPheno_df, beta_all, nRandomProbes){
  
  totProbes <- nrow (beta_all)
  
  random_n <- sample (1:totProbes, nRandomProbes, replace = FALSE)
  
  tRandomProbes <- data.frame(t(betas[random_n ,]))
  
  randomMvalues <- log2((tRandomProbes)/(1-tRandomProbes))
  
  final <- cbind (cpgsPheno_df, randomMvalues)
  
  return (final)
}