
############################################
# funciton to compute power, given p-values 
############################################

ComputePower <- function (file, method){
  

   rawp <-   file [, paste0("pval_", method)]

   ind <- ifelse ( rawp < 0.05, 1, 0)  
  
   nSig <- sum(ind == 1)  
   nTotal <- nrow (file)
   power <-  nSig/nTotal   
  
   oneStat <- cbind (nSig, nTotal, power)
   
   colnames (oneStat) [1] <- paste0("nSig_",   method)
   colnames (oneStat) [2] <- paste0("nTotal_", method)
   colnames (oneStat) [3] <- paste0("power_",  method)
   
 return (oneStat)  
}

