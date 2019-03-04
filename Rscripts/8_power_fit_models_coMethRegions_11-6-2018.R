
###################################
## fit models to coMethylated cpgs
##################################

results.dir <- "./results_coMeth/"


source ("../FUNCTIONS/FitModels.R")
source ("../FUNCTIONS/getCorr.R")


source ("../FUNCTIONS/CpGsInRegionV2.R")

library(coMethDMR)
library(bumphunter)

locationFile <- readRDS (paste0(data.dir.ad, "cpg.locations.RDS"))

minCorr_vec <- c(0.5, 0.8)
ncpgs_vec <- c(3, 5, 8) 
nrep <- 10

for (c in 1:2){
  
  for (n in 1:3){
    
    c1 <- minCorr_vec[c]
    c2 <- ifelse (c1 == 0.5, 0.8, 1)
    
    for (fold in 1:2){
      
      for (rep in 1:nrep){
      
        corr <- c1*10
        
        ########## all cpgs in the regions
        
        datasetName <- paste0("RegionPheno_C", corr, 
                             "_ncpg", ncpgs_vec[n], 
                             "_rep", rep, "_fold",fold, ".RDS")
        
        temp <- readRDS ( paste0("../DATA_simulated/", datasetName))
        
        # take out means and age
        temp2 <- lapply (temp, function(item)
          t(subset (item, select = -c(means, age))) )
        
        markRegions_lsdf <- lapply (temp2,  MarkComethylatedCpGs, 
                                    rDropThresh_num = 0.4)
        
        contRegions_lsdf <- lapply (markRegions_lsdf, FindComethylatedRegions, minCpGs_int = 3)
        
        predCpgs_lsvec <- lapply(contRegions_lsdf, CpgsInRegion)
        
        datContRegions_lsdf <- mapply (function (x, y){x[, c(y, "age", "means")]}, temp, predCpgs_lsvec, SIMPLIFY = FALSE )
        
        res_all <- data.frame (matrix(nrow = 0, ncol = 11))
        
        for (r in 1: length(datContRegions_lsdf)) {
        
        #for (r in 1: 2) {
        
          
          print ("################################################################################### ")
          print (datasetName)
          print (paste0("region ", r, ": ", names(datContRegions_lsdf[r]) ))
          
          one_sim <- datContRegions_lsdf[[r]]
          
          med_corr <- GetMedCorr (one_sim)
          
          print (med_corr)
          
          min_corr <- GetMinCorr (one_sim)
          
          print (min_corr)
          
          max_corr <- GetMaxCorr(one_sim)
          
          print (max_corr)
          
          tall <- get_tall_data (dat = one_sim)
          
          pval_lm_mean <- fit_lm_mean (dat = one_sim)
          
          pval_lm_median <- fit_lm_med (dat = one_sim)
          
          pval_slmm <- fit_slmm (talldat = tall)
          
          pval_gee  <- fit_gee (talldat = tall)
          
          pval_rclmm <- fit_rclmm (talldat = tall)
          
          one <- data.frame (sub(".RDS", "", datasetName), r, names(datContRegions_lsdf[r]), 
                             pval_lm_mean, pval_lm_median, pval_slmm, pval_gee, pval_rclmm, med_corr, min_corr, max_corr)
          
          res_all <- rbind (res_all, one)
          
        }
          colnames (res_all)[1:3] <- c("dataset", "region", "location")

          saveRDS (res_all, paste0(results.dir, "allCpgs_", datasetName))

          write.csv (res_all, paste0(results.dir, "allCpgs_", sub(".RDS", ".csv", datasetName)), row.names = FALSE)
      }
    }
  }
}


