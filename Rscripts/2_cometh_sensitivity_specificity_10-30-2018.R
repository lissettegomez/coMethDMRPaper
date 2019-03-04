
source ("./functions/contiguous_regions_noCpGsortByLocation.R")
source ("./functions/ComethCpgEval_SeSp.R")
source ("./functions/CpGsInRegion.R")

library(psych)
library(bumphunter)
library(coMethDMR)

# compute sensitivity and specificity for cometh selected cpg probes

locationFile <- readRDS ("./data/cpg.locations.RDS")

minCorr_vec <- c(0.5, 0.8)
ncpgs_vec <- c(3, 5, 8) 
nrep <- 10
r.drop <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9)

results_all <- data.frame (matrix(nrow = 0, ncol = 9))

for (c in 1:2){
  
  for (n in 1:3){
    
    c1 <- minCorr_vec[c]
    c2 <- ifelse (c1 == 0.5, 0.8, 1)
    
    for (fold in 1:2){
      
      for (rep in 1:nrep){
      
        ########## true positive cpgs in regions
        corr <- c1*10
        
        tp_lsdf  <- readRDS (paste0("./data/DATA_simulated/RegionPheno_C", 
                                    corr, "_ncpg", ncpgs_vec[n], "_rep", rep, ".RDS") )  
        
        tp2_lsdf <-  lapply (tp_lsdf, function(item)
                                        t(subset (item, select = -c(means, age))) ) 
          
        trueCpgs_lsvec <- lapply (tp2_lsdf, function(item) rownames (item))
        
        
        ########## all cpgs in the regions
        temp <- readRDS (paste0("./data/DATA_simulated/RegionPheno_C", 
                       corr, "_ncpg", ncpgs_vec[n], "_rep", rep, "_fold",fold, ".RDS"))
        
        # take out means and age
        temp2 <- lapply (temp, function(item)
                                 t(subset (item, select = -c(means, age))) )
        
        
        allCpgs_lsvec <- lapply (temp2, function (item) rownames (item) )
        
        ######### cpgs in coMeth regions
        
        for (r in 1:9) {
          
          rdrop <- 10*r.drop[r]
          
          markRegions_lsdf <- lapply (temp2,  MarkComethylatedCpGs, 
                                            rDropThresh_num = r.drop[r])
          
          contRegions_lsdf <- lapply (markRegions_lsdf, FindComethylatedRegions, minCpGs_int = 3)
          
          # test <- betaMatrix_ex4
          # marks_df <- MarkComethylatedCpGs(t(test), rDropThresh_num = 0.6)
          # FindComethylatedRegions(marks_df)
          
          
          if(length(contRegions_lsdf) > 0) {
              
              predCpgs_lsvec <- lapply(contRegions_lsdf, CpgsInRegion)
            
              ########## compute precision and recall for each simulation dataset
              
              result_lsdf <- mapply (ComethCpgEval_SeSp, allCpgs_lsvec, 
                                     trueCpgs_lsvec, 
                                     predCpgs_lsvec, SIMPLIFY = FALSE)
              
              results_df <- do.call (rbind, result_lsdf)
              
              ##### add region names
              strings <- strsplit(row.names(results_df), '\\.')
              
              first <- do.call (rbind,  lapply (strings, function (x) x[[1]]))
              
              results_df$region <- first
              
              ######### output 
              
              write.csv (results_df, paste0("results/resultsCoMeth_C", corr, "_ncpg", ncpgs_vec[n], 
                                            "_rep", rep, "_fold",fold, "_rdrop_", rdrop,".csv"), row.names = FALSE )
              
              results2_df <- unique(subset(results_df, select = c(region, sensitivity, specificity)))
              
              results2_df$c1 <- c1
              results2_df$c2 <- c2
              results2_df$ncpgs <- ncpgs_vec[n]
              results2_df$rep   <- rep
              results2_df$fold  <- fold
              results2_df$rdrop <- rdrop
              
              results_all <- rbind (results_all, results2_df)
          }
          
          
        }  
      }
    }
  }
}
