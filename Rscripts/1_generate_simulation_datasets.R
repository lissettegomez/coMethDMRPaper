
source ("./functions/GetCorr.R")
source ("./functions/AddPhenoOrderedByMean.R")
source ("./functions/SelectRegions.R")
source ("./functions/AddRandomProbes.R")
source ("./functions/OrderCpGsByLocation.R")

######################################################################
####  generate simulation datsets
# - take CGI regions with ncpgs, and min_corr > value
# - add phenotye (age) = ranked same way as mean of cpgs
#   - output simulation datasets
# - add random probes
#   - output simualtion datasets
####################################################################

# datasets
tmvalCgi_lsdf <- readRDS ("./data/tmvalCgi_lsdf.RDS" )

regionStats   <- readRDS ("./data/regionStats_5.RDS")

betas <- readRDS ("./data/PFCbetas.bmiq.RDS")

cpgLocations <- readRDS ("./data/cpg.locations.RDS")


#test <- selectRegions (regionInfo = regionStats, minCorr = 0.8, ncpgs = 5)

minCorr_vec <- c(0.5, 0.8)
ncpgs_vec <- c(3, 5, 8)
nrep <- 10

for (c in 1:2){

  for (n in 1:length(ncpgs_vec)){

    # 1. Select regions with minCorr & ncpgs
    # regions with minCorr = 0.5 to 0.8, and 0.8 - 1

    c1 <- minCorr_vec[c]
    c2 <- ifelse (c1 == 0.5, 0.8, 1)

    topRegions_vec <- SelectRegions (regionInfo = regionStats,
                                     minCorr = c1,
                                     minCorr_ceiling = c2,
                                     ncpgs = ncpgs_vec[n])

    topRegions_lsdf <- tmvalCgi_lsdf[topRegions_vec]

    # order cpgs by location
    topRegionsOrd_lsdf <- lapply (topRegions_lsdf, OrderCpGs, locationFile = cpgLocations)

      for (rep in 1:nrep){

      # 2. add pseudo-age ordered in the same way as means of cpgs
      set.seed (c*n*rep)
      TPRegionPheno_lsdf <- lapply (topRegionsOrd_lsdf, AddPhenoOrderedByMean)

      # these are true positive probes
      corr <- minCorr_vec[c]*10

      saveRDS (TPRegionPheno_lsdf, paste0("./data/DATA_simulated/RegionPheno_C",
                                        corr, "_ncpg", ncpgs_vec[n], "_rep", rep, ".RDS"))

      # 3. add random probes
          for (fold in 1:2){
            nRProbes <- ncpgs_vec[n]*fold

            RandomRegionPheno_lsdf <-  lapply (TPRegionPheno_lsdf, AddRandomProbes,
                                         beta_all = betas, nRandomProbes = nRProbes)

            # datasets that ends with "fold" have random probes added
            saveRDS (RandomRegionPheno_lsdf, paste0("./data/DATA_simulated/RegionPheno_C",
                                               corr, "_ncpg", ncpgs_vec[n], "_rep", rep, "_fold",fold, ".RDS")  )
          }
     }
  }
}
