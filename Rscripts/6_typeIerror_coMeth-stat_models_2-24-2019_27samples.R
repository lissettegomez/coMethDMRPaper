
#### type I error analysis
# take cometh regions
# generate pseudo age - randomly, 1000 times
# fit each model

source ("./functions/FitModelsV4.R")
source ("./functions/getCorr.R")


### 1. generate pseudo age randomly

set.seed (915)

age_mat <- matrix (0, nrow = 27, ncol = 1000)

for (i in 1:1000){

  age_mat [, i] <- rpois (n = 27, lambda = 65)

}

age_df <- as.data.frame (age_mat)

saveRDS (age_df, "age27_df.RDS")


### 2. fit each model

betas_df <- readRDS ("./data/PFCbetas.bmiq.RDS")

coMethRegions_ls  <- readRDS ("coMeth_ls.RDS")$CpGsSubregions

coMethRegions_lsdf <- lapply (coMethRegions_ls[1:10], function(x)  betas_df [rownames(betas_df) %in% unlist(x),] )

saveRDS (coMethRegions_lsdf, "coMethRegions_lsdf.RDS")



age_df <- readRDS ("age27_df.RDS")


sink("all_models_fitting_details_coMeth_27samples.txt")

res_all <- data.frame (matrix (nrow = 0, ncol = 10))

nRegions <- 10

nReps <- 1000

for (r in 1 : nRegions){

  for (j in 1:nReps) {

      set.seed (8*r)

      print ("################################################################################### ")
      print (paste0 ("r = ", r, ", region = ", names(coMethRegions_lsdf[r]), ", repetition = ", j ))

      oneMethyl_df <- t(coMethRegions_lsdf[[r]])

      age <- age_df[,j]

      one_sim <- data.frame (oneMethyl_df, age)

      ####### get corr
      med_corr <- GetMedCorr (one_sim)

      print (med_corr)

      min_corr <- GetMinCorr (one_sim)

      print (min_corr)

      max_corr <- GetMaxCorr(one_sim)

      print (max_corr)


      ####### fit model

      tall <- get_tall_data (dat = one_sim, "age")

      pval_lm_mean <- fit_lm_mean (dat = one_sim, "age")

      pval_lm_median <- fit_lm_med (dat = one_sim, "age")


      res_slmm <- fit_slmm (talldat = tall, "age")
      pval_slmm <- res_slmm [[1]]
      stat_slmm <- res_slmm [[2]]

      res_gee  <- fit_gee (talldat = tall, "age")
      pval_gee <- res_gee [[1]]
      stat_gee <- res_gee [[2]]

      res_rclmm <- fit_rclmm (talldat = tall, "age")
      pval_rclmm <- res_rclmm[[1]]
      stat_rclmm <- res_rclmm[[2]]

      one <- data.frame (r, j, pval_lm_mean, pval_lm_median,
                         pval_gee, stat_gee,
                         pval_slmm, stat_slmm,
                         pval_rclmm, stat_rclmm, med_corr, min_corr, max_corr)

     res_all <- rbind (res_all, one)

    }
}

saveRDS (res_all, "res_27samples_coMeth.RDS")

sink()

# type I error rate
res_all <- readRDS ("res_27samples_coMeth.RDS")


testSize_27samples <- data.frame (matrix (nrow = 0, ncol = 5))


for (i in 1:10){

  oneRep <- res_all[res_all$r == i ,]

  lm_mean   <- nrow(oneRep[oneRep$pval_lm_mean < 0.05 ,])/nrow(oneRep)
  lm_median <- nrow(oneRep[oneRep$pval_lm_median < 0.05 ,])/nrow(oneRep)
  gee       <- nrow(oneRep[oneRep$pval_gee < 0.05 ,])/nrow(oneRep)

  slmm      <- nrow(oneRep[oneRep$pval_slmm < 0.05 ,])/nrow(oneRep)
  rclmm     <- nrow(oneRep[oneRep$pval_rclmm < 0.05 ,])/nrow(oneRep)

  typeIerror <- data.frame (lm_mean, lm_median, gee, slmm, rclmm)

  testSize_27samples <- rbind (testSize_27samples, typeIerror)

}


write.csv (testSize_27samples, "testSize_27samples_coMeth.csv")
