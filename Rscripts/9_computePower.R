
source ("./functions/ComputePower.R")

minCorr_vec <- c(0.5, 0.8)
ncpgs_vec <- c(3, 5, 8)
nrep <- 10

AllRes <- data.frame (matrix(nrow = 0, ncol = 16))

for (c in 1:2){

  for (n in 1:3){

    c1 <- minCorr_vec[c]
    c2 <- ifelse (c1 == 0.5, 0.8, 1)

    for (fold in 1:2){



      for (rep in 1:nrep){

        corr <- c1*10

        ########## all cpgs in the regions

        oneRep <- readRDS (paste0("results/allCpgs_RegionPheno_C", corr,
                              "_ncpg", ncpgs_vec[n],
                              "_rep", rep, "_fold",fold, ".RDS"))

        res1 <- data.frame(ComputePower (file = oneRep, method = "lm_mean"))
        res2 <- data.frame(ComputePower (file = oneRep, method = "lm_median"))
        res3 <- data.frame(ComputePower (file = oneRep, method = "slmm"))
        res4 <- data.frame(ComputePower (file = oneRep, method = "gee"))
        res5 <- data.frame(ComputePower (file = oneRep, method = "rclmm"))

        resOneRep <- cbind (res1, res2, res3, res4, res5)

        resOneRep$mincorr <- minCorr_vec[c]
        resOneRep$ncpgs <- ncpgs_vec[n]
        resOneRep$fold <- fold
        resOneRep$rep <- rep

        AllRes <- rbind (AllRes, resOneRep)
      }


    }
  }
}

saveRDS (AllRes, "AllRes_coMeth.RDS")

write.csv (AllRes, "AllRes_CoMeth.csv", row.names = FALSE)

