
# 1. line plots for top regions with at least 3 CpGs

library(coMethDMR)
library(ggplot2)
library(reshape)


source ("../simulations/FUNCTIONS/PlotLinesWithDots.R")

resid_df <- readRDS ("resid_df.RDS")

pdf ("./RESULTS_otherMethods/linePlotsWithDots_allMethods.pdf")

for (m in 1:7){

  res <- readRDS (paste0("./RESULTS_otherMethods/top_", methods[m], ".RDS"))

  switch(methods[m],
         "randCoef" = {res$pvalMethod <- res$pValue},
         "simple" = {res$pvalMethod <- res$pValue},
         "IMA_median" = {res$pvalMethod <- res$P.Value},
         "IMA_mean" = {res$pvalMethod <- res$P.Value},
         "seqlm" = {res$pvalMethod <- res$p.value},
         "GEE" = {res$pvalMethod <- res$pValue},
         "combp" = {res$pvalmethod <- res$z_p}

  )

  #pdf (paste0("./RESULTS_otherMethods/linePlotsWithDots_", methods[m], ".pdf"))

  for (i in 1:nrow(res)){

    if (methods[m] != "combp"){
       title <- paste0 ( "method = ", methodsLabel[m],  "\n(", i, ") ",
                    res$regionName[i],
                    ", pVal = ", formatC(res$pval[i], format = "e", digits = 2),
                     ", fdr = ", round(res$fdr[i],4))
    } else {

      title <- paste0 ( "method = ", methodsLabel[m],  "\n(", i, ") ",
                        res$regionName[i],
                        ", pVal = ", formatC(res$pval[i], format = "e", digits = 2),
                        ", sidak_p = ", round(res$z_sidak_p[i],4))
    }

  PlotLinesWithDots (
      regionName_char = as.character(res$regionName[i]),
      mvalues_df = resid_df,
      pheno_df = pheno_df,
      title_char = title
    )


  }



}

dev.off()


