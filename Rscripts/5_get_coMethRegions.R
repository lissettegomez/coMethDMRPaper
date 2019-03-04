
library(coMethDMR)

pfc_df <- readRDS ("./data/PFCbetas.bmiq.RDS")

cgi_ls <- system.file ("extdata",
                       "ISLAND3_200.RDS",
                 package = 'coMethDMR',
                 mustWork = TRUE)

coMeth_ls <- CoMethAllRegions (
 betaMatrix = pfc_df,
 file = cgi_ls,
 fileType = "RDS",
 arrayType = "450k",
 returnAllCpGs = FALSE
)

saveRDS (coMeth_ls, "coMeth_ls.RDS")
