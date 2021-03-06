library(coMethDMR)
library(DMRcate)

###########################################################

#  analysis of PFC data form Lunnon et al. (2014)
#  - using linear mixed model

###########################################################

## 1. import dasen normalized data

pfc_df <- readRDS ("./data/pfc_df.RDS")

pfcPheno_df <- readRDS ("./data/pfcPheno_df.RDS")


## 2. load closeby regions that are in CGIs -> compute co-methylated regions

library(parallel)

clust <- makeCluster(15)
clusterEvalQ(cl = clust, library(coMethDMR))

clusterExport(cl = clust, varlist = c("result.dir" ) )


pfc_cgi_rdrop0_4_ls <- CoMethAllRegions (
              betaMatrix = pfc_df,
              regionType = "ISLAND",
              arrayType = "450k",
              rDropThresh_num = 0.4,
              returnAllCpGs = FALSE
        )

saveRDS (pfc_cgi_rdrop0_4_ls, "pfc_cgi_rdrop0_4_ls.RDS")

## 3. Test one region

pfc_cgi_ls <- readRDS ("pfc_cgi_rdrop0_4_ls.RDS")

oneRegion_df <- pfc_df [pfc_cgi_ls[[1]] ,]

identical(as.character(pfcPheno_df$Sample), colnames(pfc_df))

lmmTest (betaOne_df = oneRegion_df,
         pheno_df = pfcPheno_df,
         contPheno_char = "stage",
         covariates_char = c("age.brain", "sex", "Mplate", "prop.neuron"),
         modelType = "randCoef",
         arrayType = "450k")


## 4. test all regions

### need to make sure Mplate and sex are factors

str(pfcPheno_df)

pfc_cgi_rdrop0_4_ls <- readRDS ("pfc_cgi_rdrop0_4_ls.RDS"))

res_cgi_df1 <- lmmTestAllRegions(
                    beta_df = pfc_df,
                    region_ls = pfc_cgi_rdrop0_4_ls,
                    pheno_df = pfcPheno_df,
                    contPheno_char = "stage",
                    covariates_char = c("age.brain", "sex", "Mplate", "prop.neuron"),
                    modelType = "randCoef",
                    arrayType = "450k"
                  )


write.csv (res_cgi_df1, "res_cgi_randCoef.csv", row.names = FALSE)

res_cgi_df2 <- lmmTestAllRegions(
  beta_df = pfc_df,
  region_ls = pfc_cgi_rdrop0_4_ls,
  pheno_df = pfcPheno_df,
  contPheno_char = "stage",
  covariates_char = c("age.brain", "sex", "Mplate", "prop.neuron"),
  modelType = "simple",
  arrayType = "450k"
)

write.csv (res_cgi_df2, "res_cgi_simple.csv", row.names = FALSE)


gc()
stopCluster(clust)

