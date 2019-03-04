
## Figure showing randCoef model improves specificity

library(coMethDMR)

### 1. Read in datasets

pfc_df <- readRDS ("./data/PFCbetas.bmiq.RDS")

pheno_df <- readRDS ("./data/pfcPheno_df.RDS")

pheno_df$Mplate <- factor(pheno_df$Mplate)

str(pheno_df)

### 2. get data read for computing regionwise p-values

CpGs <- c("cg04977002", "cg06756193", "cg11179286", "cg12513911", "cg20154264")

betas_mat <- t(pfc_df [CpGs ,])

mvalues_df <- data.frame (log2 (betas_mat / (1 - betas_mat)))

mvaluesPheno_df <- merge (pheno_df, mvalues_df, by.x = "Sample", by.y = "row.names")

mvaluesPheno_df$Sample <- factor (mvaluesPheno_df$Sample)


### 4. regionwise p-value by each method

source ("./functions/FitModelsV2.R")

print ("**** 1. linear model with mean mvalues as outcome ***")
fit_lm_mean (dat = mvaluesPheno_df, contPheno_char = "stage")

print ("**** 2. linear model with median mvalues as outcome ***")
fit_lm_med (dat = mvaluesPheno_df, contPheno_char = "stage")

mvaluesTall_df <- get_tall_data(dat = mvaluesPheno_df, contPheno_char = "stage")

print ("**** 3. GEE ****")
fit_gee (talldat = mvaluesTall_df, contPheno_char = "stage")

print ("**** 4. simple linear mixed model ***")
fit_slmm (talldat = mvaluesTall_df, contPheno_char = "stage")

print ("**** 5. random coeff mixed model ***")
fit_rclmm (talldat = mvaluesTall_df, contPheno_char = "stage")

### 5. figure for line plot
source ("./functions/PlotLinesWithDots.R")

CpGsToTest_char <- CpGs

PlotLinesWithDots (CpGsToTest_char, mvalues_df = pfc_df, pheno_df)


