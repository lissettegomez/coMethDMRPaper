### 1. Import dasen normalized methylation data

library(data.table)
dasen <- fread('GSE43414_dasen_geo_all_cohorts.csv', header = T, sep = ',')
dasen_df <- data.frame (dasen)

### 2. Import phenotype data for pfc samples

library(readxl)
pheno <- read_excel( "./data/pheno_PFC.xlsx", col_names = FALSE)

pheno$barcode <- sub("barcode: ", "", pheno$...12)

pheno$Xbarcode <- paste0("X", pheno$barcode)

pheno$stage <- sub("braak.stage: ", "", pheno$...14)

pheno <- pheno[pheno$stage != "Exclude" ,]

pheno$stage <- as.numeric (pheno$stage)

pheno$Sample <- pheno$...1

pheno$subject.id <- pheno$...11

pheno$Mplate <- as.factor(substr(pheno$...12, 9,19))

pheno$sex <- as.factor (pheno$...15)

pheno$age.brain <- as.numeric(substr(pheno$...17, 11, 14))
  
pheno_df <- pheno[, c("barcode", "Xbarcode", "stage", "Sample", "subject.id", "Mplate", "sex", "age.brain")]

# saveRDS (pheno_df, "pheno_df.RDS")
# 
# write.csv (pheno_df, "pheno_df.csv", row.names = FALSE)

#### 3. subset columns in dasen dataset that are pfc samples

row.names(dasen_df) <- dasen_df$V1 

pfc_df <- dasen_df[, pheno_df$Xbarcode]

identical (colnames(pfc_df), pheno_df$Xbarcode)

colnames(pfc_df) <- pheno_df$Sample

saveRDS (pfc_df, "pfc_df.RDS")

#### 4. Estimate Cell type proportions

load("~/CETS_Image.RData");

## get reference profile from Caucasions + controls
idx <- list(controlNeuron = pdBrain$celltype == "N" & pdBrain$diag == "Control" & pdBrain$ethnicity == "Caucasian", 
            controlGlia   = pdBrain$celltype == "G" & pdBrain$diag == "Control" & pdBrain$ethnicity == "Caucasian")

refProfile <- getReference(brain, idx)
head(refProfile)

## estimate proportions of neurons in PFC samples

# limit to 10,000 cpgs in the refProfile dataset
pfc <- pfc_df

selected <- rownames(pfc) %in% rownames(refProfile)
table(selected)

pfc.refcpgs <- pfc[selected ,]

# estimate proportion of neurons
prop <- data.frame (estProportion(pfc.refcpgs, profile = refProfile))

colnames(prop) <- "prop.neuron"

saveRDS (prop, "pfc.prop.neurons.RDS")

write.csv (prop, "pfc.prop.neurons.csv")

##### 5. Make final phenotype file
pfcPheno_df <- merge ( pheno, prop, by.x="Sample", by.y="row.names")

# factors these to have correct factor levels
pfcPheno_df$Mplate <- factor(pfcPheno_df$Mplate)

pfcPheno_df$Sample <- factor(pfcPheno_df$Sample)

pfcPheno_df$subject.id <- as.factor (pfcPheno_df$subject.id)

str(pfcPheno_df)

saveRDS (pfcPheno_df, "pfcPheno_df.RDS")


