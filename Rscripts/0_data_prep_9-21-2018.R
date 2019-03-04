
source ("../FUNCTIONS/getCorr.R")

###### CGI regions in beta values (dataframe)  -------------------------------------
## rowname = cpg id
## colname = sample ids

pfc   <- readRDS ("pfc.closeCpGs.regions.final.RDS")

betas <- readRDS ("PFCbetas.bmiq.RDS") 

betaCgi_df <- lapply(pfc, function (item)
  betas [ which(rownames(betas)%in%item),])  

########### CGI regions in M values (dataframe)---------------------------------

# get Mvalues
mvalCgi_lsdf <- lapply (betaCgi_df,
                        function (item) log2((item)/(1-item)))


######### transposed Mvalue cgi regions ----------------------------------------------
tmvalCgi_lsdf <- lapply (mvalCgi_lsdf,
                         function (item)  data.frame(t(item)))


######### cgi region stats (dataframe) ---------------------------------------
# varaibles = region, min_corr, nProbes, nProbe_c, 
#             corr_range (range of pairwise correlations)

# get correlations
minCorrs <- lapply (tmvalCgi_lsdf, GetMinCorr) 

minCorrs_df <- data.frame(do.call ( rbind, minCorrs))
colnames (minCorrs_df)[1] <- "min_Corr"

minCorrs_df$region <- row.names (minCorrs_df)

# get number of cpgs
nProbes <- lapply(tmvalCgi_lsdf, ncol)

nProbes_df <- data.frame (do.call (rbind, nProbes))

colnames(nProbes_df) [1] <- "nProbes"

nProbes_df$region <- row.names(nProbes_df)

### merge 
regionStats <- merge (minCorrs_df, nProbes_df, by = "region")

# binned nprobes = 8,9,10 into 1 category
regionStats$nProbes_c  <- ifelse (regionStats$nProbes < 8, regionStats$nProbes, 8 ) 

# those with corr >= 0.5
regionStats_5 <- regionStats[regionStats$min_Corr >=0.5 ,] 
regionStats_5$corr_range <- ifelse (regionStats_5$min_Corr < 0.8, "0.5-0.8", "0.8-1.0")

table (regionStats_5$corr_range, regionStats_5$nProbes_c)

table (regionStats_5$corr_range, regionStats_5$nProbes)

saveRDS (regionStats_5, "regionStats_5.RDS")








