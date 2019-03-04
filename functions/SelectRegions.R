
# input: 
# regionInfo = data frame of region, min_Corr, nProbes
# minCorr = min pairwise corr between the cpgs
# ncpgs = nubmer of cpgs in the region

# output: 
# character vector of region names

SelectRegions <- function (regionInfo, minCorr, minCorr_ceiling, ncpgs){
  
  regionsOrdered <- regionStats[order(regionStats$min_Corr, decreasing = TRUE) ,]
  
  regionsOrdered <- regionsOrdered[regionsOrdered$nProbes_c == ncpgs ,]
  
  select <- (regionsOrdered$min_Corr >= minCorr) & (regionsOrdered$min_Corr < minCorr_ceiling)
  
  top <- regionsOrdered[select, ]
  
  return (top$region)
}

#test <- SelectRegions (regionInfo = regionStats, minCorr = 0.8, ncpgs = 5)
