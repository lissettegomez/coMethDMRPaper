######## orders CpGs for input into contiguous.region function

#### input: 
# cpg_df = df with rows = samples ids, columns = cpg ids

#### output: 
# cpgOrd_df = df with cpgs ordered by location

OrderCpGs <- function (cpg_df, locationFile){
  
  cluster.location <- locationFile [which (locationFile$cpg %in% colnames(cpg_df)),]
  
  cluster.location.ordered <- cluster.location[order(cluster.location$CHR,
                                                     as.numeric(as.character(cluster.location$MAPINFO))),]
  
  cpgOrdered_df <- cpg_df [, as.character(cluster.location.ordered$cpg)]
  
  return (cpgOrdered_df)

}


## test: cpg_df <- topRegions_lsdf[[10]]