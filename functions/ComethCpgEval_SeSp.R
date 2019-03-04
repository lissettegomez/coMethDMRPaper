### input: 
# allcpgs - vector of all cpgs in the region
# truecpgs - vector of true cpgs 
# predcpgs - vector of predicted cpgs from cometh

#### output: 
# Sens = Probability (predicted positive | actual positive)
# Specificity =   Pr (predicted negative | actual negative)



ComethCpgEval_SeSp <- function (allcpgs, truecpgs, predcpgs){
  
  cpgs <- data.frame(allcpgs)
  colnames (cpgs) <- "cpg"
  
  cpgs$actual <- ifelse (cpgs$cpg %in% truecpgs, 1, 0)
  
  cpgs$pred   <- ifelse (cpgs$cpg %in% predcpgs, 1, 0)
  
  
  cpgs$sensitivity <- sum(cpgs$actual & cpgs$pred) / sum(cpgs$actual)
  
  cpgs$specificity <- sum( (cpgs$actual == 0) & (cpgs$pred == 0)) / sum(cpgs$actual == 0)
    
    return (as.data.frame(cpgs))
  
}

# allcpgs <- allCpgs_lsvec[[10]]
# 
# predcpgs <- predCpgs_lsvec[[10]]
# 
# truecpgs <- trueCpgs_lsvec [[10]]

