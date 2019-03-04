# add pseudo age - mean of mvalues for cpgs

AddPhenoOrderedByMean <- function (methyl_mat){
  
  age <- sort( rpois (n = 27, lambda = 65) )
  
  tall <- data.frame(methyl_mat)
  
  tall$means <- rowMeans (tall)
  
  tall <- tall[order(tall$means) ,]
  
  tall2 <- cbind (tall, age)
  
  return (tall2)
}