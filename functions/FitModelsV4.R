## linear model with median as outcome

## dat = data frame, rowsnames = sample ids, col = cpgs + pheno info (include variable age)

fit_lm_med <- function (dat, contPheno_char){

  tryCatch({

  cpg.index <- grep("cg", colnames(dat))

  cpgs_df <- dat[, cpg.index]

  med <- apply(cpgs_df, 1, median)

  #test <- rowMeans(dat[, 1:4])

  modFormula <- paste("med ~", contPheno_char )

  mod <- lm (as.formula (modFormula), data = dat)

  pval_lm_median <- coef(summary(mod)) [2,4]

  #cor(means, p)

  print (summary(mod))

  return (pval_lm_median)

  }, error=function(e){})

}

###### linear model with mean as outcome

# dat = data frame, rowsnames = sample ids, col = cpgs + pheno info (include variable age)

fit_lm_mean <- function (dat, contPheno_char){

  tryCatch({

  cpg.index <- grep("cg", colnames(dat))

  cpgs_df <- dat[, cpg.index]

  means <- apply(cpgs_df, 1, mean)

  modFormula <- paste("means ~", contPheno_char )

  mod <- lm (as.formula (modFormula), data = dat)

  pval_lm_mean <- coef(summary(mod)) [2,4]

  print (summary(mod))

  return (pval_lm_mean)

  }, error = function(e){})

}

## get tall data for fitting mixed model & gee
# dat = data frame, row.names = sample ids, variables = cpgs + pheno type (include variable age)
# output = data frame, tall data, with variables sample, age, cpg, mvalue

get_tall_data <- function (dat, contPheno_char){

  require (reshape2)

  wide <- dat

  wide$sample <- row.names(wide)

  cpg.index <- grep("cg", colnames(wide))

  tall <- melt (wide,
                id.vars = c("sample", contPheno_char),
                measure.vars = colnames (wide)[cpg.index],
                value.name = "mvalue",
                variable.name = "cpg")

  tall <- tall[order(tall$sample, tall$cpg) ,]

  tall$sample <- factor(tall$sample)

  tall$cpg <- as.character(tall$cpg)

  return (tall)
}

## simple linear mixed model

fit_slmm <- function (talldat, contPheno_char){

  require (lmerTest)

  tryCatch({

    modFormula <- paste("mvalue ~", contPheno_char, "+ (1|sample)")

    f <- lmer(as.formula (modFormula), talldat)

    stat_slmm <- coef(summary(f)) [2 , 4]

    pval_slmm <- 2 * ( 1- pnorm (abs(stat_slmm)))

    print(summary(f))

    return (list (pval_slmm, stat_slmm))

  }, error=function(e){})
}

#fit_slmm (talldat = tall_simdata)


## gee

fit_gee <- function (talldat, contPheno_char){

  tryCatch({

  require(gee)
  set.seed (08232018)

  modFormula <- paste("mvalue ~", contPheno_char)

  m <- gee(formula = as.formula(modFormula), id = sample, data = talldat,
           family = gaussian, corstr="exchangeable", maxiter = 1000, silent = TRUE)
  print(summary(m))

  est <- data.frame(m$coefficients)

  se <- data.frame(sqrt(diag(m$robust.variance)))

  wald.z <- est / se

  p.value <- 2*( 1- pnorm(abs(wald.z[,1])))

  #names(p.value) <- paste0("gee_", row.names(wald.z))

  print (p.value)

  #return(p.value[2])

  #return(wald.z[2,1])

  pval_gee <- p.value[2]

  stat_gee <- wald.z[2,1]

  return (list(pval_gee, stat_gee))

  }, error=function(e){})

}

### random coefficient mixed model

fit_rclmm <- function (talldat, contPheno_char){

  require (lmerTest)

  tryCatch({

    modFormula <- paste0("mvalue ~ ", contPheno_char, " + (", contPheno_char, "|cpg) + (1|sample)")

    f <- lmer(as.formula(modFormula) , talldat)

    print(summary(f))

    stat_rclmm <- coef(summary(f)) [2, 4]

    pval_rclmm <- 2 * ( 1- pnorm (abs(stat_rclmm)))

    return (list(pval_rclmm, stat_rclmm))

  }, error=function(e){})

}

#fit_rclmm (talldat = tall_simdata)


