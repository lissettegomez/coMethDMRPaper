
# regionName_char - region in the format of chrx:xxxx-xxxx
# mvalues_df - data frame of (adjusted) mvalues, with rownames = cpgs, column = sample IDs
# pheno_df - data frame of phenotype values, with variable Sample and stage

PlotLinesWithDots <- function (regionName_char, mvalues_df, pheno_df, title_char){

  ### Extract individual CpGs in the region ###
  CpGsToTest_char <- GetCpGsInRegion(regionName_char, arrayType = "450k")

  ### Transpose betaMatrix from wide to long ###
  CpGsMvalue_df <- as.data.frame (
    t (
      mvalues_df [
        which(rownames(mvalues_df) %in% CpGsToTest_char),
        ]
    )
  )


  CpGsMvalue_df$Sample <- row.names(CpGsMvalue_df)

  #### reshaping data for plotting
  CpGMvalueTall_df <- melt (
    CpGsMvalue_df,
    id = "Sample",
    variable_name = "ProbeID"
  )

  mvaluePheno_df <- merge (
    pheno_df,
    CpGMvalueTall_df,
    by = "Sample"
  )

  ### center by mean
  cpgMeansByStage <- aggregate (
    value ~ stage + ProbeID,
    data = mvaluePheno_df,
    mean
  )

  cpgMeansByStage$beta <- 2^(cpgMeansByStage$value)/ (2^cpgMeansByStage$value + 1)

  # cpgMeans <- aggregate (
  #   value ~ ProbeID,
  #   data = cpgMeansByStage,
  #   mean
  # )
  #
  # colnames(cpgMeans)[2] <- "mean"
  #
  # cpgMeansCenter <- merge (cpgMeansByStage, cpgMeans, by="ProbeID")
  #
  # cpgMeansCenter$centered_value <- cpgMeansCenter$value - cpgMeansCenter$mean

  # p <- ggplot(cpgMeansCenter, aes(x = stage, y = centered_value, color = ProbeID) ) +
  #   geom_point(size=3, shape=21) +
  #   geom_line() +
  #   ggtitle(paste0("region = ", regionName_char)) +
  #   theme_bw() +
  #   theme(legend.position = c(0.8, 0.7)) +
  #   scale_y_continuous( limits = c(min(means.cpg.center$value.centered),
  #                                  max(means.cpg.center$value.centered)))

  p <- ggplot(cpgMeansByStage, aes(x = stage, y = value, color = ProbeID) ) +
    geom_point ( size=3, shape=21) +
    geom_smooth(method = lm, se = FALSE) +
    theme_bw() +
    ggtitle( title_char ) +
    xlab ("Braak stage") +
    ylab ("Corrected mean methylation M values")

  #  scale_y_continuous(limits = c(-0.5, 0.5))


  print (p)


  # ggplot(test, aes(x = stage, y = value) ) +
  #    geom_point ( size=3, shape=21) +
  #    geom_smooth(method = lm, se = FALSE) +
  #    theme_bw() +
  #   scale_y_continuous(limits = c(-0.5, 0.5))
  #
  #
  #    ggtitle(paste0("region = ", regionName_char))


}
