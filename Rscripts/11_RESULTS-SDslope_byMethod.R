

## Create var method and keep var probePvals in datasets
getDat <- function(method, label) {
  # Call in data set
  dat <- read.csv(paste0("./data/probePvals_", method, ".csv"))
  
  ## Get nCpGs for method seqlm data
  if (method == "seqlm"){
    dat$nCpGs <- dat$length
  }
  
  # Keep regions with at least 3 CpGs
  dat <- dat[dat$nCpGs >= 3, ]
  
  # Compute variance of slopeEstimate by regionName
  slopeEstimateSd <- aggregate(
    slopeEstimate ~ regionName,
    data = dat,
    sd
  )
  
  names(slopeEstimateSd)[2] <- "slopeEstimateSd"
  
  # Create variable method
  slopeEstimateSd$method <- paste0(label)
  
  # order by significance level
  
  top10RegionsOrdered <- unique( factor(dat$regionName)) [1:10]
  
  SdOrdered <- slopeEstimateSd [top10RegionsOrdered ,]
  
  # Return new data set
  return(SdOrdered )
}

IMA_mean <- getDat(method = "IMA_mean", label = "IMA_mean")
IMA_median <- getDat(method = "IMA_median", label = "IMA_median")

randCoef <- getDat(method = "randCoef", label = "coMethDMR_randCoef")
simple   <- getDat(method = "simple", label = "coMethDMR_simple")

seqlm <- getDat(method = "seqlm", label = "seqlm")
GEE   <- getDat(method = "GEE", label = "Aclust_GEE")

## Combine datasets for plotting
Alldat <- rbind(IMA_mean, IMA_median, randCoef, simple, seqlm, GEE)

write.csv (Alldat, "allMethods_std_slopeEstimates.csv")

## Box plot of slope p-values across methods
library(ggplot2)

pdf("slopeEstimateSd by methods.pdf")

p <- ggplot(Alldat, aes(x = method, y = slopeEstimateSd, fill = method)) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(width = 0.2) +
  geom_jitter(shape=16, position=position_jitter(0.01)) + 
  #geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5) +
  ggtitle("Top 10 regions") +
  xlab("Method") +
  ylab("Standard Deviations of Slope Estimates") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)

dev.off()