
library(ChIPpeakAnno)
library(coMethDMR)

# methods <-      c("randCoef", "IMA_mean", "IMA_median", "GEE")
# methodsLabel <- c("coMethDMR_randCoef", "IMA_mean", "IMA_median", "Aclust_GEE")
#
# methods      <- c("simple", "IMA_mean", "IMA_median", "GEE")
# methodsLabel <- c("coMethDMR_simple", "IMA_mean", "IMA_median", "Aclust_GEE")

methods      <- c("randCoef", "simple", "combp")
methodsLabel <- c("coMethDMR_randCoef", "coMethDMR_simple", "combp")


n <- length (methods)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols = gg_color_hue(n)

pdf ( paste0("./RESULTS_otherMethods/overlap DMRs - coMethDMR_combp - PFC data.pdf"))

ranges.list <- list()


for (d in 1: n){

  temp <- read.csv ( paste0("./RESULTS_otherMethods/sig_",methods[d], ".csv"))

  temp$method <- methodsLabel[d]

  ranges.list[[d]] <- RegionsToRanges(temp$regionName)

}

p <- makeVennDiagram(ranges.list,
                     NameOfPeaks = methodsLabel [1:n],
                     totalTest = 1000, by="region",
                     main = paste0("Venn Diagram for PFC data"),
                     fill=cols[1:n],
                     cat.pos = c(-20, 0, 20),
                     cat.dist=c(0.05, 0.05, 0.02)
                     )
print (p)

dev.off()
