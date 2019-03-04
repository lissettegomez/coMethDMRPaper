
library(pROC)

#### 1. datasets - put all datasets in a summary file ---------------------------

minCorr_vec <- c(0.5, 0.8)
ncpgs_vec <- c(3, 5, 8)
nrep <- 10
r.drop <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9)

all <- data.frame (matrix(nrow = 0, ncol = 7))

for (c in 1:2){

  for (n in 1:3){

    for (fold in 1:2){

      for (r in 1:9) {

        allReps <- data.frame (matrix(nrow = 0, ncol = 3))

        for (rep in 1:nrep){

                c1 <- minCorr_vec[c]
                c2 <- ifelse (c1 == 0.5, 0.8, 1)

                rdrop <- 10*r.drop[r]
                corr <- c1*10

                one <- read.csv (paste0("results/resultsCoMeth_C", corr, "_ncpg", ncpgs_vec[n],
                                              "_rep", rep, "_fold",fold, "_rdrop_", rdrop,".csv"))

                one$pred <- factor(one$pred)

                # assign auc to be 0 if only one class is predicted

                if (nlevels((one$pred)) >1) {
                  one_auc <- auc(one$pred, one$actual)
                } else one_auc = 0


                # one value per region
                oneRep <- unique (subset (one, select = -c(cpg, actual, pred)))

                allReps <- rbind (allReps, oneRep)
        }

         avSen <- mean (allReps$sensitivity)

         avSp  <- mean (allReps$specificity)

         minCorr <- minCorr_vec[c]
         ncpgs   <- ncpgs_vec[n]
         rdrop   <- r.drop[r]
         fold    <- fold

         oneResult <- cbind (minCorr, ncpgs, rdrop, fold, avSen, avSp, one_auc)

         all <- rbind (oneResult, all)
         }
      }
    }
}

vars <- c("minCorr", "ncpgs", "fold", "rdrop")

all [vars] <- lapply (all[vars], factor)


all$grp <- as.factor(paste0("minCorr", all$minCorr, "_ncpgs", all$ncpgs, "_fold", all$fold))

saveRDS (all, "all_results.RDS")

################################# plots --------------------------------------------------------

all <- readRDS ("all_results.RDS")

all <- subset (all, select = -grp)

library(lattice)
library(latticeExtra)

all <- all[order(all$ncpgs, all$minCorr, all$fold) ,]

pdf ("XYplots_sen_sp_auc.pdf", width = 10, height = 7)

ThreeXyPlots <- function (fold) {

  temp <- all[all$fold == fold,]

  p <- doubleYScale(
    xyplot (avSp ~ rdrop | factor (ncpgs,   labels = c("ncpgs = 3", "ncpgs =5", "ncpgs = 8")) +
              factor (minCorr, labels = c("minCorr = 0.5", "minCorr = 0.8")),
            ylab = "specificity",
            data = temp, main = paste0("fold = ", fold)) ,
    xyplot (avSen ~ rdrop | factor (ncpgs,   labels = c("ncpgs = 3", "ncpgs =5", "ncpgs = 8")) +
              factor (minCorr, labels = c("minCorr = 0.5", "minCorr = 0.8")),
            ylab = "sensitivity",
            data = temp, main = paste0("fold = ", fold)),
    add.ylab2 = TRUE)
  
  print(p)
  
  u <- xyplot (one_auc ~ rdrop | factor (ncpgs,   labels = c("ncpgs = 3", "ncpgs =5", "ncpgs = 8")) +
                 factor (minCorr, labels = c("minCorr = 0.5", "minCorr = 0.8")),
               ylab = "AUC",
               data = temp, main = paste0("fold = ", fold))
  
  print (u)
  
}

ThreeXyPlots (fold = 1)

ThreeXyPlots (fold = 2)

dev.off()

###################### output value of rdrop that achived max ----------------------------------------

max_rdrop <- aggregate(avSen ~ minCorr + ncpgs + fold, data = all, max)
maxSen <- merge (max_rdrop, all, by = c("minCorr", "ncpgs", "fold", "avSen"))
write.csv (maxSen, "rdrop_for_maxSen.csv", row.names = FALSE)

max_rdrop <- aggregate(avSp ~ minCorr + ncpgs + fold, data = all, max)
maxSp <- merge (max_rdrop, all, by = c("minCorr", "ncpgs", "fold", "avSp"))
write.csv (maxSp, "rdrop_for_maxSp.csv", row.names = FALSE)

max_rdrop <- aggregate(one_auc ~ minCorr + ncpgs + fold, data = all, max)
maxAuc <- merge (max_rdrop, all, by = c("minCorr", "ncpgs", "fold", "one_auc"))
write.csv (maxAuc, "rdrop_for_auc.csv", row.names = FALSE)

write.csv (all, "rdrop_for_all.csv", row.names = FALSE)



