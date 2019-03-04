# plot power for different methods

library(tidyr)

AllRes <- readRDS ("AllRes.RDS")

power <- subset (AllRes, select = c(power_slmm, power_rclmm,
                                    mincorr, ncpgs, rep, fold))

AllResCoMeth <- readRDS ("AllRes_coMeth.RDS")

powerCoMeth <- subset (AllResCoMeth, select = c(power_slmm, power_rclmm,  mincorr, ncpgs, rep, fold))

colnames (powerCoMeth) [1:2] <- c("power_coMeth_slmm", "power_coMeth_rclmm")


both <- merge (power, powerCoMeth, by = c("mincorr", "ncpgs", "rep", "fold"))

# reshape data

tall <- gather(data = both, key = method, value = power, -c(mincorr, ncpgs, rep, fold))

tall$method <- sub("power_", "", tall$method)

av <- aggregate ( power ~ method + mincorr + ncpgs + fold, data = tall, mean)

av$fold    <- factor(paste0("fold = ", av$fold))
av$ncpgs   <- factor(av$ncpgs)
av$mincorr <- factor(paste0("minCorr = ", av$mincorr))

str(av)

av$method <- factor(av$method, levels = c("slmm", "rclmm", "coMeth_slmm", "coMeth_rclmm"))

av <- av [order(av$fold, av$ncpgs, av$mincorr, av$method) ,]

av$method <- as.character(av$method)

library(ggplot2)

pdf ("power_both.pdf")

ggplot (data = av, aes (x = ncpgs, y = power, fill = method)) +
 geom_bar(stat="identity", position= position_dodge()) + theme_bw() +
 facet_grid(fold ~ mincorr )

dev.off ()




