#!/opt/conda/envs/R_env/bin/Rscript --vanilla
# Sanja Priselac
# 26/04/2021

#########################################################
library(dplyr)
#########################################################

###########################################################################################################
tab.original <- read.csv(snakemake@input[[1]],  sep='')
tab.qn <- read.csv(snakemake@input[[2]])
tab.dp05 <- read.csv(snakemake@input[[3]])
tab.dp1 <- read.csv(snakemake@input[[4]])
tab.dp25 <- read.csv(snakemake@input[[5]])
tab.dp5 <- read.csv(snakemake@input[[6]])
tab.dp10 <- read.csv(snakemake@input[[7]])
tab.nn <- read.csv(snakemake@input[[8]])


tab.original$results <- rep(0, nrow(tab.original))
tab.original[which(tab.original$Log2.Fold_Change.>1 & tab.original$FDR < 0.05), 'results'] <- 1
tab.original[which(tab.original$Log2.Fold_Change.< (-1) & tab.original$FDR < 0.05), 'results'] <- (-1)
tab.original$Genes <- tab.original$Gene
tab.original$logFC <- tab.original$Log2.Fold_Change.

###########################################################################################################


test.orders <- function(tab.main, tab.sub, n, res = 1) {
  if (res > 0) {
    tab.main <- tab.main[which(tab.main$results == 1), ]
    tab.main <- (tab.main[order(tab.main$FDR), ])[1:n, c("Genes", "FDR")]
  } else {
    tab.main <- tab.main[which(tab.main$results == (-1)), ]
    tab.main <- (tab.main[order(tab.main$FDR), ] )[1:n, c("Genes", "FDR")]
  }
  tab.test <- merge(tab.main, tab.sub[, c("Genes", "FDR")], by="Genes", all.x=TRUE)
  t <- friedman.test(as.matrix(tab.test[, -1]))
  num <- length(tab.sub[which(tab.sub$Genes %in% tab.test$Genes), "results"] == res)
  tab.test$FriedmanP <- c(t$p.value, rep(NA, nrow(tab.test)-1))
  tab.test$NumExprsses <- c(num, rep(NA, nrow(tab.test)-1))
  return(tab.test)
}


tab1up <- test.orders(tab.original, tab.qn, n=10)
tab1up <- rename(tab1up, FDR.original = FDR.x, FDR.qn =FDR.y)
write.csv(tab1up, file=snakemake@output[[1]], row.names = FALSE)

tab1down <- test.orders(tab.original, tab.qn, n=10, res=(-1))
tab1down <- rename(tab1down, FDR.original = FDR.x, FDR.qn =FDR.y)
write.csv(tab1down, file=snakemake@output[[8]], row.names = FALSE)

tab2up <- test.orders(tab.qn, tab.dp05, n=10)
tab2up <- rename(tab2up, FDR.qn = FDR.x, FDR.dp05 =FDR.y)
write.csv(tab2up, file=snakemake@output[[2]], row.names = FALSE)

tab2down <- test.orders(tab.qn, tab.dp05, n=10, res=(-1))
tab2down <- rename(tab2down, FDR.qn = FDR.x, FDR.dp05 =FDR.y)
write.csv(tab2down, file=snakemake@output[[9]], row.names = FALSE)

tab3up <- test.orders(tab.qn, tab.dp1, n=10)
tab3up <- rename(tab3up, FDR.qn = FDR.x, FDR.dp1 =FDR.y)
write.csv(tab3up, file=snakemake@output[[3]], row.names = FALSE)

tab3down <- test.orders(tab.qn, tab.dp1, n=10, res=(-1))
tab3down <- rename(tab3down, FDR.qn = FDR.x, FDR.dp1 =FDR.y)
write.csv(tab3down, file=snakemake@output[[10]], row.names = FALSE)

tab4up <- test.orders(tab.qn, tab.dp25, n=10)
tab4up <- rename(tab4up, FDR.qn = FDR.x, FDR.dp25 =FDR.y)
write.csv(tab4up, file=snakemake@output[[4]], row.names = FALSE)

tab4down <- test.orders(tab.qn, tab.dp25, n=10, res=-1)
tab4down <- rename(tab4down, FDR.qn = FDR.x, FDR.dp25 =FDR.y)
write.csv(tab4down, file=snakemake@output[[11]], row.names = FALSE)

tab5up <- test.orders(tab.qn, tab.dp5, n=10)
tab5up <- rename(tab5up, FDR.qn = FDR.x, FDR.dp5 =FDR.y)
write.csv(tab5up, file=snakemake@output[[5]], row.names = FALSE)

tab5down <- test.orders(tab.qn, tab.dp5, n=10, res=-1)
tab5down <- rename(tab5down, FDR.qn = FDR.x, FDR.dp5 =FDR.y)
write.csv(tab5down, file=snakemake@output[[12]], row.names = FALSE)

tab6up <- test.orders(tab.qn, tab.dp10, n=10)
tab6up <- rename(tab6up, FDR.qn = FDR.x, FDR.dp10 =FDR.y)
write.csv(tab6up, file=snakemake@output[[6]], row.names = FALSE)

tab6down <- test.orders(tab.qn, tab.dp10, n=10, res=-1)
tab6down <- rename(tab6down, FDR.qn = FDR.x, FDR.dp10 =FDR.y)
write.csv(tab6down, file=snakemake@output[[13]], row.names = FALSE)

tab7up <- test.orders(tab.original, tab.nn, n=10, res=1)
tab7up <- rename(tab7up, FDR.original = FDR.x, FDR.nn =FDR.y)
write.csv(tab7up, file=snakemake@output[[7]], row.names = FALSE)

tab7down <- test.orders(tab.original, tab.nn, n=10, res=-1)
tab7down <- rename(tab7down, FDR.original = FDR.x, FDR.nn =FDR.y)
write.csv(tab7down, file=snakemake@output[[14]], row.names = FALSE)
