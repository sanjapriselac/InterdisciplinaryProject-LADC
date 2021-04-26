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

tab.original$adj.P.Val <- tab.original$pValue
tab.original$ID <- tab.original$Gene
tab.original$logFC <- tab.original$Log2.Fold_Change.

###########################################################################################################

## tab main is the table where the first n p values are taken
test.orders <- function(tab.main, tab.sub, lf = 1) {
  if (lf > 0) {
    tab.main <- tab.main[which(tab.main$logFC > lf & tab.main$adj.P.Val < 0.05), ]
    tab.sub <- tab.sub[which(tab.sub$logFC>lf & tab.sub$adj.P.Val < 0.05), ]
  } else {
    tab.main <- tab.main[which(tab.main$logFC < lf & tab.main$adj.P.Val < 0.05), ]
    tab.sub <- tab.sub[which(tab.sub$logFC < lf & tab.sub$adj.P.Val < 0.05), ]
  }
  tab.main <- tab.main[order(tab.main$adj.P.Val), c('ID', 'adj.P.Val')]
  deparse(substitute(tab.main))
  tab.main <- rename(tab.main, adjPval.main = adj.P.Val)
  tab.test <- merge(tab.main, tab.sub[,  c('ID', 'adj.P.Val')], by='ID', all=TRUE)
  tab.test <- rename(tab.test, adjPval.sub = adj.P.Val)
  t <-  wilcox.test(tab.test[, 2], tab.test[, 3], paired=TRUE)
  tab.test$WicoxonP <- c(t$p.value, rep(NA, nrow(tab.test)-1))
  return(tab.test)
}


tab1 <- test.orders(tab.original, tab.qn)
tab1 <- rename(tab1, adjP.original = adjPval.main, adjP.qn =adjPval.sub)
write.csv(tab1, file=snakemake@output[[1]])

tab2 <- test.orders(tab.qn, tab.dp05)
tab2 <- rename(tab2, adjP.qn = adjPval.main, adjP.dp05 =adjPval.sub)
write.csv(tab2, file=snakemake@output[[2]])

tab3 <- test.orders(tab.qn, tab.dp1)
tab3 <- rename(tab3, adjP.qn = adjPval.main, adjP.dp1 =adjPval.sub)
write.csv(tab3, file=snakemake@output[[3]])

tab4 <- test.orders(tab.qn, tab.dp25)
tab4 <- rename(tab4, adjP.qn = adjPval.main, adjP.dp25 =adjPval.sub)
write.csv(tab4, file=snakemake@output[[4]])

tab5 <- test.orders(tab.qn, tab.dp5)
tab5 <- rename(tab5, adjP.qn = adjPval.main, adjP.dp5 =adjPval.sub)
write.csv(tab5, file=snakemake@output[[5]])

tab6 <- test.orders(tab.qn, tab.dp10)
tab6 <- rename(tab6, adjP.qn = adjPval.main, adjP.dp10 =adjPval.sub)
write.csv(tab6, file=snakemake@output[[6]])
