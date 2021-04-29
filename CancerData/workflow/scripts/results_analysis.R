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


test.orders <- function(tab.main, tab.sub, lf = 1) {
  if (lf > 0) {
    tab.main <- tab.main[which(tab.main$logFC > lf & tab.main$adj.P.Val < 0.05), ]
    tab.sub <- tab.sub[which(tab.sub$logFC>lf & tab.sub$adj.P.Val < 0.05), ]
  } else {
    tab.main <- tab.main[which(tab.main$logFC < lf & tab.main$adj.P.Val < 0.05), ]
    tab.sub <- tab.sub[which(tab.sub$logFC < lf & tab.sub$adj.P.Val < 0.05), ]
  }
  tab.main <- tab.main[order(tab.main$adj.P.Val), c('ID', 'adj.P.Val')]
  #deparse(substitute(tab.main))
  tab.main <- rename(tab.main, adjPval.main = adj.P.Val)
  tab.test <- merge(tab.main, tab.sub[,  c('ID', 'adj.P.Val')], by='ID', all=TRUE)
  tab.test <- rename(tab.test, adjPval.sub = adj.P.Val)
  t <-  wilcox.test(tab.test[, 2], tab.test[, 3], paired=TRUE)
  tab.test$WicoxonP <- c(t$p.value, rep(NA, nrow(tab.test)-1))
  return(tab.test)
}


tab1up <- test.orders(tab.original, tab.qn)
tab1up <- rename(tab1up, adjP.original = adjPval.main, adjP.qn =adjPval.sub)
write.csv(tab1up, file=snakemake@output[[1]], row.names = FALSE)

tab1down <- test.orders(tab.original, tab.qn, lf=-1)
tab1down <- rename(tab1down, adjP.original = adjPval.main, adjP.qn =adjPval.sub)
write.csv(tab1down, file=snakemake@output[[7]], row.names = FALSE)

tab2up <- test.orders(tab.qn, tab.dp05)
tab2up <- rename(tab2up, adjP.qn = adjPval.main, adjP.dp05 =adjPval.sub)
write.csv(tab2up, file=snakemake@output[[2]], row.names = FALSE)

tab2down <- test.orders(tab.qn, tab.dp05, lf=-1)
tab2down <- rename(tab2down, adjP.qn = adjPval.main, adjP.dp05 =adjPval.sub)
write.csv(tab2down, file=snakemake@output[[8]], row.names = FALSE)

tab3up <- test.orders(tab.qn, tab.dp1)
tab3up <- rename(tab3up, adjP.qn = adjPval.main, adjP.dp1 =adjPval.sub)
write.csv(tab3up, file=snakemake@output[[3]], row.names = FALSE)

tab3down <- test.orders(tab.qn, tab.dp1, lf=-1)
tab3down <- rename(tab3down, adjP.qn = adjPval.main, adjP.dp1 =adjPval.sub)
write.csv(tab3down, file=snakemake@output[[9]], row.names = FALSE)

tab4up <- test.orders(tab.qn, tab.dp25)
tab4up <- rename(tab4up, adjP.qn = adjPval.main, adjP.dp25 =adjPval.sub)
write.csv(tab4up, file=snakemake@output[[4]], row.names = FALSE)

tab4down <- test.orders(tab.qn, tab.dp25, lf=-1)
tab4down <- rename(tab4down, adjP.qn = adjPval.main, adjP.dp25 =adjPval.sub)
write.csv(tab4down, file=snakemake@output[[10]], row.names = FALSE)

tab5up <- test.orders(tab.qn, tab.dp5)
tab5up <- rename(tab5up, adjP.qn = adjPval.main, adjP.dp5 =adjPval.sub)
write.csv(tab5up, file=snakemake@output[[5]], row.names = FALSE)

tab5down <- test.orders(tab.qn, tab.dp5, lf=-1)
tab5down <- rename(tab5down, adjP.qn = adjPval.main, adjP.dp5 =adjPval.sub)
write.csv(tab5down, file=snakemake@output[[11]], row.names = FALSE)

tab6up <- test.orders(tab.qn, tab.dp10)
tab6up <- rename(tab6up, adjP.qn = adjPval.main, adjP.dp10 =adjPval.sub)
write.csv(tab6up, file=snakemake@output[[6]], row.names = FALSE)

tab6down <- test.orders(tab.qn, tab.dp10, lf=-1)
tab6down <- rename(tab6down, adjP.qn = adjPval.main, adjP.dp10 =adjPval.sub)
write.csv(tab6down, file=snakemake@output[[12]], row.names = FALSE)
