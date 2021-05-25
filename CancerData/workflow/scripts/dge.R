#!/opt/conda/envs/R_env/bin/Rscript --vanilla
# Sanja Priselac
# 11/05/2021

##################
library(ggplot2)
library(limma)
library(edgeR)
library(VennDiagram)
library(svglite)
library(scran)
##################

for(i in 1:7) {
  load(snakemake@input[[i]])
}

##################################################################################################################

DEGtable <- function(DEGlist) {
  res <- findMarkers(DEGlist$counts, DEGlist$samples$group, test.type = 'wilcox', add.summary=TRUE,  pval.type = "all", full.stats=TRUE)
  FDR <- exp(res@listData$NonNormal@listData$stats.Normal@listData$log.FDR)
  FDR <-  data.frame(Genes = res@listData$NonNormal@rownames, FDR = FDR)
  pValue <- exp(res@listData$NonNormal@listData$stats.Normal@listData$log.p.value)
  pValue <-  data.frame(Genes = res@listData$NonNormal@rownames, pValue = pValue)
  adj.pValue <- data.frame(Genes = res@listData$NonNormal@rownames, adjPvalue = p.adjust(pValue$pValue, "BH"))
  logFC <- apply(DEGlist$counts[, DEGlist$samples$group== "NonNormal"], 1, mean) - apply(DEGlist$counts[, DEGlist$samples$group== "Normal"], 1, mean) 
  logFC <- data.frame(Genes = attributes(logFC)$names, logFC = logFC)
  resultTable <- merge(FDR, pValue)
  resultTable <- merge(resultTable, adj.pValue)
  resultTable <- merge(resultTable, logFC)
  results <- rep(0, nrow(resultTable))
  results[which(resultTable$logFC>1 & resultTable$FDR < 0.05)] <- 1
  results[which(resultTable$logFC<(-1) & resultTable$FDR < 0.05)] <- -1
  resultTable$results <- as.factor(results)
  return(resultTable)
}

volcano_plot <- function(resultTable) {
  p <- ggplot(resultTable, aes(x=logFC, y=-log10(FDR), col=results)) + 
    geom_point(size = 1.9) +
    #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
    theme_minimal() +
    xlab(paste("Log2(Fold Change)")) +
    ylab("-log10(FDR)") +
    theme(text = element_text(size=20)) + 
    scale_colour_manual(values = c('#33FF00', '#000000', '#FF3300'))
    #geom_text(aes(label=ifelse(results!=0,as.character(Genes),'')),hjust=0, vjust=0)
  return(p)
}


##################################################################################################################
resultTableNN <- DEGtable(log2.y)
resultTable <- DEGtable(log2.y.qn)
resultTable05 <- DEGtable(y.qn.05)
resultTable1 <- DEGtable(y.qn.1)
resultTable25 <- DEGtable(y.qn.25)
resultTable5 <- DEGtable(y.qn.5)
resultTable10 <- DEGtable(y.qn.10)

vpNN <- volcano_plot(resultTableNN)
ggsave("VolcanoPlotNN.jpeg", plot=vpNN, path = "/tmp/repo/DGE_LADC/results/plots")
vp <- volcano_plot(resultTable)
ggsave("VolcanoPlot.jpeg", plot=vp, path = "/tmp/repo/DGE_LADC/results/plots")
vp05 <- volcano_plot(resultTable05)
ggsave("VolcanoPlot05.jpeg", plot=vp05, path = "/tmp/repo/DGE_LADC/results/plots")
vp1 <- volcano_plot(resultTable1)
ggsave("VolcanoPlot1.jpeg", plot=vp1, path = "/tmp/repo/DGE_LADC/results/plots")
vp25 <- volcano_plot(resultTable25)
ggsave("VolcanoPlot25.jpeg", plot=vp25, path = "/tmp/repo/DGE_LADC/results/plots")
vp5 <- volcano_plot(resultTable5)
ggsave("VolcanoPlot5.jpeg", plot=vp5, path = "/tmp/repo/DGE_LADC/results/plots")
vp10 <- volcano_plot(resultTable10)
ggsave("VolcanoPlot10.jpeg", plot=vp10, path = "/tmp/repo/DGE_LADC/results/plots")

##################################################################################################################
# Extraction of differentially expressed genes - comparison with the paper
##################################################################################################################

#### TME-related genes
length(which(resultTable$results == 1))
length(which(resultTable$results == -1))

diff.ex.genes <- read.csv(snakemake@input[[8]], sep='')
## add the results columns, 0, 1, -1
diff.ex.genes$results <- rep(0, nrow(diff.ex.genes))
diff.ex.genes[which(diff.ex.genes$Log2.Fold_Change.>1 & diff.ex.genes$pValue < 0.05), 'results'] <- 1
diff.ex.genes[which(diff.ex.genes$Log2.Fold_Change.< (-1) & diff.ex.genes$pValue < 0.05), 'results'] <- -1

## how many of paper's DE Genes are in my DE Genes?
sum(diff.ex.genes$Gene %in%  resultTable$Genes[which(resultTable$results != 0)]) 
## 77 

sum(diff.ex.genes$Gene %in%  resultTable$Genes) 
## 91 in general 
length(diff.ex.genes$Gene) 
## 2 gene not in the data set
diff.ex.genes$Gene[which(!(diff.ex.genes$Gene %in% resultTable$Genes))] ## genes "TENM4" "WISP1" not in my data set


##################################################################################################################
# Extraction of differentially expressed genes - Venn Diagramms 
##################################################################################################################

#genes.up = list(original.results =diff.ex.genes$Gene[which(diff.ex.genes$results == 1)], 
#                dge = names(which(results@.Data[, 1]==1)), 
#                dge.05 = names(which(results.05@.Data[, 1]==1)), 
#                dge.1 = names(which(results.1@.Data[, 1]==1)), 
#                dge.25 = names(which(results.25@.Data[, 1]==1)), 
#                dge.5 = names(which(results.5@.Data[, 1]==1)), 
#                dge.10 = names(which(results.10@.Data[, 1]==1)))
#genes.down = list(original.results =diff.ex.genes$Gene[which(diff.ex.genes$results == -1)], 
#                  dge = names(which(results@.Data[, 1]==-1)), 
#                  dge.05 = names(which(results.05@.Data[, 1]==-1)), 
#                  dge.1 = names(which(results.1@.Data[, 1]==-1)), 
#                  dge.25 = names(which(results.25@.Data[, 1]==-1)), 
#                  dge.5 = names(which(results.5@.Data[, 1]==-1)), 
#                  dge.10 = names(which(results.10@.Data[, 1]==-1)))

#venn.diagram(genes.up[3:7], filename = "/tmp/repo/DGE_LADC/results/plots/vennUp1.tiff", 
#             fill = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#008E60"), imagetype="tiff", main = 'Venn Diagram comparing differentially private results')
            
#venn.diagram(genes.up[c(1,2,5,6,7)], filename = "/tmp/repo/DGE_LADC/results/plots/vennUp2.tiff", 
#             fill = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#008E60"), imagetype="tiff", main = 'Venn Diagram comparing the original results with the results obtained by DGEA')


#venn.diagram(genes.down[3:7], filename = "/tmp/repo/DGE_LADC/results/plots/vennDown1.tiff", 
#             fill = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#008E60"), imagetype="tiff", main = 'Venn Diagram comparing the differentially private results') 
            
#venn.diagram(genes.down[c(1,2,5,6,7)], filename = "/tmp/repo/DGE_LADC/results/plots/vennDown2.tiff", 
#             fill = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#008E60"), imagetype="tiff", main = 'Venn Diagram comparing the original results with the results obtained by DGEA')


##################################################################################################################
# Extraction of differentially expressed genes - save
##################################################################################################################


write.csv(resultTable, file=snakemake@output[[1]])
write.csv(resultTable05, file=snakemake@output[[2]])
write.csv(resultTable1, file=snakemake@output[[3]])
write.csv(resultTable25, file=snakemake@output[[4]])
write.csv(resultTable5, file=snakemake@output[[5]])
write.csv(resultTable10, file=snakemake@output[[6]])
write.csv(resultTableNN, file=snakemake@output[[7]])

