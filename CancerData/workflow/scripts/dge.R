#!/opt/conda/envs/R_env/bin/Rscript --vanilla
# Sanja Priselac
# 20/04/2021

##################
library(ggplot2)
library(limma)
library(edgeR)
library(VennDiagram)
library(svglite)
##################

for(i in 1:6) {
  load(snakemake@input[[i]])
}


##################################################################################################################

condition <- log2.y.qn$samples$group 
design <- model.matrix(~ 0 + condition) 

fit <- lmFit(log2.y.qn$counts, design)
fit.05 <- lmFit(y.qn.05$counts, design)
fit.1 <- lmFit(y.qn.1$counts, design)
fit.25 <- lmFit(y.qn.25$counts, design)
fit.5 <- lmFit(y.qn.5$counts, design)
fit.10 <- lmFit(y.qn.10$counts, design)

cont.matrix <- makeContrasts(
  condition = conditionNonNormal - conditionNormal,
  levels=design)

#cont.matrix

##################################################################################################################

#fit1 <- treat(fit, lfc=log(1.2), trend=TRUE)
## all 0 in fit1

fit2 <- contrasts.fit(fit, cont.matrix)
fit2.05 <- contrasts.fit(fit.05, cont.matrix)
fit2.1 <- contrasts.fit(fit.1, cont.matrix)
fit2.25 <- contrasts.fit(fit.25, cont.matrix)
fit2.5 <- contrasts.fit(fit.5, cont.matrix)
fit2.10 <- contrasts.fit(fit.10, cont.matrix)

fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)
fit2.05 <- eBayes(fit2.05, trend=TRUE, robust=TRUE)
fit2.1 <- eBayes(fit2.1, trend=TRUE, robust=TRUE)
fit2.25 <- eBayes(fit2.25, trend=TRUE, robust=TRUE)
fit2.5 <- eBayes(fit2.5, trend=TRUE, robust=TRUE)
fit2.10 <- eBayes(fit2.10, trend=TRUE, robust=TRUE)


results <- decideTests(fit2, adjust.method = "fdr", p.value=0.05, lfc = 1) 
results.05 <- decideTests(fit2.05, adjust.method = "fdr", p.value=0.05, lfc = 1) 
results.1 <- decideTests(fit2.1, adjust.method = "fdr", p.value=0.05, lfc = 1) 
results.25 <- decideTests(fit2.25, adjust.method = "fdr", p.value=0.05, lfc = 1) 
results.5 <- decideTests(fit2.5, adjust.method = "fdr", p.value=0.05, lfc = 1) 
results.10 <- decideTests(fit2.10, adjust.method = "fdr", p.value=0.05, lfc = 1) 

volcano_plot <- function (fit, results) {
  table <- topTable(fit, coef=1, number=Inf, adjust.method='BH', sort.by="P")
  #gene_names <- table$ID   does not work with the subset
  gene_names <- attributes(table)$row.names
  xx <- table$logFC
  name <- attributes(fit$contrasts)$dimnames$Contrasts
  yy <- table$adj.P.Val
  yy <- -log10(yy)
  result <- as.factor(unname(results@.Data[gene_names,1]))
  data <- data.frame(xx, yy, result)
  p <- ggplot(data, aes(x=xx, y=yy, col=result)) + 
    geom_point(size = 1.9) +
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
    #theme_minimal() +
    xlab(paste("Log folds for", name)) +
    ylab("-log10(adj.p.value)") +
    ggtitle(paste0("Volcano plot for ", name))
  return(p)  
}

vp <- volcano_plot(fit2, results)
ggsave("VolcanoPlot.pdf", plot=vp, path = "/tmp/repo/DGE_LADC/results/plots")
vp05 <- volcano_plot(fit2.05, results.05)
ggsave("VolcanoPlot05.pdf", plot=vp05, path = "/tmp/repo/DGE_LADC/results/plots")
vp1 <- volcano_plot(fit2.1, results.1)
ggsave("VolcanoPlot1.pdf", plot=vp1, path = "/tmp/repo/DGE_LADC/results/plots")
vp25 <- volcano_plot(fit2.25, results.25)
ggsave("VolcanoPlot25.pdf", plot=vp25, path = "/tmp/repo/DGE_LADC/results/plots")
vp5 <- volcano_plot(fit2.5, results.5)
ggsave("VolcanoPlot5.pdf", plot=vp5, path = "/tmp/repo/DGE_LADC/results/plots")
vp10 <- volcano_plot(fit2.10, results.10)
ggsave("VolcanoPlot10.pdf", plot=vp10, path = "/tmp/repo/DGE_LADC/results/plots")

##################################################################################################################
# Extraction of differentially expressed genes - comparison with the paper
##################################################################################################################

#### TME-related genes
length(which(results@.Data == 1))
length(which(results@.Data == -1))

diff.ex.genes <- read.csv(snakemake@input[[7]], sep='')
## add the results columns, 0, 1, -1
diff.ex.genes$results <- rep(0, nrow(diff.ex.genes))
diff.ex.genes[which(diff.ex.genes$Log2.Fold_Change.>1 & diff.ex.genes$pValue < 0.05), 'results'] <- 1
diff.ex.genes[which(diff.ex.genes$Log2.Fold_Change.< (-1) & diff.ex.genes$pValue < 0.05), 'results'] <- -1

my.diff.ex.genes <- (attributes(results@.Data)$dimnames[[1]])[which(results@.Data != 0)]

## how many of paper's DE Genes are in my DE Genes?
sum(diff.ex.genes$Gene %in% my.diff.ex.genes) 
## logged2 (0.05 added) data, limma-trend 48  
## logged2 (0.1 added) data, limma-trend 47

sum(diff.ex.genes$Gene %in% attributes(results@.Data)$dimnames[[1]]) 
## 91 in general 
length(diff.ex.genes$Gene) 
## 2 gene not in the data set
diff.ex.genes$Gene[which(!(diff.ex.genes$Gene %in% attributes(results@.Data)$dimnames[[1]]))] ## genes "TENM4" "WISP1" not in my data set


##################################################################################################################
# Extraction of differentially expressed genes - Venn Diagramms 
##################################################################################################################

## Condition:
DEGenes <- topTable(fit2, coef=1, number=Inf, adjust.method='BH', sort.by="P", genelist = attributes(log2.y.qn$counts)[[2]][[1]])
DEGenes.05 <- topTable(fit2.05, coef=1, number=Inf, adjust.method='BH', sort.by="P", genelist = attributes(log2.y.qn$counts)[[2]][[1]])
DEGenes.1 <- topTable(fit2.1, coef=1, number=Inf, adjust.method='BH', sort.by="P", genelist = attributes(log2.y.qn$counts)[[2]][[1]])
DEGenes.25 <- topTable(fit2.25, coef=1, number=Inf, adjust.method='BH', sort.by="P", genelist = attributes(log2.y.qn$counts)[[2]][[1]])
DEGenes.5 <- topTable(fit2.5, coef=1, number=Inf, adjust.method='BH', sort.by="P", genelist = attributes(log2.y.qn$counts)[[2]][[1]])
DEGenes.10 <- topTable(fit2.10, coef=1, number=Inf, adjust.method='BH', sort.by="P", genelist = attributes(log2.y.qn$counts)[[2]][[1]])


genes.up = list(original.results =diff.ex.genes$Gene[which(diff.ex.genes$results == 1)], 
                dge = names(which(results@.Data[, 1]==1)), 
                dge.05 = names(which(results.05@.Data[, 1]==1)), 
                dge.1 = names(which(results.1@.Data[, 1]==1)), 
                dge.25 = names(which(results.25@.Data[, 1]==1)), 
                dge.5 = names(which(results.5@.Data[, 1]==1)), 
                dge.10 = names(which(results.10@.Data[, 1]==1)))
genes.down = list(original.results =diff.ex.genes$Gene[which(diff.ex.genes$results == -1)], 
                  dge = names(which(results@.Data[, 1]==-1)), 
                  dge.05 = names(which(results.05@.Data[, 1]==-1)), 
                  dge.1 = names(which(results.1@.Data[, 1]==-1)), 
                  dge.25 = names(which(results.25@.Data[, 1]==-1)), 
                  dge.5 = names(which(results.5@.Data[, 1]==-1)), 
                  dge.10 = names(which(results.10@.Data[, 1]==-1)))

venn.diagram(genes.up[3:7], filename = "/tmp/repo/DGE_LADC/results/plots/vennUp1.tiff", 
             fill = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#008E60"), imagetype="tiff", main = 'Venn Diagram comparing differentially private results')
            
venn.diagram(genes.up[c(1,2,5,6,7)], filename = "/tmp/repo/DGE_LADC/results/plots/vennUp2.tiff", 
             fill = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#008E60"), imagetype="tiff", main = 'Venn Diagram comparing the original results with the results obtained by DGEA')


venn.diagram(genes.down[3:7], filename = "/tmp/repo/DGE_LADC/results/plots/vennDown1.tiff", 
             fill = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#008E60"), imagetype="tiff", main = 'Venn Diagram comparing the differentially private results') 
            
venn.diagram(genes.down[c(1,2,5,6,7)], filename = "/tmp/repo/DGE_LADC/results/plots/vennDown2.tiff", 
             fill = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#008E60"), imagetype="tiff", main = 'Venn Diagram comparing the original results with the results obtained by DGEA')


##################################################################################################################
# Extraction of differentially expressed genes - save
##################################################################################################################


save(DEGenes, file=snakemake@output[[1]])
save(DEGenes.05, file=snakemake@output[[2]])
save(DEGenes.1, file=snakemake@output[[3]])
save(DEGenes.25, file=snakemake@output[[4]])
save(DEGenes.5, file=snakemake@output[[5]])
save(DEGenes.10, file=snakemake@output[[6]])

