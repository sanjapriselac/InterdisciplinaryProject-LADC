#!/opt/conda/envs/R_env/bin/Rscript --vanilla
# Sanja Priselac
# 22/04/2021

########################
library(edgeR)
library(limma)
library(ggplot2)
library(gplots)
library(rmutil)
library(RColorBrewer)
########################

load(snakemake@input[[1]])

############################################################################################
# Normalisation
############################################################################################

y <- y.tme

############################################################################################
# Quantile normalisation function
############################################################################################

quantileNormalisation <- function(y, filepath, save=FALSE) {
  ## y is DGEList object
  # save the indices of the sort
  orders <- apply(y$counts, 2, order, decreasing = TRUE)
  x.sort <- apply(y$counts, 2, sort, decreasing = TRUE)
  r.means <- rowMeans(x.sort)
  if (save) {
    write.csv(r.means, file=filepath, row.names=FALSE)
  }
  for (i in 1:ncol(x.sort)) {
    y$counts[orders[, i], i] <- r.means
  }
  return(y)
} 

############################################################################################
# Laplace Quantile normalisation function 
############################################################################################

quantileNormalisation.dp <- function(y, b, filepath, seed=100) {
  ## y is DGEList object
  # save the indices of the sort
  orders <- apply(y$counts, 2, order, decreasing = TRUE)
  x.sort <- apply(y$counts, 2, sort, decreasing = TRUE)
  r.means <- rowMeans(x.sort)
  
  ## laplace distributed random variables
  set.seed(seed)
  lap <- rlaplace(nrow(y$counts), m=0, s=b)
  ## add the noise and sort the vector 
  r.means <- sort(r.means + lap, decreasing = TRUE)
  negatives <- which(r.means < 0)
  if (length(negatives) > 0) {
    warning(cat(length(negatives), 'r.means values negative and set to 0'))
    r.means[negatives] <- 0
  }
  write.csv(r.means, file=filepath, row.names=FALSE)
  for (i in 1:ncol(x.sort)) {
    y$counts[orders[, i], i] <- r.means
  }
  return(y)
} 

sensitivity <- function(y) {
  ## y is DGEList
  x.sort <- apply(y$counts, 2, sort, decreasing = TRUE)
  return(sum(apply(x.sort, 1, max) - apply(x.sort, 1, min))/ncol(x.sort))
}

############################################################################################
# Density plot, Sensitivity evaluation
############################################################################################

sensitivity(y)  
## TME Genes: 98.47
## all coding genes: 

## density plot
y.plot <- as.vector(t(y$counts[1:10, ]))
y.plot <- data.frame(values=y.plot, Gene= rep(rownames(y.tme$counts)[1:10], each=594))
p <- ggplot(y.plot, aes(values, fill = Gene)) + 
  geom_density(alpha = 0.2) + 
  xlab("Expression Values") +
  xlim(0, 300) + 
  theme(text = element_text(size=15))

y.qn <- quantileNormalisation(y, "/tmp/repo/DGE_LADC/results/mean_vector.csv", TRUE) 

log2.y <- y.tme
log2.y$counts <- t(apply(log2.y$counts + 1, 1, log2)) ## why does apply change the matrix?

## heatmap
heatmap.2(log2.y$counts[, c(which(color == 'green '), which(color != 'green '))], trace='none', scale = "none", dendrogram = "none", 
          ColSideColors =color[c(which(color == 'green '), which(color != 'green '))],  Colv=FALSE, sepcolor="white", margins=c(1, 4), cexCol = 0.4) 
legend("topright", title = "Tissues",legend=c("Normal","Not Normal"), fill=c("green", "blue"), cex=0.8, box.lty=0)

log2.y.qn <- quantileNormalisation(log2.y)

#plotDensities(log2.y, legend=FALSE)
## TME genes: better than original, but not necesarly needed
## all coding genes: 4 peaks for 4 combinations in a microarray: apears only in alive, only in dead, in both or none


############################################################################################
# Differentially Private Quantile normalisation
############################################################################################

eps <- c(0.5, 1, 2.5, 5, 10)
b <- floor(sensitivity(y)/eps)

## calculate laplacian 
y.qn.05 <- quantileNormalisation.dp(y.tme, b[1], "/tmp/repo/DGE_LADC/results/mean_vector_05.csv") ## 136 negateives
y.qn.1 <- quantileNormalisation.dp(y.tme, b[2], "/tmp/repo/DGE_LADC/results/mean_vector_1.csv") ## 136 negatives
y.qn.25 <- quantileNormalisation.dp(y.tme, b[3], "/tmp/repo/DGE_LADC/results/mean_vector_25.csv") ## 131 negatives
y.qn.5 <- quantileNormalisation.dp(y.tme, b[4], "/tmp/repo/DGE_LADC/results/mean_vector_5.csv")  ## 113 negatives
y.qn.10 <- quantileNormalisation.dp(y.tme, b[5], "/tmp/repo/DGE_LADC/results/mean_vector_10.csv") ## 95 negatives

#log2.y <- y.tme
y.qn.05$counts <- t(apply(y.qn.05$counts + 1, 1, log2)) ## why does apply change the matrix?
y.qn.1$counts <- t(apply(y.qn.1$counts + 1, 1, log2))
y.qn.25$counts <- t(apply(y.qn.25$counts + 1, 1, log2))
y.qn.5$counts <- t(apply(y.qn.5$counts + 1, 1, log2))
y.qn.10$counts <- t(apply(y.qn.10$counts + 1, 1, log2))

#plotDensities(y.qn.25, col=myPalette[1:16], legend=FALSE) 

## save the  
save(log2.y.qn, file=snakemake@output[[1]])
save(y.qn.05, file=snakemake@output[[2]])
save(y.qn.1, file=snakemake@output[[3]])
save(y.qn.25, file=snakemake@output[[4]])
save(y.qn.5, file=snakemake@output[[5]])
save(y.qn.10, file=snakemake@output[[6]])
save(log2.y, file=snakemake@output[[7]])

############################################################################################
# Data clustering (only for quantile-normalised data)
############################################################################################

pca <- prcomp(t(log2.y$counts[-which(apply(log2.y$counts, 1, var) == 0), ]), scale = T)
#plot(pca)
#summary(pca)

loadings <- pca$rotation
scores <- pca$x

data <- data.frame(scores, Condition = log2.y$samples$group)
plotPCs <- function(i, j) {
  x_lab <- paste("PC", i,": ", round(summary(pca)$importance[2,i],3)*100, "% variance", sep="")
  y_lab <- paste0("PC", j,": ", round(summary(pca)$importance[2,j],3)*100, "% variance")
  
  p <- ggplot(data, aes(x=data[,i], y=data[,j], color=Condition)) +
    geom_point() +
    xlab(x_lab) +
    ylab(y_lab) 
    labs(color="Condition")
  return(p)
}

pc12 <- plotPCs(1, 2)
ggsave("PCA12.pdf", plot=pc12, path = "/tmp/repo/DGE_LADC/results/plots")

pc23 <- plotPCs(2, 3)
ggsave("PCA23.pdf", plot=pc23, path = "/tmp/repo/DGE_LADC/results/plots")

pc34 <- plotPCs(3, 4)  
ggsave("PCA34.pdf", plot=pc34, path = "/tmp/repo/DGE_LADC/results/plots")

#############################################################################################################
