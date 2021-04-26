#!/opt/conda/envs/R_env/bin/Rscript --vanilla
# Sanja Priselac
# 22/04/2021

## This script contains the commands to produce the expression matrix

#########################################
library(stringr) 
library(readr) 
#########################################

path <- snakemake@input[[1]]

############################################################################################################################
## read the raw data - create the counts matrix 
############################################################################################################################

files <- dir(path)

c.matrix <- matrix(nrow=60483, ncol=594)


file <- dir(file.path(path, files[1]))
filepath <- file.path(path, files[1], file[which(grepl('\\gz$', file))])

P <- read.table(gzfile(filepath),sep="\t")
# order for genes
vec <- P$V2[order(P$V1)]
genes.names <- P$V1[order(P$V1)]
rownames(c.matrix) <- genes.names
c.matrix[, 1] <- vec

for (i in 2:length(files)) {
  print(i)
  file <- dir(file.path(path, files[i]))
  filepath <- file.path(path, files[i], file[which(grepl('\\gz$', file))])
  
  P <- read.table(gzfile(filepath),sep="\t")
  # order for genes
  vec <- P$V2[order(P$V1)]
  ## check if the gene names are the same
  if (sum(genes.names == P$V1[order(P$V1)]) != 60483) {
    print(i)
  } 
  c.matrix[, i] <- vec
}

## remove the version from the Ensamble Genes' names
rownames(c.matrix) <- str_replace(rownames(c.matrix), pattern = ".[0-9]+$", replacement = "")

expression.data <- list(expression.matrix = c.matrix, file_id = files)
save(expression.data, file = snakemake@output[[1]])


