#!/opt/conda/envs/R_env/bin/Rscript --vanilla
# Sanja Priselac
# 21/04/2021

#####################################
library(biomaRt)
library(edgeR)
library(data.table) ## should be added 
#####################################

#####################################################################################################################################
## Load the expression matrix 
#####################################################################################################################################

load(snakemake@input[[1]])

expression.matrix <- expression.data[[1]]
file_id <- expression.data[[2]]

############################################################################################################################
## map the ensamble ids to the gene names - annotation
############################################################################################################################

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

ens_ids <- rownames(expression.matrix)
mapping <- getBM(attributes=c('ensembl_gene_id',
                              'external_gene_name'),
                 #values = ens_ids,
                 values=c('protein_coding'),
                 mart = ensembl)

#str(mapping)

## just protein - coding
genes <- merge(data.frame(ensembl_gene_id = ens_ids), mapping,  by="ensembl_gene_id")

expression.matrix <- expression.matrix[which(rownames(expression.matrix) %in% genes$ensembl_gene_id == TRUE), ]

sum(order(rownames(expression.matrix)) == order(genes$ensembl_gene_id))
# same order ==> take new names for the genes 

rownames(expression.matrix) <- genes$external_gene_name

expression.data.annotated <- list(expression.matrix, file_id)
save(expression.data.annotated, file = "/tmp/repo/DGE_LADC/results/expression_matrix_annotated.RData")

#####################################################################################################################################
## samples - add the condition 
#####################################################################################################################################

samples <- as.data.frame(fread(snakemake@input[[2]]))

## order in the same way as file_id
samples <- samples[order(samples$`File ID`), ]
order(samples$`File ID`) == order(file_id)

samples$Type <- ifelse(samples$`Sample Type` == "Solid Tissue Normal", 'Normal', 'NonNormal')
samples$Type <- as.factor(samples$Type)

#####################################################################################################################################
## read the TME genes   
#####################################################################################################################################

tme.genes.table <- read.csv(snakemake@input[[3]], header=TRUE)
tme.genes <- c(unname(unlist(tme.genes.table[1, -c(1,2)])), unname(unlist(tme.genes.table[2, -c(1,2)])))

## 276 / 282 genes in the data set
subset.genes <- which(rownames(expression.matrix) %in% tme.genes == TRUE)

#####################################################################################################################################
## create DGE objects 
##  * all coding genes 
##  * TME genes
#####################################################################################################################################

iszero <- function(vec) { return(sum(vec == 0)/ length(vec)) }
zeros <- apply(expression.matrix, 1, iszero)

y <- DGEList(expression.matrix[-which(zeros == 1), ], group=samples$Type)

# tme.genes[which(tme.genes %in% rownames(expression.matrix) == FALSE)]
## 6 genes from TME missing : "WISP1"  "GPR124" "LPPR4"  "TXNDC3" "ODZ4" "FYB"   
y.tme <- DGEList(expression.matrix[subset.genes, ], group = samples$Type)


save(y, file=snakemake@output[[1]])
save(y.tme, file=snakemake@output[[2]])
