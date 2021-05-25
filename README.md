# InterdisciplinaryProject
Interdisciplinary Project in Data Science - Bioinformatics

Author: Sanja Priselac 

### What is this repository for ? ###

Project in Interdisciplinary Project in Data Science 

### Folder structure ###
*taken from https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html*

*data*: data divided in the folders, here only source URL provided	
 * cancerData: data taken from TCGA database
 * Metadata: metadata about the project

*workflow*: files needed for executing the workflow, Snakemake file, the R scripts in *scripts* folder, .yaml environment specifications in *envs* folder

*results*: partial results (no *expression_matrix.RData*, *expression_matrix_annotated.RData* and *y_coding.RData* due to big file sizes), RData files in *R* folder and plots in *plots* folder

*config*: empty for now


### How to run ? ###
1) Activate the snakemake environment (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) 
2) Go to the *workflow* folder containing the Snakemake file
3) execute *snakemake --cores [] --use-conda* command 

