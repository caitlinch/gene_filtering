# empirical_treelikeness/code/0_Whelan2015_data_formatting.R
# Caitlin Cherryh, 2021
## This script recreates the 100kb genomic window alignments used to estimate gene trees in Pease et al (2016)

## Whelan et al (2015) paper: 
#       Nathan V. Whelan, Kevin M. Kocot, Leonid L. Moroz, Kenneth M. Halanych 2015. Ctenophora is sister to all other animals.
#           Proceedings of the National Academy of Sciences May 2015, 112 (18) 5773-5778; DOI: 10.1073/pnas.1503453112

## Whelan et al (2015) data:
#       Whelan, Nathan; M. Kocot, Kevin; Moroz, Leonid L.; Kenneth M. Halanych 2016. 
#           Error, signal, and the placement of Ctenophora sister to all other animals. figshare. Dataset. 
#           https://doi.org/10.6084/m9.figshare.1334306.v3 
#

## This script:
# 1. Reads in the supermatrix and the gene locations from dataset "10"
# 2. Separates each gene into a different file
# 3. Outputs some summary information about the genes

## Open required libraries
library(phylotools)
library(phangorn)

## Parameters
# dataset_dir <- unzipped "Final_Datasets.tar.gz" folder (downloaded from figshare)
# output_dir <- file where alignments for individual genes will be saved

dataset_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Whelan2015/"
output_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Whelan2015/10_genes/"

## Open supermatrix file and file containing gene locations
supermat_file <- paste0(dataset_dir, "10/Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.phy")
location_file <- paste0(dataset_dir, "10/Dataset10_GeneList_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.txt")

# Open supermatrix as an alignment using Phylotools and write it out as a fasta file
fasta_supermat_file <- gsub("\\.phy", "\\.fa", supermat_file)
if (file.exists(fasta_supermat_file) == FALSE){
  p <- read.phylip(supermat_file)
  dat2fasta(p, outfile = fasta_supermat_file)
  # Open fasta file
  aa_mat <- as.matrix(read.aa(file = fasta_supermat_file, format = "fasta"))
} else {
  # Open fasta file
  aa_mat <- as.matrix(read.aa(file = fasta_supermat_file, format = "fasta"))
}

# Open location file as lines
locations <- readLines(location_file)

## Save each gene
print('Saving individual genes:')
all_gene_lengths <- c()
all_gene_names <- c()
for (i in 1:length(locations)){
  # Get the ith gene from the list of genes
  gene_line <- locations[i]
  # Split the line to determine the gene name
  gene_line <- gsub(";", "", gene_line)
  gene_line_split <- strsplit(gene_line, "=")[[1]]
  gene_name <- gsub(" ", "", gene_line_split[1])
  print(gene_name)
  # Split the line to determine the gene start and end position
  gene_range <- gsub(" ", "", gene_line_split[2])
  gene_range_split <- strsplit(gene_range, "-")[[1]]
  gene_start <- as.numeric(gene_range_split[[1]])
  gene_end <- as.numeric(gene_range_split[[2]])
  
  # Save name and length of gene
  gene_length <- length(gene_start:gene_end)
  all_gene_lengths <- c(all_gene_lengths, gene_length)
  all_gene_names <- c(all_gene_names, gene_name)
  
  # Subset the supermatrix to get the sites (columns) for this gene
  # To subset matrix: matrix[row, col]
  gene_mat <- aa_mat[, gene_start:gene_end]
  
  # Assemble file name for gene
  gene_file <- paste0(output_dir, gene_name, ".fa")
  # Write gene to file
  write.FASTA(gene_mat, file = gene_file)
}

# Create dataframe with gene name and lengths
df <- data.frame(dataset = "Whelan2015_10", gene_name = all_gene_names, gene_length = all_gene_lengths)
write.csv(df, file = paste0(dataset_dir, "Whelan2015_10_gene_lengths.csv"), row.names = FALSE)



