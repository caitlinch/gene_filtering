# empirical_treelikeness/code/0_Whelan2015_data_formatting.R
## This script separates the supermatrix from Whelan et al (2017) into fasta alignments for individual genes
# Caitlin Cherryh, 2022

## Whelan et al (2017) paper: 
#     Whelan, N.V., Kocot, K.M., Moroz, T.P. et al. Ctenophore relationships and their placement as the sister group to all other animals. 
#         Nat Ecol Evol 1, 1737–1746 (2017). https://doi.org/10.1038/s41559-017-0331-3

## Whelan et al (2017) data:
#     Whelan, Nathan; M. Kocot, Kevin; Moroz, Tatiana P.; Mukherjee, Krishanu; Williams, Peter; Paulay, Gustav; et al. (2017): 
#         Ctenophora Phylogeny Datasets and Core Orthologs. figshare. Dataset. https://doi.org/10.6084/m9.figshare.4484138.v1 

## This script:
# 1. Reads in the supermatrix and the gene locations from dataset "10" (the "Metazoa_Choano_RCFV_strict" dataset from the original paper, used to estimate the tree in Fig. 2.)
# 2. Separates each gene into a different file
# 3. Outputs some summary information about the genes



#### 1. Parameters ####
## Specify directories:
# dataset_dir   <- directory that contains the Metazoa_Choano_RCFV_strict.phy file from figshare and the partioning scheme for individual genes
# output_dir    <- file where alignments for individual genes will be saved
## Specify file names for supermatrix and gene partition file (in RAxML format):
# supermat_file <- name for RAxML supermatrix file (.phy file)
# location_file <- name for RAxML gene partition file (.txt file)

### Caitlin's paths ###
dataset_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Whelan2017/raw_data/"
output_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Whelan2017/genes/"
supermat_file <- paste0(dataset_dir, "Metazoa_Choano_RCFV_strict.phy")
location_file <- paste0(dataset_dir, "Whelan2017_gene_partitions.txt")
### End of Caitlin's paths ###



#### 2. Open required libraries ####
library(phylotools)
library(phangorn)



#### 3. Code body ####
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
df <- data.frame(dataset = "Whelan2017", gene_name = all_gene_names, gene_length = all_gene_lengths)
write.csv(df, file = paste0(dataset_dir, "Whelan2017_gene_lengths.csv"), row.names = FALSE)


