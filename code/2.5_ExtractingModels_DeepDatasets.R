### gene_filtering/code/2.5_ExtractingModels_DeepDataset.R
## R program to estimate trees from treelike or non-treelike loci
# Caitlin Cherryh 2023

## This script:
# 1. Determines the best model for each locus that is NOT a free rate model (because for the deep datasets, using free rate models made tree estimation computationally intractable)



##### Step 1: Set file paths and run variables #####
# input_names       <- set name(s) for the dataset(s) - make sure input_names is in same order as alignment_dir 
#                      (e.g. for 2 datasets, put same dataset first and same dataset last)
# gene_tree_dir     <- directory containing the results from estimating the gene tree for each locus in IQ-Tree
# output_dir        <- where the output csv file will be stored
# maindir           <- "gene_filtering" repository location (github.com/caitlinch/gene_filtering)

### Caitlin's paths ###
run_location = "server"
if (run_location == "local"){
  input_names <- c("1KP")
  gene_tree_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/"
  output_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/"
  maindir <- "~/Documents/Repositories/empirical_treelikeness/"
} else if (run_location=="server"){
  input_names <- c("1KP")
  gene_tree_dir <- "/data/caitlin/empirical_treelikeness/Output/"
  output_dir <- "/data/caitlin/empirical_treelikeness/Output/"
  maindir <- "/data/caitlin/empirical_treelikeness/code/"
  
}
### End of Caitlin's paths ###



##### Step 2: Extract models for each locus #####
# Source the file containing the function required for extracting ModelFinder details from .iqtree files
source(paste0(maindir, "func_empirical.R"))

# Iterate through each locus in each dataset and identify the best model and the best model without FreeRate model parameters (+R)
# Save the results as a dataframe per dataset
for (dataset in input_names){
  print(paste0("Dataset: ", dataset))
  # Assemble name of folder containing gene tree results for this dataset
  dataset_dir <- paste0(gene_tree_dir, dataset, "/")
  # Get the names of each loci
  loci_names <- list.files(dataset_dir)
  # Iterate through each loci name and collect the models and AIC/BIC values from ModelFinder
  loci_models_list <- lapply(loci_names, get.ModelFinder.models, dataset = dataset, 
                             dataset_directory = dataset_dir, return.best.model.without.free.rates = TRUE)
  loci_models_df <- as.data.frame(do.call(rbind, loci_models_list))
  # Save as csv
  loci_models_file <- paste0(output_dir, "02.5_", dataset, "_loci_models_noFreeRates.csv")
  write.csv(loci_models_df, loci_models_file)
}


