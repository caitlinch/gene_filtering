### Small script to apply Divvier and trimAl to single gene alignments from Strassert 2021 dataset
# Same process they followed in the original paper

# Set path to each executable and to the folder containing the filtered single gene alignments from Strassert 2021
divvier_path <- "/Users/caitlincherryh/Documents/Executables/Divvier-1.01/divvier"
trimal_path <- "/Users/caitlincherryh/Documents/Executables/trimal-1.4.1/source/trimal"
gene_folder <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Strassert2021/01_filtered_genes_from_supplement/"
edited_gene_folder <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Strassert2021/02_trimAL_Divvier_filtered_genes_only/"

# Get list of all alignments in the gene folder
all_genes <- paste0(gene_folder, list.files(gene_folder))

# Apply divvier using the -partial option
# command will be ./divvier –partial -mincol 4 -divvygap myfile.fas
run.divvier <- function(alignment_path, divvier_path){
  command <- paste0(divvier_path, " –partial -mincol 4 -divvygap ", alignment_path)
  system(command)
}
# Run divvier on all alignments
lapply(all_genes, run.divvier, divvier_path)

# List all alignments in gene folder
all_genes <- paste0(gene_folder, list.files(gene_folder))
# Extract all genes that were processed in divvier (will have ".divvy." in name)
divvy_genes <- grep("divvy", all_genes, value = TRUE)

# Apply trimal using the command "-gt 0.05"
# Command will be trimal -in <inputfile> -out <outputfile> -(other options)
run.trimal <- function(alignment_path, trimal_path){
  trimal_gene <- gsub("divvy.fas", "divvy.trimal.fas", alignment_path)
  command <- paste0(trimal_path, " -in ", alignment_path, " -out ", trimal_gene, " -gt 0.05")
  system(command)
}
# Run trimal on all alignments
lapply(divvy_genes, run.trimal, trimal_path)

# Get a list of all the genes that have been run through both divvier and trimal
all_genes <- paste0(list.files(gene_folder))
complete_genes <- grep("divvy.trimal.fas", all_genes, value = TRUE)
original_genes <- paste0(gene_folder, complete_genes)
moved_genes <- paste0(edited_gene_folder, complete_genes)
# Copy the version of the gene that was run through both programs into a fresh folder
for (i in 1:length(complete_genes)){
  file.copy(from = original_genes[i], to = moved_genes[i], overwrite = TRUE, recursive = FALSE, copy.mode = TRUE, copy.date = FALSE)
}
