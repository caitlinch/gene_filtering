### empirical_treelikeness/code/func_analysis.R
## R functions to analyse empirical treelikeness data
# Caitlin Cherryh 2021


# Given an alignment, this function will check the IQ-Tree .log file and .iqtree file for warnings and return all warnings and the corresponding alignment
check.for.IQTree.warnings <- function(alignment_path){
  # Open the .iqtree and .log files
  dotiqtree_path <- paste0(alignment_path, ".iqtree")
  dotiqtree_file <- readLines(dotiqtree_path)
  dotlog_path <- paste0(alignment_path, ".log")
  dotlog_file <- readLines(dotlog_path)
  
  # Check for any warnings
  dotiqtree_warnings <- dotiqtree_file[grep("WARNING", dotiqtree_file)]
  dotlog_warnings <- dotlog_file[grep("WARNING", dotlog_file)]
  
  # Extract loci name from alignment
  alignment_path_base <- basename(alignment_path)
  alignment_path_split <- strsplit(alignment_path_base, "\\.")[[1]]
  loci_name <- paste0(alignment_path_split[1:(length(alignment_path_split) - 1)], collapse = ".")
  
  # If one or more warnings was found, output a dataframe of the warnings
  if ((length(dotiqtree_warnings) + length(dotlog_warnings)) > 0){
    # Output warnings
    warning_df <- data.frame(loci_path = alignment_path, loci = loci_name, file = c(rep(".iqtree", length(dotiqtree_warnings)), rep(".log", length(dotlog_warnings))),
                             warnings = c(dotiqtree_warnings, dotlog_warnings))
    return(warning_df)
  }
}



# This function looks in a single IQ-Tree folder and checks which loci are present
get.loci.from.analysis <- function(folder, output_folder){
  # Create a filename for the output file
  op_filename <- paste0(output_folder, gsub("_IQ-Tree_partition","",basename(folder)), "_loci.csv")
  # Get all files in the folder
  folder_files <- list.files(folder)
  # Remove any partition files from the folder
  nonpartition_files <- grep("partition", folder_files, invert = TRUE, value = TRUE)
  # Get the location of each alignment
  loci_files <- paste0(folder, nonpartition_files)
  # Get the loci name from each alignment
  loci_names <- gsub(".fa", "", nonpartition_files)
  # Combine the names and locations of each alignment
  loci_df <- data.frame(loci_name = loci_names, loci_file = loci_files)
  # Output this information about this analysis
  write.csv(loci_df, op_filename, row.names = FALSE)
}





  

