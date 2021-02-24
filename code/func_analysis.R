### empirical_treelikeness/code/func_analysis.R
## R functions to analyse empirical treelikeness data
# Caitlin Cherryh 2021

all_alignments <- paste0("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/",
                         c("Vanderpool2020/ORTHOMCL11779/ORTHOMCL11779.nex", "Vanderpool2020/ORTHOMCL11783/ORTHOMCL11783.nex", "Vanderpool2020/ORTHOMCL11784/ORTHOMCL11784.nex",
                           "Vanderpool2020/ORTHOMCL11785/ORTHOMCL11785.nex", "Vanderpool2020/ORTHOMCL14552/ORTHOMCL14552.nex"))

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
