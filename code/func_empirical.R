### empirical_treelikeness/code/func_empirical.R
## R functions to process and modify empirical sequence alignments, specifically Rob Lanfear's "BenchmarkAlignments"
## BenchmarkAlignments and metadata available here: https://github.com/roblanf/BenchmarkAlignments and have CC0 or CCBY licenses
# Caitlin Cherryh 2019



cutSpecies <- function(alignment_path, keep, output_path){
  # open nexus file and get species names
  n <- read.nexus.data(alignment_path)
  # check which of the keep species are in the input names and adjust accordingly
  in_names <- names(n)
  keep <- keep[keep %in% in_names]
  # trim the nexus file to contain only the keep species
  n_trim <- n[keep]
  # output the nexus file under the output_path
  write.nexus.data(n_trim, file = output_path, format = "dna",interleaved = TRUE, datablock = FALSE) # write the output as a nexus file
  # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
  nexus <- readLines(output_path) # open the new nexus file
  ind <- grep("BEGIN CHARACTERS",nexus)+2 # find which line
  nexus[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE;" # replace the line
  writeLines(nexus,output_path) # output the edited nexus file
}
