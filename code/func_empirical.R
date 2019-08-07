### empirical_treelikeness/code/func_empirical.R
## R functions to process and modify empirical sequence alignments, specifically Rob Lanfear's "BenchmarkAlignments"
## BenchmarkAlignments and metadata available here: https://github.com/roblanf/BenchmarkAlignments and have CC0 or CCBY licenses
# Caitlin Cherryh 2019


# Given a nexus file, this function removes all species except those specified and outputs a nexus file of the new (smaller) alignment
cutSpecies <- function(alignment_path, keep, output_path_provided = "FALSE", output_path){
  # open nexus file and get species names
  n <- read.nexus.data(alignment_path)
  # check which of the keep species are in the input names and adjust accordingly
  in_names <- names(n)
  keep <- keep[keep %in% in_names]
  # trim the nexus file to contain only the keep species
  n_trim <- n[keep]
  # create the output path variable
  if (output_path_provided == FALSE){
    out_path <- alignment_path
  } else if (output_path_provided == TRUE){
    out_path <- output_path
  }
  # output the nexus file under the out_path
  write.nexus.data(n_trim, file = out_path, format = "dna",interleaved = TRUE, datablock = FALSE) # write the output as a nexus file
  # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
  nexus <- readLines(out_path) # open the new nexus file
  ind <- grep("BEGIN CHARACTERS",nexus)+2 # find which line
  nexus[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE;" # replace the line
  writeLines(nexus,out_path) # output the edited nexus file
}

# Function to call IQ-tree, and estimate a maximum likelihood tree and corresponding the site concordance factors
# Site concordance factors (sCF) are the fraction of decisive alignmen sites supporting that branch
# sCF Citation: Minh B.Q., Hahn M., Lanfear R. (2018) New methods to calculate concordance factors for phylogenomic datasets. https://doi.org/10.1101/487801
sCF <- function(iqtree_path,alignment_path, num_cores = "AUTO"){
  # Check if the tree file already exists and if it doesn't, run IQ-tree and create it
  if (file.exists(paste0(alignment_path,".treefile")) == FALSE){
    # Given an alignment, estimate the maximum likelihood tree
    # to estimate: iqtree -s ALN_FILE -p PARTITION_FILE --prefix concat -bb 1000 -nt AUTO
    call <- paste0(iqtree_path," -s ",alignment_path," -nt ",num_cores," -redo -safe")
    system(call)
    # Create the command and call it in the system
    # for sCF: iqtree -t concat.treefile -s ALN_FILE --scf 100 --prefix concord -nt 10
    treefile <- paste0(alignment_path,".treefile")
    call <- paste0(iqtree_path," -t ",treefile," -s ",alignment_path," -- scf 100"," -nt ",num_cores," -redo -safe")
    system(call) # call IQ-tree!
  }
}