# Code to read through nexus files and remove the empty sequences (all gaps, all ? or all either)
BA_dir <- "/Users/caitlincherryh/Documents/Repositories/BenchmarkAlignments_DataSubSet/"

files <- list.files(BA_dir,recursive = TRUE) # list all files
als <- paste0(BA_dir,files[grep(".nex",files)]) # get all the nexus files
als <- als[!als %in% als[grep(".nex.",als)]] # remove all non alignment files to leave only alignments
als <- als[!als %in% als[grep("bootstrapReplicate",als)]] # remove all bootstrap alignments (if any present)
n_files <- als[!als %in% als[grep("alignment.nex",als)]] # remove full alignments (only want to run per loci)

# Remove N or n or other ambiguous characters in the files
for (n_file in n_files){
  print(n_file)
  # for each alignment:
  #     - copy the sequence out from the new alignment: temp <- as.numeric(new_aln$X)
  #     - replace the non 18s with the generated sequence of the right length: temp[which(new_aln$X !=18)] <- new_seq
  #     - replace the new seq into the new aln: new_aln$X <- temp
  # Open the new alignment as a nexus file
  n <- read.nexus.data(n_file) # open file
  n_new <- n # copy file
  # Get the names of all the sequences
  seq_names <- names(n_new)
  print(paste0("number of names: ",length(seq_names)))
  print(paste0("number of unique names: ",length(unique(seq_names))))
  # Iterate through the names
  for (seq_name in seq_names){
    original_seq <- n[[seq_name]] # get the original empirical sequence
    new_seq <- n_new[[seq_name]] # get the new simulated sequence that has the same name
    n_inds <- which(original_seq == "n") # find out which sites are A/G/T/C in the original alignment
    new_seq[n_inds] <- "?" # switch the Ns for unknowns
    N_inds <- which(original_seq == "N") # find out which sites are A/G/T/C in the original alignment
    new_seq[N_inds] <- "?" # switch the Ns for unknowns
    x_inds <- which(original_seq == "R") # find out which sites are A/G in the original alignment
    new_seq[x_inds] <- "?" # switch the x's for unknowns
    x_inds <- which(original_seq == "r") # find out which sites are A/G in the original alignment
    new_seq[x_inds] <- "?" # switch the x's for unknowns
    x_inds <- which(original_seq == "M") # find out which sites are A/C in the original alignment
    new_seq[x_inds] <- "?" # switch the x's for unknowns
    x_inds <- which(original_seq == "m") # find out which sites are A/C in the original alignment
    new_seq[x_inds] <- "?" # switch the x's for unknowns
    x_inds <- which(original_seq == "W") # find out which sites are A/T in the original alignment
    new_seq[x_inds] <- "?" # switch the x's for unknowns
    x_inds <- which(original_seq == "w") # find out which sites are A/T in the original alignment
    new_seq[x_inds] <- "?" # switch the x's for unknowns
    x_inds <- which(original_seq == "S") # find out which sites are G/C in the original alignment
    new_seq[x_inds] <- "?" # switch the x's for unknowns
    x_inds <- which(original_seq == "s") # find out which sites are G/C in the original alignment
    new_seq[x_inds] <- "?" # switch the x's for unknowns
    x_inds <- which(original_seq == "K") # find out which sites are G/T in the original alignment
    new_seq[x_inds] <- "?" # switch the x's for unknowns
    x_inds <- which(original_seq == "k") # find out which sites are G/T in the original alignment
    new_seq[x_inds] <- "?" # switch the x's for unknowns
    x_inds <- which(original_seq == "Y") # find out which sites are C/T in the original alignment
    new_seq[x_inds] <- "?" # switch the x's for unknowns
    x_inds <- which(original_seq == "y") # find out which sites are C/T in the original alignment
    new_seq[x_inds] <- "?" # switch the x's for unknowns
    x_inds <- which(original_seq == "V") # find out which sites are A/G/C in the original alignment
    new_seq[x_inds] <- "?" # switch the x's for unknowns
    x_inds <- which(original_seq == "v") # find out which sites are A/G/C in the original alignment
    new_seq[x_inds] <- "?" # switch the x's for unknowns
    x_inds <- which(original_seq == "H") # find out which sites are A/C/T in the original alignment
    new_seq[x_inds] <- "?" # switch the x's for unknowns
    x_inds <- which(original_seq == "h") # find out which sites are A/C/T in the original alignment
    new_seq[x_inds] <- "?" # switch the x's for unknowns
    x_inds <- which(original_seq == "D") # find out which sites are A/G/T in the original alignment
    new_seq[x_inds] <- "?" # switch the x's for unknowns
    x_inds <- which(original_seq == "d") # find out which sites are A/G/T in the original alignment
    new_seq[x_inds] <- "?" # switch the x's for unknowns
    x_inds <- which(original_seq == "B") # find out which sites are G/C/T in the original alignment
    new_seq[x_inds] <- "?" # switch the x's for unknowns
    x_inds <- which(original_seq == "b") # find out which sites are G/C/T in the original alignment
    new_seq[x_inds] <- "?" # switch the x's for unknowns
    n_new[[seq_name]] <- new_seq
  }
  # Output the final alignment (same parameters and gaps as input alignment) as a nexus file
  print("output nexus file")
  write.nexus.data(n_new,file = n_file, format = "dna", interleaved = FALSE)
  # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
  nexus_edit <- readLines(n_file) # open the new nexus file
  ind <- grep("BEGIN DATA",nexus_edit)+2 # find which line
  nexus_edit[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE;" # replace the line
  writeLines(nexus_edit,n_file) # output the edited nexus file
}

# Remove empty sequences in the files
for (n_file in n_files){
  print(n_file)
  n <- read.nexus.data(n_file)
  # Initialise a new sequence
  copy_names <- c()
  seq_names <- names(n)
  # Iterate through the names and add non-missing sequences to the new sequence
  for (seq_name in seq_names){
    seq <- n[[seq_name]] # get the original empirical sequence
    chars <- unique(seq)
    if (setequal(chars,c("?")) || setequal(chars,c("-")) || setequal(chars,c("-","?"))){
      # If the only characters are the empty character or the question mark character, ignore this sequence
      next
    } else {
      # If the sequence contains genetic information, include it in the new sequence
      copy_names <- c(copy_names,seq_name)
    }
  }
  n_new <- n[copy_names]
  print(str(n_new))
  write.nexus.data(n_new,file = n_file, format = "dna", interleaved = FALSE)
  # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
  nexus_edit <- readLines(n_file) # open the new nexus file
  ind <- grep("BEGIN DATA",nexus_edit)+2 # find which line
  nexus_edit[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE;" # replace the line
  writeLines(nexus_edit,n_file) # output the edited nexus file
}

