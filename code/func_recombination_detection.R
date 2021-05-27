### empirical_treelikeness/code/func_recombination_detection.R
## R functions to process loci from empirical phylogenetic datasets and apply recombination detection methods
# Caitlin Cherryh 2021

library(seqinr)
library(ape)

# Wrapper to take a dataframe row and feed it into the apply.recombination.detection.methods function
recombination.detection.wrapper <- function(index, df, executable_paths, iqtree_num_threads){
  loci_row <- df[index,]
  apply.recombination.detection.methods(loci_row = loci_row, executable_paths, iqtree_num_threads)
}



# Apply the recombination detection methods for a loci and then estimate the tree using IQ-Tree2
apply.recombination.detection.methods <- function(loci_row, executable_paths, iqtree_num_threads){
  # Create a new folder for the results from this locus
  alignment_folder <- paste0(loci_row$output_folder, loci_row$loci_name, "/")
  if (dir.exists(alignment_folder) == FALSE){
    dir.create(alignment_folder)
  }
  
  # Create the names for the output files for this locus
  recombination_results_file <- paste0(alignment_folder, loci_row$dataset, "_", loci_row$loci_name, 
                                       "_RecombinationDetection_results.csv")
  # Check if this output file exists
  # If it doesn't, run the recombination detection methods
  if (file.exists(recombination_results_file) == FALSE){
    # Set this directory as the current working directory (so output files from executables will get saved here)
    setwd(alignment_folder)
    # Copy the alignment into the output folder
    loci_row_path <- loci_row$loci_path
    file_extension_list <- strsplit(loci_row_path, "\\.")[[1]]
    file_extension <- file_extension_list[length(file_extension_list)]
    alignment_path <- paste0(alignment_folder, loci_row$loci_name, ".", file_extension)
    if (file.exists(alignment_path) == FALSE){
      file.copy(loci_row_path, alignment_path, overwrite = TRUE)
    }
    
    # Remove empty taxa from the alignment
    remove.empty.taxa(alignment_path, loci_row$alphabet)
    # Remove any sequences from the alignment that have more than the allowable proportion of ambiguous/missing sites or gaps
    if (is.na(loci_row$allowable_proportion_missing_sites) == FALSE){
      # If the allowable proportion of missing sites is NOT NA, run the pruning function
      # The pruning function will remove any sequences where the proportion of missing sites is GREATER than the allowable proportion of missing sites
      prune.taxa.by.length(alignment_path, loci_row$allowable_proportion_missing_sites, loci_row$alphabet, write_output_text = TRUE, alignment_folder)
    }
    
    # Run 3seq
    threeseq_results <- run.3seq(alignment_path, alignment_folder, threeseq_path = executable_paths[["3seq"]])
    # Run PhiPack
    phipack_results <- run.phipack(alignment_path, alignment_folder, phipack_path = executable_paths[["PHIPack"]], loci_row$alphabet)
    # Run GeneConv
    geneconv_results <- run.geneconv(alignment_path, alignment_folder, geneconv_path = executable_paths[["GeneConv"]], loci_row$alphabet)
    # Run IQ-Tree
    call.IQTREE.empirical.UFB(alignment_path, executable_paths[["IQTree"]], loci_row$best_model, num_bootstraps = 1000, num_threads = iqtree_num_threads)
    
    # Extract basic information about the alignment for the output csv file
    # Check file extension of alignment
    if (file_extension == "fasta" | file_extension == "fa" | file_extension == "fna" | file_extension == "ffn" | 
        file_extension == "faa" | file_extension == "frn" | file_extension == "fas"){
      # open alignment
      if (loci_row$alphabet == "dna"){
        seqtype = "DNA"
      } else if (loci_row$alphabet == "protein"){
        seqtype = "AA"
      }
      f <- read.fasta(file = alignment_path, seqtype = seqtype)
      n_taxa <- length(f)
      n_char <- length(f[[1]])
    } else if (file_extension == "nex" | file_extension == "nexus"){
      # open alignment
      n <- read.nexus.data(alignment_path)
      n_taxa <- length(n)
      n_char <- length(unlist(n[1]))
    }
    
    # Calculate extra variables for output
    # Calculate proportion of recombinant sequences
    prop_recomb_seq <- round(as.numeric(threeseq_results["num_distinct_recombinant_sequences"])/n_taxa,3)
    # Extract tree from iqtree file
    alignment_files <- list.files(alignment_folder)
    tree_file <- paste0(alignment_folder, grep("treefile", alignment_files, value = TRUE))
    newick_tree <- readLines(tree_file)
    
    # Combine results into a nice row for output csv file
    results_vec <- c(loci_row$dataset, loci_row$loci_name, loci_row$alphabet, loci_row$best_model,
                     n_taxa, n_char, threeseq_results, prop_recomb_seq, phipack_results,
                     geneconv_results, newick_tree)
    results_df <- data.frame(as.list(results_vec))
    result_names <- c("dataset", "loci_name", "alphabet", "best_model", "n_taxa", "n_bp",
                      names(threeseq_results), "proportion_recombinant_sequences", names(phipack_results),
                      names(geneconv_results), "tree")
    names(results_vec) <- result_names
    names(results_df) <- result_names
    write.csv(results_df, file = recombination_results_file, row.names = FALSE)
  } else if (file.exists(recombination_results_file) == TRUE){
    results_df <- read.csv(recombination_results_file)
    result_names <- names(results_df)
    results_vec <- as.character(results_df[1,])
    names(results_vec) <- result_names
  }
  return(results_vec)
}



# Run the 3SEQ program on one alignment and return the key variables calculated by the program
run.3seq <- function(alignment_path, alignment_folder, threeseq_path){
  # Change to the log (storage for log files) folder for this alignment - means that 3seq and Phi files will be saved into a unique folder
  setwd(alignment_folder)
  
  # Determine file extension
  file_extension <- tail(strsplit(alignment_path,"\\.")[[1]],1)
  
  # Only run 3SEQ and collect sCF values if the path is an alignment (don't need bootstrap replicates for 3SEQ/sCF)
  # Check if 3SEQ has already been run. Only run 3SEQ if the log file doesn't exist
  if (file.exists(paste0(alignment_folder,"3s.log")) == FALSE){
    if (file_extension == "fasta" | file_extension == "fa" | file_extension == "fna" | file_extension == "ffn" | 
        file_extension == "faa" | file_extension == "frn" | file_extension == "fas"){
      # 3SEQ only reads in phylip or fasta format - if the alignment is already in that format, call 3SEQ on the alignment
      seq_command <- paste0(threeseq_path," -f ", alignment_path)
      system(seq_command) #call 3SEQ
    } else if (file_extension == "nex" | file_extension == "nexus"){
      # 3SEQ only reads in Phylip or fasta format - need to convert if the alignment is a nexus file (using the nexus data opened above)
      # Assemble a name for a copy of the alignment as a fasta file
      fasta.name <- gsub(file_extension, ".fasta", alignment_path)
      # Check if the fasta file already exists
      if (file.exists(fasta.name) == FALSE){
        # If the fasta version doesn't exist, write out the nexus sequence in fasta formt
        n <- read.nexus.data(alignment_path)
        write.fasta(sequences = n, names = names(n), file.out = fasta.name) # output alignment as a fasta format
      }
      # There is now a definitely a fasta format version of this alignment
      # Assemble the 3SEQ command using the new fasta alignment
      seq_command <- paste0(threeseq_path, " -f ", fasta.name)
      # Call 3SEQ
      system(seq_command)
    }
  }
  # Now, collect the results from the 3SEQ log file
  seq_file <- paste0(alignment_folder,"3s.log")
  seq_log <- readLines(seq_file) # open file
  ind      <- grep("Number of recombinant triplets",seq_log) # find the number of recombinant triplets line index
  num_trips <- seq_log[ind]
  num_trips <- strsplit(num_trips,":")[[1]][2] # extract the number of recombinant triplets
  num_trips <- trimws(num_trips) # trim the whitespace from the number of triplets
  ind      <- grep("Number of distinct recombinant sequences",seq_log) # find the number of distinct recombinant sequences line index
  num_dis <- seq_log[ind]
  num_dis <- strsplit(num_dis,":")[[1]][2] # extract the number of distinct recombinant sequences
  num_dis <- trimws(num_dis) # trim the whitespace from the number of distinct recombinant sequences
  # null hypothesis is of clonal evolution - need significant p-value to accept the alternative hypothesis
  ind      <- grep("Rejection of the null hypothesis of clonal evolution",seq_log) # find the p value line index
  seq_sig <- seq_log[ind]
  seq_sig <- strsplit(seq_sig,"=")[[1]][2] # extract the p value
  seq_sig <- trimws(seq_sig) # trim the whitespace from the number of distinct recombinant sequences
  
  # Collect variables into a nice named vector and output them
  threeseq_results <- c(num_trips, num_dis, seq_sig)
  names(threeseq_results) <- c("num_distinct_recombinant_triplets","num_distinct_recombinant_sequences","3SEQ_p_value")
  return(threeseq_results)
}



# Run the PhiPack program on one alignment and return the key variables calculated by the program
run.phipack <- function(alignment_path, alignment_folder, phipack_path, seqtype){
  # Set working directory to alignment folder, so that PHI output files are saved correctly
  setwd(alignment_folder)
  
  # Check whether PHIPack has already been run
  if (file.exists(paste0(alignment_folder,"Phi.log")) == FALSE){
    
    # Determine file extension
    file_extension <- tail(strsplit(alignment_path,"\\.")[[1]],1)
    
    # Check file extension of alignment
    if (file_extension == "fasta" | file_extension == "fa" | file_extension == "fna" | file_extension == "ffn" | 
        file_extension == "faa" | file_extension == "frn" | file_extension == "fas"){
      # 3SEQ only reads in phylip or fasta format - if the alignment is already in that format, call 3SEQ on the alignment
      if (seqtype == "dna"){
        phi_command <- paste0(phipack_path, " -f ", alignment_path, " -v -o -p") # assemble system command
      } else if (seqtype == "protein"){
        phi_command <- paste0(phipack_path, " -f ", alignment_path, " -v -o -p -t A") # assemble system command
      }
      system(phi_command) #call PHIPack
    } else if (file_extension == "nex" | file_extension == "nexus"){
      # PHIPack only reads in Phylip or fasta format - need to convert if the alignment is a nexus file (using the nexus data opened above)
      # Assemble a name for a copy of the alignment as a fasta file
      fasta.name <- gsub(file_extension, ".fasta", alignment_path)
      # Check if the fasta file already exists
      if (file.exists(fasta.name) == FALSE){
        # If the fasta version doesn't exist, write out the nexus sequence in fasta formt
        n <- read.nexus.data(alignment_path)
        write.fasta(sequences = n, names = names(n), file.out = fasta.name) # output alignment as a fasta format
      }
      # There is now a definitely a fasta format version of this alignment
      # Assemble the 3SEQ command using the new fasta alignment
      if (seqtype == "dna"){
        phi_command <- paste0(phipack_path, " -f ", alignment_path, " -v -o -p") # assemble system command
      } else if (seqtype == "protein"){
        phi_command <- paste0(phipack_path, " -f ", alignment_path, " -v -o -p -t A") # assemble system command
      }
      # Call 3SEQ
      system(phi_command)
    }
  }
  
  # Extract significance from Phi Pack output
  phi_file <- paste0(alignment_folder,"Phi.log")
  phi_file <- readLines(phi_file)
  # Collect p-values
  ind      <- grep("p-Value",phi_file)
  NSS_ind <- ind + 3
  NSS_sig <- as.numeric(strsplit(strsplit(phi_file[NSS_ind], ":")[[1]][2], "\\(")[[1]][1])
  maxchi_ind <- ind + 4
  maxchi_sig <- as.numeric(strsplit(strsplit(phi_file[maxchi_ind], ":")[[1]][2], "\\(")[[1]][1])
  phi_permutation_ind <- ind + 5
  phi_permutation_sig <- as.numeric(strsplit(strsplit(phi_file[phi_permutation_ind], ":")[[1]][2], "\\(")[[1]][1])
  phi_normal_ind <- ind + 6
  phi_normal_sig <- as.numeric(strsplit(phi_file[phi_normal_ind], ":")[[1]][2])
  # Collect PHI values
  ind <- grep("PHI Values", phi_file)
  mean_ind <- ind + 4
  mean_vals <- strsplit(strsplit(phi_file[mean_ind],":")[[1]][2], "      ")[[1]]
  mean_vals <- as.numeric(mean_vals[mean_vals != ""])
  var_ind <- ind + 5
  var_vals <- strsplit(strsplit(phi_file[var_ind],":")[[1]][2], "      ")[[1]]
  var_vals <- as.numeric(var_vals[var_vals != ""])
  obs_ind <- ind + 6
  obs_vals <- strsplit(strsplit(phi_file[obs_ind],":")[[1]][2], "      ")[[1]]
  obs_vals <- as.numeric(obs_vals[obs_vals != ""])
  
  # Collate output 
  phi_results <- c(mean_vals, var_vals, obs_vals, NSS_sig, maxchi_sig, phi_permutation_sig, phi_normal_sig)
  names(phi_results) <- c("analytical_PHI_mean_value", "permutation_PHI_mean_value",
                          "analytical_PHI_variance_value", "permutation_PHI_variance_value",
                          "analytical_PHI_observed_value", "permutation_PHI_observed_value",
                          "NSS_p_value", "max_chi_squared_p_value", 
                          "PHI_permutation_p_value", "PHI_normal_p_value")
  
  # Return PHIPack results
  return(phi_results)
}



# Run the GeneConv program on one alignment and return the key variables calculated by the program
run.geneconv <- function(alignment_path, alignment_folder, geneconv_path, seqtype){
  # Set working directory to alignment folder, so that PHI output files are saved correctly
  setwd(alignment_folder)
  
  # Determine file extension
  file_extension <- tail(strsplit(alignment_path,"\\.")[[1]],1)
  
  # Check whether GeneConv has already been run
  if (file.exists(paste0(gsub(file_extension, "", alignment_path), "frags")) == FALSE){
    #generate seed by taking the numbers from today's date
    seed <- gsub("-", "", Sys.Date())
    # Assemble and call the system command
    if (seqtype == "dna"){
      geneconv_command <- paste0(geneconv_path, " ", basename(alignment_path), " /w", seed, " /lp /sp /sb -ExpFormat -Nolog > ",
                                 gsub(paste0(".", file_extension), "_geneconv_terminal.txt", basename(alignment_path)))
    } else if (seqtype == "protein"){
      geneconv_command <- paste0(geneconv_path, " -Seqfile=", basename(alignment_path), " /w", seed, " /p /lp /sp /sb -ExpFormat -Nolog > ",
                                 gsub(paste0(".", file_extension), "_geneconv_terminal.txt", basename(alignment_path)))
    }
    system(geneconv_command)
  }
  
  # Open geneconv outputted terminal text file to extract results
  geneconv_file_path <- paste0(alignment_folder, grep("_terminal_text.fas", list.files(alignment_folder), value = TRUE))
  geneconv_file <- readLines(geneconv_file_path)
  # Check whether geneconv ran successfully
  check_ind <- grep("Only one polymorphism: Too few to analyze!", geneconv_file)
  # If check_ind returns as integer(0), that means geneconv could not run (due to too few polymorphisms)
  # If check_ind returns as a number, that means geneconv ran successfully
  if (identical(integer(0), check_ind) == TRUE){
    # Open geneconv .frag file to extract results
    geneconv_file_path <- paste0(alignment_folder, grep(".frag", list.files(alignment_folder), value = TRUE))
    geneconv_file <- readLines(geneconv_file_path)
    # Extract seed
    ind <- grep("# The starting random number seed is", geneconv_file)
    seed_line <- geneconv_file[ind]
    starting_seed <- gsub(" ", "", gsub("\\.", "", gsub("# The starting random number seed is", "", seed_line)))
    # Extract number of significant fragments
    ind <- grep("Global lists:", geneconv_file)
    global_line <- geneconv_file[ind]
    global_line <- strsplit(gsub(" significant fragments", "", gsub("# Global lists:", "", global_line)), " and ")[[1]]
    global_inner_pair <- gsub(" ", "", gsub(" I", "", global_line[1]))
    if (global_inner_pair == "no"){global_inner_pair = 0}
    global_outer_seq <- gsub(" ", "", gsub(" O", "", global_line[2]))
    if (global_outer_seq == "no"){global_outer_seq = 0}
    pairwise_line <- geneconv_file[ind+1]
    pairwise_line <- strsplit(gsub(" significant fragments", "", gsub("# Pairwise lists:", "", pairwise_line)), " and ")[[1]]
    pairwise_inner_pair <- gsub(" ", "", gsub(" I", "", pairwise_line[1]))
    if (pairwise_inner_pair == "no"){pairwise_inner_pair = 0}
    pairwise_outer_seq <- gsub(" ", "", gsub(" O", "", pairwise_line[2]))
    if (pairwise_outer_seq == "no"){pairwise_outer_seq = 0}
    # Extract minimum p-values
    ind <- grep("# Simulated P-values are based on ", geneconv_file)
    inner_line <- geneconv_file[ind+4]
    inner_line_split <- strsplit(inner_line, "    ")[[1]]
    inner_line_split <- inner_line_split[inner_line_split != ""]
    inner_line_vars <- gsub(" ", "", inner_line_split[2:5])
    inner_line_var_names <- c("geneconv_inner_fragment_maximum_blast-like_score", "geneconv_inner_fragment_simulated_p_value", "geneconv_inner_fragment_sd_above_sim_mean", "geneconv_inner_fragment_sd_of_sim")
    outer_line <- geneconv_file[ind+6]
    outer_line_split <- strsplit(outer_line, "    ")[[1]]
    outer_line_split <- outer_line_split[outer_line_split != ""]
    outer_line_vars <- gsub(" ", "", outer_line_split[2:5])
    outer_line_var_names <- c("geneconv_outer_fragment_maximum_blast-like_score", "geneconv_outer_fragment_simulated_p_value", "geneconv_outer_fragment_sd_above_sim_mean", "geneconv_outer_fragment_sd_of_sim")
    # Assemble geneconv results
    geneconv_results <- c(starting_seed, global_inner_pair, global_outer_seq, pairwise_inner_pair, pairwise_outer_seq, inner_line_vars, outer_line_vars)
    names(geneconv_results) <- c("geneconv_seed", "geneconv_num_global_inner_fragments", "geneconv_num_global_outer-sequence_fragments", "geneconv_num_pairwise_inner_fragments", 
                                 "geneconv_num_pairwise_outer-sequence_fragments", inner_line_var_names, outer_line_var_names) 
  } else if (identical(integer(0), check_ind) == FALSE){
    # Geneconv could not run - return NA for all variables
    geneconv_results <- rep(NA, 13)
    # Assemble names
    inner_line_var_names <- c("geneconv_inner_fragment_maximum_blast-like_score", "geneconv_inner_fragment_simulated_p_value", "geneconv_inner_fragment_sd_above_sim_mean", "geneconv_inner_fragment_sd_of_sim")
    outer_line_var_names <- c("geneconv_outer_fragment_maximum_blast-like_score", "geneconv_outer_fragment_simulated_p_value", "geneconv_outer_fragment_sd_above_sim_mean", "geneconv_outer_fragment_sd_of_sim")
    names(geneconv_results) <- c("geneconv_seed", "geneconv_num_global_inner_fragments", "geneconv_num_global_outer-sequence_fragments", "geneconv_num_pairwise_inner_fragments", 
                                 "geneconv_num_pairwise_outer-sequence_fragments", inner_line_var_names, outer_line_var_names)
  }
  # Return geneconv results
  return(geneconv_results)
}



