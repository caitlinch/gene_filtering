### empirical_treelikeness/code/func_empirical.R
## R functions to process and modify empirical sequence alignments, specifically Rob Lanfear's "BenchmarkAlignments"
## BenchmarkAlignments and metadata available here: https://github.com/roblanf/BenchmarkAlignments and have CC0 or CCBY licenses
# Caitlin Cherryh 2021

# Packages required 
library(phytools)
library(ape)
library(phangorn)
library(parallel)
library(seqinr)

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





# Functions to run test statistics on empirical datasets
empirical.runTS <- function(alignment_path, program_paths, bootstrap_id, iqtree.num_threads, iqtree.num_quartets, iqtree.model = "MFP", alphabet, dataset_name){
  # Want to only do the tree proportion, 3SEQ and sCF
  print("in empirical.runTS")
  print(alignment_path)
  # extract the alignment folder from the alignment path
  alignment_folder <- paste0(dirname(alignment_path),"/")
  output_id <- gsub(".nex","",basename(alignment_path))
  output_id <- gsub(".fasta","",output_id)
  
  # Make the name for the testStatistics output file for this alignment
  results_file <- paste0(alignment_folder,output_id,"_testStatistics.csv")
  # Check if the testStatistics.csv file for this alignment already exists
  if (file.exists(results_file) == FALSE){
    # If the testStatistics.csv doesn't exist, run the test statistics and collate the results
    
    # Create some folder and filenames
    if (bootstrap_id == "alignment"){
      # Get the alignment name and remove the extension to get the loci name
      loci_name <- output_id
      # Extract the dataset name (basename of alignment folder: element after last "/" in alignment_folder)
      dataset <- dataset_name
    } else {
      # If the alignment is a bootstrap replicate, need to remove the bootstrap rep number to get the loci name
      loci_name <- output_id # get the basis of the loci name
      loci_list <- unlist(strsplit(loci_name, "_")) # break the alignment name into bits
      max_ind <- grep("bootstrapReplicate",loci_list) - 1 # which ind is the bootstrapReplicate at?
      loci_list <- loci_list[1:max_ind] # get only parts of alignment name
      loci_name <- paste(loci_list,collapse="_") # squash the alignment name together
      # Extract the dataset name (basename of alignment folder: element after second-last "/" in alignment_folder)
      dataset <- dataset_name
    }
    
    if (bootstrap_id == "alignment"){
      # If this is not a bootstrap replicate, you need to create a folder to store the program logs in
      # Otherwise they will get overwritten for each locus you run - this means you can keep them for late
      log_folder <- paste0(alignment_folder,loci_name,"/")
      # make a rep id - e.g. an id to store in output df so you can identify what's a simulation and what's from empirical data
      rep_id <- loci_name
    } else {
      # If the run is a bootstrap replicate, you can just save the information in its folder (only 1 alignment per folder in bootstrap replicate folders)
      log_folder <- alignment_folder
      # make a rep id - e.g. an id to store in output df so you can identify what's a simulation and what's from empirical data
      rep_id <- bootstrap_id
    }
    
    
    # If the log file doesn't exist, create it 
    if (dir.exists(log_folder) == FALSE){
      dir.create(log_folder) # create a new folder to store the log files from the executables for this loci in
    }
    
    # Open the nexus file and get the number of taxa and the number of characters 
    file_name_list <- strsplit(basename(alignment_path),"\\.")[[1]]
    file_extension <- file_name_list[length(file_name_list)]
    
    if (file_extension == "nex" | file_extension == "nexus"){
      n <- read.nexus.data(alignment_path)
      n_taxa <- length(n)
      n_char <- length(unlist(n[1]))
    } else if (file_extension == "fasta" | file_extension == "fa" | file_extension == "fna" | file_extension == "ffn" | file_extension == "faa" | file_extension == "frn"){
      f <- read.fasta(file = alignment_path)
      n_taxa <- length(f)
      n_char <- length(unlist(f[1]))
    }
    
    
    # Run IQ-tree on the alignment (if it hasn't already been run)
    if (file.exists(paste0(alignment_path, ".iqtree")) == FALSE){
      call.IQTREE.empirical(alignment_path, iqtree_path = program_paths["IQTree"], iqtree.model, num_threads = iqtree.num_threads) 
    }
    initial_iqtree_tree <- paste0(alignment_path,".treefile")
    
    # Change to the log (storage for log files) folder for this alignment - means that 3seq and Phi files will be saved into a unique folder
    setwd(log_folder)
    
    # Only run 3SEQ and collect sCF values if the path is an alignment (don't need bootstrap replicates for 3SEQ/sCF)
    if (bootstrap_id == "alignment"){
      # Check if 3SEQ has already been run. Only run 3SEQ if the log file doesn't exist
      if (file.exists(paste0(log_folder,"3s.log")) == FALSE){
        if (file_extension == "fasta" | file_extension == "fa" | file_extension == "fna" | file_extension == "ffn" | file_extension == "faa" | file_extension == "frn"){
          # 3SEQ only reads in phylip or fasta format - if the alignment is already in that format, call 3SEQ on the alignment
          seq_command <- paste0(program_paths[["3seq"]]," -f ", alignment_path)
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
          seq_command <- paste0(program_paths[["3seq"]], " -f ", fasta.name)
          # Call 3SEQ
          system(seq_command)
        }
      }
      # Now, collect the results from the 3SEQ log file
      seq_file <- paste0(log_folder,"3s.log")
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
      # record proportion of recombinant sequences
      prop_recomb_seq <- as.numeric(num_dis)/as.numeric(n_taxa)
      
      # Calculate the site concordance factor results
      # Run IQ-tree on the alignment (if it hasn't already been run), and get the site concordance factor results
      sCF <- calculate.empirical.sCF(alignment_path, iqtree_path = program_paths["IQTree"], iqtree.model, num_threads = iqtree.num_threads, num_scf_quartets = iqtree.num_quartets)
      scf_mean_value <- sCF$mean_scf
      scf_median_value <- sCF$median_scf
    } else {
      # If this is a bootstrap replicate, don't need either the sCF or the 3SEQ values
      # These are recorded from the original alignment
      num_trips <- "NA"
      num_dis <- "NA"
      seq_sig <- "NA"
      prop_recomb_seq <- "NA"
      scf_mean_value <- "NA"
      scf_median_value <- "NA"
    }
    
    # Change back to directory containing alignments and iqtree files
    setwd(alignment_folder)
    
    # Run trimmed version of the NeighborNet tree proportion
    # Call the test statistic functions
    initial_iqtree_tree <- paste0(alignment_path,".treefile")
    nn_trimmed <- tree.proportion(iqpath = program_paths[["IQTree"]], splitstree_path = program_paths[["SplitsTree"]],
                                  path = alignment_path, network_algorithm = "neighbournet", trimmed = TRUE,
                                  tree_path = initial_iqtree_tree, run_IQTREE = FALSE, seq_type = alphabet)
    
    
    # Name the test statistics file using the output id (this way if it's a  bootstrap replicate, it adds the replicate number!)
    print(paste0("output results for ",output_id))
    results_file <- paste0(alignment_folder,output_id,"_testStatistics.csv")
    # Make somewhere to store the results
    df_names <- c("dataset", "loci", "bootstrap_replicate_id", "n_taxa", "n_sites", "alignment_file",
                  "3SEQ_num_recombinant_triplets", "3SEQ_num_distinct_recombinant_sequences", "3SEQ_prop_recombinant_sequences", "3SEQ_p_value",
                  "tree_proportion","sCF_mean", "sCF_median")
    df <- data.frame(matrix(nrow=0,ncol=length(df_names))) # create an empty dataframe of the correct size
    op_row <- c(dataset, loci_name, rep_id, n_taxa, n_char, alignment_path,
                num_trips, num_dis, prop_recomb_seq, seq_sig,
                nn_trimmed, scf_mean_value, scf_median_value) # collect all the information
    df <- rbind(df,op_row,stringsAsFactors = FALSE) # place row in dataframe
    names(df) <- df_names # add names to the df so you know what's what
    write.csv(df,file = results_file, row.names = FALSE)
    
    if (bootstrap_id == "alignment"){
      # Repeat the above to create an output folder for the sCF values
      print(paste0("output sCF values for ",output_id))
      results_file <- paste0(alignment_folder,output_id,"_sCF_branch.csv")
      # Make somewhere to store the results
      op_row <- c(dataset,loci_name,rep_id,n_taxa,n_char,alignment_path,sCF$all_scfs) # collect all the information
      df <- data.frame(matrix(nrow=0,ncol=length(op_row))) # create an empty dataframe of the correct size
      df_names <- c("dataset","loci","bootstrap_replicate_id","n_taxa","n_sites","alignment_file", sCF$branch_ids)
      df <- rbind(df,op_row,stringsAsFactors = FALSE) # place row in dataframe
      names(df) <- df_names # add names to the df so you know what's what
      write.csv(df,file = results_file, row.names = FALSE)
    }
  }
}





# Function to calculate one empirical bootstrap for an empirical alignment
do1.empirical.parametric.bootstrap <- function(bootstrap_id, empirical_alignment_path, empirical_alignment_row, alignment_params, program_paths, iqtree.num_threads, iqtree.num_quartets){
  print("in do1.empirical.parametric.bootstrap")
  print(empirical_alignment_path)
  print(bootstrap_id)
  # Create the folder for this replicate, gather and create filenames
  loci_name <- empirical_alignment_row$loci_name
  bootstrap_name <- paste0(loci_name,"_",bootstrap_id) # this will be the name of the alignment
  bootstrap_folder <- paste0(dirname(empirical_alignment_path),"/",bootstrap_name,"/") # folder to store results from this bootstrap in
  # If the bootstrap folder doesn't exist, create it
  if (dir.exists(bootstrap_folder) == FALSE){
    dir.create(bootstrap_folder)
  }
  
  # Get the file extension of the empirical alignment path
  file_name_list <- strsplit(basename(empirical_alignment_path),"\\.")[[1]]
  file_extension <- file_name_list[length(file_name_list)]
  
  # Make the file names for the files created by this function
  shuffled_alignment_path <- paste0(bootstrap_folder,bootstrap_name,"_shuffled_noGaps.",file_extension)
  bootstrap_alignment_path <- paste0(bootstrap_folder,bootstrap_name,".",file_extension)
  empirical_alignment_tree_path <- paste0(empirical_alignment_path,".treefile")
  empirical_alignment_tree <- read.tree(empirical_alignment_tree_path)
  
  # First generate completely new DNA using the params
  # Create an alignment for this replicate using the alignment params - name will be loci_bootstrapReplicateXXXX
  if (file.exists(bootstrap_alignment_path) == FALSE) {
    if (empirical_alignment_row$alphabet == "dna"){
      new_aln <- generate.DNA.alignment(alignment_params, empirical_alignment_tree)
    } else if (empirical_alignment_row$alphabet == "protein"){
      new_aln <- generate.AA.alignment(alignment_params, empirical_alignment_tree)
    }
    
    # Second, randomise sites in new_aln (otherwise masking will mean the gamma categories disproportionately get affected)
    print("shuffling alignment")
    #   Turn the phyDat into a dataframe and shuffle the rows using sample()
    #   Rows are shuffled as columns represent sequences - shuffling the rows keeps the relationships between species (shown by rows) 
    #   but rearranges the order the relationships occur in
    new_aln_df <- as.data.frame(new_aln) # turn the new alignment into a dataframe
    new_aln_df <- new_aln_df[sample(nrow(new_aln_df)),] # sample the rows randomly (this will randomly distribute the gamma categories throughout)
    rownames(new_aln_df) <- 1:nrow(new_aln_df) # reset the row names as 1:nrows (they got mixed up when the sampling occurred)
    # Turn the dataframe back into an alignment:
    if (file_extension == "nex" | file_extension == "nexus"){
      if (empirical_alignment_row$alphabet == "dna"){
        new_aln_shuffled <- as.phyDat(new_aln_df, type = "DNA")
        #write the shuffled alignment to disk
        write.phyDat(new_aln_shuffled, file = shuffled_alignment_path, format = "nexus", interleave = TRUE)
        # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
        nexus_edit <- readLines(shuffled_alignment_path) # open the new nexus file
        ind <- grep("BEGIN DATA",nexus_edit)+2 # find which line
        nexus_edit[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE;" # replace the line
        writeLines(nexus_edit,shuffled_alignment_path) # output the edited nexus file
      } else if (empirical_alignment_row$alphabet == "protein"){
        new_aln_shuffled <- as.phyDat(new_aln_df, type = "AA")
        #write the shuffled alignment to disk
        write.phyDat(new_aln_shuffled, file = shuffled_alignment_path, format = "nexus", interleave = TRUE)
        # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
        nexus_edit <- readLines(shuffled_alignment_path) # open the new nexus file
        ind <- grep("BEGIN DATA",nexus_edit)+2 # find which line
        nexus_edit[ind] <- "  FORMAT DATATYPE=PROTEIN MISSING=? GAP=- INTERLEAVE;" # replace the line
        writeLines(nexus_edit,shuffled_alignment_path) # output the edited nexus file
      }
      
    } else if (file_extension == "fasta" | file_extension == "fa" | file_extension == "fna" | 
                file_extension == "ffn" | file_extension == "faa" | file_extension == "frn"){
      if (empirical_alignment_row$alphabet == "dna"){
        new_aln_shuffled <- as.phyDat(new_aln_df, type = "DNA")
        #write the shuffled alignment to disk
        write.phyDat(new_aln_shuffled, file = shuffled_alignment_path, format = "fasta", colsep = "")
      } else if (empirical_alignment_row$alphabet == "protein"){
        new_aln_shuffled <- as.phyDat(new_aln_df, type = "AA")
        #write the shuffled alignment to disk
        write.phyDat(new_aln_shuffled, file = shuffled_alignment_path, format = "fasta", colsep = "")
      }
    }
    
    # Third, mask each alignment with the gaps and unknown characters from the original sequence 
    print("masking alignment")
    # for each alignment:
    #     - copy the sequence out from the new alignment: temp <- as.numeric(new_aln$X)
    #     - replace the non 18s with the generated sequence of the right length: temp[which(new_aln$X !=18)] <- new_seq
    #     - replace the new seq into the new aln: new_aln$X <- temp
    # If nexus file used for tree estimation, open the nexus file and recreate the gaps
    if (file_extension == "nex" | file_extension == "nexus"){
      n <- read.nexus.data(empirical_alignment_path)
      n_taxa <- length(n)
      n_char <- length(unlist(n[1]))
      # Open the new alignment as a nexus file
      n_new <- read.nexus.data(shuffled_alignment_path)
      # Get the names of all the sequences
      seq_names <- names(n_new)
      print(paste0("number of names: ",length(seq_names)))
      print(paste0("number of unique names: ",length(unique(seq_names))))
      # Iterate through the names and add the gaps from the original sequence into the bootstrap sequence
      for (seq_name in seq_names){
        original_seq <- n[[seq_name]] # get the original empirical sequence
        new_seq <- n_new[[seq_name]] # get the new simulated sequence that has the same name
        gap_inds <- which(original_seq == "-") # find out which sites are a gap in the original alignment
        unknown_inds <- which(original_seq == "?") # find out which sites are unknown in the original alignment
        new_seq[gap_inds] <- "-" # add the gaps into the simulated alignment
        new_seq[unknown_inds] <- "?" # add the unknowns into the simulated alignment
        n_new[[seq_name]] <- new_seq
      }
      # Output the final alignment (same parameters and gaps as input alignment) as a nexus file
      print("output nexus file")
      if (empirical_alignment_row$alphabet == "dna"){
        # write out the nexus with the gaps added
        write.nexus.data(n_new,file = bootstrap_alignment_path, format = "dna", interleaved = TRUE)
        # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
        nexus_edit <- readLines(bootstrap_alignment_path) # open the new nexus file
        ind <- grep("BEGIN DATA",nexus_edit)+2 # find which line
        nexus_edit[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE;" # replace the line
        writeLines(nexus_edit,bootstrap_alignment_path) # output the edited nexus file
      } else if (empirical_alignment_row$alphabet == "protein"){
        # write out the nexus with the gaps added
        write.nexus.data(n_new,file = bootstrap_alignment_path, format = "protein", interleaved = TRUE)
        # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
        nexus_edit <- readLines(bootstrap_alignment_path) # open the new nexus file
        ind <- grep("BEGIN DATA",nexus_edit)+2 # find which line
        nexus_edit[ind] <- "  FORMAT DATATYPE=PROTEIN MISSING=? GAP=- INTERLEAVE;" # replace the line
        writeLines(nexus_edit,bootstrap_alignment_path) # output the edited nexus file
      }
      
    } else if (file_extension == "fasta" | file_extension == "fa" | file_extension == "fna" |
               file_extension == "ffn" | file_extension == "faa" | file_extension == "frn"){
      # If the fasta file was used for the gene tree estimation, open the fasta file and recreate the gaps
      if (empirical_alignment_row$alphabet == "dna"){
        fasta_type = "DNA"
      } else if (empirical_alignment_row$alphabet == "protein"){
        fasta_type = "AA"
      }
      f <- read.fasta(empirical_alignment_path, seqtype = fasta_type)
      n_taxa <- length(f)
      n_char <- length(unlist(f[1]))
      # Open the new alignment as a nexus file
      f_new <- read.fasta(shuffled_alignment_path, seqtype = fasta_type)
      # Get the names of all the sequences
      seq_names <- names(f_new)
      print(paste0("number of names: ",length(seq_names)))
      print(paste0("number of unique names: ",length(unique(seq_names))))
      # Iterate through the names
      for (seq_name in seq_names){
        original_seq <- f[[seq_name]] # get the original empirical sequence
        new_seq <- f_new[[seq_name]] # get the new simulated sequence that has the same name
        gap_inds <- which(original_seq == "-") # find out which sites are a gap in the original alignment
        unknown_inds <- which(original_seq == "?") # find out which sites are unknown in the original alignment
        new_seq[gap_inds] <- "-" # add the gaps into the simulated alignment
        new_seq[unknown_inds] <- "?" # add the unknowns into the simulated alignment
        f_new[[seq_name]] <- new_seq
      }
      # Save the bootstrap alignment with the gaps added in the same place as in the original alignment
      if (empirical_alignment_row$alphabet == "dna"){
        #write the alignment with gaps to disk
        write.fasta(f_new, seq_names, file.out = bootstrap_alignment_path, )
      } else if (empirical_alignment_row$alphabet == "protein"){
        #write the alignment with gaps to disk
        write.fasta(f_new, seq_names, file.out = bootstrap_alignment_path)
      }
    }
    
  }
  # Run all the test statistics
  # bootstrap_id will be "bootstrapReplicateXXXX" where XXXX is a number
  print("run test statistics")
  
  # If the alignment comes from Vanderpool 2020, there is no specified best model (the original best model was determined using MFP in IQ-Tree)
  # Need to pass the model from the empirical alignment into the bootstrap replicates (this will be from the .iqtree file)
  # Alternatively, if the alignment comes from Misof 2014 or 1KP, we know the best model and we should use that model for both the original 
  #   alignment and the parametric bootstrap
  if (empirical_alignment_row$best_model == "MFP"){
    MFP_model <- alignment_params$parameters[9,2]
    empirical.runTS(alignment_path = bootstrap_alignment_path, program_paths = program_paths, bootstrap_id = bootstrap_id, 
                    iqtree.num_threads, iqtree.num_quartets, iqtree.model = MFP_model, 
                    alphabet = empirical_alignment_row$alphabet, dataset_name = empirical_alignment_row$dataset)
  } else {
    empirical.runTS(alignment_path = bootstrap_alignment_path, program_paths = program_paths, bootstrap_id = bootstrap_id, 
                    iqtree.num_threads, iqtree.num_quartets, iqtree.model = empirical_alignment_row$best_model, 
                    alphabet = empirical_alignment_row$alphabet, dataset_name = empirical_alignment_row$dataset)
    
    
  }
}





# Function to act as outer shell for all the layers in performing test statistics and parametric bootstraps
empirical.bootstraps.wrapper <- function(loci_number, loci_df, program_paths, number_of_replicates, iqtree.num_threads, iqtree.num_quartets, num_of_cores){
  print("in empirical.bootstraps.wrapper")
  # Extract loci name and path from loci_df using the row number
  loci_row <- loci_df[loci_number,]
  loci_name <- loci_row$loci_name
  empirical_alignment_path <- loci_row$loci_path
  print(paste0("Number = ",loci_number," (of ",nrow(loci_df),"), Dataset = ",loci_row$dataset,", loci = ",loci_name))
  # Create a folder to store the files for this loci in
  alignment_folder <- paste0(loci_row$output_folder,loci_name,"/")
  # Check if this folder exists - and if it doesn't, create it!
  if (dir.exists(alignment_folder) == FALSE){
    dir.create(alignment_folder)
  }
  # Set this as the working directory
  setwd(alignment_folder)
  # Create output file names, the name of the loci and the file path of the loci location
  collated_ts_file <- paste0(alignment_folder,loci_name,"_testStatistics_collatedBSReplicates.csv")
  collated_sCF_file <- paste0(alignment_folder,loci_name,"_branchSCF_collatedBSReplicates.csv")
  p_value_file  <- paste0(alignment_folder,loci_name,"_pValues.csv")
  # Create file names for saving the IQTree parameters (e.g. rate matrix, model of sequence evolution, gamma categories)
  parameters_file <- paste0(alignment_folder,loci_name,"_parameterValues.csv")
  gamma_categories_file <- paste0(alignment_folder,loci_name,"_gammaCategories.csv")
  aa_frequency_file <- paste0(alignment_folder,loci_name,"_amino_acid_frequencies.csv")
  rate_matrix_file <- paste0(alignment_folder,loci_name,"_QRateMatrix.csv")
  ts_file <- paste0(alignment_folder,loci_name,"_testStatistics.csv")
  
  # Only run this section if the p-value csv has not been created yet (skip reruns)
  if (file.exists(p_value_file) == FALSE){
    # If alignment is a nexus, copy it into the alignment folder
    # If it's not, rewrite it into a nexus and copy to the alignment folder
    # From now on, use the copy for analysis (leaving the original untouched)
    empirical_alignment_path <- copy.alignment.as.nexus(loci_row$loci_path, alignment_folder, loci_name, loci_row)
    
    # Remove empty taxa from the alignment
    remove.empty.taxa(empirical_alignment_path, loci_row$alphabet)
    # Remove any sequences from the alignment that have more than the allowable proportion of ambiguous/missing sites or gaps
    if (!is.na(loci_row$allowable_proportion_missing_sites)){
      # If the allowable proportion of missing sites is NOT NA, run the pruning function
      # The pruning function will remove any sequences where the proportion of missing sites is GREATER than the allowable proportion of missing sites
      prune.taxa.by.length(empirical_alignment_path, loci_row$allowable_proportion_missing_sites, loci_row$alphabet, write_output_text = TRUE, alignment_folder)
    }
    
    # If the alignment file is a nexus file, make sure that it doesn't contain any characters that aren't allowed in the format (e.g. R and Y for nexus DNA files)
    # If it is a fasta file, ignore this setion and just fun the fasta file
    file_type <- tail(strsplit(empirical_alignment_path,"\\.")[[1]],1)
    if (file_type == "nex" || file_type == "nexus" || file_type == "nxs"){
      invalid_character_check <- check.invalid.nexus.characters(empirical_alignment_path, loci_row$alphabet)
    } else if (file_type == "fasta" || file_type == "faa" || file_type == "fna" || 
               file_type == "fa" || file_type == "ffn" || file_type == "frn" || file_type == "fas") {
      invalid_character_check <- "use_FASTA"
    }
    
    # Check whether IQ-Tree has been run
    if (file.exists(paste0(empirical_alignment_path,".treefile.cf.stat")) == FALSE){
      print("need to run IQ-Tree")
      # Run IQ-tree on the alignment (if it hasn't already been run), and get the sCF results
      print("run IQTree and estimate sCFs")
      if (invalid_character_check == "use_FASTA"){
        if (file_type == "fasta" || file_type == "faa" || file_type == "fna" || 
            file_type == "fa" || file_type == "ffn" || file_type == "frn" || file_type == "fas"){
          # If the file is already a fasta, call the function using the fasta file
          scfs <- calculate.empirical.sCF(alignment_path = empirical_alignment_path, iqtree_path = program_paths[["IQTree"]], 
                                          alignment_model = loci_row$best_model, num_threads = iqtree.num_threads, 
                                          num_scf_quartets = iqtree.num_quartets)
        } else if (file_type == "nex" || file_type == "nexus" || file_type == "nxs"){
          # If the nexus file contains ambiguous characters that aren't permitted by the format, run IQ-Tree using the FASTA file
          # Change the nexus file extension to a fasta file extension
          empirical_alignment_path <- gsub(".nex",".fasta",empirical_alignment_path)
          # Call the function using the newly selected fasta file
          scfs <- calculate.empirical.sCF(alignment_path = empirical_alignment_path, iqtree_path = program_paths[["IQTree"]], 
                                          alignment_model = loci_row$best_model, num_threads = iqtree.num_threads, 
                                          num_scf_quartets = iqtree.num_quartets)
        }
      } else {
        # If the nexus file doesn't contain any ambiguous characters, use the nexus file to run IQ-Tree
        scfs <- calculate.empirical.sCF(alignment_path = empirical_alignment_path, iqtree_path = program_paths[["IQTree"]], 
                                        alignment_model = loci_row$best_model, num_threads = iqtree.num_threads, 
                                        num_scf_quartets = iqtree.num_quartets)
      }
    }
    
    # Calculate the test statistics if it hasn't already been done
    if (file.exists(ts_file) == FALSE){
      print("run test statistics")
      empirical.runTS(empirical_alignment_path, program_paths, bootstrap_id = "alignment", iqtree.num_threads, iqtree.num_quartets, 
                      iqtree.model = loci_row$best_model, alphabet = loci_row$alphabet, dataset_name = loci_row$dataset)
    }
    
    #Check that the test statistic file ran ok 
    if (file.exists(ts_file) == FALSE){
      print("need to rerun test statistics")
      empirical.runTS(empirical_alignment_path, program_paths, bootstrap_id = "alignment", iqtree.num_threads, iqtree.num_quartets, 
                      iqtree.model = loci_row$best_model, alphabet = loci_row$alphabet, dataset_name = loci_row$dataset)
    }
    
    #Extract the parameters from the .iqtree log file.
    print("get simulation params")
    params <- get.simulation.parameters(paste0(empirical_alignment_path,".iqtree"))
    write.csv(params$parameters,file = parameters_file, row.names = TRUE)
    write.csv(params$gamma_categories,file = gamma_categories_file, row.names = TRUE)
    if (loci_row$alphabet == "dna") {
      write.csv(params$Q_rate_matrix,file = rate_matrix_file, row.names = TRUE)
    } else if (loci_row$alphabet == "protein"){
      write.csv(params$frequency,file = aa_frequency_file, row.names = TRUE)
    }
    
    
    # Create the bootstrap ids (pad out to 4 digits) - should be "bootstrapReplicateXXXX" where XXXX is a number
    bootstrap_ids <- paste0("bootstrapReplicate",sprintf("%04d",1:number_of_replicates))
    
    # Run all the bootstrap ids that HAVEN'T already been run (e.g. in previous attempts) using lapply (feed info into do1.empirical.parametric.bootstrap)
    # If the alignment doesn't exist OR the test statistic csv doesnt exist, this indicates a complete run has not previously been done 
    # These bootstrap replicates will thus be calculates
    # This should save A BUNCH of time because it means if the test statistic file exists, you don't have to run any calculations
    print(paste0("Prepare ",number_of_replicates," bootstrap replicates"))
    bs_als <- paste0(alignment_folder,loci_name,"_",bootstrap_ids,"/",loci_name,"_",bootstrap_ids,".nex")
    ts_csvs <- paste0(alignment_folder,loci_name,"_",bootstrap_ids,"/",loci_name,"_",bootstrap_ids,"_testStatistics.csv")
    missing_als <- bs_als[!file.exists(bs_als)]
    missing_testStatistics <- bs_als[!file.exists(ts_csvs)]
    # Collate the missing files and identify the alignments to rerun
    all_to_run <- unique(c(missing_als,missing_testStatistics))
    ids_to_run <- bootstrap_ids[which((bs_als %in% all_to_run))]
    print(paste0("Run ",length(ids_to_run)," (of ", number_of_replicates, ") bootstrap replicates"))
    if(length(ids_to_run)>0){
      mclapply(ids_to_run, do1.empirical.parametric.bootstrap, empirical_alignment_path = empirical_alignment_path, empirical_alignment_row = loci_row,
               alignment_params = params, program_paths = program_paths, iqtree.num_threads = iqtree.num_threads, iqtree.num_quartets = iqtree.num_quartets, 
               mc.cores = num_of_cores)
    }
    
    # Before you can collate all the bootstrap files, you need to check every bootstrap ran and rerun the failed ones
    # Generate the names of each alignment, the test statistics csvs, the .iqtree files, the treefiles, the likelihood mapping files
    # Check which of these files are missing
    print("check for missing alignments")
    bs_als <- paste0(alignment_folder,"/",loci_name,"_",bootstrap_ids,"/",loci_name,"_",bootstrap_ids,".nex")
    ts_csvs <- paste0(alignment_folder,"/",loci_name,"_",bootstrap_ids,"/",loci_name,"_",bootstrap_ids,"_testStatistics.csv")
    missing_als <- bs_als[!file.exists(bs_als)]
    missing_iqtree <- bs_als[!file.exists(paste0(bs_als,".iqtree"))]
    missing_tree <- bs_als[!file.exists(paste0(bs_als,".treefile"))]
    missing_testStatistics <- bs_als[!file.exists(ts_csvs)]
    # Collate the missing files and identify the alignments to rerun
    all_missing <- unique(c(missing_als,missing_iqtree,missing_tree,missing_testStatistics))
    als_to_rerun <- bootstrap_ids[which((bs_als %in% all_missing))]
    print(paste0("Number of missing alignments to rerun = ",length(als_to_rerun)))
    # Rerun the missing als
    if (length(als_to_rerun)>0){
      mclapply(ids_to_run, do1.empirical.parametric.bootstrap, empirical_alignment_path = empirical_alignment_path, empirical_alignment_row = loci_row,
               alignment_params = params, program_paths = program_paths, iqtree.num_threads = iqtree.num_threads, iqtree.num_quartets = iqtree.num_quartets, 
               mc.cores = num_of_cores)
    }
    
    # collate the sCF by branch distributions bootstrap info into 1 file and write it to disk
    print("collate sCF from bootstraps")
    collated_scf_branch_df <- collate.bootstraps(directory = alignment_folder, file.name = "sCF_branch", id = loci_name, output.file.name = collated_sCF_file)
    write.csv(collated_scf_branch_df, file = collated_sCF_file, row.names = FALSE)
    
    # collate the test statistics bootstrap info into 1 file
    print("collate test statistics from bootstraps")
    p_value_df <- collate.bootstraps(directory = alignment_folder, file.name = "testStatistics", id = loci_name, output.file.name = collated_ts_file)
    # add the column with the bootstrap replicates and "alignment"
    new_bootstrap_ids <- p_value_df$bootstrap_replicate_id # copy col
    aln_id <- grep(new_bootstrap_ids[!grepl("bootstrapReplicate",new_bootstrap_ids)],new_bootstrap_ids) # get which element of col is the alignment
    new_bootstrap_ids[aln_id] <- "alignment"
    p_value_df$bootstrap_replicate_id <- new_bootstrap_ids
    # Get the single row that's just the alignment
    op_p_value_df <- p_value_df[p_value_df$bootstrap_replicate_id == "alignment",]
    
    # Calculate the p-values and add them to the original test statistic dataframe
    print("calculate p values")
    # Calculate the p_values of the variables of interest
    # Calculate the p-values for each test statistic
    op_p_value_df$tree_proportion_p_value         <- calculate.p_value(p_value_df$tree_proportion, p_value_df$bootstrap_replicate_id)
    op_p_value_df <- op_p_value_df[,c('dataset','loci','bootstrap_replicate_id','n_taxa','n_sites','alignment_file','X3SEQ_num_recombinant_triplets',
                                      'X3SEQ_num_distinct_recombinant_sequences','X3SEQ_prop_recombinant_sequences','X3SEQ_p_value','tree_proportion',
                                      'tree_proportion_p_value','sCF_mean','sCF_median')]
    # Output the p-values file
    write.csv(op_p_value_df,file = p_value_file, row.names = FALSE)
  }
}





# Function to call IQ-tree and estimate a maximum likelihood tree 
call.IQTREE.empirical <- function(alignment_path, iqtree_path, alignment_model, num_threads = "AUTO"){
  # Check if the tree file already exists and if it doesn't, run IQ-tree and create it
  if (file.exists(paste0(alignment_path,".treefile")) == FALSE){
    # Given an alignment, estimate the maximum likelihood tree
    # to estimate: iqtree -s ALN_FILE -p PARTITION_FILE --prefix concat -bb 1000 -nt AUTO
    call <- paste0(iqtree_path," -s ",alignment_path," -m ",alignment_model," -nt ", num_threads," -redo -safe")
    system(call)
  }
}





# Function to call IQ-tree, and estimate a maximum likelihood tree and corresponding the site concordance factors
# Site concordance factors (sCF) are the fraction of decisive alignmen sites supporting that branch
# sCF Citation: Minh B.Q., Hahn M., Lanfear R. (2018) New methods to calculate concordance factors for phylogenomic datasets. https://doi.org/10.1101/487801
calculate.empirical.sCF <- function(alignment_path, iqtree_path, alignment_model, num_threads = "AUTO", num_scf_quartets = 100){
  # Check if the tree file already exists and if it doesn't, run IQ-tree and create it
  if (file.exists(paste0(alignment_path,".treefile")) == FALSE){
    call.IQTREE.empirical(alignment_path, iqtree_path, alignment_model, num_threads)
  }
  if (file.exists(paste0(alignment_path,".treefile.cf.stat")) == FALSE){
    # Create the command and call it in the system
    # for sCF: iqtree -t concat.treefile -s ALN_FILE --scf 100 --prefix concord -nt 10
    treefile <- paste0(alignment_path,".treefile")
    call <- paste0(iqtree_path," -t ",treefile," -s ",alignment_path," --scf ",num_scf_quartets," -nt ","1"," -redo -safe")
    system(call) # call IQ-tree!
  }
  # retrieve the sCF from the output
  scf_results <- extract.sCF.results(alignment_path)
  return(scf_results)
}





# This function takes a list of names and a file to a partitioning scheme and returns a vector of models
# (one per name, in the same order as the input name vector)
# This small function calls get.one.model below
model.from.partition.scheme <- function(names,model_path,dataset){
  # Open the file
  model_file <- readLines(model_path)
  # Use lapply to get the model of evolution for each loci
  if (dataset == "Misof2014"){
    all_m <- unlist(lapply(names, get.Misof.model, char_lines = model_file))
  } else if (dataset == "1KP"){
    all_m <- unlist(lapply(names, get.OKP.model, char_lines = model_file))
  }
  
  # Return the list of models of evolution
  return(all_m)
}

# This function looks up a single loci name in a partition file and returns the associated model of evolution
# This small function is called by model.from.partition.scheme above (using lapply)
get.Misof.model <- function(name,char_lines){
  # name is a loci name
  # char_lines is a .txt file opened using readLines
  name_ind <- grep(name,char_lines)
  if (identical(name_ind,integer(0))){
    # If there's no matches in the file for this name, return NA
    m = NA
  } else {
    # If there is a match for this name in the file, find the match
    line <- char_lines[name_ind]
    # Break the line into columns using the dividers ("|")
    split_line <- strsplit(line,"\\|")[[1]]
    # Take the second element - this will be the model
    m <- split_line[2]
    # Trim the white space
    m <- gsub(" ","",m)
  }
  # Return the model
  return(m)
}

# This function looks up a single loci name in a partition file and returns the associated model of evolution
# This small function is called by model.from.partition.scheme above (using lapply)
get.OKP.model <- function(name,char_lines){
  # name is a loci name
  # char_lines is a .txt file opened using readLines
  name_ind <- grep(name,char_lines)
  if (identical(name_ind,integer(0))){
    # If there's no matches in the file for this name, return NA
    m = NA
  } else {
    # If there is a match for this name in the file, find the match
    line <- char_lines[name_ind]
    # Break the line into columns using the dividers ("|")
    split_line <- strsplit(line,":")[[1]]
    # Take the second element - this will be the model
    m <- split_line[3]
    # Trim the white space
    m <- gsub(" ","",m)
  }
  # Return the model
  return(m)
}




# Function to copy alignment to new folder in output folder and convert to nexus if necessary
copy.alignment.as.nexus <- function(alignment_path, alignment_folder, loci_name, loci_row){
  file_type <- tail(strsplit(alignment_path,"\\.")[[1]],1)
  new_path <- paste0(alignment_folder,loci_name,".nex")
  if (file_type == "nex" || file_type == "nexus" || file_type == "nxs"){
    if (file.exists(new_path) == FALSE){
      file.copy(from = alignment_path, to = new_path)
    }
    alignment_path <- new_path
  } else if (file_type == "fasta" || file_type == "faa" || file_type == "fna" || 
             file_type == "fa" || file_type == "ffn" || file_type == "frn" || file_type == "fas") {
    new_fasta_path <- paste0(alignment_folder,loci_name,".fasta")
    # Move and save fasta file
    if (file.exists(new_fasta_path) == FALSE){
      file.copy(from = alignment_path, to = new_fasta_path)
    }
    # Convert fasta to nexus and save nexus in new folder
    if (file.exists(new_path) == FALSE){
      # Set details for fasta data
      if (loci_row$alphabet == "dna") {
        seq_type = "DNA" 
      } else if (loci_row$alphabet == "protein") {
        seq_type = "AA"
      }
      if (loci_row$dataset == "Vanderpool2020"){
        # read in the fasta data
        f_data <- read.fasta(alignment_path, seqtype = seq_type)
        f_new <- list()
        # Remove any " " in the alignment
        for (s in names(f_data)){
          f_new[[s]] <- c(f_data[[s]][grep(" ", f_data[[s]], invert = TRUE)])
        }
        # write the output as a nexus file to the output folder for this alignment
        write.nexus.data(f_new, file = new_path, format = loci_row$alphabet, interleaved = FALSE, datablock = FALSE)
      } else {
        f_data <- read.fasta(alignment_path, seqtype = seq_type)
        write.nexus.data(f_data, file = new_path, format = loci_row$alphabet, interleaved = FALSE, datablock = FALSE)
      }
      # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
      nexus <- readLines(new_path) # open the new nexus file
      ind <- grep("BEGIN CHARACTERS",nexus)+2 # find which line
      if (loci_row$alphabet == "dna"){
        nexus[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=-;" # replace the line
      } else if (loci_row$alphabet == "protein"){
        nexus[ind] <- "  FORMAT MISSING=? GAP=- DATATYPE=PROTEIN;" # replace the line
      }
      writeLines(nexus,new_path) # output the edited nexus file
      alignment_path <- new_path
    }
  }
  # Depending on which dataset you are using, return either the nexus path or the fasta path
  if (loci_row$dataset == "Vanderpool2020"){
    # Use nexus path
    return_path = new_path
  } else if (loci_row$dataset == "1KP"){
    # Use fasta path
    return_path = new_fasta_path
  } else if (loci_row$dataset == "Strassert2021"){
    # use fasta path
    return_path = new_fasta_path
  }
  return(return_path)
}





# Function to copy alignment to new folder in output folder and convert to nexus if necessary
copy.alignment.as.nexus.tpts <- function(alignment_path, alignment_folder, loci_name, loci_alphabet){
  file_type <- tail(strsplit(alignment_path,"\\.")[[1]],1)
  new_path <- paste0(alignment_folder,loci_name,".nex")
  if (file_type == "nex" || file_type == "nexus" || file_type == "nxs"){
    if (file.exists(new_path) == FALSE){
      file.copy(from = alignment_path, to = new_path)
    }
    alignment_path <- new_path
  } else if (file_type == "fasta" || file_type == "faa" || file_type == "fna" || 
             file_type == "fa" || file_type == "ffn" || file_type == "frn" || file_type == "fas") {
    new_fasta_path <- paste0(alignment_folder,loci_name,".fasta")
    # Move and save fasta file
    if (file.exists(new_fasta_path) == FALSE){
      file.copy(from = alignment_path, to = new_fasta_path)
    }
    # Convert fasta to nexus and save nexus in new folder
    if (file.exists(new_path) == FALSE){
      # Set details for fasta data
      if (loci_alphabet == "dna") {
        seq_type = "DNA" 
      } else if (loci_alphabet == "protein") {
        seq_type = "AA"
      }
      # read in the fasta data
      f_data <- read.fasta(alignment_path, seqtype = seq_type)
      # write the output as a nexus file to the output folder for this alignment
      write.nexus.data(f_data, file = new_path, format = loci_alphabet, interleaved = TRUE, datablock = FALSE)
      # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
      nexus <- readLines(new_path) # open the new nexus file
      ind <- grep("BEGIN CHARACTERS",nexus)+2 # find which line
      if (loci_alphabet == "dna"){
        nexus[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE;" # replace the line
      } else if (loci_alphabet == "protein"){
        nexus[ind] <- "  FORMAT MISSING=? GAP=- DATATYPE=PROTEIN INTERLEAVE;" # replace the line
      }
      writeLines(nexus,new_path) # output the edited nexus file
      alignment_path <- new_path
    }
  }
  return(new_path)
}





# Function to remove empty sequences from a nexus file (either AA or DNA)
remove.empty.taxa <- function(alignment_path, seq_type){
  file_type <- tail(strsplit(alignment_path,"\\.")[[1]],1)
  # If file is a nexus, open nexus and remove empty sequences
  if (file_type == "nex" || file_type == "nexus" || file_type == "nxs"){
    n <- read.nexus.data(alignment_path)
    # Initialise a new sequence
    copy_names <- c()
    seq_names <- names(n)
    # Iterate through the names and add non-missing sequences to the new sequence
    for (seq_name in seq_names){
      seq <- n[[seq_name]] # get the original empirical sequence
      chars <- toupper(unique(seq))
      if (setequal(chars,c("?")) || setequal(chars,c("-")) || setequal(chars,c("-","?")) || setequal(chars,c("X")) ||
          setequal(chars,c("-", "X")) || setequal(chars,c("-", "N")) || setequal(chars,c("N"))){
        # If the only characters are the empty character or the question mark character, ignore this sequence
        next
      } else {
        # If the sequence contains genetic information, include it in the new sequence
        copy_names <- c(copy_names,seq_name)
      }
    }
    n_new <- n[copy_names]
    write.nexus.data(n_new,file = alignment_path, format = seq_type, interleaved = TRUE)
    # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
    nexus_edit <- readLines(alignment_path) # open the new nexus file
    ind <- grep("BEGIN DATA",nexus_edit)+2 # find which line
    if (seq_type == "dna"){
      nexus_edit[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE;" # replace the line
    } else if (seq_type == "protein"){
      nexus_edit[ind] <- "  FORMAT MISSING=? GAP=- DATATYPE=PROTEIN INTERLEAVE;" # replace the line
    }
    writeLines(nexus_edit,alignment_path) # output the edited nexus file
    
  } else if (file_type == "fasta" || file_type == "faa" || file_type == "fna" || 
             file_type == "fa" || file_type == "ffn" || file_type == "frn" || file_type == "fas") {
    # Change seq_type to be the right parameters for opening a fasta file
    if (seq_type == "dna"){
      fasta_type = "DNA"
    } else if (seq_type == "protein"){
      fasta_type = "AA"
    }
    # If file is a fasta, open fasta and remove empty sequences
    f <- read.fasta(alignment_path, fasta_type)
    # Initialise a new sequence
    copy_names <- c()
    seq_names <- names(f)
    # Iterate through the names and add non-missing sequences to the new sequence
    for (seq_name in seq_names){
      seq <- f[[seq_name]] # get the original empirical sequence
      chars <- toupper(unique(seq))
      if (setequal(chars,c("?")) || setequal(chars,c("-")) || setequal(chars,c("-","?")) || setequal(chars,c("X")) ||
          setequal(chars,c("-", "X")) || setequal(chars,c("-", "N")) || setequal(chars,c("N"))){
        # If the only characters are the empty character or the question mark character, ignore this sequence
        next
      } else {
        # If the sequence contains genetic information, include it in the new sequence
        copy_names <- c(copy_names,seq_name)
      }
    }
    # Copy across the non-empty sequences
    f_new <- f[copy_names]
    # Write the new fasta file
    write.fasta(f_new, copy_names, file.out = alignment_path, open = "w")
  }
}




# Function to remove empty sequences from a nexus file (either AA or DNA)
prune.taxa.by.length <- function(alignment_path, proportion_allowed_missing, seq_type, write_output_text = FALSE, output_folder){
  # Set which characters you don't want in your sequences based on sequence type
  if (seq_type == "dna"){
    missing_chars <- "Z|O|N|X|-|\\.|\\~|\\*|\\?"
  } else if (seq_type == "protein"){
    missing_chars <- "X|-|\\.|\\~|\\*|\\?"
  }
  # Get file type
  file_type <- tail(strsplit(alignment_path,"\\.")[[1]],1)
  # If file is a nexus
  if (file_type == "nex" || file_type == "nexus" || file_type == "nxs"){
    # Read in nexus file
    n <- read.nexus.data(alignment_path)
    # Initialise a new sequence
    copy_names <- c()
    seq_names <- names(n)
    output_text <- c("Taxa, proportion_missing")
    # Iterate through the names. If a sequence has more missing/ambiguous characters than the allowed proportion, do not include it
    for (seq_name in seq_names){
      seq <- toupper(n[[seq_name]]) # get the original empirical sequence
      proportion_missing <- sum(lengths(regmatches(seq, gregexpr(missing_chars, seq))))/length(seq)
      output_text <- c(output_text,paste0(seq_name," , ",proportion_missing))
      if (proportion_missing<=proportion_allowed_missing){
        copy_names <- c(copy_names,seq_name)
      }
    }
    n_new <- n[copy_names]
    write.nexus.data(n_new,file = alignment_path, format = seq_type, interleaved = FALSE)
    # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
    nexus_edit <- readLines(alignment_path) # open the new nexus file
    ind <- grep("BEGIN DATA",nexus_edit)+2 # find which line
    if (seq_type == "dna"){
      nexus_edit[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=-;" # replace the line
    } else if (seq_type == "protein"){
      nexus_edit[ind] <- "  FORMAT MISSING=? GAP=- DATATYPE=PROTEIN;" # replace the line
    }
    writeLines(nexus_edit,alignment_path) # output the edited nexus file
  } else if (file_type == "fasta" || file_type == "faa" || file_type == "fna" || 
             file_type == "fa" || file_type == "ffn" || file_type == "frn" || file_type == "fas") {
    # If file is a fasta:
    # Change the type into one readable by read.fasta function
    if (seq_type == "dna"){
      fasta_type = "DNA"
    } else if (seq_type == "protein"){
      fasta_type = "AA"
    }
    # Open fasta
    f <- read.fasta(alignment_path, fasta_type)
    # Initialise a new sequence
    copy_names <- c()
    seq_names <- names(f)
    output_text <- c("Taxa, proportion_missing")
    # Iterate through the names. If a sequence has more missing/ambiguous characters than the allowed proportion, do not include it
    for (seq_name in seq_names){
      seq <- toupper(f[[seq_name]]) # get the original empirical sequence
      proportion_missing <- sum(lengths(regmatches(seq, gregexpr(missing_chars, seq))))/length(seq)
      output_text <- c(output_text,paste0(seq_name," , ",proportion_missing))
      if (proportion_missing<=proportion_allowed_missing){
        copy_names <- c(copy_names,seq_name)
      }
    }
    # Copy across the non-empty sequences
    f_new <- f[copy_names]
    # Write the new fasta file
    write.fasta(f_new, copy_names, file.out = alignment_path, open = "w")
  }
  
  # If required, output all the proportions
  if (write_output_text == TRUE){
    write(output_text, file = paste0(output_folder,"proportion_of_ambiguous_sequence.txt"))
  }
}




# Function to remove invalid characters from a nexus file (either AA or DNA)
check.invalid.nexus.characters <- function(alignment_path, seq_type){
  if (seq_type == "dna"){
    n <- read.nexus.data(alignment_path)
    seq_names <- names(n)
    # Find which sequences contain ambiguous characters
    to_edit_inds <- grep("W|w|S|s|M|m|K|k|R|r|Y|y|B|b|D|d|H|h|V|v",n)
    to_edit_seqs <- seq_names[to_edit_inds]
    # If there are more than one of these characters, convert the file into a fasta file
    # While IQ-Tree can handle ambiguous characters, the NEXUS format has different language (e.g. "[CT]" for "Y")
    # By changing to a fasta file, you can sidestep this
    if (length(to_edit_seqs) > 0){
      fasta.name <- gsub(".nex",".fasta",alignment_path)
      if (!file.exists(fasta.name)){
        write.dna(n, file = fasta.name, format = "fasta")
      }
      output_indicator = "use_FASTA"
    } else if (length(to_edit_seqs) == 0){
      output_indicator = "no_ambiguous_characters"
    }
  } else if (seq_type == "protein"){
    n <- read.nexus.data(alignment_path)
    seq_names <- names(n)
    # Find which sequences contain ambiguous characters
    to_edit_inds <- grep("X|x|B|b|Z|z|K|j|Φ|Ω|Ψ|π|ζ|+",n)
    to_edit_inds <- grep("X|x|B|b|Z|z|K|j|\\Φ|\\Ω|\\Ψ|\\π|\\ζ|\\+",n)
    to_edit_seqs <- seq_names[to_edit_inds]
    if (length(to_edit_seqs) > 0){
      fasta.name <- gsub(".nex",".fasta",alignment_path)
      if (!file.exists(fasta.name)){
        write.dna(n, file = fasta.name, format = "fasta")
      }
      output_indicator = "use_FASTA"
    } else if (length(to_edit_seqs) == 0){
      output_indicator = "no_ambiguous_characters"
    }
  }
  return(output_indicator)
}





# Quick function to randomly add/subtract to make sure that the number of sites in each gamma category is correct
fix.gammaCategory.siteNums <- function(df,num){
  print("fix gamma categories")
  if (sum(df$cat_sites)==num){
    # If the sum of the num of sites = num of sites, return the df (all is good)
    return(df)
  }
  else {
    while (sum(df$cat_sites)!=num){
      # randomly sample a row
      row <- df[sample(nrow(df),1),]
      row_id <- which(df$category==row$category)
      if (sum(df$cat_sites)<num){
        row$cat_sites <- row$cat_sites + 1
      }
      if (sum(df$cat_sites)>num){
        row$cat_sites <- row$cat_sites - 1
      }
      # replace the row
      df[row_id,]<-row
    }
    # outside the while loop: therefor return the gamma categories matrix
    return(df)
  }
}





# Quick function to estimate a phylogenetic network in SplitsTree and open it using splits format in phangorn
estimateNetwork <- function(alignment_path, splitstree_path, network_algorithm = "neighbournet",seq_type){
  call.SplitsTree(splitstree_path,alignment_path,network_algorithm,seq_type)
  # Retrieve the file name for the splits output file
  splits.filepath <- splits.filename(alignment_path)
  
  # Get the number of splits by opening the splits block and reading the number after "nsplits-"
  splits_text <- readLines(splits.filepath)
  splits_row <- splits_text[grep("nsplits",splits_text)]
  nsplits <- str_extract_all(splits_row,"[0-9]+")[[1]][2]
  
  # Check whether there are valid splits in the splits block
  # If there are, open the splits
  if (as.numeric(nsplits) > 0){
    # Otherwise, if the network does contain more than 0 splits:
    # Extract the splits 
    splits <- read.nexus.splits(splits.filepath) # Open the splits from SplitsTree
  }
  return(splits)
}




# Function to open an iqtree file and extract the tree depth
extract.total.tree.length <- function(alignment_path){
  # Open .iqtree file associated with that alignment
  iqtree_filepath <- paste0(alignment_path,".iqtree")
  if (file.exists(iqtree_filepath) == FALSE){
    return("MISSING_FILE")
  } else if (file.exists(iqtree_filepath) == TRUE){
    iq_file <- readLines(iqtree_filepath)
    ind <- grep("Total tree length",iq_file)
    str <- iq_file[ind]
    vec <- strsplit(str," ")[[1]]
    tree_length <- vec[str_detect(vec,"([0-9])")]
    return(tree_length)
  }
}




# Function to open an iqtree file and extract the sum of internal branches
extract.sum.internal.branch.length <- function(alignment_path){
  # Open .iqtree file associated with that alignment
  iqtree_filepath <- paste0(alignment_path,".iqtree")
  if (file.exists(iqtree_filepath) == FALSE){
    return("MISSING_FILE")
  } else if (file.exists(iqtree_filepath) == TRUE){
    iq_file <- readLines(iqtree_filepath)
    ind <- grep("Sum of internal branch lengths",iq_file)
    str <- iq_file[ind]
    vec <- strsplit(str," ")[[1]]
    sibl <- vec[str_detect(vec,"([0-9])")]
    sibl <- gsub("\\(","",sibl)
    sibl <- gsub("%","",sibl)
    return(sibl)
  }
}




# Function to open an iqtree file and extract the number of parsimony informative sites
extract.num.informative.sites <- function(alignment_path){
  # Open .iqtree file associated with that alignment
  iqtree_filepath <- paste0(alignment_path,".iqtree")
  if (file.exists(iqtree_filepath) == FALSE){
    return("MISSING_FILE")
  } else if (file.exists(iqtree_filepath) == TRUE){
    iq_file <- readLines(iqtree_filepath)
    ind <- grep("Number of parsimony informative sites",iq_file)
    str <- iq_file[ind]
    vec <- strsplit(str," ")[[1]]
    num_sites <- vec[str_detect(vec,"([0-9])")]
    return(num_sites)
  }
}





# Function to open an iqtree file and extract the number of parsimony informative sites
extract.treefile.tree <- function(alignment_path){
  # Open .iqtree file associated with that alignment
  iqtree_filepath <- paste0(alignment_path,".treefile")
  if (file.exists(iqtree_filepath) == FALSE){
    return("MISSING_FILE")
  } else if (file.exists(iqtree_filepath) == TRUE){
    tree <- readLines(iqtree_filepath)
    return(tree)
  }
}





# Function to extract the mean and variance of the GC content of an alignment 
# (by checking sequence by sequence then taking the mean, rather than whole alignment - allows you to get sd and var as well as mean)
calculate.GC.content <- function(alignment_file){
  al <- read.nexus.data(alignment_file)
  al <- as.DNAbin(al)
  vals <- unlist(lapply(names(al), function(x){GC.content(al[x])}))
  mean <- round(mean(vals),digits = 3)
  var <- signif(var(vals), digits = 3)
  sd <- signif(sd(vals), digits = 3)
  return(c(mean,var,sd))
}




# Function to extract the mean and variance of the GC content of an alignment 
# (by checking sequence by sequence then taking the mean, rather than whole alignment - allows you to get sd and var as well as mean)
save.tree <- function(alignment_folder,trees_folder){
  # Get all files in folder
  alignment_files <- list.files(alignment_folder)
  # Find treefiles
  treefile_files <- grep(".treefile",alignment_files, value = TRUE)
  # Extract the real tree - the only ".treefile" that isn't a ".treefile."
  treefile <- grep(".treefile.", treefile_files, invert = TRUE, value = TRUE)
  # Create variables for full name of file in current and new location
  old_treefile <- paste0(alignment_folder,"/", treefile)
  new_treefile <- paste0(trees_folder, treefile)
  # Save the file into the trees_folder
  file.copy(from = old_treefile, to = new_treefile, overwrite = TRUE, recursive = FALSE)
}




# Given a tree and a set of parameters extracted from IQTree, this function generates a DNA multiple sequence alignment
generate.DNA.alignment <- function(alignment_params, empirical_alignment_tree){
  print("creating alignment")
  # Sample code for generating a parametric DNA sequence if you have a tree
  # s1 = simSeq(t1, l = 500, type="DNA", bf=c(.25,.25,.25,.25), Q=c(1,1,1,1,1,1), rate=1)
  # s2 = simSeq(t2, l = 500, type="DNA", bf=c(.25,.25,.25,.25), Q=c(1,1,1,1,1,1), rate=1)
  # aln = c(s1, s2) # concatenate the alignments
  
  # Extract the parameters you need to enter into simSeq
  n_bp = as.numeric(alignment_params$parameters[4,2])  # sequence length should be the same as in the original sequence
  # Extract the vector form of the rate matrix 
  m <- alignment_params$Q_rate_matrix[,2:5] # extract the square block with the rates and not the header column
  Q_vec <- c(m[2,1],m[3,1],m[3,2],m[4,1],m[4,2],m[4,3]) # extract the rates 
  # Extract the base frequencies in the following order: A, C, G, T
  base_freqs <- c(as.numeric(alignment_params$parameters[[16,2]]), as.numeric(alignment_params$parameters[[17,2]]),
                  as.numeric(alignment_params$parameters[[18,2]]), as.numeric(alignment_params$parameters[[19,2]]))
  seq_type <- "DNA" # generate DNA sequence
  
  # need to generate the number of sites for each gamma category:
  g_cat <- alignment_params$gamma_categories
  # If uniform gamma categories, can just simulate a sequence straight off
  if (class(g_cat) == "character"){
    new_aln <- simSeq(x = empirical_alignment_tree, l = n_bp, type = seq_type, bf = base_freqs, Q = Q_vec)
  } else if (class(g_cat) == "data.frame") {
    # If there are gamma categories, create multiple alignments and concatenate them: one for each gamma category
    cat_sites <- round(n_bp*g_cat$proportion,digits = 0)
    g_cat$cat_sites <- cat_sites
    g_cat <- fix.gammaCategory.siteNums(g_cat,n_bp) # make sure the sum of the gamma sites category is correct
    # Get the number of gamma rate categories
    num_g_cat <- length(g_cat$category)
    # Initialise an empty sequence
    aln_exists <- FALSE
    
    # Iterate through the gamma categories and create an alignment with each of the relative rates
    for (i in 1:num_g_cat){
      row <- g_cat[i,]
      # Generate the DNA sequence using the informtion from the parameter list and the rate and number of sites for THIS GAMMA CATEGORY
      dna_sim <- simSeq(x = empirical_alignment_tree, l = row$cat_sites, type = seq_type, bf = base_freqs, Q = Q_vec, rate = row$relative_rate)
      #concatenate the sequence with the new alignment
      if (aln_exists == FALSE){
        # If the alignment doesn't exist, create it
        new_aln <- dna_sim
        aln_exists <- TRUE
      } else if (aln_exists == TRUE){
        # If the alignment does exist, concatenate it
        new_aln <- c(new_aln,dna_sim)
      }
    }
  }
  return(new_aln)
}

# Given a tree and a set of parameters extracted from IQTree, this function generates an amino acid multiple sequence alignment
generate.AA.alignment <- function(alignment_params, empirical_alignment_tree){
  print("creating alignment")
  # Extract the parameters you need to enter into simSeq
  n_bp = as.numeric(alignment_params$parameters[4,2])  # sequence length should be the same as in the original sequence
  # extract explanation of state frequencies:
  freqs_info <- gsub(" ","",alignment_params$parameters[13,2])
  # extract the amino acid frequencies
  if (freqs_info == "(empiricalcountsfromalignment)"){
    aa_order <- alignment_params$frequency$amino_acid
    aa_freq <- as.numeric(alignment_params$frequency$frequency)
  } else if (freqs_info == "(model)") {
    aa_freq <- NULL
  }
  seq_type <- "AA" # generate DNA sequence
  # Extract the AA model of sequence evolution
  model <- alignment_params$parameters[9,2]
  # check if the model contains +G, +I, or +F
  if (grepl("\\+",model) == TRUE){
    # if the model contains these additions, remove them and get the base model
    # split the string at the "+" characters
    model_split <- strsplit(model,"\\+")
    # the base model will be the first element (as the structure is Model+FreqType+RateType e.g. JTT+F+I+G)
    model <- model_split[[1]][1]
  }
  # some model names are different in IQTree and simSeq function in R - change those model names to be accepted in simSeq
  if (model == "JTTDCMut"){
    model = "JTT_DCMut"
  } else if (model == "DCMut") {
    model = "Dayhoff_DCMut"
  }
  # if there are multiple gamma categories, we will generate multiple alignments (one for each category)
  # and concatenate them together to create one alignment
  g_cat <- alignment_params$gamma_categories
  # If there are multiple gamma categories, format them nicely so they can be used to generate sequences
  if (class(g_cat) == "data.frame") {
    # need to generate the number of sites for each gamma category:
    cat_sites <- round(n_bp*g_cat$proportion,digits = 0)
    g_cat$cat_sites <- cat_sites
    g_cat <- fix.gammaCategory.siteNums(g_cat,n_bp) # make sure the sum of the gamma sites category is correct
    # Get the number of gamma rate categories
    num_g_cat <- length(g_cat$category)
  }
  # Only certain models can be fed into the simSeq function - if the model is one of these, feed it in
  # Only models needed for the Misof 2014, 1Kp and Vanderpool 2020 datasets have been added as of 21/10/2020
  if (class(g_cat) == "character"){
    # If uniform gamma categories, can just simulate a sequence straight off
    new_aln <- simSeq(x = empirical_alignment_tree, l = n_bp, type = seq_type, bf = aa_freq, model = model)
  } else if (class(g_cat) == "data.frame"){
    # If there are gamma categories, create multiple alignments and concatenate them: one for each gamma category
    # Initialise an empty sequence
    aln_exists <- FALSE
    
    # Iterate through the gamma categories and create an alignment with each of the relative rates
    for (i in 1:num_g_cat){
      row <- g_cat[i,]
      # Generate the DNA sequence using the information from the parameter list and the rate and number of sites for THIS GAMMA CATEGORY
      aa_sim <- simSeq(x = empirical_alignment_tree, l = row$cat_sites, type = seq_type, bf = aa_freq, model = model, rate = row$relative_rate)
      #concatenate the sequence with the new alignment
      if (aln_exists == FALSE){
        # If the alignment doesn't exist, create it
        new_aln <- aa_sim
        aln_exists <- TRUE
      } else if (aln_exists == TRUE){
        # If the alignment does exist, concatenate it
        new_aln <- c(new_aln,aa_sim)
      }
    }
  }
  return(new_aln)
}




# Function to add information about the alignment from the parameter file to the p-value file
add.alignment.information <- function(loci_folder){
  # Find and open the p-value, parameter, and tree files
  loci_files <- list.files(loci_folder, full.name = TRUE)
  p_value_file <- grep("_pValues.csv", loci_files, value = TRUE)
  p_value_df <- read.csv(p_value_file)
  params_file <- grep("_parameterValues.csv", loci_files, value = TRUE)
  params <- read.csv(params_file)
  all_treefiles <- grep("treefile",loci_files, value = TRUE)
  newick_treefile <- grep(".treefile.", all_treefiles, invert = TRUE, value = TRUE)[1] # Extract the only ".treefile" file that isn't ".treefile."
  alignment_file <- gsub(".treefile","",newick_treefile)
  
  # Extract the relevant variables
  p_value_df$sequence_type <- params$value[which(params$parameter == "sequence_type")]
  p_value_df$num_constant_sites <- params$value[which(params$parameter == "Number of constant sites")]
  p_value_df$num_invariant_sites <- params$value[which(params$parameter == "Number of invariant (constant or ambiguous constant) sites")]
  p_value_df$num_parsimony_informative_sites <- params$value[which(params$parameter == "Number of parsimony informative sites")]
  p_value_df$num_site_patterns <- params$value[which(params$parameter == "Number of distinct site patterns")]
  p_value_df$substitution_model <- params$value[which(params$parameter == "substitution_model")]
  p_value_df$AC_rate <- params$value[which(params$parameter == "A-C_rate")]
  p_value_df$AG_rate <- params$value[which(params$parameter == "A-G_rate")]
  p_value_df$AT_rate <- params$value[which(params$parameter == "A-T_rate")]
  p_value_df$CG_rate <- params$value[which(params$parameter == "C-G_rate")]
  p_value_df$CT_rate <- params$value[which(params$parameter == "C-T_rate")]
  p_value_df$GT_rate <- params$value[which(params$parameter == "G-T_rate")]
  p_value_df$A_freq <- params$value[which(params$parameter == "A_freq")]
  p_value_df$C_freq <- params$value[which(params$parameter == "C_freq")]
  p_value_df$G_freq <- params$value[which(params$parameter == "G_freq")]
  p_value_df$T_freq <- params$value[which(params$parameter == "T_freq")]
  p_value_df$GC_content_mean <- calculate.GC.content(alignment_file)[1]
  p_value_df$GC_content_variance <- calculate.GC.content(alignment_file)[2]
  p_value_df$GC_content_sd <- calculate.GC.content(alignment_file)[3]
  p_value_df$model_of_rate_heterogeneity <- params$value[which(params$parameter == "model_of_rate_heterogeneity")]
  p_value_df$model_of_rate_heterogeneity_line2_name <- params$value[which(params$parameter == "model_of_rate_heterogeneity_line2_name")]
  p_value_df$model_of_rate_heterogeneity_line2_value <- params$value[which(params$parameter == "model_of_rate_heterogeneity_line2_value")]
  p_value_df$total_tree_length <- extract.total.tree.length(alignment_file)
  p_value_df$sum_of_internal_branch_lengths <- extract.sum.internal.branch.length(alignment_file)[1]
  p_value_df$proportion_internal_branches <- extract.sum.internal.branch.length(alignment_file)[2]
  
  # Open the tree file
  newick_tree <- readLines(newick_treefile)
  p_value_df$tree <- newick_tree
  
  # Save new, extended p-value file
  out_file <- gsub("pValues","results",p_value_file)
  write.csv(p_value_df, out_file, row.names = FALSE)
}




# Function to take a list of loci and copy the trees for those loci into either a new folder, a new text file, or both
copy.loci.trees <- function(loci_names,loci_trees,output_folder,output_name,copy.all.individually = FALSE, copy.and.collate = TRUE){
  # Make a little dataframe of the tree information
  tree_df <- data.frame(names = loci_names, trees = loci_trees)
  # Create output file name and folder name
  op_dir <- paste0(output_folder, output_name, "/")
  op_file <- paste0(output_folder, output_name, ".txt")
  # If required, save each tree into a separate file within the same folder
  if (copy.all.individually == TRUE){
    # Create the folder to save trees in (if it doesn't already exist)
    if (!dir.exists(op_dir)){
      dir.create(op_dir)
    }
    # Iterate through the tree_df to save each tree
    lapply(1:nrow(tree_df), save.one.tree, tree_df, op_dir)
  } 
  # If required, save each tree into the same file
  if (copy.and.collate == TRUE){
    # Collate trees together
    all_tree_text <- paste0(tree_df$trees)
    # Output all trees
    write(x = all_tree_text, file = op_file)
  }
}

# Function to take a single tree from a tree dataframe and copy it into an output folder
save.one.tree <- function(row_number,tree_df,op_dir){
  row <- tree_df[row_number,]
  op_name <- paste0(op_dir, row$names, ".treefile" )
  write(x = row$trees, file = op_name)
}




# Function to estimate species tree using ASTRAL
estimate.ASTRAL.species.tree <- function(gene_tree_file, species_tree_file, log_file, ASTRAL_path){
  astral_command <- paste0("java -jar ", ASTRAL_path, " -i ", gene_tree_file, " -o ", species_tree_file, " 2> ", log_file)
  system(astral_command)
}




ASTRAL.wrapper <- function(text_file, ASTRAL_path){
  species_tree_file <- gsub(".txt", "_species.tre", text_file)
  log_file <- gsub(".txt", "_species.log", text_file)
  estimate.ASTRAL.species.tree(text_file, species_tree_file, log_file, ASTRAL_path)
}




# Function to estimate species tree of a folder full of alignments using IQ-Tree
estimate.IQTREE.species.tree <- function(gene_tree_folder, IQTREE_path){
  iqtree_command <- paste0(IQTREE_path, " -p ", gene_tree_folder, " -m MFP+MERGE -bb 1000 -nt AUTO")
  system(iqtree_command)
}



# Function to estimate species tree of a folder full of alignments using IQ-Tree
estimate.partitioned.IQTREE.species.tree <- function(partition_file, IQTREE_path){
  iqtree_command <- paste0(IQTREE_path, " -p ", partition_file, " -m MFP+MERGE -bb 1000 -nt AUTO")
  system(iqtree_command)
}



# Function to take a list of loci and copy the alignments for those loci into either a new folder
copy.loci.alignment <- function(loci_name, dataset_loci_folder, new_alignment_location){
  # Extract the list of loci from the dataset_loci_folder
  all_loci_alignments <- list.files(dataset_loci_folder)
  # Extract the loci you are looking for 
  loci_path <- paste0(dataset_loci_folder, grep(loci_name, all_loci_alignments, value = TRUE))
  # Check if the new directory exists
  if (!dir.exists(new_alignment_location)){
    dir.create(new_alignment_location)
  }
  # Extract the suffix from the loci_path - ensure you are keeping file format constant!
  suffix <- tail(strsplit(loci_path,"\\.")[[1]], 1)
  # Create the new path using the new location, loci name and the file extension
  new_loci_path <- paste0(new_alignment_location, loci_name, ".", suffix)
  # Check if the file exists in the new place
  if (!file.exists(new_loci_path)){
    file.copy(loci_path, new_loci_path, overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
  }
}



make.partition.file <- function(directory, add.charpartition = FALSE){
  # get all loci in directory
  all_loci <- list.files(directory)
  # remove any partition files
  all_loci <- grep("partition", all_loci, invert = TRUE, value = TRUE)
  # make start and end of partition files
  header <- c("#nexus","begin sets;","[loci]")
  footer <- c("end;")
  # make charset for each loci - split by codon position
  all_charsets <- unlist(lapply(all_loci, make.codon.position.charset, directory))
  # make the charpartition
  all_loci_labels <- unlist(lapply(all_loci, make.charpartition.labels))
  all_inds <- 1:length(all_loci_labels)
  loci_order <- paste0(all_inds,":",all_loci_labels,", ")
  loci_order[length(loci_order)] <- gsub(", ",";",loci_order[length(loci_order)])
  if (add.charpartition == TRUE){
    # generate charpartition
    charpartitions_vec <- c("\tcharpartition loci = ",loci_order)
    charpartitions <- paste(charpartitions_vec, collapse = "")
    # attach all the information together
    op_lines <- c(header, all_charsets, "", charpartitions, footer)
  } else
  {
    # attach all the information together
    op_lines <- c(header, all_charsets, footer)
  }
  # output partition file
  output_file <- paste0(directory, "partitions.nex")
  writeLines(op_lines, output_file)
}



make.codon.position.charset <- function(filename, directory){
  loci_name <- remove.suffix(filename)
  firstpos <- paste0("\tcharset ", loci_name, "_1stpos = ", directory, filename, ": 1-.\\3;")
  secondpos <- paste0("\tcharset ", loci_name, "_2ndpos = ", directory, filename, ": 2-.\\3;")
  thirdpos <- paste0("\tcharset ", loci_name, "_3rdpos = ", directory, filename, ": 3-.\\3;")
  file_charset <- c(firstpos, secondpos, thirdpos)
  return(file_charset)
}



make.charpartition.labels <- function(filename){
  loci_name <- remove.suffix(filename)
  labels <- c(paste0(loci_name, "_1stpos"),paste0(loci_name, "_2ndpos"),paste0(loci_name, "_3rdpos"))
  return(labels)
}



remove.suffix <- function(full_filename){
  split_filename <- strsplit(full_filename,"\\.")[[1]]
  just_filename <- paste0(split_filename[1:(length(split_filename) -1)])
  return(just_filename)  
}


##### Functions for calculating statistical test for tree proportion #####
# The functions from here onwards are to calculate the statistical test for tree proportion
# These have been modified and adapted from the functions used to perform the analyses to process just one alignment at a time
# These functions call functions found higher up in this script (e.g. call.IQTREE.empirical) and in other scripts 
# (e.g. tree.proportion from https://github.com/caitlinch/treelikeness/code/func_test_statistic.R)

# Function to perform parametric bootstrap and calculate statistical test for tree proportion
tree.proportion.statistical.test <- function(loci_path, loci_name, loci_alphabet, loci_model, loci_dataset, loci_output_folder, iqtree_path, splitstree_path, number_of_replicates, 
                                             allowable_proportion_missing_sites, iqtree.num_threads, num_of_cores){
  print("tree.proportion.test.statistic")
  # Specify path to alignment
  empirical_alignment_path <- loci_path
  print(paste0("Loci name: ",loci_name))
  # Check if the loci name already has a slash at the end
  loci_folder <- paste0(normalizePath(loci_output_folder),"/")
  # Create a folder to store the files for this loci in
  alignment_folder <- paste0(loci_folder,loci_name,"_tree_proportion/")
  # Check if this folder exists - and if it doesn't, create it!
  if (dir.exists(alignment_folder) == FALSE){
    dir.create(alignment_folder)
  }
  # Set this as the working directory
  setwd(alignment_folder)
  # Create output file names, the name of the loci and the file path of the loci location
  collated_ts_file <- paste0(alignment_folder,loci_name,"_testStatistics_collatedBSReplicates.csv")
  p_value_file  <- paste0(alignment_folder,loci_name,"_pValues.csv")
  # Create file names for saving the IQTree parameters (e.g. rate matrix, model of sequence evolution, gamma categories)
  parameters_file <- paste0(alignment_folder,loci_name,"_parameterValues.csv")
  gamma_categories_file <- paste0(alignment_folder,loci_name,"_gammaCategories.csv")
  aa_frequency_file <- paste0(alignment_folder,loci_name,"_amino_acid_frequencies.csv")
  rate_matrix_file <- paste0(alignment_folder,loci_name,"_QRateMatrix.csv")
  ts_file <- paste0(dirname(empirical_alignment_path),"/",loci_name,"_testStatistics.csv")
  
  # Only run this section if the p-value csv has not been created yet (skip reruns)
  if (file.exists(p_value_file) == FALSE){
    # If alignment is a nexus, copy it into the alignment folder
    # If it's not, rewrite it into a nexus and copy to the alignment folder
    # From now on, use the copy for analysis (leaving the original untouched)
    empirical_alignment_path <- copy.alignment.as.nexus.tpts(loci_path, alignment_folder, loci_name, loci_alphabet)
    
    # Remove empty taxa from the alignment
    remove.empty.taxa(empirical_alignment_path, loci_alphabet)
    # Remove any sequences from the alignment that have more than the allowable proportion of ambiguous/missing sites or gaps
    if (!is.na(allowable_proportion_missing_sites)){
      # If the allowable proportion of missing sites is NOT NA, run the pruning function
      # The pruning function will remove any sequences where the proportion of missing sites is GREATER than the allowable proportion of missing sites
      prune.taxa.by.length(empirical_alignment_path, allowable_proportion_missing_sites, loci_alphabet, write_output_text = TRUE, alignment_folder)
    }
    
    # Remove characters that IQ-Tree won't accept from the alignment
    # Leave only A,C,G,N,T,-
    if (loci_alphabet == "dna"){
      invalid_character_check <- check.invalid.nexus.characters(empirical_alignment_path, loci_alphabet)
    } else if (loci_alphabet == "protein"){
      invalid_character_check <- "no_ambiguous_characters"
    }
    
    # Check whether IQ-Tree has been run
    if (file.exists(paste0(empirical_alignment_path,".treefile")) == FALSE){
      print("need to run IQ-Tree")
      # Run IQ-tree on the alignment (if it hasn't already been run), and get the sCF results
      print("run IQTree and estimate sCFs")
      if (invalid_character_check == "contains_ambiguous_characters-use_FASTA"){
        # If the nexus file contains ambiguous characters that aren't permitted by the format, run IQ-Tree using the FASTA file
        empirical_alignment_path <- gsub(".nex",".fasta",empirical_alignment_path)
        call.IQTREE.empirical(empirical_alignment_path, iqtree_path, loci_model, iqtree.num_threads)
      } else {
        # If the nexus file doesn't contain any ambiguous characters, use the nexus file to run IQ-Tree
        call.IQTREE.empirical(empirical_alignment_path, iqtree_path, loci_model, iqtree.num_threads)
      }
    }
    
    # Calculate the test statistics if it hasn't already been done
    if (file.exists(ts_file) == FALSE){
      print("run test statistics")
      output.tree.proportion.csv(empirical_alignment_path, iqtree_path, splitstree_path, bootstrap_id = "alignment", iqtree.num_threads, 
                                 iqtree.model = loci_model, alphabet = loci_alphabet, dataset_name = loci_dataset)
    }
    
    #Check that the test statistic file ran ok 
    if (file.exists(ts_file) == FALSE){
      print("need to rerun test statistics")
      output.tree.proportion.csv(empirical_alignment_path, iqtree_path, splitstree_path, bootstrap_id = "alignment", iqtree.num_threads, 
                                 iqtree.model = loci_model, alphabet = loci_alphabet, dataset_name = loci_dataset)
    }
    
    #Extract the parameters from the .iqtree log file.
    print("get simulation params")
    params <- get.empirical.simulation.parameters(paste0(empirical_alignment_path,".iqtree"))
    write.csv(params$parameters,file = parameters_file, row.names = TRUE)
    write.csv(params$gamma_categories,file = gamma_categories_file, row.names = TRUE)
    if (loci_alphabet == "dna") {
      write.csv(params$Q_rate_matrix,file = rate_matrix_file, row.names = TRUE)
    } else if (loci_alphabet == "protein"){
      write.csv(params$frequency,file = aa_frequency_file, row.names = TRUE)
    }
    
    
    # Create the bootstrap ids (pad out to 4 digits) - should be "bootstrapReplicateXXXX" where XXXX is a number
    bootstrap_ids <- paste0("bootstrapReplicate",sprintf("%04d",1:number_of_replicates))
    
    # Run all the bootstrap ids that HAVEN'T already been run (e.g. in previous attempts) using lapply (feed info into do1.empirical.parametric.bootstrap)
    # If the alignment doesn't exist OR the test statistic csv doesnt exist, this indicates a complete run has not previously been done 
    # These bootstrap replicates will thus be calculates
    # This should save A BUNCH of time because it means if the test statistic file exists, you don't have to run any calculations
    print(paste0("Prepare ",number_of_replicates," bootstrap replicates"))
    bs_als <- paste0(alignment_folder,loci_name,"_",bootstrap_ids,"/",loci_name,"_",bootstrap_ids,".nex")
    ts_csvs <- paste0(alignment_folder,loci_name,"_",bootstrap_ids,"/",loci_name,"_",bootstrap_ids,"_testStatistics.csv")
    missing_als <- bs_als[!file.exists(bs_als)]
    missing_testStatistics <- bs_als[!file.exists(ts_csvs)]
    # Collate the missing files and identify the alignments to rerun
    all_to_run <- unique(c(missing_als,missing_testStatistics))
    ids_to_run <- bootstrap_ids[which((bs_als %in% all_to_run))]
    print(paste0("Run ",length(ids_to_run)," (of ", number_of_replicates, ") bootstrap replicates"))
    if(length(ids_to_run)>0){
      mclapply(ids_to_run, run.one.tpts.bootstrap, empirical_alignment_path = empirical_alignment_path, empirical_alignment_loci_name = loci_name, 
               empirical_alignment_alphabet = loci_alphabet, empirical_alignment_model = loci_model, empirical_alignment_dataset = loci_dataset,
               alignment_params = params, iqtree_path = iqtree_path, splitstree_path = splitstree_path, iqtree.num_threads = iqtree.num_threads, 
               mc.cores = num_of_cores)
    }
    
    # Before you can collate all the bootstrap files, you need to check every bootstrap ran and rerun the failed ones
    # Generate the names of each alignment, the test statistics csvs, the .iqtree files, the treefiles, the likelihood mapping files
    # Check which of these files are missing
    print("check for missing alignments")
    bs_als <- paste0(alignment_folder,"/",loci_name,"_",bootstrap_ids,"/",loci_name,"_",bootstrap_ids,".nex")
    ts_csvs <- paste0(alignment_folder,"/",loci_name,"_",bootstrap_ids,"/",loci_name,"_",bootstrap_ids,"_testStatistics.csv")
    missing_als <- bs_als[!file.exists(bs_als)]
    missing_iqtree <- bs_als[!file.exists(paste0(bs_als,".iqtree"))]
    missing_tree <- bs_als[!file.exists(paste0(bs_als,".treefile"))]
    missing_testStatistics <- bs_als[!file.exists(ts_csvs)]
    # Collate the missing files and identify the alignments to rerun
    all_missing <- unique(c(missing_als,missing_iqtree,missing_tree,missing_testStatistics))
    als_to_rerun <- bootstrap_ids[which((bs_als %in% all_missing))]
    print(paste0("Number of missing alignments to rerun = ",length(als_to_rerun)))
    # Rerun the missing als
    if (length(als_to_rerun)>0){
      mclapply(ids_to_run, do1.empirical.parametric.bootstrap, empirical_alignment_path = empirical_alignment_path, empirical_alignment_row = loci_row,
               alignment_params = params, program_paths = program_paths, iqtree.num_threads = iqtree.num_threads, iqtree.num_quartets = iqtree.num_quartets, 
               mc.cores = num_of_cores)
    }
    
    # collate the test statistics bootstrap info into 1 file
    print("collate test statistics from bootstraps")
    p_value_df <- collate.empirical.bootstraps(directory = alignment_folder, file.name = "testStatistics", id = loci_name, output.file.name = collated_ts_file)
    # add the column with the bootstrap replicates and "alignment"
    new_bootstrap_ids <- p_value_df$bootstrap_replicate_id # copy col
    aln_id <- grep(new_bootstrap_ids[!grepl("bootstrapReplicate",new_bootstrap_ids)],new_bootstrap_ids) # get which element of col is the alignment
    new_bootstrap_ids[aln_id] <- "alignment"
    p_value_df$bootstrap_replicate_id <- new_bootstrap_ids
    # Get the single row that's just the alignment
    op_p_value_df <- p_value_df[p_value_df$bootstrap_replicate_id == "alignment",]
    
    # Calculate the p-values and add them to the original test statistic dataframe
    print("calculate p values")
    # Calculate the p_values of the variables of interest
    # Calculate the p-values for each test statistic
    op_p_value_df$tree_proportion_p_value         <- calculate.2tailed.p_value(p_value_df$tree_proportion, p_value_df$bootstrap_replicate_id)
    op_p_value_df <- op_p_value_df[,c('dataset','loci','bootstrap_replicate_id','n_taxa','n_sites','alignment_file','tree_proportion','tree_proportion_p_value')]
    # Output the p-values file
    write.csv(op_p_value_df,file = p_value_file, row.names = FALSE)
  }
}





output.tree.proportion.csv <- function(alignment_path, iqtree_path, splitstree_path, bootstrap_id, iqtree.num_threads, iqtree.model = "MFP", alphabet, dataset_name){
  # Calculate tree proportion and output as csv file
  print("output.tree.proportion.csv")
  print(alignment_path)
  # extract the alignment folder from the alignment path
  alignment_folder <- paste0(dirname(alignment_path),"/")
  output_id <- gsub(".nex","",basename(alignment_path))
  # Create some folder and filenames
  if (bootstrap_id == "alignment"){
    # Get the alignment name and remove the extension to get the loci name
    loci_name <- gsub(".nex","",basename(alignment_path))
    # Extract the dataset name (basename of alignment folder: element after last "/" in alignment_folder)
    dataset <- dataset_name
  } else {
    # If the alignment is a bootstrap replicate, need to remove the bootstrap rep number to get the loci name
    loci_name <- gsub(".nex","",basename(alignment_path)) # get the basis of the loci name
    loci_list <- unlist(strsplit(loci_name, "_")) # break the alignment name into bits
    max_ind <- grep("bootstrapReplicate",loci_list) - 1 # which ind is the bootstrapReplicate at?
    loci_list <- loci_list[1:max_ind] # get only parts of alignment name
    loci_name <- paste(loci_list,collapse="_") # squash the alignment name together
    # Extract the dataset name (basename of alignment folder: element after second-last "/" in alignment_folder)
    dataset <- dataset_name
  }
  
  if (bootstrap_id == "alignment"){
    # If this is not a bootstrap replicate, you need to create a folder to store the program logs in
    # Otherwise they will get overwritten for each locus you run - this means you can keep them for late
    log_folder <- paste0(alignment_folder,loci_name,"/")
    # make a rep id - e.g. an id to store in output df so you can identify what's a simulation and what's from empirical data
    rep_id <- loci_name
  } else {
    # If the run is a bootstrap replicate, you can just save the information in its folder (only 1 alignment per folder in bootstrap replicate folders)
    log_folder <- alignment_folder
    # make a rep id - e.g. an id to store in output df so you can identify what's a simulation and what's from empirical data
    rep_id <- bootstrap_id
  }
  
  
  # If the log file doesn't exist, create it 
  if (dir.exists(log_folder) == FALSE){
    dir.create(log_folder) # create a new folder to store the log files from the executables for this loci in
  }
  
  # Open the nexus file and get the number of taxa and the number of characters 
  n <- read.nexus.data(alignment_path)
  n_taxa <- length(n)
  n_char <- length(unlist(n[1]))
  
  # Run IQ-tree on the alignment (if it hasn't already been run)
  call.IQTREE.empirical(alignment_path, iqtree_path = iqtree_path, iqtree.model, num_threads = iqtree.num_threads)
  initial_iqtree_tree <- paste0(alignment_path,".treefile")
  
  setwd(alignment_folder)
  
  # Run trimmed version of the NeighborNet tree proportion
  # Call the test statistic functions
  initial_iqtree_tree <- paste0(alignment_path,".treefile")
  tp <- tree.proportion(iqpath = iqtree_path, splitstree_path = splitstree_path,
                        path = alignment_path, network_algorithm = "neighbournet", trimmed = TRUE,
                        tree_path = initial_iqtree_tree, run_IQTREE = FALSE, seq_type = alphabet)
  
  
  # Name the test statistics file using the output id (this way if it's a  bootstrap replicate, it adds the replicate number!)
  print(paste0("output results for ",output_id))
  results_file <- paste0(alignment_folder,output_id,"_testStatistics.csv")
  # Make somewhere to store the results
  df_names <- c("dataset", "loci", "bootstrap_replicate_id", "n_taxa", "n_sites", "alignment_file","tree_proportion")
  df <- data.frame(matrix(nrow=0,ncol=length(df_names))) # create an empty dataframe of the correct size
  op_row <- c(dataset, loci_name, rep_id, n_taxa, n_char, alignment_path, tp) # collect all the information
  df <- rbind(df,op_row,stringsAsFactors = FALSE) # place row in dataframe
  names(df) <- df_names # add names to the df so you know what's what
  write.csv(df,file = results_file, row.names = FALSE)
}





# Given a .iqtree file, this function will extract the relevant parameters
get.empirical.simulation.parameters <- function(dotiqtree_file){
  # read in the IQ-TREE file to get substitution model and parameters
  iq_file <- readLines(dotiqtree_file)
  # extract the file name
  ind      <- grep("Input file name:",iq_file)
  op1      <- substr(iq_file[[ind]],18,nchar(iq_file[[ind]]))
  # extract the number of taxa and extract the length of the alignment
  ind         <- grep("Input data:",iq_file)
  input_str   <- iq_file[[ind]] # get the line that contains this info
  input_ls    <- strsplit(input_str," ")
  op2         <- input_ls[[1]][3] # extract number of sequences (number of taxa)
  op3         <- input_ls[[1]][6] # extract number of sites 
  # Extract the model of substitution (same for amino-acid and nucleotide files)
  ind         <- grep("Model of substitution: ",iq_file)
  sub_str     <- iq_file[[ind]] # get the line that contains this info
  sub_ls      <- strsplit(sub_str," ")
  op4         <- sub_ls[[1]][4]
  # Extract information about the sequence alignment
  ind <- grep("Number of constant sites:",iq_file)
  num_lines <- iq_file[c(ind:(ind+3))] # take the four lines listing the number of different kinds of sites
  num_split <- unlist(strsplit(num_lines,":")) # Split the lines at the colon
  num_names <- num_split[c(TRUE,FALSE)] # before the colon = name
  num_vals <- num_split[c(FALSE,TRUE)] # after the colon = value
  num_vals_regx <- regmatches(num_vals, gregexpr("[[:digit:]]+", num_vals)) # extract all the numbers after the colon
  # four strings = four lists of numbers: take first value from each list to get number of sites
  num_vals <- c(num_vals_regx[[1]][1],num_vals_regx[[2]][1],num_vals_regx[[3]][1],num_vals_regx[[4]][1]) 
  
  # check the type of sites - amino acid or DNA
  # if the input is DNA (nucleotide sites), gather that information
  if (input_ls[[1]][7]=="nucleotide"){
    # Extract the rate parameters
    rate1 <- as.numeric(strsplit(iq_file[[grep("A-C",iq_file)]],":")[[1]][2]) # A-C rate (same as code above, but combined 4 lines into 1 line)
    rate2 <- as.numeric(strsplit(iq_file[[grep("A-G",iq_file)]],":")[[1]][2]) # A-G rate
    rate3 <- as.numeric(strsplit(iq_file[[grep("A-T",iq_file)]],":")[[1]][2]) # A-T rate
    rate4 <- as.numeric(strsplit(iq_file[[grep("C-G",iq_file)]],":")[[1]][2]) # C-G rate
    rate5 <- as.numeric(strsplit(iq_file[[grep("C-T",iq_file)]],":")[[1]][2]) # C-T rate
    rate6 <- as.numeric(strsplit(iq_file[[grep("G-T",iq_file)]],":")[[1]][2]) # G-T rate
    
    # Extract the state frequencies
    state_freq_line <- iq_file[[grep("State frequencies",iq_file)]]
    if (state_freq_line == "State frequencies: (equal frequencies)"){
      # If the state frequencies are all equal, assign them all to 0.25 (1/4)
      sf1 <- 0.25 # pi(A) - A freq.
      sf2 <- 0.25 # pi(C) - C freq.
      sf3 <- 0.25 # pi(G) - G freq.
      sf4 <- 0.25 # pi(T) - T freq.
    } else {
      # If the state frequencies are not all equal, extract what they are
      sf1 <- as.numeric(strsplit(iq_file[[grep("pi\\(A\\)",iq_file)]],"=")[[1]][2]) # pi(A) - A freq. Remember to double backslash to escape before brackets
      sf2 <- as.numeric(strsplit(iq_file[[grep("pi\\(C\\)",iq_file)]],"=")[[1]][2]) # pi(C) - C freq.
      sf3 <- as.numeric(strsplit(iq_file[[grep("pi\\(G\\)",iq_file)]],"=")[[1]][2]) # pi(G) - G freq.
      sf4 <- as.numeric(strsplit(iq_file[[grep("pi\\(T\\)",iq_file)]],"=")[[1]][2]) # pi(T) - T freq.
    }
    
    # Extract model of rate heterogeneity
    mrh1      <- strsplit(iq_file[[grep("Model of rate heterogeneity:",iq_file)]],":")[[1]][2] # Extract model of rate heterogeneity 
    mrh2      <- as.numeric(strsplit(iq_file[[(grep("Model of rate heterogeneity:",iq_file)+1)]],":")[[1]][2]) # Line after the "model of rate heterogeneity" varies - extract it regardless of what it is 
    mrh2_name <- strsplit(iq_file[[(grep("Model of rate heterogeneity:",iq_file)+1)]],":")[[1]][1] # As the line varies, extract the name for the output dataframe
    mrh2_name <- gsub(" ","_",mrh2_name) # change the name to be easy to parse
    
    # make a list of the rows for the output dataframe
    names <- c("file_name","sequence_type","n_taxa","n_sites",num_names,"substitution_model", "A-C_rate","A-G_rate","A-T_rate","C-G_rate","C-T_rate","G-T_rate","A_freq","C_freq",
               "G_freq","T_freq","model_of_rate_heterogeneity","model_of_rate_heterogeneity_line2_name","model_of_rate_heterogeneity_line2_value")
    # Make a list of the output rows for the output dataframe
    op <- c(op1,"DNA",op2,op3,num_vals,op4,rate1,rate2,rate3,rate4,rate5,rate6,sf1,sf2,sf3,sf4,mrh1,mrh2_name,mrh2)
    # Create the output dataframe
    op_df <- data.frame(names,op, stringsAsFactors = FALSE)
    # Name the columns
    names(op_df) <- c("parameter","value")
    
    # Create the rate matrix Q
    Q_start <- grep("Rate matrix Q:",iq_file)+2
    Q_end   <- Q_start+3
    # Create the columns
    c1 <- c("A","C","G","T")
    c2 <- c()
    c3 <- c()
    c4 <- c()
    c5 <- c()
    # For each row in the iqtree file rate matrix
    for (i in Q_start:Q_end){
      # Split the row
      row <- strsplit(iq_file[[i]]," ")[[1]]
      row <- row[str_detect(row,"([0-9])")] # take only the numeric elements of the vector
      # Add the resulting values to the relevant columns
      c2 <- c(c2,as.numeric(row[1])) # convert to numeric so can use the numbers more easily later
      c3 <- c(c3,as.numeric(row[2]))
      c4 <- c(c4,as.numeric(row[3]))
      c5 <- c(c5,as.numeric(row[4]))
    }
    # Create a dataframe of the rate matrix Q
    q_df <- data.frame(c1,c2,c3,c4,c5, stringsAsFactors = FALSE)
    #Rename the columns
    names(q_df) <- c("nucleotide","A","C","G","T")
    
    # Check if the model for rate heterogeneity is uniform
    mrh1_check <- gsub(" ","",mrh1)
    if (mrh1_check=="Uniform"){
      # If the model for rate heterogeneity is uniform, don't need to create a matrix for discrete gamma rate categories
      g_df <- "Uniform"
    } else {
      # If the model isn't uniform, need to create a matrix to collect and store the gamme category information
      #Create the matrix for discrete gamma categories
      g_start <- grep(" Category",iq_file)+1 # get the index for the first line of the gamma categories matrix
      empty   <- which(iq_file=="") # get indexes of all empty lines
      empty   <- empty[empty>g_start] # get empty lines above gamma categories matrix
      g_end   <- empty[1]-1 # get end index for gamma categories matrix (one less than next empty line)
      end_line <- iq_file[g_end]
      # if the end isn't an empty line, subtract one from the end count 
      # to exclude lines like "Relative rates are computed as MEAN of the portion of the Gamma distribution falling in the category."
      # to see if this is what's happening, check whether the line starts with a numeric section (i.e. a category for the gamma rate)
      check_line <- length(strsplit(strsplit(end_line, "        " )[[1]][1]," ")[[1]])
      if (check_line > 3){
        # If the check_line is longer than 3 characters, it won't be a group for the gamma categories but an instruction
        # Instructions can be excluded from the gamma matrix (but categories can't)
        g_end = g_end - 1
      }
      # Start collecting info for the matrix
      g1 <- c() # initialise columns to store data in
      g2 <- c()
      g3 <- c()
      # Iterate through rows in gamma matrix
      for (i in g_start:g_end){
        row <- strsplit(iq_file[[i]]," ")[[1]] # split the rows on the (large amount of) " "'s in the middle
        row <- row[str_detect(row,"([0-9])")] # take only the numeric elements of the vector
        g1 <- c(g1,as.numeric(row[1])) # add the values to the columns
        g2 <- c(g2,as.numeric(row[2]))
        g3 <- c(g3,as.numeric(row[3]))
      }
      g_df <- data.frame(g1,g2,g3, stringsAsFactors = FALSE) # create a dataframe of the information
      names(g_df) <- c("category","relative_rate","proportion") # name the columns
    }
    
    # Create a list of the three dataframes
    # This will be the output 
    params <- list(op_df,g_df,q_df)
    # Name the parameters so they're easy to access once you've outputted the data
    names(params) <- c("parameters","gamma_categories","Q_rate_matrix")
    
  } else if (input_ls[[1]][7]=="amino-acid"){ # alternatively if the data is amino acid sites
    # Extract model of rate heterogeneity
    mrh1      <- strsplit(iq_file[[grep("Model of rate heterogeneity:",iq_file)]],":")[[1]][2] # Extract model of rate heterogeneity 
    mrh2      <- as.numeric(strsplit(iq_file[[(grep("Model of rate heterogeneity:",iq_file)+1)]],":")[[1]][2]) # Line after the "model of rate heterogeneity" varies - extract it regardless of what it is 
    mrh2_name <- strsplit(iq_file[[(grep("Model of rate heterogeneity:",iq_file)+1)]],":")[[1]][1] # As the line varies, extract the name for the output dataframe
    # Extract state frequencies
    sf1      <- strsplit(iq_file[[grep("State frequencies:",iq_file)]],":")[[1]][2]
    
    # make a list of the rows for the output dataframe
    names <- c("file_name","sequence_type","n_taxa","n_sites",num_names,"substitution_model","model_of_rate_heterogeneity","model_of_rate_heterogeneity_line2_name",
               "model_of_rate_heterogeneity_line2_value","state_frequencies")
    # Make a list of the output rows for the first output dataframe
    op <- c(op1,"amino-acid",op2,op3,num_vals,op4,mrh1,mrh2_name,mrh2,sf1)
    op_df <- data.frame(names,op, stringsAsFactors = FALSE)
    names(op_df) <- c("parameter","value")
    
    # Check whether a gamma matrix is needed
    mrh1_check <- gsub(" ","",mrh1)
    if (mrh1_check=="Uniform"){
      # If the model for rate heterogeneity is uniform, don't need to create a matrix for discrete gamma rate categories
      g_df <- "Uniform"
    } else {
      #Create the matrix for discrete gamma categories
      g_start <- grep(" Category",iq_file)+1 # get the index for the first line of the gamma categories matrix
      empty   <- which(iq_file=="") # get indexes of all empty lines
      empty   <- empty[empty>g_start] # get empty lines above gamma categories matrix
      g_end   <- empty[1]-1 # get end index for gamma categories matrix (one less than next empty line)
      end_line <- iq_file[g_end]
      # if the end isn't an empty line, subtract one from the end count 
      # to exclude lines like "Relative rates are computed as MEAN of the portion of the Gamma distribution falling in the category."
      # to see if this is what's happening, check whether the line starts with a numeric section (i.e. a category for the gamma rate)
      check_line <- length(strsplit(strsplit(end_line, "        " )[[1]][1]," ")[[1]])
      if (check_line > 3){
        # If the check_line is longer than 3 objects, it won't be a group for the gamma categories but an instruction
        # Eg. the relative rates line has 17 objects, because each word in the sentence is an object
        # Instructions can be excluded from the gamma matrix (but categories can't)
        g_end = g_end - 1
      }
      # initialise columns to store data in
      g1 <- c()
      g2 <- c()
      g3 <- c()
      # Iterate through rows in gamma matrix
      for (i in g_start:g_end){
        row <- strsplit(iq_file[[i]],"        ") # split the rows on the long strong of 0's in the middle
        g1 <- c(g1,as.numeric(row[[1]][1])) # add the values to the columns
        g2 <- c(g2,as.numeric(row[[1]][2]))
        g3 <- c(g3,as.numeric(row[[1]][3]))
      }
      g_df <- data.frame(g1,g2,g3, stringsAsFactors = FALSE) # create a dataframe of the information
      names(g_df) <- c("category","relative_rate","proportion") # name the columns 
    }
    
    # Check whether state frequencies are needed
    sf1_squashed <- gsub(" ","",sf1)
    if (sf1_squashed == "(empiricalcountsfromalignment)"){
      # Get starting line for frequencies
      start_ind <- grep("State frequencies:",iq_file) + 2
      # Take the 20 lines containing AA frequencies
      freq_lines <- iq_file[start_ind:(start_ind+19)]
      # Split up the frequency lines into the label and the frequency
      freq_split <- unlist(strsplit(freq_lines,"="))
      # Get the frequency
      freq_nums <- freq_split[c(FALSE,TRUE)]
      # Remove any spaces (from IQTree formatting)
      freq_nums <- gsub(" ","",freq_nums)
      # Get corresponding AA letter
      freq_names <- freq_split[c(TRUE,FALSE)]
      # Remove IQTree formatting
      freq_names <- gsub("pi\\(","",freq_names)
      freq_names <- gsub("\\)","",freq_names)
      freq_names <- gsub(" ","", freq_names)
      # Create a nice dataframe
      f_df <- data.frame("amino_acid" = freq_names,
                         "frequency" = freq_nums,
                         stringsAsFactors = FALSE)
    } else {
      f_df <- "State frequencies from model"
    }
    
    # Create a list of the dataframes - this will be the output
    params <- list(op_df,g_df,f_df)
    # Name the parameters so they're easy to access once you've outputted the data
    names(params) <- c("parameters","gamma_categories", "frequency")
  }
  
  # Now the information has been collected, create an output dataframe
  # make the output dataframe
  return(params)
}





# Function to calculate one empirical bootstrap for an empirical alignment
run.one.tpts.bootstrap <- function(bootstrap_id, empirical_alignment_path, empirical_alignment_loci_name, empirical_alignment_alphabet, empirical_alignment_model, empirical_alignment_dataset,
                                   alignment_params, iqtree_path, splitstree_path, iqtree.num_threads){
  # Remove all references to empirical_alignment_row
  
  print("in do1.empirical.parametric.bootstrap")
  print(empirical_alignment_path)
  print(bootstrap_id)
  # Create the folder for this replicate, gather and create filenames
  loci_name <- empirical_alignment_loci_name
  bootstrap_name <- paste0(loci_name,"_",bootstrap_id) # this will be the name of the alignment
  bootstrap_folder <- paste0(dirname(empirical_alignment_path),"/",bootstrap_name,"/") # folder to store results from this bootstrap in
  # If the bootstrap folder doesn't exist, create it
  if (dir.exists(bootstrap_folder) == FALSE){
    dir.create(bootstrap_folder)
  }
  shuffled_alignment_path <- paste0(bootstrap_folder,bootstrap_name,"_shuffled_noGaps.nex")
  bootstrap_alignment_path <- paste0(bootstrap_folder,bootstrap_name,".nex")
  empirical_alignment_tree_path <- paste0(empirical_alignment_path,".treefile")
  empirical_alignment_tree <- read.tree(empirical_alignment_tree_path)
  
  # First generate completely new DNA using the params
  # Create an alignment for this replicate using the alignment params - name will be loci_bootstrapReplicateXXXX
  if (file.exists(bootstrap_alignment_path) == FALSE) {
    if (empirical_alignment_alphabet == "dna"){
      new_aln <- generate.DNA.alignment(alignment_params, empirical_alignment_tree)
    } else if (empirical_alignment_alphabet == "protein"){
      new_aln <- generate.AA.alignment(alignment_params, empirical_alignment_tree)
    }
    
    # Second, randomise sites in new_aln (otherwise masking will mean the gamma categories disproportionately get affected)
    print("shuffling alignment")
    #   Turn the phyDat into a dataframe and shuffle the rows using sample()
    #   Rows are shuffled as columns represent sequences - shuffling the rows keeps the relationships between species (shown by rows) 
    #   but rearranges the order the relationships occur in
    new_aln_df <- as.data.frame(new_aln) # turn the new alignment into a dataframe
    new_aln_df <- new_aln_df[sample(nrow(new_aln_df)),] # sample the rows randomly (this will randomly distribute the gamma categories throughout)
    rownames(new_aln_df) <- 1:nrow(new_aln_df) # reset the row names as 1:nrows (they got mixed up when the sampling occurred)
    # Turn the dataframe back into an alignment:
    if (empirical_alignment_alphabet == "dna"){
      new_aln_shuffled <- as.phyDat(new_aln_df, type = "DNA")
      #write the shuffled alignment to disk
      write.phyDat(new_aln_shuffled, file = shuffled_alignment_path, format = "nexus", interleave = TRUE)
      # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
      nexus_edit <- readLines(shuffled_alignment_path) # open the new nexus file
      ind <- grep("BEGIN DATA",nexus_edit)+2 # find which line
      nexus_edit[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE;" # replace the line
      writeLines(nexus_edit,shuffled_alignment_path) # output the edited nexus file
    } else if (empirical_alignment_alphabet == "protein"){
      new_aln_shuffled <- as.phyDat(new_aln_df, type = "AA")
      #write the shuffled alignment to disk
      write.phyDat(new_aln_shuffled, file = shuffled_alignment_path, format = "nexus", interleave = TRUE)
      # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
      nexus_edit <- readLines(shuffled_alignment_path) # open the new nexus file
      ind <- grep("BEGIN DATA",nexus_edit)+2 # find which line
      nexus_edit[ind] <- "  FORMAT DATATYPE=PROTEIN MISSING=? GAP=- INTERLEAVE;" # replace the line
      writeLines(nexus_edit,shuffled_alignment_path) # output the edited nexus file
    }
    
    # Third, mask each alignment with the gaps and unknown characters from the original sequence 
    print("masking alignment")
    # for each alignment:
    #     - copy the sequence out from the new alignment: temp <- as.numeric(new_aln$X)
    #     - replace the non 18s with the generated sequence of the right length: temp[which(new_aln$X !=18)] <- new_seq
    #     - replace the new seq into the new aln: new_aln$X <- temp
    # Open the new alignment as a nexus file
    n <- read.nexus.data(empirical_alignment_path)
    n_new <- read.nexus.data(shuffled_alignment_path)
    # Get the names of all the sequences
    seq_names <- names(n_new)
    print(paste0("number of names: ",length(seq_names)))
    print(paste0("number of unique names: ",length(unique(seq_names))))
    # Iterate through the names
    for (seq_name in seq_names){
      original_seq <- n[[seq_name]] # get the original empirical sequence
      new_seq <- n_new[[seq_name]] # get the new simulated sequence that has the same name
      gap_inds <- which(original_seq == "-") # find out which sites are a gap in the original alignment
      unknown_inds <- which(original_seq == "?") # find out which sites are unknown in the original alignment
      new_seq[gap_inds] <- "-" # add the gaps into the simulated alignment
      new_seq[unknown_inds] <- "?" # add the unknowns into the simulated alignment
      n_new[[seq_name]] <- new_seq
    }
    # Output the final alignment (same parameters and gaps as input alignment) as a nexus file
    print("output nexus file")
    if (empirical_alignment_alphabet == "dna"){
      # write out the nexus with the gaps added
      write.nexus.data(n_new,file = bootstrap_alignment_path, format = "dna", interleaved = TRUE)
      # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
      nexus_edit <- readLines(bootstrap_alignment_path) # open the new nexus file
      ind <- grep("BEGIN DATA",nexus_edit)+2 # find which line
      nexus_edit[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE;" # replace the line
      writeLines(nexus_edit,bootstrap_alignment_path) # output the edited nexus file
    } else if (empirical_alignment_alphabet == "protein"){
      # write out the nexus with the gaps added
      write.nexus.data(n_new,file = bootstrap_alignment_path, format = "protein", interleaved = TRUE)
      # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
      nexus_edit <- readLines(bootstrap_alignment_path) # open the new nexus file
      ind <- grep("BEGIN DATA",nexus_edit)+2 # find which line
      nexus_edit[ind] <- "  FORMAT DATATYPE=PROTEIN MISSING=? GAP=- INTERLEAVE;" # replace the line
      writeLines(nexus_edit,bootstrap_alignment_path) # output the edited nexus file
    }
  }
  # Run all the test statistics
  # bootstrap_id will be "bootstrapReplicateXXXX" where XXXX is a number
  print("run test statistics")
  
  # If the alignment comes from Vanderpool 2020, there is no specified best model (the original best model was determined using MFP in IQ-Tree)
  # Need to pass the model from the empirical alignment into the bootstrap replicates (this will be from the .iqtree file)
  # Alternatively, if the alignment comes from Misof 2014 or 1KP, we know the best model and we should use that model for both the original 
  #   alignment and the parametric bootstrap
  if (empirical_alignment_model == "MFP"){
    MFP_model <- alignment_params$parameters[9,2]
    output.tree.proportion.csv(alignment_path = bootstrap_alignment_path, iqtree_path = iqtree_path, splitstree_path = splitstree_path, bootstrap_id = bootstrap_id, 
                               iqtree.num_threads, iqtree.model = MFP_model, alphabet = empirical_alignment_alphabet, dataset_name = empirical_alignment_dataset)
  } else {
    output.tree.proportion.csv(alignment_path = bootstrap_alignment_path, iqtree_path = iqtree_path, splitstree_path = splitstree_path, bootstrap_id = bootstrap_id, 
                               iqtree.num_threads, iqtree.model = empirical_alignment_model, alphabet = empirical_alignment_alphabet, dataset_name = empirical_alignment_dataset)
  }
}





# Function to collect bootstraps from empirical run
collate.empirical.bootstraps <- function(directory, file.name, id, output.file.name){
  # Collect all the folders within the directory
  all_files <- list.files(directory, full.names = FALSE, recursive = TRUE)
  # Now reduce that to only get folders for the particular id of interest
  id_files <- all_files[grep(id,all_files)]
  # Get the files of interest and their full file paths
  csv_paths <- id_files[grep(file.name,id_files)]
  # Remove collated test statistic file from the list of csvs
  collated_name <- basename(output.file.name)
  csv_paths <- csv_paths[!grepl(collated_name,csv_paths)]
  if (length(csv_paths) > 0){
    csv_paths <- paste0(dirname(directory),"/",basename(directory),"/",csv_paths)
    # Set the number of rows in the dataframe (will equal the number of csv files: one per simulation + one for the original empirical alignment)
    num_rows <- length(csv_paths)
    # Open all the csv files, store the results as a list
    output_list <- lapply(csv_paths, read.csv, stringsAsFactors = FALSE)
    # Reduce the dataframe from a list into a matrix
    output_df <- Reduce(rbind, output_list)
    # save the output dataframe
    write.csv(output_df, file = output.file.name, row.names = FALSE)
    return(output_df)
  } else {
    return(NA)
  }
}




# Given two vectors (one of test statistic values, and one of ids), calculates the p-value for that alignment
calculate.2tailed.p_value <- function(value_vector,id_vector){
  p_value_df <- data.frame(value_vector,id_vector, stringsAsFactors = FALSE)
  names(p_value_df) <- c("value","id")
  alignment_value <- p_value_df[which(p_value_df$id == "alignment"),1]
  if (is.na(alignment_value) == TRUE){
    # If the alignment value is NA, can't calculate a score
    # If it wasn't possible to calculate a score (will usually be a PHI score) for the alignment, output an NA
    p_value_2tail <- NA
  } else {
    # If it's possible to calculate a p-value, calculate one
    # Exclude NA rows
    p_value_df <-  p_value_df[!is.na(p_value_df$value),]
    # Find the number of bootstrap replicates and where the actual alignment value is located
    num_rows <- nrow(p_value_df) # number of bootstrap replicates + alignment value
    p_value_df <- p_value_df[order(p_value_df$value),] # order values from smallest to largest
    alignment_row <- which(p_value_df$id == "alignment") # find the ranking of the alignment value
    alignment_value <- p_value_df[alignment_row,1] # find the alignment's test statistic value
    # check whether there are other values that are the same as the alignment value
    identical_df <- subset(p_value_df,value == alignment_value)
    # if there are identical values, you don't know where the alignment actually falls within that list
    if (nrow(identical_df)>1){
      # get all the indexes of identical values
      identical_inds <- grep(alignment_value,p_value_df$value)
      # pick an ind at random
      random_identical_row <- sample(identical_inds,1)
      # For left tail probability: want to find the number of observations less than or equal to the alignment value, then divide by the number of bootstrap observations
      # If there are 8 rows and the alignment is the 5th row, then there will be 5 alignments less than or equal to the alignment value
      p_value_left <- random_identical_row/num_rows
      # For right tail probability: want to find the number of observations greater than or equal to the alignment value, then divide by the number of bootstrap observations
      # If there are 8 rows and the alignment is the 5th row, then there will be (8-5)=3 alignments greater than or equal to the alignment value
      p_value_right <- (num_rows - random_identical_row)/num_rows
    } else if (nrow(identical_df) == 1){
      # else, simply calculate the p value using the formula 
      # For left tail probability: want to find the number of observations less than or equal to the alignment value, then divide by the number of bootstrap observations
      # If there are 8 rows and the alignment is the 5th row, then there will be 5/8 alignments less than or equal to the alignment value
      p_value_left <- alignment_row/num_rows
      # For right tail probability: want to find the number of observations greater than or equal to the alignment value, then divide by the number of bootstrap observations
      # If there are 8 rows and the alignment is the 5th row, then there will be (8-5)=3 alignments greater than or equal to the alignment value
      p_value_right <- (num_rows-alignment_row)/num_rows
    }
    # To find two tailed probability, multiply the lower of those values by 2
    p_value_2tail <- 2*min(p_value_left,p_value_right) 
  }
  
  # return the p-value
  return(p_value_2tail)
}






