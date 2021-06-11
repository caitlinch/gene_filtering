# empirical_treelikeness/code/0_Pease2016_data_formatting.R
# Caitlin Cherryh, 2021
## This script recreates the 100kb genomic window alignments used to estimate gene trees in Pease et al (2016)

## Pease et al (2016) paper: 
#       Pease, J. B., Haak, D. C., Hahn, M. W., Moyle, L. 2016, Phylogenomics reveals three sources of adaptive variation during a rapid radiation, PLOS Biology, 14(2):e1002379
#       https://doi.org/10.1371/journal.pbio.1002379

## This script:
# 1. Uses mvftools to read in the HQ alignment mvf.gz file from the Pease et al (2016) datadryad and output 100kb genomic windows as .phy alignments 
# 2. Collect information about the RAxML run and the models used for each window to estimate an ML gene tree
# 3. Identify which temporary .phy alignment matches to which window
# 4. Copy the 2745 alignments used for the Pease et al (2016) ASTRAL gene tree analysis into a fresh folder, labelled with the corresponding genomic window 
# 5. Output a summary of which temporary .phy alignment matches to which window

## This script requires:
#     - RAxML: https://github.com/stamatak/standard-RAxML
#     - mvftools: https://github.com/peaselab/mvftools
#     - the Pease_etal_Tomato29acc_HQ.mvf.gz file from the datadryad for Pease et al (2016): https://datadryad.org/stash/dataset/doi:10.5061/dryad.182dv
#     - the Pease_etal_TomatoPhylo_100kbTrees.txt file from the datadryad for Pease et al (2016): https://datadryad.org/stash/dataset/doi:10.5061/dryad.182dv

#### Specify parameters ####
# raxml_path                    <- path to RAxML executable
# mvftools_path                 <- path to mvftools/mvftools.py program
# Pease2016_100kb_windows_path  <- path to the Pease_etal_TomatoPhylo_100kbTrees.txt file from the datadryad for Pease et al (2016)
# Pease2016_HQ_mvf.gz           <- path to the the Pease_etal_Tomato29acc_HQ.mvf.gz file from the datadryad for Pease et al (2016)
# Pease2016_alignments_folder   <- location of the alignment mvf.gz file, where mvftools will be run
# alignment_output_folder       <- output folder for the identified, copied and renamed genomic window alignments
# summary_output_folder         <- output folder for .csv files 
# run.mvf.tools                 <- whether to run mvftools InferTree to create alignments for genomic windows. 
#                                  TRUE if want to run, FALSE if have already run/do not want to run

raxml_path <- "/Users/caitlincherryh/Documents/Executables/standard-RAxML/raxmlHPC-AVX2"
mvftools_path <- "/Users/caitlincherryh/Documents/Executables/mvftools/mvftools.py"
Pease2016_100kb_windows_path <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/Pease_etal_TomatoPhylo_100kbTrees.txt"
Pease2016_HQ_mvf.gz <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/Pease_alignments/Pease_etal_Tomato29acc_HQ.mvf.gz"
Pease2016_alignments_folder <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/Pease_alignments/"
alignment_output_folder <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/all_window_alignments/"
summary_output_folder <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/"
run.mvf.tools = FALSE

#### Open packages ####
library(ape)
library(phangorn)

#### Functions ####
# For one temp .phy file, open it and extract information about the alignment and corresponding tree
Pease.extract.info <- function(datetime, mvftools_output_files){
  # Extract files using datetime
  all_dt_files <- grep(datetime, mvftools_output_files, value = TRUE)
  alignment_file <- grep("_temp.phy",all_dt_files, value = TRUE)
  raxml_info <- grep("RAxML_info", all_dt_files, value = TRUE)
  best_tree_file <- grep("RAxML_bestTree", all_dt_files, value = TRUE)
  # Extract information about alignment
  p <- read.dna(alignment_file, format = "sequential")
  aligndepth <- dim(p)[1]
  alignlength <- dim(p)[2]
  # Extract trees
  best_tree_text <- readLines(best_tree_file)
  # Open the RAxML info file and extract information about the model
  lines <- readLines(raxml_info)
  ind <- grep("Substitution Matrix", lines)
  substitution_matrix <- gsub(" ", "", strsplit(lines[ind],":")[[1]][2])
  ind <- grep("RAxML was called as follows:", lines)
  system_call <- lines[ind+2]
  system_call_vector <- strsplit(system_call, " ")[[1]]
  system_model <- system_call_vector[grep("-m", system_call_vector) + 1]
  # Return a row with information about the alignment
  op_row <- c(datetime, alignlength, aligndepth, substitution_matrix, system_model, best_tree_text, 
              alignment_file, raxml_info, best_tree_file)
  names(op_row) <- c("alignment_datetime", "alignlength", "aligndepth",  "substitution_matrix", "RAxML_model_input", "best_ML_tree", 
                     "alignment_file", "raxml_info_file", "treefile")
  return(op_row)
}

# Match up one window used in Pease 2016 for the ASTRAL study with one temporary .phy alignment
Pease.get.astral.window <- function(index, complete_windows_df, infertree_df, alignment_info_df, copy.alignment = FALSE, output_directory = NA){
  # Get one window from the Pease2016 run from the set of windows with all 29 species
  complete_windows_row <- complete_windows_df[index,]
  print(paste0("Window: Contig = ", complete_windows_row$`#contig`, ", Window start = ", complete_windows_row$windowstart))
  # Get the equivalent row from the local run
  local_run_row <- infertree_df[((infertree_df$`#contig` == as.integer(complete_windows_row$`#contig`)) &
                                   (as.integer(infertree_df$windowstart) == as.integer(complete_windows_row$windowstart)) &
                                   (as.integer(infertree_df$windowsize) == as.integer(complete_windows_row$windowsize))),]
  loci_name <- paste0("c", local_run_row$`#contig`, "_s", local_run_row$windowstart, "_100kb_windows")
  local_run_tree <- read.tree(text = local_run_row$tree)
  # Compare trees until you get the right tree
  match_tree_ind <- NA
  for (i in 1:nrow(alignment_info_df)){
    # Pick one row from the alignment info and read the tree in
    aln_row <- alignment_info_df[i,]
    temp_tree <- read.tree(text = aln_row$best_ML_tree)
    # If the trees have the same number of tips, compare them
    if (Ntip(local_run_tree) == Ntip(temp_tree)){
      tree_distances <- treedist(local_run_tree, temp_tree)
      spr <- SPR.dist(local_run_tree, temp_tree)
      compare_tree <- all(c(tree_distances["symmetric.difference"], tree_distances["path.difference"], spr) == 0)
      # If the trees match, store this location
      if (compare_tree == TRUE){
        match_tree_ind <- i
        print(paste0("Found tree match in alignment_info_df row = ", i))
      } 
    }
  }
  
  # Now each alignment has been tested (by checking to identify the correct tree)
  if (is.na(match_tree_ind) == FALSE){
    # If there was a match, that means that the tree was found and therefore the alignment was found!
    # Associate the window and the datetime together
    w_aln_row <- alignment_info_df[match_tree_ind, ]
    
    # Output alignment if required
    if (copy.alignment == TRUE){
      # Specify current location and location to copy alignment to
      old_alignment_path <- w_aln_row$alignment_file
      new_alignment_path <- paste0(output_directory, loci_name, ".fasta")
      copied = TRUE
      # Open alignment as phylip
      p <- read.dna(old_alignment_path, format = "sequential")
      # Save alignment as fasta
      write.dna(p, file = new_alignment_path, format = "fasta", colsep = "", nbcol = 10, colw = 10)
      
    } else {
      copied <- FALSE
      new_alignment_path <- NA
    }
    
    op_row <- c("Pease2016", w_aln_row$alignment_datetime, loci_name, local_run_row$`#contig`, local_run_row$windowstart, local_run_row$windowsize, 
                local_run_row$alignlength, local_run_row$aligndepth, w_aln_row$substitution_matrix, w_aln_row$RAxML_model_input, 
                w_aln_row$alignment_file, w_aln_row$raxml_info_file, w_aln_row$treefile, match_tree_ind, Ntips(local_run_tree), TRUE,
                copied, new_alignment_path)
    
    
  } else if (is.na(match_tree_ind) == TRUE){
    # If match_tree_ind is na, it means there was no match for this tree
    # Return "NA" where the match information would be
    op_row <- c("Pease2016", NA, loci_name, local_run_row$`#contig`, local_run_row$windowstart, local_run_row$windowsize, 
                local_run_row$alignlength, local_run_row$aligndepth, NA, NA, NA, 
                NA, NA, NA, match_tree_ind, FALSE,
                FALSE, NA)
  }
  # Give the output vector some nice labels
  names(op_row) <- c("dataset", "datetime", "loci_name", "#contig", "windowstart", "windowsize",
                     "alignlength", "aligndepth", "substitution_matrix", "RAxML_model_input", 
                     "alignment_file", "raxml_info_file", "treefile", "InferTree_info_row_number", "N_tips_InferTree_tree", "matched",
                     "alignment_copied", "alignment_copy_location")
  return(op_row)
}

#### Code body ####
# Set working directory to Pease2016 alignment folder so that the windows will be unfolded there
setwd(Pease2016_alignments_folder) 
mvftools_command <- paste0("python3 ", mvftools_path,
                           " InferTree --mvf ",Pease2016_HQ_mvf.gz,
                           " --windowsize 100000 --out trees100k.txt --contig-ids 1", 
                           "--sample-indices 0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29",
                           " --raxml-path ", raxml_path)
if (run.mvf.tools == TRUE){
  system(mvftools_command)
}

# Construct the name of the trees100k.txt file that the mvftools InferTree run created 
my_100kb_windows_path <- paste0(Pease2016_alignments_folder, "trees100k.txt")

# Open the table with information about the windows
w_df <- read.table(Pease2016_100kb_windows_path)
names(w_df)  <- c("#contig", "windowstart", "windowsize", "tree", "topology", "topoid", "alignlength", "aligndepth", "status")
# Get the list of windows used for ASTRAL in Pease 2016 (all alignments with 29 taxa)
complete_windows_df <- w_df[w_df$aligndepth == 29,]

# Extract information about the alignment
# Iterate through the alignment files and extract information
mvftools_output_files <- paste0(Pease2016_alignments_folder, list.files(paste0(Pease2016_alignments_folder)))
alignment_files <- paste0(Pease2016_alignments_folder, sort(grep(".phy", mvftools_output_files, value = TRUE)))
all_datetimes <- gsub("_temp.phy","",gsub("mvftree.","",basename(alignment_files)))
# Run the function using lapply and compile information
aln_info_list <- lapply(all_datetimes, Pease.extract.info, mvftools_output_files)
aln_info_df <- as.data.frame(do.call(rbind, aln_info_list))
alignment_info_df <- aln_info_df
write.csv(aln_info_df, file = paste0(summary_output_folder, "Pease2016_mvftools_InferTree_alignment_information.csv"))

## Match up the windows with the temporary .phy alignmment files
match_list <- lapply(1:nrow(complete_windows_df), Pease.get.astral.window, complete_windows_df = complete_windows_df, infertree_df = infertree_df,
                     alignment_info_df = alignment_info_df, copy.alignment = FALSE, output_directory = alignment_output_folder)
match_df <- as.data.frame(do.call(rbind, match_list))
write.csv(aln_info_df, file = paste0(summary_output_folder, "Pease2016_data_recreation_100kb_windows.csv"))




