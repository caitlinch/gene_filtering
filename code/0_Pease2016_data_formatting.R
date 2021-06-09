# empirical_treelikeness/code/0_Pease2016_data_formatting.R
# Code to take the tomato chromosome alignments from Pease et al (2016), create non-overlapping genomic windows of 100 kb, and save each window as a new alignment

# Caitlin Cherryh, 2021

## Open packages
library(ape)
library(phangorn)

## Specify parameters
raxml_path <- "/Users/caitlincherryh/Documents/Executables/standard-RAxML/raxmlHPC-AVX2"
mvftools_path <- "/Users/caitlincherryh/Documents/Executables/mvftools/mvftools.py"
Pease2016_100kb_windows_path <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/Pease_etal_TomatoPhylo_100kbTrees.txt"
Pease2016_HQ_mvf.gz <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/Pease_alignments/Pease_etal_Tomato29acc_HQ.mvf.gz"
Pease2016_alignments_folder <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/Pease_alignments/ch01/" #location of the alignment mvf.gz file, where mvftools will be run
my_100kb_windows_path <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/Pease_alignments/ch01/trees100k.txt"
alignment_output_folder <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/all_window_alignments/"
summary_output_folder <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/"

# Set working directory to Pease2016 alignment folder so that the windows will be unfolded there
setwd(Pease2016_alignments_folder) 
mvftools_command <- paste0("python3 ", mvftools_path,
                           " InferTree --mvf ",Pease2016_HQ_mvf.gz,
                           " --windowsize 100000 --out trees100k.txt --contig-ids 1", 
                           "--sample-indices 0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29",
                           " --raxml-path ", raxml_path)
# # Uncomment the following line to run InferTrees from mvftools, generating the 100kb window alignments
# system(mvftools_command)

# Open the table with information about the windows
tree_df <- read.table(paste0(Pease2016_alignments_folder, "trees100k.txt"))
names(tree_df) <- c("#contig", "windowstart", "windowsize", "tree", "topology", "topoid", "alignlength", "aligndepth", "status")

# Open the table with information about the windows from your run
my_tree_df <- read.table(my_100kb_windows_path)
names(my_tree_df) <- c("#contig", "windowstart", "windowsize", "tree", "topology", "topoid", "alignlength", "aligndepth", "status")

# Get a list of all the temporary files of the 100kb windows
#all_temp_files <- list.files(paste0(Pease2016_alignments_folder, "raxmltemp/"))
all_temp_files <- list.files(paste0(Pease2016_alignments_folder))
tree_files <- paste0(Pease2016_alignments_folder, sort(grep("RAxML_bestTree", all_temp_files, value = TRUE)))
phy_files <- paste0(Pease2016_alignments_folder, sort(grep(".phy", all_temp_files, value = TRUE)))
raxml_info_files <- paste0(Pease2016_alignments_folder, sort(grep("RAxML_info", all_temp_files, value = TRUE)))

# For each window:
# 1. determine which file contains the window information. Check this file has the right parameters (cross-check against tree_df)
# 2. Name this loci based on the window details
# 3. Rename and copy the file to the output folder
# 4. Open the RAxML info file and extract information about the model 
# 5. Return a row with the window location, window information, whether the cross-check was correct, loci location, model info
process.one.Pease2016.window <- function(index, df, alignment_files, raxml_files, tree_files, output_directory, copy.alignment = TRUE){
  # Use index to select correct files
  row <- df[index,]
  temp_al_file <- alignment_files[index]
  temp_raxml_file <- raxml_files[index]
  temp_tree_file <- tree_files[index]
  ## 1. Check whether this alignment is correct for this window
  # Check whether the filepaths match up
  al_datetime <- gsub("_temp.phy", "", gsub("mvftree.", "", basename(temp_al_file)))
  raxml_datetime <- gsub("RAxML_info.mvftree.", "", basename(temp_raxml_file))
  does.datetime.match <- as.character(identical(al_datetime, raxml_datetime))
  print("Next alignment")
  print(index)
  print(al_datetime)
  print(raxml_datetime)
  print(temp_al_file)
  print(temp_raxml_file)
  print(temp_tree_file)
  # Open the alignment
  p <- read.dna(temp_al_file, format = "sequential")
  # Check whether the alignment information matches up
  does.n_taxa.match <- identical(row$aligndepth, dim(p)[1])
  does.n_sites.match <- identical(row$alignlength, dim(p)[2])
  print(dim(p))
  did.raxml.work <- row$status
  # Check whether the trees are identical
  t1 <- read.tree(text = row$tree) # tree from 100kb.text
  t2 <- read.tree(file = temp_tree_file) # tree from RAxML_bestTree.mvftree.datetime
  rf <- RF.dist(t1,t2)
  wrf <- wRF.dist(t1,t2)
  SPR <- SPR.dist(t1,t2)
  KF <- KF.dist(t1,t2)
  path <- path.dist(t1,t2)
  is.tree.identical <- as.character(all.equal(rf, 0))
  ## 2. Name this loci based on the window details
  # New name = contig_windowStart_windowSize.fasta
  loci_name <-  paste0(gsub("SL2.50","",row$`#contig`), "_s", row$windowstart, "_", "100kb_windows")
  new_al_file <- paste0(output_directory, loci_name, ".fasta")
  ## 3. Output alignment as a fasta file
  if (copy.alignment == TRUE){
    write.dna(p, file = new_al_file, format = "fasta", colsep = "", nbcol = 10, colw = 10)
  }
  ## 4. Open the RAxML info file and extract information about the model
  lines <- readLines(temp_raxml_file)
  ind <- grep("Substitution Matrix", lines)
  substitution_matrix <- gsub(" ", "", strsplit(lines[ind],":")[[1]][2])
  ind <- grep("RAxML was called as follows:", lines)
  system_call <- lines[ind+2]
  system_call_vector <- strsplit(system_call, " ")[[1]]
  system_model <- system_call_vector[grep("-m", system_call_vector) + 1]
  ## 5. Return a row with information about the alignment
  op_row <- c(loci_name, "Pease2016", as.numeric(gsub("SL2.50ch", "", row$`#contig`)), row$`#contig`, row$windowstart, row$windowsize, row$alignlength, row$aligndepth, 
              substitution_matrix, system_model, al_datetime, does.datetime.match, does.n_taxa.match, does.n_sites.match, is.tree.identical, did.raxml.work, 
              rf, wrf, SPR, KF, path, new_al_file)
  names(op_row) <- c("loci_name", "dataset", "contig", "chromosome", "windowstart", "windowsize", "alignlength", "aligndepth", 
                     "substitution_matrix", "RAxML_model_input", "alignment_datetime", "datetime_match", "n_taxa_match","n_sites_match", "tree_match", "raxml_status",
                     "RF_distance", "wRF_distance", "SPR_distance", "KR_distance", "path_distance", "alignment_path")
  return(op_row)
}

# Iterate through the rows of the tree_df and apply the function, then convert results to a dataframe
# To run one alignment: process.one.Pease2016.window(1, tree_df, phy_files, raxml_info_files, tree_files, alignment_output_folder, copy.alignment = TRUE)
list <- lapply(1:nrow(my_tree_df), process.one.Pease2016.window, df = my_tree_df, alignment_files = phy_files, raxml_files = raxml_info_files, tree_files = tree_files, output_directory = alignment_output_folder, copy.alignment = TRUE)
df <- as.data.frame(do.call(rbind,list))
# Output dataframe
df_file <- paste0(summary_output_folder, "Pease2016_DataSummary_100kb_windows_all.csv")
write.csv(df, file = df_file)

# Identify which windows were included in the original study
Pease2016_windows_df <- read.table(Pease2016_100kb_windows_path)
names(Pease2016_windows_df) <- c("#contig", "windowstart", "windowsize", "tree", "topology", "topoid", "alignlength", "aligndepth", "status")
Pease2016_windows_df$loci_name <- paste0(Pease2016_windows_df$`#contig`, "_s", Pease2016_windows_df$windowstart, "_", "100kb_windows")
Pease2016_windows_df$compare <- paste0(Pease2016_windows_df$`#contig`, "_", Pease2016_windows_df$windowstart)
df$compare <- paste0(df$contig, "_", df$windowstart)
# limit to 1st chromosome
Pease_c1_df <- Pease2016_windows_df[Pease2016_windows_df$`#contig` == "1",] 
c1_df <- df[df$contig == "1",]
# Find out which windows are missing
missing_in_Pease <- setdiff(c1_df$compare, Pease_c1_df$compare)
missing_rows <- c1_df[(c1_df$compare %in% missing_in_Pease),]
kept_rows <- c1_df[!(c1_df$compare %in% missing_in_Pease),]
nrow(c1_df) == nrow(missing_rows) + nrow(kept_rows)







