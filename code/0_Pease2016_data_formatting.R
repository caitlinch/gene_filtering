# empirical_treelikeness/code/0_Pease2016_data_formatting.R
# Code to take the tomato chromosome alignments from Pease et al (2016), create non-overlapping genomic windows of 100 kb, and save each window as a new alignment

# Caitlin Cherryh, 2021

# Before you can run this code, you must have used mvftools to convert the HQ alignment from the DataDryad into fasta format
# 1. Download the DataDryad (https://datadryad.org/stash/dataset/doi:10.5061/dryad.182dv)
# 2. Unzip the HQ Alignment in Multisample Variant Format file: Pease_etal_Tomato29acc_HQ.mvf.gz
# 3. Download mvftools (https://github.com/jbpease/mvftools)
# 4. Use mvftools to convert the HQ alignment MVF file to a fasta file for each chromosome by using the following commands in the terminal:
#    cd  ~/Documents/Executables/mvftools/Pease_alignments
#    python3 mvftools.py ConvertMVF2FastaGene --mvf Pease_alignments/Pease_etal_Tomato29acc_HQ.mvf  --output-dir Pease_alignments/Tomato29acc_HQ/ --output-data "dna" 
# There are more instructions in the mvftools documentation if you need more help with this step

# Open libraries
library("ape")

tsv_file <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/DataDryad_data/Pease_etal_TomatoPhylo_GeneTrees/Pease_etal_TomatoPhylo_100kbTrees.txt"
tsv_100kb <- read.delim(tsv_file)
tsv_file <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/DataDryad_data/Pease_etal_TomatoPhylo_GeneTrees/Pease_etal_TomatoPhylo_1MbTrees.txt"
tsv_1Mb <- read.delim(tsv_file)

# Specify the lineages that are allowed in the alignments
# These are the lineages that appear in Supplementary Figure 2E of Pease et al (2016) - the lineages included in the window trees
allowed_lineages <- c('LA3475','SL2.50','LA3124','LA0429','LA0436','LA3909','LA2933','LA1269','LA1589','LA1028','LA1316','LA2172',
                      'LA1322','LA2133','LA1782','LA4117','LA1364','LA2744','LA0107','LA2964','LA0444','LA1358','LA0407','LA1777',
                      'LA0716','LA3778','LA4116','LA2951','LA4126')

# Specify the input folder - the location of the fasta format chromosome files
input_folder <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/Pease_alignments/Tomato29acc_HQ/"
# Specify the output folder - where to store the window alignments
output_folder <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/Pease_alignments/Tomato29acc_HQ_100kb_windows/"

# Collect all the chromosome files
c_fastas <- list.files(input_folder)
# Remove the Sl2.50ch00.fa file - the unplaced scaffolds were not used for the phylogenetic analysis
c_fastas <- grep("ch00", c_fastas, invert = TRUE, value = TRUE)

# Write a function to save the alignment for each window
# Iterate through the dataframe and save one window at a time
write.one.window <- function(index, window_df, chromosome, allowed_lineages){
  # Select row with information for this window using the index
  row <- window_df[index,]
  # Get the portion of the DNAbin chromosome contained by the window
  chromosome_portion <- as.list(as.matrix(chromosome)[,row$window_start:row$window_end])
  # Extract only the lineages you want to keep
  chromosome_portion <- chromosome_portion[c(allowed_lineages)]
  # Create the output file name
  output_fasta_name <- paste0(row$output_folder, row$window_name, ".fa")
  # Save chromosome portion as a fasta file
  write.dna(chromosome_portion, file = output_fasta_name, format = "fasta", colsep = "", nbcol = 10, colw = 10, append = FALSE)
}

# Iterate through the chromosome files and save non-overlapping 100kb windows as alignments in fasta format
for (c in c_fastas){
  # Open the fasta chromosome file
  chromosome <- read.FASTA(file = paste0(input_folder, c), type = "DNA")
  # Extract the name of this chromosome
  c_name <- gsub("\\.fa","",gsub("SL2.50","",c))
  # Create a dataframe of information about the windows
  chromosome_length <- length(chromosome[[1]])
  n_windows <- length(seq(0, chromosome_length, 100000))
  c_df <- data.frame(window_name = paste0(c_name, "_", sprintf("%04d", 1:n_windows)),
                     window_start = seq(0, chromosome_length, 100000), 
                     window_end = c(seq(100000, chromosome_length, 100000), chromosome_length),
                     chromosome_location = paste0(input_folder, c),
                     output_folder = output_folder
                     )
  # Iterate through one window at a time and save the genetic information from that window as a fasta file
  #lapply(1:nrow(c_df), write.one.window, window_df = c_df, chromosome = chromosome, allowed_lineages = allowed_lineages)
  # Write the window df out as a csv
  window_df_output_name <- paste0(dirname(output_folder),"/", "Pease2016_genomicWindows_100kb_", c_name, ".csv")
  write.csv(c_df, file = window_df_output_name)
}


# Investigate how many windows
all_files <- list.files(paste0(dirname(output_folder),"/"))
csv_files <- grep("\\.csv", all_files, value = TRUE)
csv_files <- paste0(dirname(output_folder), "/", csv_files)
csv_list <- lapply(csv_files, read.csv)
windows_df <- as.data.frame(do.call(rbind, csv_list))

# Check how many sites specified in mvf for chromosome 1
mvf_file <- "/Users/caitlincherryh/Downloads/doi_10.5061_dryad.182dv__v1/Pease_etal_Tomato29acc_HQ_ch01.mvf"
mvf <- readLines(mvf_file)
mvf_trimmed <- mvf[grep("1:", mvf)]
mvf_nums <- unlist(strsplit(mvf_trimmed, " "))[c(TRUE,FALSE)]
mvf_nums <- gsub("1:","",mvf_nums)
mvf_nums <- as.numeric(mvf_nums)
# Compare to how many sites there should be based on the alignment length for chromosome 1
all_nums <- seq(924, 98543114, 1)
all_nums[!(all_nums %in% mvf_nums)]

# Try and recreate the windows based on the lengths of the alignments in the 100kb tsv
ch01_tsv_100kb <- tsv_100kb[tsv_100kb$X.contig == 1,]

# This is how long all the alignments from chromosome 1 added up are: 
sum(ch01_tsv_100kb$alignlength) # 19396432

# This is how long chromosome 1 is
chromosome_length # 21918490

# There's an inconsistency between the length of the chromosome and the length of all the alignments
# Can't just take 100kb windows because there isn't every base present - windows will be incorrect lengths
# Need to map back onto reference I suppose...
# Leave this here for now - this is a complex problem to solve

windows <- tsv_100kb$alignlength
windows <- sort(windows, decreasing = TRUE)
windows[2744:2746] # 24629 24603 24564
length(which(windows > 24600)) # 2745
# So if the threshold was 24,600 then we would end up with 2745 gene trees
length(which(windows > 20000)) # 3118
# So if the threshold was 20,000 then we would end up with 3118 gene trees


csv <- read.csv("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/DataDryad_data/Pease_etal_TomatoPhylo_GeneTrees/Pease_etal_TomatoPhylo_100kbTrees.txt")
