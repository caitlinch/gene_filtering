# Script to check the quality of alignments from empirical data sets in the BenchmarkAlignments database

print("opening packages")
library(ape)

print("initialising namespace")
run_location = "mac"
# run_location = "soma"

if (run_location == "mac"){
  BA_dir <- "/Users/caitlincherryh/Documents/Repositories/BenchmarkAlignments_DataSubSet/"
  output_dir <- "/Users/caitlincherryh/Documents/Repositories/BenchmarkAlignments_DataSubSet/"
  treelikeness_dir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # where the code for treelikeness statistics and processing is
  empirical_treelikeness_dir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/" # where the code for empirical data and alignment scoring is
  exec_folder <- "/Users/caitlincherryh/Documents/Honours/Executables/"
  # Create a vector with all of the executable file paths
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("3seq","iqtree","Phi","SimBac","SplitsTree.app/Contents/MacOS/JavaApplicationStub","Aliscore_v.2.0/Aliscore.02.2.pl")
  exec_paths <- paste0(exec_folder,exec_paths)
  names(exec_paths) <- c("3seq","IQTree","Phi","SimBac","SplitsTree", "ALISCORE")
  source(paste0(treelikeness_dir,"code/func_BA.R"))
} else if (run_location=="soma"){
  BA_dir <- "/data/caitlin/treelikeness/BenchmarkAlignments_DataSubSet/"
  output_dir <- "/data/caitlin/treelikeness/BenchmarkAlignments_DataSubSet_Results/"
  treelikeness_dir <- "/data/caitlin/treelikeness/" # where the code for treelikeness statistics and processing is
  empirical_treelikeness_dir <- "/data/caitlin/empirical_treelikeness/" # where the code for empirical data and alignment scoring is
  # Create a vector with all of the executable file paths
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("/data/caitlin/linux_executables/3seq/3seq","/data/caitlin/linux_executables/iqtree/bin/iqtree","/data/caitlin/linux_executables/PhiPack/Phi",
                  "/data/caitlin/linux_executables/SimBac/SimBac","/data/caitlin/splitstree4/SplitsTree", "/data/caitlin/linux_executables/Aliscore_v.2.0/Aliscore.02.2.pl")
  names(exec_paths) <- c("3seq","IQTree","Phi","SimBac","SplitsTree","ALISCORE")
  source(paste0(treelikeness_dir,"code/func_BA_parallel.R")) # run code parallel
} 

# Source files for functions
source(paste0(treelikeness_dir,"code/func_process_data.R"))
source(paste0(treelikeness_dir,"code/func_parametric_bootstrap.R"))

# Extract the file names of the alignments
# Accepted values for order are none (as ordered in directory), from smallest to largest number of taxa ("ntaxa") or numer of partitions ("npartitions")
# If order = "user-specified", the user can select which order the datasets will be run in - providing a subset of the dataset names will run a subset of the alignments
print("Extract alignment file paths")
#als <- extract.BA.files(dir = BA_dir, order_by = "ntaxa")
als <- extract.BA.files(dir = "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/test_01_aliscoretest/", order_by = "ntaxa")

print("Testing alignment quality")
# Sample aliscore function call: aliscore(al, gaps = "5char", w = 6, aliscore_path = exec_paths[["ALISCORE"]], quality_threshold = 0.5, redo = FALSE)
# Run aliscore on each alignment in the als
lapply(als, aliscore, gaps = "5char", w = 6, aliscore_path = exec_paths[["ALISCORE"]], quality_threshold = 0.5, redo = FALSE)

alipath <- "/Users/caitlincherryh/Documents/Honours/Executables/Aliscore_v.2.0/Aliscore.02.2.pl"
#alipath <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/test_01_aliscoretest/Anderson_2013/Aliscore.02.2.pl"
p <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/test_01_aliscoretest/Anderson_2013/COI.nex"
op <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/test_01_aliscoretest/Anderson_2013/16S.nex.fasta_List_random.txt"
n <- read.nexus.data(p)
d <- as.DNAbin(n)
