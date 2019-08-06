# Script to check the quality of alignments from empirical data sets in the BenchmarkAlignments database

print("opening packages")
library(ape)

print("initialising namespace")
run_location = "mac"
# run_location = "soma"

if (run_location == "mac"){
  BA_dir <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/run1_BA_AlQualTest/"
  output_dir <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/run1_BA_AlQualTest_results/"
  treelikeness_dir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # where the code for treelikeness statistics and processing is
  empirical_treelikeness_dir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/" # where the code for empirical data and alignment scoring is
  exec_folder <- "/Users/caitlincherryh/Documents/Honours/Executables/"
  # Create a vector with all of the executable file paths
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("3seq","iqtree","Phi","SimBac","SplitsTree.app/Contents/MacOS/JavaApplicationStub","Aliscore_v.2.0/Aliscore.02.2.pl","Aliscore_v.2.0/Aliscore_module.pm")
  exec_paths <- paste0(exec_folder,exec_paths)
  names(exec_paths) <- c("3seq","IQTree","Phi","SimBac","SplitsTree", "ALISCORE","ALISCORE_module")
  source(paste0(treelikeness_dir,"code/func_BA.R"))
} else if (run_location=="soma"){
  BA_dir <- "/data/caitlin/treelikeness/BenchmarkAlignments_DataSubSet/"
  output_dir <- "/data/caitlin/treelikeness/BenchmarkAlignments_DataSubSet_Results/"
  treelikeness_dir <- "/data/caitlin/treelikeness/" # where the code for treelikeness statistics and processing is
  empirical_treelikeness_dir <- "/data/caitlin/empirical_treelikeness/" # where the code for empirical data and alignment scoring is
  # Create a vector with all of the executable file paths
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("/data/caitlin/linux_executables/3seq/3seq","/data/caitlin/linux_executables/iqtree/bin/iqtree","/data/caitlin/linux_executables/PhiPack/Phi",
                  "/data/caitlin/linux_executables/SimBac/SimBac","/data/caitlin/splitstree4/SplitsTree", "/data/caitlin/linux_executables/Aliscore_v.2.0/Aliscore.02.2.pl",
                  "/data/caitlin/linux_executables/Aliscore_v.2.0/Aliscore_module.pm")
  names(exec_paths) <- c("3seq","IQTree","Phi","SimBac","SplitsTree","ALISCORE","ALISCORE_module")
  source(paste0(treelikeness_dir,"code/func_BA_parallel.R")) # run code parallel
} 

# Source files for functions
source(paste0(treelikeness_dir,"code/func_process_data.R"))
source(paste0(treelikeness_dir,"code/func_parametric_bootstrap.R"))
source(paste0(empirical_treelikeness_dir,"code/func_ALISCORE.R"))

# Extract the file names of the alignments
# Accepted values for order are none (as ordered in directory), from smallest to largest number of taxa ("ntaxa") or numer of partitions ("npartitions")
# If order = "user-specified", the user can select which order the datasets will be run in - providing a subset of the dataset names will run a subset of the alignments
print("Extract alignment file paths")
#als <- extract.BA.files(dir = BA_dir, order_by = "ntaxa")
test_order <- c('Richart_2015','Smith_2014','Crawford_2012','Leache_2015','Meiklejohn_2016','Faircloth_2013',
                'McCormack_2013','Wood_2012','Bergsten_2013','Ran_2018_dna','Brown_2012','Cognato_2001',
                'Dornburg_2012','Prebus_2017','Sauquet_2011','Broughton_2013','Siler_2013','Devitt_2013',
                'Kawahara_2013','Cannon_2016_dna','Lartillot_2012','Oaks_2011','Worobey_2014c',
                'Rightmyer_2013','Seago_2011','Moyle_2016','Fong_2012','Unmack_2013','Sharanowski_2011',
                'Anderson_2013','Worobey_2014a','Day_2013','Branstetter_2017','Wainwright_2012','Horn_2014',
                'Tolley_2013','Reddy_2017','Murray_2013','Near_2013','Looney_2016','Pyron_2011') 
                # This is ntaxa order, excluding Worobey because some of the taxa there have names with "*" in and that will crash ALISCORE
als <- extract.BA.files(dir = BA_dir, order_by = "user-specified", user_ordered_list = test_order)
missing_als <- unlist(lapply(als,ALISCORE.output.exists)) # extract only the alignments that _haven't_ been run <- allows you to test why these don't work!

# Run aliscore on each alignment in the als
# Sample aliscore function call: aliscore(al, gaps = "5char", w = 6, aliscore_path = exec_paths[["ALISCORE"]], quality_threshold = 0.5, redo = FALSE)
# Putting in a blank r value means r isn't missing and so it will test 4*N pairs, and will not test the nodes on the tree
print("Testing alignment quality")
lapply(als, aliscore, gaps = "5char", w = 6, r = " ", aliscore_paths = c(exec_paths[["ALISCORE"]],exec_paths[["ALISCORE_module"]]), quality_threshold = 0.5, redo = FALSE)

# Collate alignment quality scoring metrics for each locus into one big old csv file
results_file <- paste0(output_dir,basename(BA_dir),"_collatedAlignmentQuality.csv")
df <- collate.bootstraps(directory = BA_dir, file.name = "locusAlignmentQuality", id = "", output.file.name = results_file) 
