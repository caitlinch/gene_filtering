### empirical_treelikeness/code/1_TestStatistics_EmpiricalData.R
## R program to apply treelikeness statistics to transcriptomes from empirical data
## Specifically designed to work with Rob Lanfear's "BenchmarkAlignments"
## BenchmarkAlignments and metadata have CC0 or CCBY licenses and are available here: https://github.com/roblanf/BenchmarkAlignments
## A number of additional software packages are required, specifically:
##     - IQTREE (Nguyen et al 2015) (http://www.iqtree.org/)
##     - 3SEQ (Lam et al 2018) (http://mol.ax/software/3seq/)
##     - PHI test (Bruen et al 2006) (https://www.maths.otago.ac.nz/~dbryant/software.html)
##     - Splitstree (Huson and Bryant 2006) (http://www.splitstree.org/)
# Caitlin Cherryh 2019



##### Step 1: Open packages #####
print("opening packages")
library(ape) # analyses of phylogenetics and evolution
library(parallel) # support for parallel computation
library(phangorn) # phylogenetic reconstruction and analysis
library(phytools) # tools for comparative biology and phylogenetics
library(seqinr) # data analysis and visualisation for biological sequence data
library(stringr) # wrappers for string operations
library(TreeSim) # simulating phylogenetic trees



##### Step 2: Set file paths and run variables #####
print("initialising namespace")
# input_dir    <- the folder containing the empirical data
#              <- where simulated alignments and output from analysis (e.g. IQ-Tree output files, 3seq output files) will be placed
# output_dir   <- for collated output and results
# treedir      <- "treelikeness" repository location (github.com/caitlinch/treelikeness)
# maindir      <- "empirical_treelikeness" repository location (github.com/caitlinch/empirical_treelikeness)
# cores_to_use <- the number of cores to use for parametric bootstrap. 1 for a single core (wholly sequential), or higher if using parallelisation.
# reps_to_do   <- the number of of bootstrap replicates to perform
# exec_paths   <- location to each executable within the folder. Attach the names of the executables so the paths can be accessed by name

# The SplitsTree executable path can be tricky to find: 
#       - in MacOS, the path is "SplitsTree.app/Contents/MacOS/JavaApplicationStub" (assuming you are in the same directory as the application)
#       - in Linux, after installing and navigating into the folder it's simply "SplitsTree"

# input_dir <- ""
# output_dir <- ""
# treedir <- ""
# maindir <- ""
# cores_to_use <- 1
# reps_to_do <- 199
# sCF_replicates <- 1000

# Create a vector with all of the executable file paths
# To access a path: exec_paths[["name"]]
# exec_folder <- "/path/to/executables/folder/"
# exec_paths <- c("3seq_executable","IQ-Tree_executable","PHIpack_executable","SplitsTree_executable")
# exec_paths <- paste0(exec_folder,exec_paths)
# names(exec_paths) <- c("3seq","IQTree","Phi","SplitsTree")

###

# run_location = "mac"
run_location = "soma"

if (run_location == "mac"){
  input_dir <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/empiricalData/Wu_2018_dna_loci/"
  output_dir <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/output/output_20190806/"
  treedir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # where the treelikeness code is
  maindir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/" # where the empirical treelikeness code is
  
  exec_folder <- "/Users/caitlincherryh/Documents/Honours/Executables/"
  # Create a vector with all of the executable file paths
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("3seq","/Users/caitlincherryh/Documents/Honours/Executables/iqtree-2.0-rc1-MacOSX/bin/iqtree",
                  "SplitsTree.app/Contents/MacOS/JavaApplicationStub")
  exec_paths <- paste0(exec_folder,exec_paths)
  names(exec_paths) <- c("3seq","IQTree","SplitsTree")
  
  # set number of cores and reps for bootstraps
  cores_to_use = 1
  reps_to_do = 9
  sCF_replicates = 1000
  
} else if (run_location=="soma"){
  input_dir <- "/data/caitlin/empirical_treelikeness/Wu_2018_dnaLoci_Primates/"
  output_dir <- "/data/caitlin/empirical_treelikeness/Wu_2018_dnaLoci_Primates_results/"
  treedir <- "/data/caitlin/treelikeness/" # where the treelikeness code is
  maindir <- "/data/caitlin/empirical_treelikeness/" # where the empirical treelikeness code is
  
  # Create a vector with all of the executable file paths
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("/data/caitlin/linux_executables/3seq/3seq","/data/caitlin/linux_executables/iqtree-2.0-rc1-Linux/bin/iqtree","/data/caitlin/splitstree4/SplitsTree")
  names(exec_paths) <- c("3seq","IQTree","SplitsTree")
  
  # set number of cores and reps for bootstraps
  cores_to_use = 30
  reps_to_do= 199
  sCF_replicates = 1000
}



##### Step 3: Source files for functions #####
# Open functions using filepaths from Step 2
source(paste0(treedir,"code/func_test_statistic.R"))
source(paste0(treedir,"code/func_process_data.R"))
source(paste0(treedir,"code/func_parametric_bootstrap.R"))
source(paste0(maindir,"code/func_empirical.R"))



##### Step 4: Extract desired loci #####
# Extract the loci from the folder of interest
full_loci <- list.files(input_dir)
full_loci <- paste0(input_dir, full_loci)
# Remove the three alignment files from the folder
full_loci <- full_loci[grep("alignment",full_loci,invert=TRUE)]
full_loci <- full_loci[grep("README",full_loci,invert=TRUE)]
# Exclude unwanted files (from previous output) from analysis
full_loci <- full_loci[grep(".nex",full_loci)]
full_loci <- full_loci[grep(".nex.",full_loci,invert=TRUE)]
# Remove first, second and third codon positions from the folder (run only whole transcriptomes)
loci <- full_loci[grep("1st",full_loci,invert=TRUE)]
loci <- loci[grep("2nd",loci,invert=TRUE)]
loci <- loci[grep("3rd",loci,invert=TRUE)]

# Create path for whole alignment file
whole_alignment <- paste0(input_dir, "alignment.nex")



##### Step 5: Trim unwanted species (optional) #####
# Initialise list of species of interest 
mammals <- c("CALLI_JAC", "MACAC_FAS", "MACAC_MUL", "PAPIO_ANU", "CHLOR_SAB", "DAUBE_MAD", "GORIL_GOR", "HOMO_SAP", "PAN_PAN", "PAN_TRO", "PONGO_ABE",
             "MICRO_MUR", "NOMAS_LEU", "OTOLE_GAR", "SAIMI_BOL", "TARSI_SYR")
mammals_commonNames <- c("Marmoset", "Crab-eating macaque", "Rhesus macaque", "Olive Baboon", "Green monkey", "Aye-aye", "Gorilla"," Human", "Bonobo",
                         "Chimpanzee", "Orangutan", "Mouse Lemur", "Gibbon", "Galago (Bushbaby)"," Squirrel monkey", "Tarsier")
mammals_speciesNames <- c("Callithrix jacchus", "Macaca fascicularis", "Macaca mulatta", "Papio anubis", "Chlorocebus sabaeus", "Daubentonia madagascariensis",
                          "Gorilla gorilla", "Homo sapiens", "Pan paniscus", "Pan troglodytes", "Pongo abelii", "Microcebus murinus", "Nomascus leucogenys",
                          "Otolemur garnettiiz", "Saimiri boliviensis", "Tarsius syrichta"  )
# Initialise list of outgroups in case rooting is needed later on
outgroups <- c("GALLU_GAL", "MELEA_GAL", "ANOLI_CAR", "PELOD_SIN", "XENOP_TRO", "LATIM_CHA", "DANIO_RER", "GASTE_ACU")
outgroups_commonName <- c("Chicken", "Turkey", "Anole lizard", "Chinese softshell turtle", "Frog", "Coelacanth", "Zebrafish", "Stickleback fish")
outgroups_speciesName <- c("Gallus gallus", "Meleagris gallopavo", "Anolis carolinensis", "Pelodiscus sinensis", "Xenopus tropicalis", "Latimeria chalumnae",
                           "Danio rerio", "Gasterosteus aculeatus")

# Call the trim function to remove all the species you want to remove from the alignment (unless you want to keep ALL the species in the alignment): 
print("trimming alignments")
#lapply(loci, cutSpecies, keep = mammals, output_path_provided = FALSE) # trim all the unwanted species from the collected loci
#cutSpecies(alignment_path = whole_alignment, keep = mammals, output_path_provided = FALSE) # trim all the unwanted species from the whole alignment file



##### Step 6: Calculate the test statistics and run the parametric bootstraps  #####
print("starting analysis")
print("apply treelikeness test statistics")
# To run locally for one alignment: empirical.bootstraps.wrapper(empirical_alignment_path = empirical_alignment_path, program_paths = program_paths,
#                                                        number_of_replicates = 9, iqtree.num_threads = AUTO, iqtree.num_quartets = 100,
#                                                        num_of_cores = 1)
# Parameter choice explanation:
#       ~ iqtree.num_thread = 1: allows parallelisation higher up in the workflow 
#                              - i.e. program will run parametric bootstrap for test statistics with multiple threads
#                              - (number of simulatenous bootstraps set by choice of cores_to_use value)
#       ~ iqtree.num_quartets = 1000: greater than 100 quartets necessary for stable sCF values
# lapply(loci,empirical.bootstraps.wrapper, program_paths = exec_paths, number_of_replicates = reps_to_do, iqtree.num_threads = 1, iqtree.num_quartets = sCF_replicates, num_of_cores = cores_to_use) 



##### Step 7: Collate test statistic results #####
print("collate results")
# Collate all the dataframes together
results_file <- paste0(output_dir,basename(input_dir),"_testStatisticResults.csv")
results_df <- collate.bootstraps(directory = input_dir, file.name = "pValues", id = "", output.file.name = results_file)



##### Step 8: Extract additional information about the alignments and trees #####
# Extracting more information for statistical tests and plots
# Extract total tree depth from each alignment's iqtree file and add to dataframe
results_df$total_tree_depth <- unlist(lapply(results_df$alignment_file, extract.total.tree.length))
# Extract the number of parsimony informative sites from each alignment's iqtree file and add to dataframe
results_df$num_parsimony_informative_sites <- unlist(lapply(results_df$alignment_file, extract.num.informative.sites))
# Extract the GC content for each species in the alignment and add summary statistics to dataframe
GC_vector <- unlist(lapply(results_df$alignment_file,calculate.GC.content))
results_df$GC_content_mean <- GC_vector[c(TRUE,FALSE,FALSE)]
results_df$GC_content_variance <- GC_vector[c(FALSE,TRUE,FALSE)]
results_df$GC_content_sd <-  GC_vector[c(FALSE,FALSE,TRUE)]
# Extract the gene trees (.treefile for each locus) from each alignment's IQ-Tree run and add them to the dataframe
results_df$newick_tree <- unlist(lapply(results_df$alignment_file, extract.treefile.tree))

# Save the newly extended dataframe
results_file <- paste0(output_dir,basename(input_dir),"_completeResults.csv")
write.csv(results_df, file = results_file, row.names = FALSE)


