# Script to apply test statistics and parametric bootstrap to empirical data sets in the BenchmarkAlignments database

print("opening packages")
library(ape)
library(parallel)
library(phangorn)
library(phytools)
library(seqinr)
library(stringr)
library(TreeSim)

print("initialising namespace")

run_location = "mac"
# run_location = "soma"


if (run_location == "mac"){
  input_dir <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/empiricalData/Wu_2018_dna_loci/"
  output_dir <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/output/output_20190806/"
  treedir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # where the code is
  maindir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/"
  exec_folder <- "/Users/caitlincherryh/Documents/Honours/Executables/"
  # Create a vector with all of the executable file paths
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("3seq","iqtree","Phi","SimBac","SplitsTree.app/Contents/MacOS/JavaApplicationStub")
  exec_paths <- paste0(exec_folder,exec_paths)
  names(exec_paths) <- c("3seq","IQTree","Phi","SimBac","SplitsTree")
  source(paste0(treedir,"code/func_BA.R"))
} else if (run_location=="soma"){
  input_dir <- "/data/caitlin/treelikeness/BenchmarkAlignments_DataSubSet/"
  output_dir <- "/data/caitlin/treelikeness/BenchmarkAlignments_DataSubSet_Results/"
  treedir <- "/data/caitlin/treelikeness/" # where the code is
  # Create a vector with all of the executable file paths
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("/data/caitlin/linux_executables/3seq/3seq","/data/caitlin/linux_executables/iqtree/bin/iqtree","/data/caitlin/linux_executables/PhiPack/Phi",
                  "/data/caitlin/linux_executables/SimBac/SimBac","/data/caitlin/splitstree4/SplitsTree")
  names(exec_paths) <- c("3seq","IQTree","Phi","SimBac","SplitsTree")
  source(paste0(treedir,"code/func_BA_parallel.R")) # run code parallel
}

# Source files for functions
source(paste0(treedir,"code/func_test_statistic.R"))
source(paste0(treedir,"code/func_process_data.R"))
source(paste0(treedir,"code/func_parametric_bootstrap.R"))

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

# Extract the loci from the folder
loci <- list.files(input_dir)




# Calculate the test statistics and run the bootstraps
print("starting analysis")
# To run for one alignment: empirical.bootstraps.wrapper(empirical_alignment_path = empirical_alignment_path, program_paths = program_paths, number_of_replicates = 9)
if (run_location=="soma"){
  #lapply(als,empirical.bootstraps.wrapper, program_paths = exec_paths, number_of_replicates = 99) 
} else if (run_location=="mac"){
  #lapply(als,empirical.bootstraps.wrapper, program_paths = exec_paths, number_of_replicates = 9) 
}

# Collate all the results
#results_file <- paste0(output_dir,basename(BA_dir),"_completeResults.csv")
#df <- collate.bootstraps(directory = BA_dir, file.name = "pValues", id = "", output.file.name = results_file)

