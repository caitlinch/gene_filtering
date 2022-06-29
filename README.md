# Gene Filtering
#### The impact of filtering loci with evidence of recombination on species tree topology for four empirical datasets

Caitlin Cherryh

Friday May 13, 2022

***
### Summary
This github repository contains scripts used to investigate the impact of filtering loci with evidence of recombination on species tree topology. 

An underlying assumption in phylogenetics is that each site in a loci shares an identical evolutionary history that fits a single bifurcating tree. However, this assumption is broken by biological processes such as introgression or recombination. We selected four empirical datasets and investigated whether removing loci identified as putatively recombinant impacted species tree topology. To do so, we selected three tests for recombination detection (PHI, MaxChi, and GeneConv). We applied each test to each loci in each dataset . Then we used the results to break the loci into subsets. For each test, the set of loci was broken into a subset of loci that passed that test and a subset of loci that failed each test (i.e. loci that were identified as putatively recombinant). We then estimated species trees from each subset with both summary coalescent (ASTRAL-III) and maximum likelihood (IQ-Tree2) tree estimation methods. Finally, we compared the goodness of fit and topology of each tree.

If you replicate any part of these analyses or use functions from these scripts, please cite this repository. Thank you!

#### Contents
+ Scripts
    + All scripts necessary to completely replicate this analysis are included in the `code/` folder
    + Each script contains a detailed description and a list of software necessary to run that script
+ Species and subset trees
    + Species trees are included in the `species_trees/` folder. 
    + Each tree is labelled by the dataset, the test for recombination estimation, and the tree estimation method (e.g. `Tomatoes_PHI_pass_ASTRAL_species.tre` means the tree is estimated from the Tomato dataset loci that passed the PHI test using the summary coalescent method ASTRAL-III).
+ Comparison trees
    + Trees used to investigate the Primates dataset included in the `primate_tree_topologies/` folder
    + These trees are used with script `code/4_DataAnalysis.R` to perform an AU test for each loci to identify the best topology (out of a set of three) around two branches within the Primates tree
+ Tomato gene trees
    + Text file containing all gene trees for the Tomatoes dataset included in the `tomato_cloudogram/` folder
    + Used with the script `code/5_Plots.R` to generate a cloudogram
+ Output
    + `.csv` files containing output from each stage of the analysis can be found in the `output/` folder
+ Instructions for replication
    + Instructions for replicating these analyses, along with details about the datasets and software used, are in this `README.md` file.

***
### Instructions to reproduce the analyses:
To fully replicate the analyses, follow these steps:

1. Download the datasets and software specified below
2. Create the conda environment `gene_filtering` using the `environment.yaml` file
3. Prepare the alignments for analysis
    + The Primates and Plants datasets do not require any preparation
    + To prepare the Tomatoes dataset, run the script `code/0_Pease2016_data_formatting.R`. This script replicates the 2745 100kbp windows used for phylogenetic analysis in Pease et. al. (2016)
    + To prepare the Metazoans dataset, run the script `code/0_Whelan2017_data_formatting.R`. This script separates the supermatrix into alignments of individual loci. 
4. Estimate gene trees and apply tests for recombination detection
    + Run the script `code/1_RecombinationDetection.R` to apply three tests for recombination detection (PHI, MaxChi, and GENECONV) to each locus and estimate a gene tree for each locus
5. Estimate species trees
    + Based on the results of the tests for recombination detection, the loci for each dataset will be subsetted into categories that pass and fail each test. For each subset we estimate a summary coalescent tree in ASTRAL-III and a maximum likelihood tree in IQ-Tree2. We also estimate a tree from the unfiltered datasets (i.e. a species tree using all loci) with both ASTRAL-III and IQ-Tree2
    + For each test, the tree estimated from the loci that passed the test will be compared to the tree estimated from loci that failed the test and the tree estimated from the whole unfiltered dataset
        + For trees estimated in ASTRAL-III, the absolute goodness of fit will be compared using the QuartetNetworkGoodnessFit.jl package
        + For trees estimated in IQ-Tree2, the absolute goodness of fit will be compared using the AU test
        + We will compare the tree topology using the Robinson-Foulds and weighted Robinson-Foulds distances
    + To estimate and compare trees, run the script `code/3_Species_Tree_Comparison.R`
    + If you are using the Plants dataset, estimating a partitioned maximum likelihood tree where the model for each loci is selected by ModelFinder in IQ-Tree becomes intractable when free-rate models are included. In this case, you can run script `code/2.5_ExtractingModels_DeepDatasets.R` to select the best model for each loci that is not a free-rate model
    + If you are using the Plants dataset, estimating a maximum likelihood tree in IQ-Tree is computationally intractable. To manage this, we estimated maximum likelihood trees for the Plants dataset in RAxML-ng with free-rate models excluded and no bootstraps
6. Data analysis and plotting
    + The script `code/4_DataAnalysis.R` performs data analysis and plots, including:
        + Identifying branches that are present in one tree and not the other (and vice versa)
        + Plotting the distribution of branch lengths
        + Plotting the distribution of support values (either Posterior Probability for trees estimated in ASTRAL-III or UltraFast Bootstraps for trees estimated in IQ-Tree)

***
### Empirical datasets
For these analyses I used four empirical phylogenetic datasets:

- **Primates**
  - **Paper**: Vanderpool D., Minh B.Q., Lanfear R., Hughes D., Murali S. et al. 2020. Primate phylogenomics uncovers multiple rapid radiations and ancient interspecific introgression. *PLOS Biology* 18(12):e3000954. https://doi.org/10.1371/journal.pbio.3000954
  - **Data**: Vanderpool, Dan et al. 2020. Supplementary data for: Primate phylogenomics uncovers multiple rapid radiations and ancient interspecific introgression. Dryad. Dataset. https://doi.org/10.5061/dryad.rfj6q577d
  - **Data analysed**: 1730_Alignments_FINAL.tar.gz (Used for Figure 1 in Vanderpool et. al. 2020)
- **Tomatoes**
  - **Paper**: Pease, J.B., Haak, D.C., Hahn, M. W., Moyle, L. 2016. Phylogenomics reveals three sources of adaptive variation during a rapid radiation, *PLOS Biology*, 14(2):e1002379. https://doi.org/10.1371/journal.pbio.1002379
  - **Data**: Pease, James B., Haak, D.C., Hahn, M.W., Moyle, L.C. 2016. Data from: Phylogenomics reveals three sources of adaptive variation during a rapid radiation. Dryad. Dataset. https://doi.org/10.5061/dryad.182dv
  - **Data analysed**: Pease_etal_Tomato29acc_HQ.mvf.gz (Used for Figure 2 in Pease et. al. (2016))
- **Metazoans**
  - **Paper**: Whelan, N.V., Kocot, K.M., Moroz, T.P. et al. 2017. Ctenophore relationships and their placement as the sister group to all other animals. *Nature Ecology and Evolution* 1:1737–1746 https://doi.org/10.1038/s41559-017-0331-3
  - **Data**: Whelan, N.V., Kocot, K.M., Moroz, T.P., Mukherjee, K. et al. 2017. Ctenophora Phylogeny Datasets and Core Orthologs. Figshare. Dataset. https://doi.org/10.6084/m9.figshare.4484138.v1 
  - **Data analysed**: Metazoa_Choano_RCFV_strict.phy (Used for Figure 2 in Whelan et. al. (2017))
- **Plants**
  - **Paper**: Leebens-Mack, J.H., Barker, M.S., Carpenter, E.J., Deyholos, M.K. 2019. One thousand plant transcriptomes and the phylogenomics of green plants. *Nature* 574:679–685. https://doi.org/10.1038/s41586-019-1693-2
  - **Data**: Siavash, M., & Sayyari, E. 2019. smirarab/1kp: zenodo (zenodo). Zenodo. https://doi.org/10.5281/zenodo.3255100
  - **Data analysed**: alignments-FAA-masked.tar.bz (Used for Figure 2 of Leebens-Mack et. al. (2019))

***
### Software

- **mvftools** 
  - Pease, J.B. and Rosenzweig, B.K. 2018. Encoding data using biological principles: the Multisample Variant Format for phylogenomics and population genomics. *IEEE/ACM Transactions on Computational Biology and Bioinformatics*. 15(4):1231-1238. http://www.dx.doi.org/10.1109/tcbb.2015.2509997
  - https://github.com/peaselab/mvftools
  - Used to process the Tomatoes alignment from MVF to fasta file format.
- **PHIPack**
  - Bruen, T.C, Philippe, H., Bryant, D. 2006. A simple and robust statistical test for detecting the presence of recombination. *Genetics*. 172(4):2665–2681. https://doi.org/10.1534/genetics.105.048975
  - https://www.maths.otago.ac.nz/~dbryant/software.html
  - Used both PHI and MaxChi tests for recombination detection
- **GeneConv**
  - Sawyer, S.A. 1989. Statistical tests for detecting gene conversion. *Molecular Biology and Evolution* 6:526-538.
  - https://www.math.wustl.edu/~sawyer/geneconv/
  - Used for recombination detection
- **IQTREE2**
  - Minh, B.Q., Schmidt, H.A., Chernomor, O., Schrempf, D., Woodhams, M.D., von Haeseler, A., Lanfear, R. 2020 IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic rra, *Molecular Biology and Evolution*, 37(5):1530–1534, https://doi.org/10.1093/molbev/msaa015
  - http://www.iqtree.org/
  - Used to estimate gene trees for each loci, to estimate maximum likelihood species trees, and to apply the AU test to trees estimated with maximum likelihood methods
- **ASTRAL-III**
  - Zhang, C., Rabiee, M., Sayyari, E., Mirarab, S. 2018. ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees. *BMC Bioinformatics* 19:153. https://doi.org/10.1186/s12859-018-2129-y
  - https://github.com/smirarab/ASTRAL
  - Used to estimate summary coalescent species tree
- **RAxML-ng**
  - Kozlov, A.M., Darriba, D., Flouri, T., Morel, B., Stamatakis, A., 2019. RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. *Bioinformatics*. 35(21):4453–4455. https://doi.org/10.1093/bioinformatics/btz305
  - https://github.com/amkozlov/raxml-ng
  - Used to estimate maximum likelihood trees for the Plants dataset (was computationally intractable in IQ-Tree2)
- **Julia** 
  - Bezanson, J., Edelman, A., Karpinski, S. and Shah, V.B. 2017. Julia: a fresh approach to numerical computing. *SIAM Review*. 59(1):65-98. https://epubs.siam.org/doi/10.1137/141000671
  - https://julialang.org/
  - Necessary for the QuartetNetworkGoodnessFit.jl Julia package
- **PhyloNetworks.jl** 
  - Solís-Lemus, C., Bastide, P. and Ané, C. 2017. PhyloNetworks: a package for phylogenetic networks. *Molecular Biology and Evolution* 34(12):3292–3298. https://doi.org/10.1093/molbev/msx235
  - https://github.com/crsl4/PhyloNetworks.jl
  - Necessary for the QuartetNetworkGoodnessFit.jl Julia package
- **QuartetNetworkGoodnessFit.jl** 
  - Cai, R., Ané, C. 2021. Assessing the fit of the multi-species network coalescent to multi-locus data. *Bioinformatics*. 37(5):634–641. https://doi.org/10.1093/bioinformatics/btaa863
  - https://github.com/cecileane/QuartetNetworkGoodnessFit.jl
  - Used to apply the QuartetNetworkoodnessFit test to trees estimated under the coalescent. 

***
### Citation information
If you use these scripts, please cite this github repository:

Cherryh, C. 2022. Gene Filtering. [Electronic resource: R code]. Available at: https://github.com/caitlinch/gene_filtering (Accessed: 13 May 2022)



