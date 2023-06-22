### gene_filtering/code/func_plots.R
## R functions to facilitate nice plotting
# Caitlin Cherryh 2023

library(phytools) # Functions: nodeHeights
library(dplyr)
library(glue)
library(ggtree)



get.DNA.site.info <- function(loci_name, alignment_dir){
  # Open an alignment and return the number of phylogenetically informative sites and number of segregating sites in that alignment
  all_als <- list.files(alignment_dir)
  aln <- paste0(alignment_dir, grep(loci_name, all_als, value = TRUE))
  f <- read.dna(aln, format = "fasta")
  ss <- seg.sites(f)
  n_ss <- length(ss)
  n_pis <- ips::pis(f, what = "absolute")
  n_sites <- dim(f)[2]
  op <- c(n_ss, n_pis, n_sites)
  names(op) <- c("n_seg_sites", "n_pis", "n_sites")
  return(op)
}



rescale.tree.length <- function(tree, scaled_length){
  # Rescale the total length (i.e. height) of a tree
  
  tree$edge.length <- tree$edge.length / max(nodeHeights(tree)[,2]) * scaled_length
  return(tree)
}



rescale.multiphylo <- function(trees, scaled_length){
  # Wrapper for feeding list of trees into rescale.tree.length using lapply
  
  for (i in 1:length(trees)){
    trees[[i]] <- rescale.tree.length(trees[[i]], scaled_length)
  }
  return(trees)
}



color.code.primate.clades <- function(p_tree, color = TRUE, concatenated = TRUE){
  # Quick function to colour code clades in primates dataset based on tree estimation method
  
  # Set which clade differs and which taxa remain identical
  if (color == TRUE){
    if (concatenated == TRUE){
      variable_species <- c("Cercocebus atys","Mandrillus leucophaeus","Papio anubis","Theropithecus gelada")
      congruent_species <- c("Aotus nancymaae","Saimiri boliviensis","Cebus capucinus imitator","Callithrix jacchus",
                             "Carlito syrichta","Otolemur garnettii","Microcebus murinus","Propithecus coquereli",
                             "Galeopterus variegatus","Tupaia chinensis","Mus musculus","Nomascus leucogenys",
                             "Pongo abelii","Gorilla gorilla","Homo sapiens","Pan paniscus","Pan troglodytes",
                             "Chlorocebus sabaeus","Macaca nemestrina","Macaca fascicularis","Macaca mulatta",
                             "Colobus angolensis palliatus","Piliocolobus tephrosceles","Rhinopithecus bieti",
                             "Rhinopithecus roxellana")
    } else if (concatenated == FALSE){
      variable_species <- c("Aotus nancymaae","Saimiri boliviensis","Cebus capucinus imitator","Callithrix jacchus")
      congruent_species <- c("Carlito syrichta","Otolemur garnettii","Microcebus murinus","Propithecus coquereli",
                             "Galeopterus variegatus","Tupaia chinensis","Mus musculus","Nomascus leucogenys",
                             "Pongo abelii","Gorilla gorilla","Homo sapiens","Pan paniscus","Pan troglodytes",
                             "Chlorocebus sabaeus","Macaca nemestrina","Macaca fascicularis","Macaca mulatta",
                             "Cercocebus atys","Mandrillus leucophaeus","Papio anubis","Theropithecus gelada",
                             "Colobus angolensis palliatus","Piliocolobus tephrosceles","Rhinopithecus bieti",
                             "Rhinopithecus roxellana")
    }
  } else if (color == FALSE){
    variable_species <- c()
    congruent_species <- c("Aotus nancymaae","Saimiri boliviensis","Cebus capucinus imitator","Callithrix jacchus",
                           "Carlito syrichta","Otolemur garnettii","Microcebus murinus","Propithecus coquereli",
                           "Galeopterus variegatus","Tupaia chinensis","Mus musculus","Nomascus leucogenys",
                           "Pongo abelii","Gorilla gorilla","Homo sapiens","Pan paniscus","Pan troglodytes",
                           "Chlorocebus sabaeus","Macaca nemestrina","Macaca fascicularis","Macaca mulatta",
                           "Cercocebus atys","Mandrillus leucophaeus","Papio anubis","Theropithecus gelada",
                           "Colobus angolensis palliatus","Piliocolobus tephrosceles","Rhinopithecus bieti",
                           "Rhinopithecus roxellana")
  }
  
  # Create dataframe with tip information
  tip_df <- data.frame(taxa = c(variable_species, congruent_species),
                       clade = c(rep("Variable", length(variable_species)), rep("Congruent", length(congruent_species)) ),
                       color = c(rep("gray60", length(variable_species)), rep("black", length(congruent_species)) ) 
  )
  tip_lab_df <- dplyr::mutate(tip_df, 
                              lab = glue('italic("{taxa}")'),
                              name = glue("<i style='color:{color}'>{taxa}</i>") )
  
  # Return the tip label dataframe
  return(tip_lab_df)
}



color.code.comparison.clades <- function(p_tree, variable = "Cebidae"){
  # Quick function to colour code clades in primates dataset for supplementary figures (based on likelihood for each tree topology)
  
  if (variable == "Cebidae"){
    clade_a <- c("Callithrix_jacchus")
    clade_b <- c("Aotus_nancymaae")
    clade_c <- c("Saimiri_boliviensis","Cebus_capucinus_imitator")
    clade_d <- c("Carlito_syrichta", "Galeopterus_variegatus", "Mus_musculus", "Tupaia_chinensis", "Microcebus_murinus", "Propithecus_coquereli",
                 "Otolemur_garnettii")
    clade_e <- c("Cercocebus_atys", "Mandrillus_leucophaeus", "Papio_anubis", "Theropithecus_gelada", "Macaca_fascicularis", "Macaca_mulatta",
                 "Macaca_nemestrina", "Chlorocebus_sabaeus", "Colobus_angolensis_palliatus", "Piliocolobus_tephrosceles", "Rhinopithecus_bieti",
                 "Rhinopithecus_roxellana", "Gorilla_gorilla", "Homo_sapiens", "Pan_paniscus", "Pan_troglodytes", "Pongo_abelii", "Nomascus_leucogenys")
    taxa_names <- c(clade_a, clade_b, clade_c, clade_d, clade_e)
    clean_taxa_names <- gsub("_", " ", taxa_names)
    
    # Create dataframe with tip information
    tip_df <- data.frame(taxa = taxa_names,
                         clean_taxa_names = clean_taxa_names,
                         clade = c(rep("clade_a", length(clade_a)), rep("clade_b", length(clade_b)),
                                   rep("clade_c", length(clade_c)), rep("clade_d", length(clade_d)),
                                   rep("clade_e", length(clade_e))),
                         color = c(rep("Sky blue", length(clade_a)), rep("Orange", length(clade_b)),
                                   rep("Bluish green", length(clade_c)), rep("Black", length(clade_d)),
                                   rep("Gray60", length(clade_e))) ) 
    
  } else if (variable == "Comparison"){
    clade_a <- c("Carlito_syrichta")
    clade_b <- c("Callithrix_jacchus", "Aotus_nancymaae", "Cebus_capucinus_imitator", "Saimiri_boliviensis")
    clade_c <- c("Galeopterus_variegatus", "Mus_musculus", "Tupaia_chinensis", "Microcebus_murinus", "Propithecus_coquereli", "Otolemur_garnettii")
    clade_d <- c("Cercocebus_atys", "Mandrillus_leucophaeus", "Papio_anubis", "Theropithecus_gelada", "Macaca_fascicularis", "Macaca_mulatta",
                 "Macaca_nemestrina", "Chlorocebus_sabaeus", "Colobus_angolensis_palliatus", "Piliocolobus_tephrosceles", "Rhinopithecus_bieti",
                 "Rhinopithecus_roxellana", "Gorilla_gorilla", "Homo_sapiens", "Pan_paniscus", "Pan_troglodytes", "Pongo_abelii", "Nomascus_leucogenys")
    taxa_names <- c(clade_a, clade_b, clade_c, clade_d)
    clean_taxa_names <- gsub("_", " ", taxa_names)
    
    # Create dataframe with tip information
    tip_df <- data.frame(taxa = taxa_names,
                         clean_taxa_names = clean_taxa_names,
                         clade = c(rep("clade_a", length(clade_a)), rep("clade_b", length(clade_b)),
                                   rep("clade_c", length(clade_c)), rep("clade_d", length(clade_d))),
                         color = c(rep("Sky blue", length(clade_a)), rep("Orange", length(clade_b)),
                                   rep("Bluish green", length(clade_c)), rep("Gray60", length(clade_d))) )
  }
  
  tip_lab_df <- dplyr::mutate(tip_df, 
                              lab = glue('italic("{clean_taxa_names}")'),
                              name = glue("<i style='color:{color}'>{clean_taxa_names}</i>") )
  
  # Return the tip label dataframe
  return(tip_lab_df)
}



comparison.clade.tip.labels <- function(clade_tree){
  # Quick function to nicely label clades from deep split supplementary figure
  taxa <- clade_tree$tip.label
  clean_taxa_names <- gsub("_", " ", taxa)
  tip_df <- data.frame(taxa = taxa,
                       clean_taxa_names = clean_taxa_names)
  tip_lab_df <- dplyr::mutate(tip_df, 
                              lab = glue('italic("{clean_taxa_names}")'))
  return(tip_lab_df)
}



rename.tomato.tips <- function(tip_vector){
  # Quick function to take in vector of tomatoes species numbers and return 
  
  new_tip_names <- unlist(lapply(tip_vector, find.one.tomato.tip, include.genus = FALSE))
  return(new_tip_names)
}



find.one.tomato.tip <- function(single_tip, include.genus = FALSE){
  # Quick function to take in a single species number from the tomatoes tree and return full scientific name for that species
  
  tip_list <- list("SL2.50" = c("Solanum lyco. 'Heinz 1706'", "S. lyco. 'Heinz 1706'", "Solanum lycopersicum 'Heinz 1706'", "S. lycopersicum 'Heinz 1706'"),
                   "LA3475" = c("Solanum lycopersicum 3475", "S. lycopersicum 3475"),
                   "LA3124" = c("Solanum cheesmaniae 3124", "S. cheesmaniae 3124"),
                   "LA0429" = c("Solanum galapagense 0429", "S. galapagense 0429"),
                   "LA0436" = c("Solanum cheesmaniae 0436", "S. cheesmaniae 0436"),
                   "LA3909" = c("Solanum galapagense 3909", "S. galapagense 3909"),
                   "LA2933" = c("Solanum lycopersicum 2933", "S. lycopersicum 2933", "Solanum lycopersicum var. cerasiforme 2933", "S. lycopersicum var. cerasiforme 2933"), 
                   "LA1269" = c("Solanum pimpinellifolium 1269", "S. pimpinellifolium 1269"),
                   "LA1589" = c("Solanum pimpinellifolium 1589", "S. pimpinellifolium 1589"),
                   "LA1028" = c("Solanum chmielewskii 1028", "S. chmielewskii 1028"),
                   "LA1316" = c("Solanum chmielewskii 1316", "S. chmielewskii 1316"),
                   "LA2172" = c("Solanum arcanum 2172", "S. arcanum 2172"),
                   "LA1322" = c("Solanum neorickii 1322", "S. neorickii 1322"), 
                   "LA2133" = c("Solanum neorickii 2133", "S. neorickii 2133"),
                   "LA4116" = c("Solanum sitiens 4116", "S. sitiens 4116"),
                   "LA2951" = c("Solanum lycopersicoides 2951", "S. lycopersicoides 2951"),
                   "LA4126" = c("Solanum lycopersicoides 4126", "S. lycopersicoides 4126"),
                   "LA0716" = c("Solanum pennellii 0716", "S. pennellii 0716"),
                   "LA3778" = c("Solanum pennellii 3778", "S. pennellii 3778"),
                   "LA0407" = c("Solanum habrochaites 0407", "S. habrochaites 0407"),
                   "LA1777" = c("Solanum habrochaites 1777", "S. habrochaites 1777"),
                   "LA1364" = c("Solanum huaylasense 1364", "S. huaylasense 1364"),
                   "LA1782" = c("Solanum chilense 1782", "S. chilense 1782"),
                   "LA4117" = c("Solanum chilense 4117A", "S. chilense 4117A"),
                   "LA2744" = c("Solanum peruvianum 2744", "S. peruvianum 2744"),
                   "LA2964" = c("Solanum peruvianum 2964", "S. peruvianum 2964"),
                   "LA0107" = c("Solanum corneliomuelleri 0107", "S. corneliomuelleri 0107"),
                   "LA0444" = c("Solanum corneliomuelleri 0444", "S. corneliomuelleri 0444"), 
                   "LA1358" = c("Solanum huaylasense 1358", "S. huaylasense 1358"),
                   "LA1360" = c("Solanum huaylasense 1360", "S. huaylasense 1360"))
  if (single_tip %in% names(tip_list)){
    # If this taxa is included in the tip list, use the taxa number to look up the taxa species name
    if (include.genus == FALSE){
      new_tip_name <- tip_list[[single_tip]][2]
    } else if (include.genus == TRUE){
      new_tip_name <- tip_list[[single_tip]][1]
    }
  } else {
    # If this taxa is not included in the tip list it may be because you have already renamed it. Return the name as is.
    new_tip_name <- single_tip
  }
  return(new_tip_name)
}



reformat.congruent.tomato.clades <- function(tom_tree){
  # Quick function to remove the sections of the tomato tree that remain identical
  
  # Want to keep all Peruvianum species:
  Peruvianum = c("LA1364", "LA2744", "LA1358", "LA0107", "LA0444", "LA2964", "LA1782", "LA4117")
  # Want to keep one tip from each other clade:
  Esculentum = "SL2.50"
  Arcanum = "LA1322"
  Hirsutum = "LA0407"
  Clade_Outgroup = c("LA4116")
  root_outgroup = c("LA4116", "LA4126", "LA2951")
  # Collate the list of taxa to keep
  taxa_to_keep <- c(Peruvianum, Esculentum, Arcanum, Hirsutum, Clade_Outgroup)
  # Root at outgroup
  tom_tree <- root(tom_tree, root_outgroup)
  # Drop all other taxa
  tom_tree <- keep.tip(tom_tree, taxa_to_keep)
  # Rename the remaining taxa into the clade name
  tom_tree$tip.label[which(tom_tree$tip.label == "SL2.50")] <- "Esculentum clade"
  tom_tree$tip.label[which(tom_tree$tip.label == "LA1322")] <- "Arcanum clade"
  tom_tree$tip.label[which(tom_tree$tip.label == "LA0407")] <- "Hirsutum clade"
  tom_tree$tip.label[which(tom_tree$tip.label == "LA4116")] <- "Outgroup"
  # Return the reformatted tree
  return(tom_tree)
}



color.code.tomato.clades <- function(tom_tree, taxa.numbers = FALSE, trimmed = TRUE, color = TRUE){
  # Quick function to colour code clades in tomatoes dataset
  
  # Which tips are in which clade:
  if (trimmed == FALSE){
    # Include all taxa
    if (taxa.numbers == TRUE){
      # Build data frame from taxa numbers
      Esculentum = c("SL2.50", "LA0436", "LA3909", "LA0429", "LA3124", "LA3475", "LA2933", "LA1269", "LA1589")
      Arcanum = c("LA1322", "LA2133", "LA2172", "LA1028", "LA1316")
      Peruvianum = c("LA1364", "LA2744", "LA1358", "LA0107", "LA0444", "LA2964", "LA1782", "LA4117")
      Hirsutum = c("LA0407", "LA1777", "LA0716", "LA3778")
      Clade_Outgroup = c("LA4116", "LA4126", "LA2951")
    } else if (taxa.numbers == FALSE){
      # Build dataframe from taxa names
      Esculentum = c("S. lyco. 'Heinz 1706'", "S. cheesmaniae 0436", "S. galapagense 3909", "S. galapagense 0429", 
                     "S. cheesmaniae 3124", "S. lycopersicum 3475", "S. lycopersicum 2933", "S. pimpinellifolium 1269",
                     "S. pimpinellifolium 1589")
      Arcanum = c("S. neorickii 1322", "S. neorickii 2133", "S. arcanum 2172", "S. chmielewskii 1028", "S. chmielewskii 1316")
      Peruvianum = c("S. huaylasense 1364", "S. peruvianum 2744", "S. huaylasense 1358", "S. corneliomuelleri 0107", 
                     "S. corneliomuelleri 0444", "S. peruvianum 2964", "S. chilense 1782", "S. chilense 4117A")
      Hirsutum = c("S. habrochaites 0407", "S. habrochaites 1777", "S. pennellii 0716", "S. pennellii 3778")
      Clade_Outgroup = c("S. sitiens 4116", "S. lycopersicoides 4126", "S. lycopersicoides 2951")
    }
  } else if (trimmed == TRUE){
    # Taxa have been trimmed to one species per each group, except for Peruvianum (includes all Peruvianum species)
    if (taxa.numbers == TRUE){
      # Build data frame from taxa numbers
      Esculentum = "Esculentum clade"
      Arcanum = "Arcanum clade"
      Peruvianum = c("LA1364", "LA2744", "LA1358", "LA0107", "LA0444", "LA2964", "LA1782", "LA4117")
      Hirsutum = "Hirsutum clade"
      Clade_Outgroup = "Outgroup"
    } else if (taxa.numbers == FALSE){
      # Build dataframe from taxa names
      Esculentum = "Esculentum clade"
      Arcanum = "Arcanum clade"
      Peruvianum = c("S. huaylasense 1364", "S. peruvianum 2744", "S. huaylasense 1358", "S. corneliomuelleri 0107", 
                     "S. corneliomuelleri 0444", "S. peruvianum 2964", "S. chilense 1782", "S. chilense 4117A")
      Hirsutum = "Hirsutum clade"
      Clade_Outgroup = "Outgroup"
    }
  }
  
  if (color == TRUE){
    # Create dataframe with tip information
    tip_df <- data.frame(taxa = c(Esculentum, Arcanum, Peruvianum, Hirsutum, Clade_Outgroup),
                         clade = c(rep("Esculentum", length(Esculentum)), rep("Arcanum", length(Arcanum)), rep("Peruvianum", length(Peruvianum)), 
                                   rep("Hirsutum", length(Hirsutum)), rep("Outgroup", length(Clade_Outgroup))),
                         color = c(rep("firebrick3", length(Esculentum)), rep("goldenrod3", length(Arcanum)), rep("darkgreen", length(Peruvianum)), 
                                   rep("navy", length(Hirsutum)), rep("black", length(Clade_Outgroup))))
  } else if (color == FALSE){
    # Create dataframe with tip information
    tip_df <- data.frame(taxa = c(Esculentum, Arcanum, Peruvianum, Hirsutum, Clade_Outgroup),
                         clade = c(rep("Esculentum", length(Esculentum)), rep("Arcanum", length(Arcanum)), rep("Peruvianum", length(Peruvianum)), 
                                   rep("Hirsutum", length(Hirsutum)), rep("Outgroup", length(Clade_Outgroup))),
                         color = c("black") )
  }
  
  # Mutate to add italic formatting and style for parsing with expression() function
  tip_lab_df <- dplyr::mutate(tip_df, 
                              lab = glue('italic("{taxa}")'),
                              name = glue("<i style='color:{color}'>{taxa}</i>") ) 
  
  # Return the tip label dataframe
  return(tip_lab_df)
}



reformat.congruent.metazoan.clades <- function(m_tree, trim = "FALSE"){
  # Quick function to remove the sections of the Metazoan tree that are not interesting
  
  if (trim == "FALSE"){
    # Keep all species in all clades
    Bilateria = c("Homo_sapiens", "Strongylocentrotus_purpatus", "Hemithris_psittacea", "Capitella_teleta", "Drosophila_melanogaster",
                  "Daphnia_pulex")
    Cnidaria = c("Hydra_vulgaris", "Bolocera_tuediae", "Aiptasia_pallida", "Hormathia_digitata", "Nematostella_vectensis", "Acropora_digitifera",
                 "Eunicella_verrucosa", "Hydra_viridissima", "Hydra_oligactis", "Physalia_physalia", "Abylopsis_tetragona",
                 "Craseo_lathetica", "Nanomia_bijuga", "Agalma_elegans", "Periphyla_periphyla")
    Placozoa = c("Trichoplax_adhaerens")
    Porifera = c("Cliona_varians", "Sycon_coactum", "Sycon_ciliatum", "Corticium_candelabrum", "Oscarella_carmela", "Hyalonema_populiferum",
                 "Aphrocallistes_vastus", "Rossella_fibulata", "Sympagella_nux", "Ircinia_fasciculata", "Chondrilla_nucula", "Amphimedon_queenslandica",
                 "Petrosia_ficiformis", "Spongilla_lacustris", "Pseudospongosorites_suberitoides", "Mycale_phylophylla", "Latrunculia_apicalis", "Crella_elegans",
                 "Kirkpatrickia_variolosa")
    Ctenophora = c("Euplokamis_dunlapae", "Vallicula_sp", "Coeloplana_astericola", "Hormiphora_californica", "Hormiphora_palmata",
                   "Pleurobrachia_pileus", "Pleurobrachia_bachei", "Pleurobrachia_sp_South_Carolina_USA", "Cydippida_sp_Maryland_USA",
                   "Callianira_Antarctica", "Mertensiidae_sp_Antarctica", "Mertensiidae_sp_Washington_USA", "Cydippida_sp",
                   "Dryodora_glandiformis", "Lobatolampea_tetragona", "Beroe_abyssicola", "Beroe_sp_Antarctica", "Beroe_ovata",
                   "Beroe_sp_Queensland_Australia", "Beroe_forskalii", "Ocyropsis_sp_Bimini_Bahamas", "Ocyropsis_crystallina",
                   "Ocyropsis_sp_Florida_USA", "Bolinopsis_infundibulum", "Mnemiopsis_leidyi", "Bolinopsis_ashleyi", 
                   "Lobata_sp_Punta_Arenas_Argentina", "Eurhamphaea_vexilligera", "Cestum_veneris", "Ctenophora_sp_Florida_USA")
    Clade_Outgroup = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis")
    # Collate the list of taxa to keep
    taxa_to_keep <- c(Bilateria, Cnidaria, Placozoa, Porifera, Ctenophora, Clade_Outgroup)
  } else if (trim == "Keep_Ctenophora"){
    # Keep all Ctenophora species and one species from each other clade
    Bilateria = "Homo_sapiens"
    Cnidaria = "Hydra_vulgaris"
    Placozoa = "Trichoplax_adhaerens"
    Porifera = "Cliona_varians"
    Ctenophora = c("Euplokamis_dunlapae", "Vallicula_sp", "Coeloplana_astericola", "Hormiphora_californica", "Hormiphora_palmata",
                   "Pleurobrachia_pileus", "Pleurobrachia_bachei", "Pleurobrachia_sp_South_Carolina_USA", "Cydippida_sp_Maryland_USA",
                   "Callianira_Antarctica", "Mertensiidae_sp_Antarctica", "Mertensiidae_sp_Washington_USA", "Cydippida_sp",
                   "Dryodora_glandiformis", "Lobatolampea_tetragona", "Beroe_abyssicola", "Beroe_sp_Antarctica", "Beroe_ovata",
                   "Beroe_sp_Queensland_Australia", "Beroe_forskalii", "Ocyropsis_sp_Bimini_Bahamas", "Ocyropsis_crystallina",
                   "Ocyropsis_sp_Florida_USA", "Bolinopsis_infundibulum", "Mnemiopsis_leidyi", "Bolinopsis_ashleyi", 
                   "Lobata_sp_Punta_Arenas_Argentina", "Eurhamphaea_vexilligera", "Cestum_veneris", "Ctenophora_sp_Florida_USA")
    Clade_Outgroup = "Salpingoeca_rosetta"
    # Collate the list of taxa to keep
    taxa_to_keep <- c(Bilateria, Cnidaria, Placozoa, Porifera, Ctenophora, Clade_Outgroup)
  } else if (trim == "Trim_all"){
    # Keep one species from each clade
    Bilateria = "Homo_sapiens"
    Cnidaria = "Hydra_vulgaris"
    Placozoa = "Trichoplax_adhaerens"
    Porifera = "Cliona_varians"
    Ctenophora = "Dryodora_glandiformis"
    Clade_Outgroup = "Salpingoeca_rosetta"
    # Collate the list of taxa to keep
    taxa_to_keep <- c(Bilateria, Cnidaria, Placozoa, Porifera, Ctenophora, Clade_Outgroup)
  }
  
  # Root at outgroup
  root_outgroup = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis")
  m_tree <- root(m_tree, root_outgroup)
  
  # If trimming, drop taxa and rename clades
  if (trim == "Keep_Ctenophora"){
    # Drop all other taxa
    m_tree <- keep.tip(m_tree, taxa_to_keep)
    # Rename the remaining taxa into the clade name
    m_tree$tip.label[which(m_tree$tip.label == "Salpingoeca_rosetta")] <- "Choanoflagellata"
    m_tree$tip.label[which(m_tree$tip.label == "Cliona_varians")] <- "Porifera"
    m_tree$tip.label[which(m_tree$tip.label == "Trichoplax_adhaerens")] <- "Placozoa"
    m_tree$tip.label[which(m_tree$tip.label == "Hydra_vulgaris")] <- "Cnidaria"
    m_tree$tip.label[which(m_tree$tip.label == "Homo_sapiens")] <- "Bilateria"
  } else if (trim == "Trim_all") {
    # Drop all other taxa
    m_tree <- keep.tip(m_tree, taxa_to_keep)
    # Rename the remaining taxa into the clade name
    m_tree$tip.label[which(m_tree$tip.label == "Salpingoeca_rosetta")] <- "Choanoflagellata"
    m_tree$tip.label[which(m_tree$tip.label == "Cliona_varians")] <- "Porifera"
    m_tree$tip.label[which(m_tree$tip.label == "Trichoplax_adhaerens")] <- "Placozoa"
    m_tree$tip.label[which(m_tree$tip.label == "Hydra_vulgaris")] <- "Cnidaria"
    m_tree$tip.label[which(m_tree$tip.label == "Dryodora_glandiformis")] <- "Ctenophora"
    m_tree$tip.label[which(m_tree$tip.label == "Homo_sapiens")] <- "Bilateria"
  }
  
  # Return the reformatted tree
  return(m_tree)
}



color.code.metazoan.clades <- function(m_tree, trimmed = "No", color = TRUE){
  # Quick function to color code Metazoan clades
  
  if (trimmed == "No" | trimmed == "FALSE"){
    # Keep all species in all clades
    Bilateria = c("Homo_sapiens", "Strongylocentrotus_purpatus", "Hemithris_psittacea", "Capitella_teleta", "Drosophila_melanogaster",
                  "Daphnia_pulex")
    Cnidaria = c("Hydra_vulgaris", "Bolocera_tuediae", "Aiptasia_pallida", "Hormathia_digitata", "Nematostella_vectensis", "Acropora_digitifera",
                 "Eunicella_verrucosa", "Hydra_viridissima", "Hydra_oligactis", "Physalia_physalia", "Abylopsis_tetragona",
                 "Craseo_lathetica", "Nanomia_bijuga", "Agalma_elegans", "Periphyla_periphyla")
    Placozoa = c("Trichoplax_adhaerens")
    Porifera = c("Cliona_varians", "Sycon_coactum", "Sycon_ciliatum", "Corticium_candelabrum", "Oscarella_carmela", "Hyalonema_populiferum",
                 "Aphrocallistes_vastus", "Rossella_fibulata", "Sympagella_nux", "Ircinia_fasciculata", "Chondrilla_nucula", "Amphimedon_queenslandica",
                 "Petrosia_ficiformis", "Spongilla_lacustris", "Pseudospongosorites_suberitoides", "Mycale_phylophylla", "Latrunculia_apicalis", "Crella_elegans",
                 "Kirkpatrickia_variolosa")
    Ctenophora = c("Euplokamis_dunlapae", "Vallicula_sp", "Coeloplana_astericola", "Hormiphora_californica", "Hormiphora_palmata",
                   "Pleurobrachia_pileus", "Pleurobrachia_bachei", "Pleurobrachia_sp_South_Carolina_USA", "Cydippida_sp_Maryland_USA",
                   "Callianira_Antarctica", "Mertensiidae_sp_Antarctica", "Mertensiidae_sp_Washington_USA", "Cydippida_sp",
                   "Dryodora_glandiformis", "Lobatolampea_tetragona", "Beroe_abyssicola", "Beroe_sp_Antarctica", "Beroe_ovata",
                   "Beroe_sp_Queensland_Australia", "Beroe_forskalii", "Ocyropsis_sp_Bimini_Bahamas", "Ocyropsis_crystallina",
                   "Ocyropsis_sp_Florida_USA", "Bolinopsis_infundibulum", "Mnemiopsis_leidyi", "Bolinopsis_ashleyi", 
                   "Lobata_sp_Punta_Arenas_Argentina", "Eurhamphaea_vexilligera", "Cestum_veneris", "Ctenophora_sp_Florida_USA")
    Clade_Outgroup = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis")
    # Prepare columns for dataframe creation
    taxa_to_keep <- c(Bilateria, Cnidaria, Placozoa, Porifera, Ctenophora, Clade_Outgroup)
    taxa_groups <- c(rep("Bilateria", length(Bilateria)), rep("Cnidaria", length(Cnidaria)), rep("Placozoa", length(Placozoa)), 
                     rep("Porifera", length(Porifera)), rep("Ctenophora", length(Ctenophora)), rep("Choanoflagellata", length(Clade_Outgroup)))
    clades_to_keep <- c()
    clade_groups <- c()
    if (color == TRUE){
      taxa_colors <- c(rep("black", length(Bilateria)), rep("red", length(Cnidaria)), rep("green", length(Placozoa)), 
                       rep("yellow", length(Porifera)), rep("blue", length(Ctenophora)), rep("gray", length(Clade_Outgroup)))
      clade_colors <- c()
    } else if (color == FALSE){
      taxa_colors <- c(rep("black", length(Bilateria)), rep("black", length(Cnidaria)), rep("black", length(Placozoa)), 
                       rep("black", length(Porifera)), rep("black", length(Ctenophora)), rep("black", length(Clade_Outgroup)))
      clade_colors <- c()
    }
  } else if (trimmed == "Keep_Ctenophora"){
    # Keep all Ctenophora species and one species from each other clade
    Bilateria = "Bilateria"
    Cnidaria = "Cnidaria"
    Placozoa = "Placozoa"
    Porifera = "Porifera"
    Ctenophora = c("Euplokamis_dunlapae", "Vallicula_sp", "Coeloplana_astericola", "Hormiphora_californica", "Hormiphora_palmata",
                   "Pleurobrachia_pileus", "Pleurobrachia_bachei", "Pleurobrachia_sp_South_Carolina_USA", "Cydippida_sp_Maryland_USA",
                   "Callianira_Antarctica", "Mertensiidae_sp_Antarctica", "Mertensiidae_sp_Washington_USA", "Cydippida_sp",
                   "Dryodora_glandiformis", "Lobatolampea_tetragona", "Beroe_abyssicola", "Beroe_sp_Antarctica", "Beroe_ovata",
                   "Beroe_sp_Queensland_Australia", "Beroe_forskalii", "Ocyropsis_sp_Bimini_Bahamas", "Ocyropsis_crystallina",
                   "Ocyropsis_sp_Florida_USA", "Bolinopsis_infundibulum", "Mnemiopsis_leidyi", "Bolinopsis_ashleyi", 
                   "Lobata_sp_Punta_Arenas_Argentina", "Eurhamphaea_vexilligera", "Cestum_veneris", "Ctenophora_sp_Florida_USA")
    Clade_Outgroup = "Choanoflagellata"
    # Prepare columns for dataframe creation
    taxa_to_keep <- Ctenophora
    taxa_groups <- rep("Ctenophora", length(Ctenophora))
    clades_to_keep <- c(Bilateria, Cnidaria, Placozoa, Porifera, Clade_Outgroup)
    clade_groups <- c(rep("Bilateria", length(Bilateria)), rep("Cnidaria", length(Cnidaria)), rep("Placozoa", length(Placozoa)), 
                      rep("Porifera", length(Porifera)), rep("Choanoflagellata", length(Clade_Outgroup)))
    if (color == TRUE){
      taxa_colors <- rep("blue", length(Ctenophora))
      clade_colors <- c(rep("black", length(Bilateria)), rep("red", length(Cnidaria)), rep("green", length(Placozoa)), 
                        rep("yellow", length(Porifera)), rep("gray", length(Clade_Outgroup)))
    } else if (color == FALSE){
      taxa_colors <- rep("black", length(Ctenophora))
      clade_colors <- c(rep("black", length(Bilateria)), rep("black", length(Cnidaria)), rep("black", length(Placozoa)), 
                        rep("black", length(Porifera)), rep("black", length(Clade_Outgroup)))
    }
  } else if (trimmed == "Trim_all"){
    # Keep one species from each clade
    Bilateria = "Bilateria"
    Cnidaria = "Cnidaria"
    Placozoa = "Placozoa"
    Porifera = "Porifera"
    Ctenophora = "Ctenophora"
    Clade_Outgroup = "Choanoflagellata"
    # Prepare columns for dataframe creation
    taxa_to_keep <- c()
    taxa_clades <- c()
    clades_to_keep <- c(Bilateria, Cnidaria, Placozoa, Porifera, Ctenophora, Clade_Outgroup)
    clade_groups <- c(rep("Bilateria", length(Bilateria)), rep("Cnidaria", length(Cnidaria)), rep("Placozoa", length(Placozoa)), 
                      rep("Porifera", length(Porifera)), rep("Ctenophora", length(Ctenophora)), rep("Choanoflagellata", length(Clade_Outgroup)))
    if (color == TRUE){
      taxa_colors <- c()
      clade_colors <- c(rep("black", length(Bilateria)), rep("red", length(Cnidaria)), rep("green", length(Placozoa)), 
                        rep("yellow", length(Porifera)), rep("blue", length(Ctenophora)), rep("gray", length(Clade_Outgroup)))
    } else if (color == FALSE){
      taxa_colors <- c()
      clade_colors <- c(rep("black", length(Bilateria)), rep("black", length(Cnidaria)), rep("black", length(Placozoa)), 
                        rep("black", length(Porifera)), rep("black", length(Ctenophora)), rep("black", length(Clade_Outgroup)))
    }
  }
  
  # Create dataframe with tip information for taxa
  if (length(taxa_to_keep) == 0){
    taxa_lab_df <- NA
  } else if (length(taxa_to_keep) > 0) {
    all_taxa <- taxa_to_keep
    all_taxa_split <- strsplit(all_taxa, "_")
    taxa_tip_df <- data.frame(taxa = all_taxa,
                              taxa_prettyprint = gsub("_" ," ", all_taxa),
                              generic_name = sapply(all_taxa_split, "[[", 1),
                              generic_initial =  as.character(sapply(sapply(all_taxa_split, "[[", 1), 
                                                                     function(x){toupper(paste(substring(x, 1, 1), collapse = ""))})), 
                              specific_name = sapply(all_taxa_split, get.specific.taxa.name, 1),
                              clade = taxa_groups,
                              color = taxa_colors)
    taxa_lab_df <- dplyr::mutate(taxa_tip_df, 
                                 short_lab = glue('italic("{generic_initial}. {specific_name}")'),
                                 short_lab_noformat = glue("{generic_initial}. {specific_name}"),
                                 long_lab = glue('italic("{taxa_prettyprint}")'),
                                 short_name = glue("<i style='color:{color}'>{generic_initial} {specific_name}</i>"),
                                 long_name = glue("<i style='color:{color}'>{taxa_prettyprint}</i>") ) 
  }
  
  # Create dataframe with tip information for taxa
  if (length(clades_to_keep) == 0){
    clade_lab_df <- NA
  } else if (length(clades_to_keep) > 0) {
    clades <- clades_to_keep
    clade_tip_df <- data.frame(taxa = clades,
                               taxa_prettyprint = clades,
                               generic_name = clades,
                               generic_initial = clades, 
                               specific_name = rep("", length(clades)),
                               clade = clade_groups,
                               color = clade_colors)
    clade_lab_df <- dplyr::mutate(clade_tip_df, 
                                  short_lab = glue('italic("{taxa_prettyprint}")'),
                                  long_lab = glue('italic("{taxa_prettyprint}")'),
                                  short_name = glue("<i style='color:{color}'>{taxa_prettyprint}</i>"),
                                  long_name = glue("<i style='color:{color}'>{taxa_prettyprint}</i>") ) 
  }
  
  # Attach dataframes (if they both exist), or else pick the one that does exist
  if ((class(clade_lab_df) == "data.frame") & (class(taxa_lab_df) == "data.frame")){
    tip_lab_df <- rbind(taxa_lab_df, clade_lab_df)
  } else if ((class(clade_lab_df) == "logical") & (class(taxa_lab_df) == "data.frame")){
    tip_lab_df <- taxa_lab_df
  } else if ((class(clade_lab_df) == "data.frame") & (class(taxa_lab_df) == "logical")){
    tip_lab_df <- clade_lab_df
  }
  
  # Return the tip label dataframe
  return(tip_lab_df)
}

# Small function to get 2nd word onwards from taxa name for Metazoan species or Plants species
get.specific.taxa.name <- function(t, n){
  if (length(t) > n) {
    new_t <- paste(t[(n + 1):length(t)], collapse = " ")
  }
  else {
    new_t <- paste(t, collapse = " ")
  }
  return(new_t)
}



determine.outlier.tree.file <- function(i, outlier_df, tree_file_paths){
  # Small function to determine which species tree is related to which outlier_df row
  
  # Extract row
  row <- outlier_df[i,]
  
  # Determine dataset
  if (row$dataset == "1KP") {
    ds <- "Plants"
  } else if (row$dataset == "Pease2016"){
    ds <- "Tomatoes"
  } else if (row$dataset == "Whelan2017"){
    ds <- "Metazoan"
  } else if (row$dataset == "Vanderpool2017"){
    ds <- "Primates"
  }
  
  # Determine tree estimation method
  if (row$support_value_type == "PP"){
    method = "ASTRAL"
  } else if (row$support_value_type == "BS"){
    method = "CONCAT"
  }
  
  # Determine test
  if (row$tree1 == "None"){
    test = "NoTest"
    test_full = "NoTest"
  } else {
    test = row$test
  }
  
  # Determine whether test was pass or fail
  if (test != "NoTest"){
    test_bool <- row$tree1
    if (test_bool == "Pass"){
      test_type = "pass"
    } else if (test_bool == "Fail"){
      test_type = "fail"
    }
    test_full = paste0(test, "_", test_type)
  }
  
  # Identify file that has the matching dataset, tree estimation method, and full_test (from vector of file paths)
  row_file <- grep(ds, grep(method, grep(test_full, tree_file_paths, value = TRUE), value = TRUE), value = TRUE)
  
  # Return row file
  return(row_file)
}



determine.outlier.plot.name <- function(i, outlier_df){
  # Small function to determine which species tree is related to which outlier_df row
  
  # Extract row
  row <- outlier_df[i,]
  
  # Determine dataset
  if (row$dataset == "1KP") {
    ds <- "Plants"
  } else if (row$dataset == "Pease2016"){
    ds <- "Tomatoes"
  } else if (row$dataset == "Whelan2017"){
    ds <- "Metazoan"
  } else if (row$dataset == "Vanderpool2017"){
    ds <- "Primates"
  }
  
  # Determine tree estimation method
  if (row$support_value_type == "PP"){
    method = "ASTRAL"
  } else if (row$support_value_type == "BS"){
    method = "CONCAT"
  }
  
  # Determine test
  if (row$tree1 == "None"){
    test = "NoTest"
    test_full = "NoTest"
  } else {
    test = row$test
  }
  
  # Determine whether test was pass or fail
  if (test != "NoTest"){
    test_bool <- row$tree1
    if (test_bool == "Pass"){
      test_type = "pass"
    } else if (test_bool == "Fail"){
      test_type = "fail"
    }
    test_full = paste0(test, "_", test_type)
  }
  
  # Create the output file name
  op_file_name <- paste0(ds, "_", test_full, "_", method, "_")
  
  # Return output file name
  return(op_file_name)
}



reformat.ASTRAL.tree.for.plotting <- function(tree, add.arbitrary.terminal.branches = FALSE, terminal.branch.length = 1, strip.nodes = FALSE,
                                              scale.tree.length = FALSE, new.tree.length = 1){
  ## Function to nicely format an ASTRAL tree for plotting
  # ASTRAL does not output terminal branch lengths, so these are arbitrarily given the length 1
  # ASTRAL includes posterior probabilities: these are stripped
  # Also have option to scale trees so they align nicely for ggdensitree
  
  # Add terminal branch lengths of arbitrary length
  if (add.arbitrary.terminal.branches == TRUE){
    # Identify which edges correspond to terminal branches
    n = length(tree$tip.label)
    missing_branch_indexes <- sapply(1:n,function(x,y)   which(y==x),y=tree$edge[,2])
    # Set length of terminal branches to 1 (or to whatever the terminal.branch.length variable is)
    tree$edge.length[missing_branch_indexes] <- terminal.branch.length
  }
  
  # Remove node support values
  if (strip.nodes == TRUE){
    # Don't want node labels - for quartnetGoFtest, need just taxa and branch lengths
    tree$node.label <- NULL
  }
  
  # Rescale tree length
  if (scale.tree.length == TRUE){
    tree$edge.length <- tree$edge.length / max(nodeHeights(tree)[,2]) * new.tree.length
  }
  
  # Return the tree
  return(tree)
}



label.plant.taxa <- function(tips, classification_df){
  # Function to take in a list of tip from a tree from the 1000 Plants transcriptome dataset and return a nice dataframe for plotting the tree in ggplot
  
  # Reduce classification_df to only contain species that are present in the list of tips
  classification_df <- classification_df[(classification_df$Code %in% tips), ]
  # Rearrange classification_df to have rows in same order as tips
  classification_df <- classification_df[match(tips, classification_df$Code), ]
  # Split all the taxa by the spaces to create the short name (e.g. `S. name` instead of `Species name`)
  all_taxa_split <- strsplit(classification_df$Species, " ")
  # Create dataframe for labelling tips
  taxa_tip_df <- data.frame(code = classification_df$Code,
                            taxa = classification_df$Species,
                            generic_name = sapply(all_taxa_split, "[[", 1),
                            generic_initial =  as.character(sapply(sapply(all_taxa_split, "[[", 1), 
                                                                   function(x){toupper(paste(substring(x, 1, 1), collapse = ""))})), 
                            specific_name = sapply(all_taxa_split, get.specific.taxa.name, 1),
                            clade = gsub("_", " ", classification_df$Very.Brief.Classification) )
  # Mutate the labels dataframe to add formatting for taxa labels
  taxa_lab_df <- dplyr::mutate(taxa_tip_df, 
                               short_lab = glue('italic("{generic_initial}. {specific_name}")'),
                               long_lab = glue('italic("{taxa}")'),
                               short_name = glue("<i style='color:black'>{generic_initial} {specific_name}</i>"),
                               long_name = glue("<i style='color:black'>{taxa}</i>") )
  # Return the formatted labels 
  return(taxa_lab_df)
}



extract.one.plant.tip <- function(tip, classification_df){
  # Function to take one tip from the 1000 Plants transcriptome dataset and return the taxa and classification
  
  # Extract the row containing information about this taxa
  row <- classification_df[which(classification_df$Code == tip),]
  # Return that row
  return(row)
}



reduce.tree.to.clades <- function(tree, labels_df, annotations_df){
  # Function to take a label df and a tree and reduce each clade to a single tip
  ## NOTE: This is not a good way to do this, because two clades (Chlorophyta and Streptophyte_algae) are not monophyletic for the Plants dataset
  
  # Use ggtree::fortify() to create a dataframe from the tree
  tree_df <- subset(fortify(tree), isTip)
  # Attach the information from the tips dataframe
  tree_df <- cbind(tree_df, labels_df[match(tree_df$label, labels_df$code), ])
  tree_df <- cbind(tree_df, annotations_df[match(tree_df$label, annotations_df[(annotations_df$Code %in% tree$tip.label), ]$Code), ])
  
  # Extract a list of the unique clades
  unique_clades <- unique(tree_df$clade)
  
  c = "Eudicots"
  # Create a new copy of the tree
  tree_new <- tree
  
  # Iterate through each unique clade
  for (c in unique_clades){
    # Reduce a copy of the dataframe to taxa from this clade only
    c_df <- tree_df[tree_df$clade == c, ]
    # Extract the nodes from this clade
    c_nodes <- c_df$node
    # Extract the list of tips from this clade
    c_tips <- c_df$label
    # Remove any ANAGrade tips (not a clade)
    if (c == "ANAGrade"){
      tree_new <- drop.tip(tree_new, c_tips, trim.internal = TRUE)
    }
    # Check if the clade is monophyletic
    mono_check <- is.monophyletic(tree, c_tips)
    print(paste0(c, " : ", mono_check))
    if (mono_check == TRUE){
      # If the clade is monophyletic, reduce it to one tip picked at random
      new_c_tip <- c_tips[1]
      drop_c_tips <- c_tips[2:length(c_tips)]
      # Drop the tips
      tree_new <- drop.tip(tree_new, drop_c_tips, trim.internal = TRUE)
    }
  } # end iterating through clades
}




