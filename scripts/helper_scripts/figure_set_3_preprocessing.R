#### specify settings #### 
random_seed <- 1
set.seed(random_seed)

# settings <- list(c(T, F, T, F, T), #experiment to check non-kronecker corrmat effects
#                  c(T, F, F, F, T),
#                  c(F, F, F, F, T))[[3]]
# settings <- list(c(T, F, F, F, T), #experiment to check trait / tissue corrmats using only intersecting genes
#                  c(F, F, F, F, T))[[2]]
settings <- c(F, F, F, F, T)
use_random_DE_genes <- settings[1]
use_random_Pred_genes <- settings[2]
use_kronecker_interactions <- settings[3]
estimate_interaction_corrmat <- settings[4]
use_all_genes_for_trait_tiss_corrmats <- settings[5]

estimate_trait_corr_mats <- F
use_tissue_specific_corr_mats <- F
estimate_composite_corrmat <- F
estimate_tissue_corrmat <- F
fit_model <- F

#### load libraries ####
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(cmdstanr)
library(posterior)
library(caret)
library(MASS)
library(ks)
library(data.table)
library(edgeR)
library(limma)
library(DESeq2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(org.Rn.eg.db)
library(biomaRt)
library(plotrix)
library(ggplot2)
library(testit)
library(circlize)
library(jpeg)
library(foreach)
library(doParallel)
library(pracma)
library(MotrpacRatTraining6mo) # v1.6.0
# also attaches MotrpacRatTraining6moData v1.8.0


#### define functions ####
source(file = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/helper_scripts/deg-trait_functions.R")

#### get ortholog map ####
gene_map <- MotrpacRatTraining6moData::RAT_TO_HUMAN_GENE
gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID <- gsub(gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID, pattern = "\\..*", replacement = "")

#### load TWAS results ####
motrpac_gtex_map = c('t30-blood-rna'='Whole_Blood',
                     't52-hippocampus'='Brain_Hippocampus',
                     't53-cortex'='Brain_Cortex',
                     't54-hypothalamus'='Brain_Hypothalamus',
                     't55-gastrocnemius'='Muscle_Skeletal',
                     't56-vastus-lateralis'='Muscle_Skeletal',
                     't58-heart'='Heart_Left_Ventricle',
                     't59-kidney'='Kidney_Cortex',
                     't60-adrenal'='Adrenal_Gland',
                     't61-colon'='Colon_Transverse',
                     't62-spleen'='Spleen',
                     't63-testes'='Testis',
                     't64-ovaries'='Ovary',
                     't66-lung'='Lung',
                     't67-small-intestine'='Small_Intestine_Terminal_Ileum',
                     't68-liver'='Liver',
                     't70-white-adipose'='Adipose_Subcutaneous')

#now let's use the TWAS results from Barbeira et al. 2021:
# obtained from <https://zenodo.org/record/3518299/files/spredixcan_eqtl.tar.gz>
twas_results_directory <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/eqtl/"
gwas_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/imputed_gwas_hg38_1.1/"
gwas_summary_files <- list.files(gwas_dir)
gwas_summary_files <- gwas_summary_files[-grep(gwas_summary_files, pattern = "README")]
gwas_names <- stringr::str_replace_all(gwas_summary_files, ".txt.gz", "")
gwas_names <- stringr::str_replace_all(gwas_names, "imputed_", "")

available_twas_traits <- list.files(path = twas_results_directory)
available_twas_tissues <- table(sapply(available_twas_traits, function(x) strsplit(x, "__PM__")[[1]][2]))
available_twas_traits <- table(sapply(gsub(x = available_twas_traits, "spredixcan_igwas_gtexmashrv8_", ""), 
                                      function(x) strsplit(x, "__P")[[1]][1]))

#load gwas / twas metadata
trait_categories <- read.csv("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/gwas_metadata.csv", header = T)
trait_categories$Category[trait_categories$new_Phenotype == "Multiple_Sclerosis"] <- "Psychiatric-neurologic" #correct original coding from "Cardiometabolic"
traitwise_partitions <- trait_categories[,c("Tag", "Category")]

#helpful objects
tissue_abbr <- MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[grep("t[0-9][0-9]", 
                                                                     x = names(MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV))]
tissue_abbr_rev <- setNames(names(tissue_abbr), tissue_abbr)

#compile the twas output into one big data frame
use_all_cats <- T
if(!exists("some.twas")){
  
  all_twas <- lapply(setNames(gwas_names, gwas_names), function(trait) data.table());
  for(trait in gwas_names){
    print(trait)
    all_tissues <- lapply(setNames(names(motrpac_gtex_map), names(motrpac_gtex_map)), function(tissue) data.table());
    for(tissue in names(motrpac_gtex_map)){
      twas <- fread(file = paste0(twas_results_directory, "spredixcan_igwas_gtexmashrv8_", trait, "__PM__", motrpac_gtex_map[names(motrpac_gtex_map) == tissue], ".csv"))
      twas <- twas[,c("gene", "gene_name", "zscore", "pvalue")]
      twas <- twas[-which(is.na(twas$zscore) & is.na(twas$pvalue)),]
      twas$tissue <- tissue
      twas$trait <- trait
      all_tissues[[tissue]] <- twas
    }
    all_twas[[trait]] <- do.call(rbind, all_tissues)
  }
  all.twas <- do.call(rbind, all_twas)
  all.twas$trait_category <- traitwise_partitions$Category[match(all.twas$trait, traitwise_partitions$Tag)]
  
  #filtering by category
  if(use_all_cats){
    salient.categories <- unique(traitwise_partitions$Category)
  } else {
    salient.categories <- c("Cardiometabolic", "Aging", "Anthropometric", 
                            "Immune", "Psychiatric-neurologic")
  }
  some.twas <- all.twas[all.twas$trait_category %in% salient.categories]
  salient_twas <- setNames(unique(some.twas$trait), unique(some.twas$trait))
  
  #cleanup memory
  rm("all.twas")
  rm("all_twas")
  
}

#load or calculate the IHW results for these output
if(use_all_cats){
  
  if(!file.exists("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ihw_results_all.twas.RData")){
    ihw_twas <- as.data.frame(cbind(trait_tissue = as.factor(paste0(some.twas$trait, "~", some.twas$tissue)), pvalue = some.twas$pvalue))
    ihw_results <- IHW::ihw(pvalue ~ trait_tissue, data = ihw_twas, alpha = 0.05)  
    save(some.twas, file = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/PrediXcan_output_all.twas.RData")
    save(ihw_results, file = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ihw_results_all.twas.RData")
  } else {
    load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ihw_results_all.twas.RData")
  }
  
} else {
  
  if(!file.exists("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ihw_results_some.twas.RData")){
    ihw_twas <- as.data.frame(cbind(trait_tissue = as.factor(paste0(some.twas$trait, "~", some.twas$tissue)), pvalue = some.twas$pvalue))
    ihw_results <- IHW::ihw(pvalue ~ trait_tissue, data = ihw_twas, alpha = 0.05)  
    save(ihw_results, file = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ihw_results_some.twas.RData")
  } else {
    load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ihw_results_some.twas.RData")
  }
  
}
ihw_pvals <- IHW::adj_pvalues(ihw_results)
rm(ihw_results)

#randomize p-values if taking that approach
if(use_random_Pred_genes){
  ihw_pvals <- sample(ihw_pvals)
}
some.twas$adj_pvalue <- ihw_pvals

#### load in the motrpac DE results ####
if(!exists("rna_dea_ensembl")){
  rna_dea_ensembl <- MotrpacRatTraining6moData::GENE_UNIVERSES$ensembl_gene$TRNSCRPT
}

#check if we've already compiled and formatted focal node information and load if so, otherwise process
genes_tested_in_transcriptome_DEA_x_tissue <- MotrpacRatTraining6moData::GENE_UNIVERSES$ensembl_gene$TRNSCRPT
all_tested_genes <- unique(unlist(genes_tested_in_transcriptome_DEA_x_tissue))

node_metadata_path <- paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/node_metadata_list", 
                             ifelse(use_random_DE_genes, ".random", ""), ".RData")
if(!file.exists(node_metadata_path)){
  
  node_sets <- MotrpacRatTraining6moData::GRAPH_COMPONENTS$node_sets
  nodes_to_look_at_list <- list(c("1w_F1_M1", "1w_F-1_M-1"),
                                c("2w_F1_M1", "2w_F-1_M-1"),
                                c("4w_F1_M1", "4w_F-1_M-1"),
                                c("8w_F1_M1", "8w_F-1_M-1"))
  setNames(nodes_to_look_at_list, paste0(2^(0:3), "w"))
  
  node_metadata_list <- lapply(nodes_to_look_at_list, function(nodes_to_look_at){
    node_metadata <- lapply(setNames(nodes_to_look_at, nodes_to_look_at), function(node_to_look_at)
      cbind(do.call(rbind, strsplit(node_sets[node_to_look_at][[node_to_look_at]][grep(
        node_sets[node_to_look_at][[node_to_look_at]], pattern = "TRNSCRPT")], ";")), node_to_look_at))
    node_metadata <- as.data.table(do.call(rbind, node_metadata))
    colnames(node_metadata) <- c("ome","tissue","ensembl_gene", "node")
    node_metadata$rat_gene_symbol <- gene_map$RAT_SYMBOL[match(node_metadata$ensembl_gene, gene_map$RAT_ENSEMBL_ID)]
    node_metadata$human_gene_symbol <- gene_map$HUMAN_ORTHOLOG_SYMBOL[match(node_metadata$ensembl_gene, gene_map$RAT_ENSEMBL_ID)]
    node_metadata$human_ensembl_gene <- gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID[match(node_metadata$ensembl_gene, gene_map$RAT_ENSEMBL_ID)]
    node_metadata$human_ensembl_gene <- gsub(node_metadata$human_ensembl_gene, pattern = "\\..*", replacement = "")
    node_metadata$cluster <- paste0(node_metadata$tissue, "-", node_metadata$node)
    node_metadata$cluster <- paste0(node_metadata$tissue, "-", "sex_homogeneous_changing")
    node_metadata
  })
  
  if(use_random_DE_genes){
    #shuffle gene indices, preserving within-tissue structure of tested genes
    #need to sample each tissue's gene universe and pick new assignments, respecting patterns of overlap
    #can do this iteratively -- first shuffle indices, reassign DE genes according to new map, and
    #if assignment is not in a tissue's gene universe, redo the assignment
    #bit hacky but much easier than the full solution 
    
    random_gene_map <- setNames(all_tested_genes, sample(all_tested_genes))
    invalid_rgm <- T
    #not very efficient and can get trapped, but works just fine here
    #can create prop_assignments outside loop and only iterate over the failed ones
    #also have an indicator reinitialize every n rounds through
    while(invalid_rgm){
      prop_assignments <- do.call(rbind, node_metadata_list)[,c("tissue", "ensembl_gene")]
      prop_assignments$rand_ensembl_gene <- random_gene_map[prop_assignments$ensembl_gene]
      tiss_x_randgene <- split(prop_assignments$rand_ensembl_gene, prop_assignments$tissue)
      bad_matches <- unique(unlist(lapply(setNames(names(tiss_x_randgene), names(tiss_x_randgene)), function(tiss){
        tiss_x_randgene[[tiss]][!(tiss_x_randgene[[tiss]] %in% genes_tested_in_transcriptome_DEA_x_tissue[[tiss]])]
      })))
      if(length(bad_matches) == 0){
        invalid_rgm <- F  
      } else {
        good_matches <- setdiff(unique(prop_assignments$rand_ensembl_gene), bad_matches)
        random_gene_map <- c(random_gene_map[random_gene_map %in% good_matches], 
                             {
                               rematch <- random_gene_map[!(random_gene_map %in% good_matches)]
                               setNames(sample(rematch), names(rematch))
                             }
        )  
      }
      
    }
    
    # now construct the randomized node_metadata_list object
    node_metadata_list <- lapply(node_metadata_list, function(node_metadata_sublist){
      new_sublist <- node_metadata_sublist
      new_sublist$ensembl_gene <- random_gene_map[new_sublist$ensembl_gene]
      new_sublist$rat_gene_symbol <- gene_map$RAT_SYMBOL[match(new_sublist$ensembl_gene, gene_map$RAT_ENSEMBL_ID)]
      new_sublist$human_gene_symbol <- gene_map$HUMAN_ORTHOLOG_SYMBOL[match(new_sublist$ensembl_gene, gene_map$RAT_ENSEMBL_ID)]
      new_sublist$human_ensembl_gene <- gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID[match(new_sublist$ensembl_gene, gene_map$RAT_ENSEMBL_ID)]
      new_sublist$human_ensembl_gene <- gsub(new_sublist$human_ensembl_gene, pattern = "\\..*", replacement = "")
      new_sublist
    })
    names(node_metadata_list) <- paste0(2^(0:3), "w")
  }
  
  save(node_metadata_list, file = node_metadata_path)  
  
} else {
  
  load(file = node_metadata_path)
  
}

#look at overlap between tissues
na.remove <- function(x) x[!is.na(x)]
DE_genes_x_tissue <- split(node_metadata_list$`8w`$human_ensembl_gene, node_metadata_list$`8w`$tissue)
DE_genes_x_tissue <- lapply(lapply(DE_genes_x_tissue, na.remove), unique)
all_genes_x_tissue <- lapply(rna_dea_ensembl, function(x) unique(na.remove(gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID[match(x, gene_map$RAT_ENSEMBL_ID)])))
if(use_all_genes_for_trait_tiss_corrmats){
  if(estimate_tissue_corrmat){
    tissue_corr_mat <- estimate_correlations(y = DE_genes_x_tissue, t = all_genes_x_tissue, print_progress = T)
    save(tissue_corr_mat, file = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/tissue_corr_mat", 
                                        ifelse(use_random_Pred_genes | use_random_DE_genes, "_random-genes", ""), 
                                        ifelse(use_all_genes_for_trait_tiss_corrmats, "", "_intersecting-genes"), 
                                        ".RData"))  
  } else {
    load(file = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/tissue_corr_mat", 
                       ifelse(use_random_Pred_genes | use_random_DE_genes, "_random-genes", ""), 
                       ifelse(use_all_genes_for_trait_tiss_corrmats, "", "_intersecting-genes"), 
                       ".RData"))
  }
  
  # pheatmap::pheatmap(tissue_corr_mat, breaks = -10:10/10, color = colorspace::diverging_hcl(21))
}

#find some useful summary stats for ortholog mapping

#expressed genes
mean(all_tested_genes %in% gene_map$RAT_ENSEMBL_ID)
range(sapply(genes_tested_in_transcriptome_DEA_x_tissue, function(tiss_i_genes){
  mean(tiss_i_genes %in% gene_map$RAT_ENSEMBL_ID)
}))

#DE genes
mean(unique(node_metadata_list$`8w`$ensembl_gene) %in% gene_map$RAT_ENSEMBL_ID)
range(sapply(split(node_metadata_list$`8w`$ensembl_gene, node_metadata_list$`8w`$tissue), 
             function(tiss_i_genes){
               mean(tiss_i_genes %in% gene_map$RAT_ENSEMBL_ID)
             }))

#find set intersect
n_deg_sigtwas_intersect <- as.data.frame(matrix(0, nrow = length(names(motrpac_gtex_map)), 
                                                ncol = length(salient_twas), 
                                                dimnames = list(names(motrpac_gtex_map), salient_twas)))
timepoint = "8w"
adj_pvalue_alpha <- 0.05
sig_twas_by_tissue <- lapply(setNames(names(motrpac_gtex_map), names(motrpac_gtex_map)), 
                             function(tissue_i) some.twas[some.twas$tissue == tissue_i & some.twas$adj_pvalue < adj_pvalue_alpha,])

tissue_code <- data.frame(tissue_name_release = names(MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV),
                          abbreviation = MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV)
for(tissue_i in names(motrpac_gtex_map)){
  
  print(tissue_i)
  
  sig_twas_by_trait <- lapply(setNames(salient_twas, salient_twas), 
                              function(trait_i) sig_twas_by_tissue[[tissue_i]][sig_twas_by_tissue[[tissue_i]]$trait == trait_i,])
  sig_twas_by_trait_genes <- lapply(setNames(salient_twas, salient_twas), 
                                    function(trait_i) sig_twas_by_trait[[trait_i]]$gene_name)
  
  
  n_intersect <- sapply(sig_twas_by_trait_genes, function(twas_genes) 
    length(intersect(twas_genes, node_metadata_list[[timepoint]]$human_gene_symbol[
      node_metadata_list[[timepoint]]$tissue == tissue_code$abbreviation[tissue_code$tissue_name_release == tissue_i]])))
  
  n_deg_sigtwas_intersect[tissue_i,names(n_intersect)] <- n_intersect
  
}

#tidy this up
n_deg_sigtwas_intersect <- n_deg_sigtwas_intersect[,order(apply(n_deg_sigtwas_intersect, 2, sum), decreasing = T)]
n_deg_sigtwas_intersect <- n_deg_sigtwas_intersect[,apply(n_deg_sigtwas_intersect, 2, sum) != 0]
rownames(n_deg_sigtwas_intersect) <- MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[rownames(n_deg_sigtwas_intersect)]
n_deg_sigtwas_intersect <- n_deg_sigtwas_intersect[order(match(rownames(n_deg_sigtwas_intersect), MotrpacRatTraining6moData::TISSUE_ORDER)),]
n_deg_sigtwas_intersect <- n_deg_sigtwas_intersect[nrow(n_deg_sigtwas_intersect):1,]
n_deg_sigtwas_intersect <- n_deg_sigtwas_intersect[!(rownames(n_deg_sigtwas_intersect) %in% c("OVARY", "TESTES")),]

#colsums
tissues <- setNames(rownames(n_deg_sigtwas_intersect), rownames(n_deg_sigtwas_intersect))
tissue_code <- MotrpacRatTraining6moData::TISSUE_ABBREV_TO_CODE

#max is all tested genes in all tissues?
genes_tested_in_transcriptome_DEA <- unique(unlist(MotrpacRatTraining6moData::GENE_UNIVERSES$ensembl_gene$TRNSCRPT))
all_orthologs_tested <- na.remove(gene_map$HUMAN_ORTHOLOG_SYMBOL[match(genes_tested_in_transcriptome_DEA, gene_map$RAT_ENSEMBL_ID)])

#or just in each tissue? this is more principled, since above doesn't correspond to *all* genes, just the multi-tissue expressed ones
genes_tested_in_transcriptome_DEA_x_tiss <- MotrpacRatTraining6moData::GENE_UNIVERSES$ensembl_gene$TRNSCRPT
all_orthologs_tested_x_tiss <- lapply(genes_tested_in_transcriptome_DEA_x_tiss, function(tiss_genes) 
  na.remove(gene_map$HUMAN_ORTHOLOG_SYMBOL[match(tiss_genes, gene_map$RAT_ENSEMBL_ID)])
)

sig_twas_by_trait_genes_matrix <- t(sapply(setNames(tissues, tissues), function(tissue_i){
  sig_twas_by_trait <- lapply(setNames(salient_twas, salient_twas), function(trait_i) 
    sig_twas_by_tissue[[tissue_code[tissue_i]]][sig_twas_by_tissue[[tissue_code[tissue_i]]]$trait == trait_i,])
  return(
    sapply(setNames(salient_twas, salient_twas), function(trait_i){
      twas_genes <- sig_twas_by_trait[[trait_i]]$gene_name
      twas_genes <- twas_genes[!is.na(twas_genes)]
      return(length(intersect(all_orthologs_tested_x_tiss[[tissue_i]], twas_genes)))
    }))
}))
sig_twas_by_trait_genes_matrix <- sig_twas_by_trait_genes_matrix[,colnames(n_deg_sigtwas_intersect)]
sig_twas_by_trait_genes_range <- apply(apply(sig_twas_by_trait_genes_matrix, 2, range), 2, paste0, collapse = " - ")
sig_twas_by_trait_genes <- apply(sig_twas_by_trait_genes_matrix, 2, mean)

prop_twas_are_degs <- t(sapply(rownames(n_deg_sigtwas_intersect), function(tissue) n_deg_sigtwas_intersect[tissue,] / sig_twas_by_trait_genes))
prop_twas_are_degs <- apply(prop_twas_are_degs, 2, unlist)
prop_twas_are_degs <- prop_twas_are_degs[,order(apply(prop_twas_are_degs, 2, mean), decreasing = T)]

#rowsums -- also need trait-specificity here

#pooling across all traits
all_twas_genes_tested <- unique(some.twas$gene_name) 
possible_genes <- intersect(all_orthologs_tested, all_twas_genes_tested)
compatible_twas_genes <- some.twas$gene_name %in% possible_genes
twas_by_tissue <- lapply(setNames(unique(some.twas$tissue), MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[unique(some.twas$tissue)]), function(tiss) {
  print(tiss)
  some.twas[some.twas$tissue == tiss & compatible_twas_genes,]
})
compatible_twas_genes <- some.twas$gene_name %in% possible_genes

if(!file.exists("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/twas_by_tissue_orig.RData")){
  twas_by_tissue_orig <- lapply(setNames(unique(some.twas$tissue), 
                                         MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[unique(some.twas$tissue)]), 
                                function(tiss) {print(tiss)
                                  some.twas[some.twas$tissue == tiss & compatible_twas_genes,]
                                })
  save(twas_by_tissue_orig, file = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/twas_by_tissue_orig.RData")
} else {
  load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/twas_by_tissue_orig.RData")
}

n_genes_in_nodes_matrix <- t(sapply(tissues, function(tiss){
  print(tiss)
  sapply(setNames(salient_twas, salient_twas), function(trait_i){
    twas_genes_tested <- twas_by_tissue_orig[[tiss]]$gene_name[twas_by_tissue[[tiss]]$trait == trait_i]
    length(intersect(twas_genes_tested, node_metadata_list[[timepoint]]$human_gene_symbol[node_metadata_list[[timepoint]]$tissue == tiss]))
  }) 
}))


#splitting into tisssue x trait pairs
all_twas_genes_tested_tiss_x_trait <- lapply(split(some.twas, some.twas$tissue), function(tissue_specific_twas){
  split(tissue_specific_twas$gene_name, tissue_specific_twas$trait)
})
possible_genes_tiss_x_trait <- lapply(tissues, function(tiss_i){
  lapply(salient_twas, function(trait_i){
    intersect(all_orthologs_tested_x_tiss[[tiss_i]], all_twas_genes_tested_tiss_x_trait[[TISSUE_ABBREV_TO_CODE[tiss_i]]][[trait_i]])
  })
})
if(!file.exists("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/twas_by_tissue.RData")){
  twas_by_tissue <- mclapply(tissues, function(tiss_i){
    do.call(rbind, lapply(salient_twas, function(trait_i){
      system(sprintf('echo "%s "', paste0(tiss_i, collapse="")))
      some.twas[some.twas$tissue == TISSUE_ABBREV_TO_CODE[tiss_i] & 
                  some.twas$trait == trait_i & 
                  some.twas$gene_name %in% possible_genes_tiss_x_trait[[tiss_i]][[trait_i]],]
    }))
  }, mc.cores = 8)
  save(twas_by_tissue, file = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/twas_by_tissue.RData")
} else {
  load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/twas_by_tissue.RData")
}

n_genes_in_nodes_matrix <- t(sapply(tissues, function(tiss){
  print(tiss)
  sapply(setNames(salient_twas, salient_twas), function(trait_i){
    twas_genes_tested <- twas_by_tissue[[tiss]]$gene_name[twas_by_tissue[[tiss]]$trait == trait_i]
    length(intersect(twas_genes_tested, node_metadata_list[[timepoint]]$human_gene_symbol[node_metadata_list[[timepoint]]$tissue == tiss]))
  }) 
}))


n_genes_in_nodes_matrix <- n_genes_in_nodes_matrix[rownames(n_deg_sigtwas_intersect),
                                                   colnames(n_deg_sigtwas_intersect)]
n_genes_in_nodes <- apply(n_genes_in_nodes_matrix, 1, mean)
n_genes_in_nodes_range <- apply(apply(n_genes_in_nodes_matrix, 1, range), 2, paste0, collapse = " - ")

prop_degs_are_twas <- (sapply(colnames(n_deg_sigtwas_intersect), function(trait) 
  n_deg_sigtwas_intersect[,trait] / n_genes_in_nodes[rownames(n_deg_sigtwas_intersect)]))
prop_degs_are_twas <- apply(prop_degs_are_twas, 2, unlist)
prop_degs_are_twas <- prop_degs_are_twas[,order(apply(prop_degs_are_twas, 2, mean), decreasing = T)]

#### estimate corrmats under bivariate probit ####

#look at overlap between traits
# if(use_tissue_specific_corr_mats){
#   if(estimate_trait_corr_mats){
#     adj_pvalue_alpha <- 0.05
#     twas_tiss <- some.twas$tissue[1]
#     tissues <- setNames(unique(some.twas$tissue), unique(some.twas$tissue))
#     trait_corr_mats <- lapply(tissues, function(twas_tiss){
#       twas_sub <- some.twas[some.twas$tissue == twas_tiss,]
#       twas_sub$gene <- gsub(twas_sub$gene, pattern = "\\..*", replacement = "")
#       twas_sub_hits <- twas_sub[twas_sub$adj_pvalue < adj_pvalue_alpha,]
#       na.remove <- function(x) x[!is.na(x)]
#       twas_genes_x_trait <- lapply(lapply(split(twas_sub_hits$gene, twas_sub_hits$trait), na.remove), unique)
#       all_genes_x_trait <- lapply(lapply(split(twas_sub$gene, twas_sub$trait), na.remove), unique)
#       trait_corr_mat <- estimate_correlations(y = twas_genes_x_trait, t = all_genes_x_trait, print_progress = T)
#       trait_corr_mat
#     })
#     #comment out, going with the composite strategy
#     #save(trait_corr_mats, file = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/trait_corr_mat", ifelse(use_random_Pred_genes, "_random-genes", ""), ".RData"))  
#   } else {
#     #load(file = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/trait_corr_mat", 
#     #ifelse(use_random_Pred_genes, "_random-genes", ""), ".RData"))
#   }
# 
#   trait_corr_mats <- trait_corr_mats[tissue_abbr_rev[rownames(n_deg_sigtwas_intersect)]]
#   
#   #obtain an average across traits?
#   traitset <- colnames(n_deg_sigtwas_intersect)
#   sapply(traitset, function(ti1) sapply(traitset, function(ti1) 1+1))
#   trait_corr_mat <- trait_corr_mats[[1]]
# }

#estimate composite trait corr mat?
if(use_all_genes_for_trait_tiss_corrmats){
  if(!use_tissue_specific_corr_mats){
    if(estimate_composite_corrmat){
      adj_pvalue_alpha <- 0.05
      tissues <- setNames(tissue_abbr_rev[rownames(n_deg_sigtwas_intersect)], tissue_abbr_rev[rownames(n_deg_sigtwas_intersect)])
      
      #get lists of genes
      trait_x_tissue_data <- lapply(tissues, function(twas_tiss){
        cat(paste0(twas_tiss, " "))
        twas_sub <- some.twas[some.twas$tissue == twas_tiss,]
        twas_sub$gene <- gsub(twas_sub$gene, pattern = "\\..*", replacement = "")
        twas_sub_hits <- twas_sub[twas_sub$adj_pvalue < adj_pvalue_alpha,]
        na.remove <- function(x) x[!is.na(x)]
        twas_genes_x_trait <- lapply(lapply(split(twas_sub_hits$gene, twas_sub_hits$trait), na.remove), unique)
        all_genes_x_trait <- lapply(lapply(split(twas_sub$gene, twas_sub$trait), na.remove), unique)
        list(hits = twas_genes_x_trait, all = all_genes_x_trait)
      })
      
      #reorder lists by trait
      traits <- setNames(colnames(n_deg_sigtwas_intersect), colnames(n_deg_sigtwas_intersect))
      trait_hits <- lapply(traits, function(trait_i) lapply(tissues, function(twas_tiss) trait_x_tissue_data[[twas_tiss]]$hits[[trait_i]]))
      trait_all <- lapply(traits, function(trait_i) lapply(tissues, function(twas_tiss) trait_x_tissue_data[[twas_tiss]]$all[[trait_i]]))
      
      trait_corr_mat_multi <- estimate_correlations_multi(y = trait_hits, t = trait_all, print_progress = T)
      
      save(trait_corr_mat_multi, file = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/trait_corr_mat_multi", 
                                               ifelse(use_random_Pred_genes | use_random_DE_genes, "_random-genes", ""), ifelse(use_all_genes_for_trait_tiss_corrmats, "_all-genes", "_intersecting-genes"), ".RData"))  
    } else {
      load(file = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/trait_corr_mat_multi", 
                         ifelse(use_random_Pred_genes | use_random_DE_genes, "_random-genes", ""), 
                         ifelse(use_all_genes_for_trait_tiss_corrmats, "", "_intersecting-genes"), 
                         ".RData"))
    }
    
    trait_corr_mat <- trait_corr_mat_multi
    # pheatmap::pheatmap(trait_corr_mat, breaks = -10:10/10, color = colorspace::diverging_hcl(21))
    
  }
}

if(estimate_interaction_corrmat){
  adj_pvalue_alpha <- 0.05
  tissues <- setNames(tissue_abbr_rev[rownames(n_deg_sigtwas_intersect)], tissue_abbr_rev[rownames(n_deg_sigtwas_intersect)])
  
  #get lists of genes
  trait_x_tissue_data <- lapply(tissues, function(twas_tiss){
    cat(paste0(twas_tiss, " "))
    de_genes <- DE_genes_x_tissue[[tissue_abbr[twas_tiss]]]
    twas_sub <- twas_by_tissue[[tissue_abbr[twas_tiss]]] #twas_by_tissue already subsetted down to mutually possible hits
    twas_sub$gene <- gsub(twas_sub$gene, pattern = "\\..*", replacement = "")
    twas_x_DE_hits <- twas_sub[twas_sub$adj_pvalue < adj_pvalue_alpha & twas_sub$gene %in% de_genes,]
    na.remove <- function(x) x[!is.na(x)]
    pairwise_hits_tiss_x_trait <- lapply(lapply(split(twas_x_DE_hits$gene, twas_x_DE_hits$trait), na.remove), unique)
    pairwise_all_tiss_x_trait <- lapply(lapply(split(twas_sub$gene, twas_sub$trait), na.remove), unique)
    list(hits = pairwise_hits_tiss_x_trait, all = pairwise_all_tiss_x_trait)
  })
  
  #reorder lists by trait
  #specifically order as traits kron tissues, 
  #so traits on major blocks containing tissues, 15 major blocks of 114 each
  traits <- setNames(colnames(n_deg_sigtwas_intersect), colnames(n_deg_sigtwas_intersect))
  hits_pairwise <- lapply(traits, function(trait_i) lapply(tissues, function(twas_tiss) 
    trait_x_tissue_data[[twas_tiss]]$hits[[trait_i]]))
  hits_pairwise <- unlist(hits_pairwise, recursive = F)
  
  all_pairwise <- lapply(traits, function(trait_i) lapply(tissues, function(twas_tiss) 
    trait_x_tissue_data[[twas_tiss]]$all[[trait_i]]))
  all_pairwise <- unlist(all_pairwise, recursive = F)
  
  #compute and save
  interaction_corr_mat <- estimate_correlations(y = hits_pairwise, t = all_pairwise, print_progress = T, ncores = 4)
  save(interaction_corr_mat, file = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/interaction_corr_mat", 
                                           ifelse(use_random_Pred_genes | use_random_DE_genes, "_random-genes", ""), ".RData"))  
} else if(!use_kronecker_interactions) {
  load(file = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/interaction_corr_mat", 
                     ifelse(use_random_Pred_genes | use_random_DE_genes, "_random-genes", ""), ".RData"))
}

#ok, here we can estimate the trait and tissue corrmats
#but not just looking at DE genes or TWAS hits, but both at the same time
#it should use the same input as the interaction corrmat above, 
#but the function estimate_correlations_multi(), like for the trait_corr_mat_multi
compute_intersecting_gene_trait_tiss_corrmats <- T
if(!use_all_genes_for_trait_tiss_corrmats){
  if(compute_intersecting_gene_trait_tiss_corrmats){
    adj_pvalue_alpha <- 0.05
    tissues <- setNames(tissue_abbr_rev[rownames(n_deg_sigtwas_intersect)], tissue_abbr_rev[rownames(n_deg_sigtwas_intersect)])
    traits <- setNames(colnames(n_deg_sigtwas_intersect), colnames(n_deg_sigtwas_intersect))
    
    #get lists of genes
    trait_x_tissue_data <- lapply(tissues, function(twas_tiss){
      cat(paste0(twas_tiss, " "))
      de_genes <- DE_genes_x_tissue[[tissue_abbr[twas_tiss]]]
      twas_sub <- twas_by_tissue[[tissue_abbr[twas_tiss]]] #twas_by_tissue already subsetted down to mutually possible hits
      twas_sub$gene <- gsub(twas_sub$gene, pattern = "\\..*", replacement = "")
      twas_x_DE_hits <- twas_sub[twas_sub$adj_pvalue < adj_pvalue_alpha & twas_sub$gene %in% de_genes,]
      na.remove <- function(x) x[!is.na(x)]
      pairwise_hits_tiss_x_trait <- lapply(lapply(split(twas_x_DE_hits$gene, twas_x_DE_hits$trait), na.remove), unique)
      pairwise_all_tiss_x_trait <- lapply(lapply(split(twas_sub$gene, twas_sub$trait), na.remove), unique)
      list(hits = pairwise_hits_tiss_x_trait, all = pairwise_all_tiss_x_trait)
    })
    
    #compute and save
    
    #traits
    trait_hits <- lapply(traits, function(trait_i) lapply(tissues, function(tissue_i) 
      trait_x_tissue_data[[tissue_i]]$hits[[trait_i]]))
    trait_all <- lapply(traits, function(trait_i) lapply(tissues, function(tissue_i) 
      trait_x_tissue_data[[tissue_i]]$all[[trait_i]]))
    trait_corr_mat <- estimate_correlations_multi(y = trait_hits, t = trait_all, print_progress = T)
    save(trait_corr_mat, file = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/trait_corr_mat", 
                                       ifelse(use_random_Pred_genes | use_random_DE_genes, "_random-genes", ""), 
                                       ifelse(use_all_genes_for_trait_tiss_corrmats, "_all-genes", "_intersecting-genes"), 
                                       ".RData"))
    
    #tissues
    tissue_hits <- lapply(tissues, function(tissue_i) lapply(traits, function(trait_i) 
      trait_x_tissue_data[[tissue_i]]$hits[[trait_i]]))
    tissue_all <- lapply(tissues, function(tissue_i) lapply(traits, function(trait_i) 
      trait_x_tissue_data[[tissue_i]]$all[[trait_i]]))
    tissue_corr_mat <- estimate_correlations_multi(y = tissue_hits, t = tissue_all, print_progress = T)
    rownames(tissue_corr_mat) <- colnames(tissue_corr_mat) <- tissue_abbr[rownames(tissue_corr_mat)]
    save(tissue_corr_mat, file = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/tissue_corr_mat", 
                                        ifelse(use_random_Pred_genes | use_random_DE_genes, "_random-genes", ""), 
                                        ifelse(use_all_genes_for_trait_tiss_corrmats, "_all-genes", "_intersecting-genes"), 
                                        ".RData"))
  } else {
    load(file = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/trait_corr_mat", 
                       ifelse(use_random_Pred_genes | use_random_DE_genes, "_random-genes", ""), 
                       ifelse(use_all_genes_for_trait_tiss_corrmats, "", "_intersecting-genes"), 
                       ".RData"))
    load(file = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/tissue_corr_mat", 
                       ifelse(use_random_Pred_genes | use_random_DE_genes, "_random-genes", ""), 
                       ifelse(use_all_genes_for_trait_tiss_corrmats, "", "_intersecting-genes"), 
                       ".RData"))
  }
}


#### obtain a bayesian estimate of the twas proportion and perform an 'enrichment' analysis ####

#get total number of possible hits
tissues <- setNames(rownames(n_deg_sigtwas_intersect), rownames(n_deg_sigtwas_intersect))
# total_number_of_possible_hits <- length(intersect(all_orthologs_tested, all_twas_genes_tested))

# load('/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/transcript_rna_seq_20211008.RData')
# rna_dea <- transcript_rna_seq$timewise_dea
# rm(transcript_rna_seq)
# rna_dea <- rna_dea[rna_dea$comparison_group == "8w",]
rna_dea <- MotrpacRatTraining6mo::combine_da_results(assays="TRNSCRPT")
rna_dea <- rna_dea[rna_dea$comparison_group == "8w",]
rna_dea$tissue_abbreviation <- rna_dea$tissue
rna_dea_genes <- lapply(setNames(tissues, tissues), function(tiss){
  genesymbols <- unique(gene_map$HUMAN_ORTHOLOG_SYMBOL[match(unique(rna_dea$feature_ID[rna_dea$tissue_abbreviation == tiss]),
                                                             gene_map$RAT_ENSEMBL_ID)])
  genesymbols <- genesymbols[!is.na(genesymbols)]
})

total_number_of_possible_hits_matrix <- t(sapply(setNames(tissues, tissues), function(tiss){
  print(tiss)
  sapply(setNames(salient_twas, salient_twas), function(trait_i){
    twas_genes_tested <- twas_by_tissue[[tiss]]$gene_name[twas_by_tissue[[tiss]]$trait == trait_i]
    motrpac_genes_tested <- rna_dea_genes[[tiss]]
    length(intersect(twas_genes_tested, motrpac_genes_tested))
  })
}))
total_number_of_possible_hits_matrix <- total_number_of_possible_hits_matrix[rownames(n_deg_sigtwas_intersect), 
                                                                             colnames(n_deg_sigtwas_intersect)]

#compile basic dataframe
data_subset <- data.frame(count = as.integer(unlist(c(n_deg_sigtwas_intersect))))
data_subset$tissue <- rep(rownames(n_deg_sigtwas_intersect), ncol(n_deg_sigtwas_intersect))
data_subset$trait <- unlist(lapply(colnames(n_deg_sigtwas_intersect), function(tri) rep(tri, nrow(n_deg_sigtwas_intersect))))
data_subset$TWAS_Hit <- "YES"
data_subset <- data_subset[data_subset$tissue != "HYPOTH",] #remove bc no DEGs observed (sampling from the prior is less efficient)
traits <- unique(data_subset$trait)
tissues <- tissues_intersect.model <- unique(data_subset$tissue)
trait_cats <- salient.categories
d <- list(cell_count = data_subset$count,
          total = sapply(1:nrow(data_subset), function(i) total_number_of_possible_hits_matrix[data_subset$tissue[i], data_subset$trait[i]]),
          row_count = sapply(1:nrow(data_subset), function(i) n_genes_in_nodes_matrix[data_subset$tissue[i], data_subset$trait[i]]),
          col_count = sapply(1:nrow(data_subset), function(i) sig_twas_by_trait_genes_matrix[data_subset$tissue[i], data_subset$trait[i]]),
          row_index = match(data_subset$tissue, tissues),
          col_index = match(data_subset$trait, traits),
          colcat_index = match(traitwise_partitions$Category[match(traits, traitwise_partitions$Tag)] , trait_cats),
          row_n = length(tissues),
          col_n = length(traits),
          colcat_n = length(trait_cats))

#some quick EDA
par(mfrow = c(1,1))
plot(d$cell_count / d$row_count, (d$col_count - d$cell_count) / (d$total - d$row_count), pch = 19, col = adjustcolor(1,0.5),
     cex = sqrt(d$row_count / max(d$row_count))); abline(0,1,lwd=2,lty=2,col=2)
plot(log10(d$cell_count / d$row_count), log10((d$col_count - d$cell_count) / (d$total - d$row_count)), pch = 19, col = adjustcolor(1,0.5),
     cex = sqrt(d$row_count / max(d$row_count))); abline(0,1,lwd=2,lty=2,col=2)
data_subset[order(d$cell_count / d$row_count, decreasing = T),]


# base = paste0("deviation_from_expected_logodds_split_the_difference_2", ifelse(use_random_DE_genes, "_randomgenes", ""))
# stan_program <- '

# '

# base = paste0("deviation_from_expected_logodds_split_the_difference_MVN_priors", ifelse(use_random_DE_genes, "_randomgenes", ""))
# stan_program <- '

# '

#data goes tiss1-trait1, tiss2-trait1, tiss3-trait1...

#more elaborate dataset accommodating non-independence among traits and tissues
if(use_kronecker_interactions){
  full_interaction_mat <- kronecker(trait_corr_mat[traits,traits], tissue_corr_mat[tissues, tissues])
  L_interaction <- t(chol(full_interaction_mat))
} else {
  colrow_names <- paste0(rep(traits, each = length(tissues)), ".", rep(tissue_abbr_rev[tissues], times = length(traits)))
  #right order, has corr = approx. 0.35ish with above
  sub_interaction_corr_mat <- interaction_corr_mat[colrow_names, colrow_names]
  L_interaction <- t(chol(sub_interaction_corr_mat))
}
d <- list(cell_count = data_subset$count,
          total = sapply(1:nrow(data_subset), function(i) 
            total_number_of_possible_hits_matrix[data_subset$tissue[i], data_subset$trait[i]]),
          row_count = sapply(1:nrow(data_subset), function(i) 
            n_genes_in_nodes_matrix[data_subset$tissue[i], data_subset$trait[i]]),
          col_count = sapply(1:nrow(data_subset), function(i) 
            sig_twas_by_trait_genes_matrix[data_subset$tissue[i], data_subset$trait[i]]),
          row_index = match(data_subset$tissue, tissues),
          col_index = match(data_subset$trait, traits),
          colcat_index = match(traitwise_partitions$Category[match(traits, traitwise_partitions$Tag)], trait_cats),
          row_n = length(tissues),
          col_n = length(traits),
          colcat_n = length(trait_cats),
          L_interaction = L_interaction,
          L_col = t(chol(trait_corr_mat[traits,traits])),
          L_row = t(chol(tissue_corr_mat[tissues, tissues]))
)

stan_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/Stan_models/"
#this incorporates pairwise interaction term
#because when pooling across traits (for trait cats), we want to allow independence
#and when pooling within traits (across tissues) we do too

base = paste0("deviation_from_expected_logodds_split_the_difference_informative-MVN-matrix-normal_priors") 

stan_program <- paste0(readLines(paste0(stan_dir, base, ".stan")), collapse = "\n")

# base = paste0("deviation_from_expected_logodds_split_the_difference", ifelse(use_random_DE_genes, "_randomgenes", ""))

#try the even simpler model out again
# data_subset <- data1
# data_subset <- data_subset[data_subset$tissue != "HYPOTH",]
# #this is cos Stan can't initialize when sampling 0 ~ binomial(x, p = 0)
# data_subset <- data_subset[sapply(1:nrow(data_subset), function(i) sig_twas_by_trait_genes_matrix[data_subset$tissue[i], data_subset$trait[i]]) != 0,]
# 
# traits <- unique(data_subset$trait)
# tissues <- tissues_intersect.model <- unique(data_subset$tissue)
# trait_cats <- salient.categories
# d <- list(cell_count = data_subset$count,
#           total = sapply(1:nrow(data_subset), function(i) total_number_of_possible_hits_matrix[data_subset$tissue[i], data_subset$trait[i]]),
#           row_count = sapply(1:nrow(data_subset), function(i) n_genes_in_nodes_matrix[data_subset$tissue[i], data_subset$trait[i]]),
#           col_count = sapply(1:nrow(data_subset), function(i) sig_twas_by_trait_genes_matrix[data_subset$tissue[i], data_subset$trait[i]]),
#           row_index = match(data_subset$tissue, tissues),
#           col_index = match(data_subset$trait, traits),
#           colcat_index = match(traitwise_partitions$Category[match(traits, traitwise_partitions$Tag)] , trait_cats),
#           row_n = length(tissues),
#           col_n = length(traits),
#           colcat_n = length(trait_cats),
#           overall_prob = sapply(1:nrow(data_subset), function(i) sig_twas_by_trait_genes_matrix[data_subset$tissue[i], data_subset$trait[i]]) /
#             sapply(1:nrow(data_subset), function(i) total_number_of_possible_hits_matrix[data_subset$tissue[i], data_subset$trait[i]]),
#           n = nrow(data_subset)
#           )
# 
# #sanity check
# all(d$cell_count[d$col_count / d$total < 1E-6] == 0)
# plot(d$cell_count / d$row_count, d$overall_prob)
# plot(logit(d$overall_prob), logit(d$cell_count / d$row_count))
# abline(0,1)
# hist(logit(d$cell_count / d$row_count) - logit(d$overall_prob))
# a <- logit(d$cell_count / d$row_count) - logit(d$overall_prob)
# mean(a[a!=Inf & a!=-Inf])
# 

# base = paste0("deviation_from_expected_logodds_fixed", ifelse(use_random_DE_genes, "_randomgenes", ""))
# stan_program <- '

# '

print(use_random_DE_genes)
print(use_random_Pred_genes)
print(use_kronecker_interactions)
print(estimate_interaction_corrmat)
print(use_all_genes_for_trait_tiss_corrmats)

if(fit_model){
  
  #compile model
  if(!exists("curr_stan_program") || stan_program != curr_stan_program){
    curr_stan_program <- stan_program
    f <- write_stan_file(stan_program)
  }
  mod <- cmdstan_model(f)
  
  #write model
  write_stan_file(stan_program, dir = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/mcmc_output/", 
                  basename = paste0(base, ifelse(use_all_cats, "_allCats", "")))
  write_stan_json(d, paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/mcmc_output/", 
                            paste0(base, ifelse(use_all_cats, "_allCats", "")), 
                            ifelse(use_random_DE_genes, "_randomgenes", ""),
                            ifelse(use_all_genes_for_trait_tiss_corrmats, "", "_intersecting-genes-marginal-corrmats"),
                            ifelse(!use_kronecker_interactions, "_pairwise-interactions", ""),
                            ".json"))
  
  #fit model
  out <- mod$sample(chains = 4, iter_sampling = 2.5E3, iter_warmup = 2.5E3, data = d, parallel_chains = 4, 
                    adapt_delta = 0.95, refresh = 50, init = 0.1, max_treedepth = 15, thin = 1)
  summ <- out$summary()
  summ[order(summ$ess_bulk),]
  summ[order(summ$rhat, decreasing = T),]
  save(out, file = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/mcmc_output/", 
                          base, 
                          ifelse(use_all_cats, "_allCats"), 
                          ifelse(use_random_DE_genes, "_randomgenes", ""),
                          ifelse(use_all_genes_for_trait_tiss_corrmats, "", "_intersecting-genes-marginal-corrmats"),
                          ifelse(!use_kronecker_interactions, "_pairwise-interactions", ""),
                          ".cmdStanR.fit"))
  save(summ, file = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/mcmc_output/", 
                           base, 
                           ifelse(use_all_cats, "_allCats"), 
                           ifelse(use_random_DE_genes, "_randomgenes", ""),
                           ifelse(use_all_genes_for_trait_tiss_corrmats, "", "_intersecting-genes-marginal-corrmats"),
                           ifelse(!use_kronecker_interactions, "_pairwise-interactions", ""),
                           ".cmdStanR.summ"))
} else {
  load(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/mcmc_output/", 
              base, 
              ifelse(use_all_cats, "_allCats"), 
              ifelse(use_random_DE_genes, "_randomgenes", ""),
              ifelse(use_all_genes_for_trait_tiss_corrmats, "", "_intersecting-genes-marginal-corrmats"),
              ifelse(!use_kronecker_interactions, "_pairwise-interactions", ""),
              ".cmdStanR.fit"))
}

#### examine MCMC output ####
samps <- data.frame(as_draws_df(out$draws()))

dev.off()
hist(samps$overall_bias, freq = F, xlab = "overall difference (shared across all counts)",
     main = latex2exp::TeX(paste0("posterior distribuion of overall difference term $\\alpha$ (", ifelse(use_random_DE_genes, "random", "observed"), " genesets)")))
mean(samps$overall_bias > 0)

prop_greater_than_0 <- function(x) mean(x>0)

cellbias <- apply(subset_samps("cell_bias", c("raw", "sd", "L"), samps = samps), 2, prop_greater_than_0)
sum(cellbias > 0.975)
sum(cellbias < 0.025)
hist(cellbias, breaks = 100)

celltotalbias <- apply(subset_samps("cell_total_prob_bias", c("raw", "sd", "logodds", "L"), samps = samps), 2, prop_greater_than_0)
sum(celltotalbias > 0.975)
sum(celltotalbias < 0.025)
hist(celltotalbias, breaks = 100, main = "Tissue x Trait Posterior Probability of Positive Enrichment",
     xlab = "Proportion Posterior Mass > 0")
head(cbind(data_subset$tissue, data_subset$trait, celltotalbias)[
  celltotalbias > 0.95 | celltotalbias < 0.05,], 20)

rowbias <- apply(subset_samps("row_bias", c("raw", "sd", "L"), samps = samps) + 
                   samps$overall_bias, 2, prop_greater_than_0)
sum(rowbias > 0.9)
sum(rowbias < 0.1)
cbind(tissues[order(rowbias)], sort(rowbias))

colbias <- apply(subset_samps("col_bias", c("raw", "sd", "L"), samps = samps) + 
                   samps$overall_bias +
                   subset_samps("colcat_bias", c("raw", "sd", "L"), samps = samps)[,
                                                                                   match(traitwise_partitions$Category[
                                                                                     match(traits, traitwise_partitions$Tag)], 
                                                                                     trait_cats)], 
                 2, prop_greater_than_0)
# colbias <- apply(subset_samps("col_bias", c("raw", "sd", "L"), samps = samps), 2, prop_greater_than_0)
hist(colbias)
sum(colbias > 0.9)
sum(colbias < 0.1)
trait_key <- setNames(trait_categories$new_Phenotype, trait_categories$Tag)
cbind(trait_key[traits[order(colbias)]], sort(colbias))

# colcatbias <- apply(subset_samps("colcat_bias", c("raw", "sd", "L"), samps = samps), 2, prop_greater_than_0)
colcatbias <- apply(subset_samps("colcat_bias", c("raw", "sd", "L"), samps = samps) + samps$overall_bias, 
                    2, prop_greater_than_0)
hist(colcatbias)
sum(colcatbias > 0.9)
sum(colcatbias < 0.1)
cbind(trait_cats[order(colcatbias)], sort(colcatbias))

#process correlation matrices
if(any(grepl("L_", colnames(samps)))){
  #columns
  L_col_bias_samps <- subset_samps("L_col_bias", c("raw"), samps = samps)
  col_bias_corrs_samps <- do.call(abind::abind, list(lapply(1:nrow(samps), function(i){L <- matrix(unlist(L_col_bias_samps[i,]), length(traits), length(traits)); L %*% t(L)}), along = 3))
  col_bias_corrs_mean <- apply(col_bias_corrs_samps, c(1,2), mean)
  col_bias_corrs_gr0.5 <- apply(col_bias_corrs_samps, c(1,2), prop_greater_than_0)
  hist(col_bias_corrs_gr0.5[upper.tri(col_bias_corrs_gr0.5)])
  
  #rows
  L_row_bias_samps <- subset_samps("L_row_bias", c("raw"), samps = samps)
  row_bias_corrs_samps <- do.call(abind::abind, list(lapply(1:nrow(samps), function(i){L <- matrix(unlist(L_row_bias_samps[i,]), length(tissues), length(tissues)); L %*% t(L)}), along = 3))
  row_bias_corrs_mean <- apply(row_bias_corrs_samps, c(1,2), mean)
  row_bias_corrs_gr0.5 <- apply(row_bias_corrs_samps, c(1,2), prop_greater_than_0)
  # hist(row_bias_corrs_gr0.5[upper.tri(row_bias_corrs_gr0.5)])
  # hist(row_bias_corrs_mean[upper.tri(row_bias_corrs_mean)])
}

# rowcatbias <- apply(subset_samps("rowcat_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)
# sum(rowcatbias > 0.95)
# sum(rowcatbias < 0.05)
# cbind(tissue_cats[order(rowcatbias)], sort(rowcatbias))


#generate a "significance" matrix
signif_threshold <- 0.05
signif_df <- data_subset[1:(nrow(data_subset)), c("tissue", "trait")]
signif_df$prob_diff_is_positive <- apply(subset_samps("cell_total_prob_bias", c("raw", "sd", "logodds"), samps = samps), 2, prop_greater_than_0)
signif_df$signif <- 0
signif_df$signif[signif_df$prob_diff_is_positive > (1-signif_threshold)] <- 1
signif_df$signif[signif_df$prob_diff_is_positive < signif_threshold] <- -1

signif_df[signif_df$signif != 0,]

signif_matrix <- reshape(signif_df[,-match("prob_diff_is_positive", colnames(signif_df))], idvar = "tissue", timevar = "trait", direction = "wide")
rownames(signif_matrix) <- signif_matrix$tissue
signif_matrix <- signif_matrix[,-match("tissue", colnames(signif_matrix))]
colnames(signif_matrix) <- gsub(colnames(signif_matrix), pattern = "signif.", replacement = "")

#add back in the HYPOTHALAMUS if we removed it
if(nrow(signif_matrix) != nrow(n_deg_sigtwas_intersect)){
  signif_matrix <- signif_matrix[rownames(n_deg_sigtwas_intersect),]
  rownames(signif_matrix) <- rownames(n_deg_sigtwas_intersect)
  signif_matrix[is.na(signif_matrix)] <- 0
}

if(!all(colnames(signif_matrix) == colnames(n_deg_sigtwas_intersect)) & all(rownames(signif_matrix) == rownames(n_deg_sigtwas_intersect))){
  stop("something's wrong with the significance matrix")
}
