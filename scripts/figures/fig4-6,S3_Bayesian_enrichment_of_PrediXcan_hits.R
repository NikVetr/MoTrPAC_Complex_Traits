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
  
  pheatmap::pheatmap(tissue_corr_mat, breaks = -10:10/10, color = colorspace::diverging_hcl(21))
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
twas_by_tissue <- mclapply(tissues, function(tiss_i){
  do.call(rbind, lapply(salient_twas, function(trait_i){
    system(sprintf('echo "%s "', paste0(tiss_i, collapse="")))
    some.twas[some.twas$tissue == TISSUE_ABBREV_TO_CODE[tiss_i] & 
                some.twas$trait == trait_i & 
                some.twas$gene_name %in% possible_genes_tiss_x_trait[[tiss_i]][[trait_i]],]
  }))
}, mc.cores = 8)
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
    pheatmap::pheatmap(trait_corr_mat, breaks = -10:10/10, color = colorspace::diverging_hcl(21))
    
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
  hist(row_bias_corrs_gr0.5[upper.tri(row_bias_corrs_gr0.5)])
  hist(row_bias_corrs_mean[upper.tri(row_bias_corrs_mean)])
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

#### figure 4 ####
plot_figure_4 <- T
if(plot_figure_4){

#graphical parameters
incl_significance <- T
incl_cell_totals <- F
trait_category_legend_below = T
use_tissue_cols_for_cols <- F
opacity_power_scaler <- 0.25
opacity_white_threshold <- 0.85
use_counts <- T
prop_TWAS <- T
order_by_counts <- F
order_by_posterior_enrichment <- T
group_by_tissue_type <- T
use_range_for_maginal_labels <- F
subset_to_traits <- T
focal_traitcats <- c("Cardiometabolic", "Immune", "Endocrine system", "Anthropometric", "Allergy", "Aging")

#for grouping by tissues
tissue_cats <- list(circulation = c("BLOOD", "HEART", "SPLEEN"),
                    skeletal_muscle = c("SKM-GN", "SKM-VL"),
                    other = rev(c("ADRNL", "KIDNEY", "LUNG", "LIVER")),
                    adipose = c("WAT-SC"),
                    brain = c("CORTEX", "HYPOTH", "HIPPOC"),
                    GI = c("SMLINT", "COLON"))
tissue_cats <- rev(tissue_cats)
disp_amount <- 0.5

#some easy graphical params
cols = list(Tissue=MotrpacRatTraining6moData::TISSUE_COLORS[names(motrpac_gtex_map)], 
            Time=MotrpacRatTraining6moData::GROUP_COLORS[paste0(c(1,2,4,8), "w")],
            Sex=MotrpacRatTraining6moData::SEX_COLORS[c('male','female')])
mars <- list(c(6,
               0 + ifelse(subset_to_traits, 2, 0),
               6 + ifelse(group_by_tissue_type, disp_amount * 4.5, 0),
               6.5 + ifelse(subset_to_traits, 1, 0)),
             c(3.5,7,3,7)+0.5, 
             c(3.5,1,3,14)+0.5,
             c(4.25,2.5,2.5,2.5), 
             c(4.25,2.5,2.5,2.5), 
             c(10,6.5,1.5,2.5))

#preliminary processing for additional graphical params
if(subset_to_traits){
  trait_subset <- colnames(n_deg_sigtwas_intersect)[
    traitwise_partitions$Category[match(colnames(n_deg_sigtwas_intersect), traitwise_partitions$Tag)] %in% 
      focal_traitcats]
  categories_represented <- unique(traitwise_partitions$Category[match(trait_subset, traitwise_partitions$Tag)])
} else {
  trait_subset <- colnames(n_deg_sigtwas_intersect)
}
categories_represented <- unique(traitwise_partitions$Category[match(trait_subset, traitwise_partitions$Tag)])

if(use_counts){
  table_to_use <- n_deg_sigtwas_intersect  
} else {
  if(prop_TWAS){
    table_to_use <- round(prop_twas_are_degs * 1000)  
  } else{
    table_to_use <- round(prop_degs_are_twas * 1000)  
  }
}

if(order_by_counts){
  table_to_use <- table_to_use[,colnames(n_deg_sigtwas_intersect)]
  signif_matrix_to_use <- signif_matrix[,colnames(n_deg_sigtwas_intersect)]
} else if(order_by_posterior_enrichment) {
  trait_order <- traits[order(colbias, decreasing = T)]
  trait_order <- trait_order[trait_order %in% colnames(table_to_use)]
  table_to_use <- table_to_use[,trait_order]
  signif_matrix_to_use <- signif_matrix[,match(trait_order, colnames(signif_matrix))]
} else {
  table_to_use <- table_to_use[,colnames(prop_twas_are_degs)]
  signif_matrix_to_use <- signif_matrix[,colnames(prop_twas_are_degs)]
}

if(subset_to_traits){
  table_to_use <- table_to_use[,trait_subset]
  signif_matrix_to_use <- signif_matrix_to_use[,trait_subset]
}


if(group_by_tissue_type){
  tissue_disps <- unlist(lapply(1:length(tissue_cats), function(tci) rep(disp_amount * (tci), length(tissue_cats[[tci]]))))
  tissue_cats_bars_ylocs <- cbind(start = (c(0, cumsum(unlist(lapply(tissue_cats, function(tc) length(tc))) + disp_amount)) + 2 * disp_amount)[-(length(tissue_cats)+1)], 
                                  end = cumsum(unlist(lapply(tissue_cats, function(tc) length(tc))) + disp_amount) + disp_amount)
  
  #get in correct order
  table_to_use <- table_to_use[unlist(tissue_cats),]
  signif_matrix_to_use <- signif_matrix_to_use[unlist(tissue_cats),]
  
}

if(use_range_for_maginal_labels){
  n_genes_in_nodes_label <- apply(apply(n_genes_in_nodes_matrix, 1, range), 2, paste0, collapse = " - ")
  sig_twas_by_trait_genes_label <- apply(apply(sig_twas_by_trait_genes_matrix, 2, range), 2, paste0, collapse = " - ")
} else {
  n_genes_in_nodes_label <- apply(n_genes_in_nodes_matrix, 1, max)
  sig_twas_by_trait_genes_label <- apply(sig_twas_by_trait_genes_matrix, 2, max)
}

if(length(salient.categories) < 9){
  category_colors <- RColorBrewer::brewer.pal(length(salient.categories), "Dark2")
  names(category_colors) <- sort(salient.categories)  
} else {
  category_colors <- RColorBrewer::brewer.pal(5, "Dark2")
  names(category_colors) <- sort(c("Cardiometabolic", "Aging", "Anthropometric", 
                                   "Immune", "Psychiatric-neurologic"))
  cats_leftover <- setdiff(salient.categories, names(category_colors))
  cols_leftover <- setdiff(RColorBrewer::brewer.pal(8, "Dark2"), category_colors)
  n_cols_needed <- length(cats_leftover) - length(cols_leftover)
  if(n_cols_needed > 0){
    extra_cols <- c(cols_leftover, RColorBrewer::brewer.pal(n_cols_needed, "Set2"))
  } else {
    extra_cols <- cols_leftover
  }
  names(extra_cols) <- cats_leftover
  category_colors <- c(category_colors, extra_cols)
}
cols$category <- category_colors

if(use_tissue_cols_for_cols){
  heatmap_cols <- sapply((1:max(table_to_use) / max(table_to_use))^opacity_power_scaler, function(opcty) 
    adjustcolor("black", opcty))
} else {
  heatmap_cols <- viridisLite::viridis(n = max(table_to_use, na.rm = T)*100+1)
  heatmap_cols <- heatmap_cols[round(log(1:max(table_to_use, na.rm = T)) / log(max(table_to_use, na.rm = T)) * max(table_to_use, na.rm = T) * 100 + 1)]
}

#load appropriate model fit
if(base != "deviation_from_expected_logodds_split_the_difference_informative-MVN-matrix-normal_priors"){
  
  base = "deviation_from_expected_logodds_split_the_difference_informative-MVN-matrix-normal_priors"
  load(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/mcmc_output/", 
              base, 
              ifelse(use_all_cats, "_allCats"), 
              ifelse(use_random_DE_genes, "_randomgenes", ""),
              ifelse(!use_kronecker_interactions, "_pairwise-interactions", ""),
              ".cmdStanR.fit"))
  
  samps <- data.frame(as_draws_df(out$draws()))
  prop_greater_than_0 <- function(x) mean(x>0)
  cellbias <- apply(subset_samps("cell_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)
  celltotalbias <- apply(subset_samps("cell_total_prob_bias", c("raw", "sd", "logodds"), samps = samps), 2, prop_greater_than_0)
  rowbias <- apply(subset_samps("row_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)
  colbias <- apply(subset_samps("col_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)
  trait_key <- setNames(trait_categories$new_Phenotype, trait_categories$Tag)
  colcatbias <- apply(subset_samps("colcat_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)
  
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
}

#start the plotting!
cairo_pdf(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/figures/fig4_intersect-enrichment_redux.pdf"), 
          width = 2000 / 72 * ncol(table_to_use) / 80, 
          height = 1200 / 72 + ifelse(group_by_tissue_type, disp_amount * 0.75, 0), 
          family="Arial Unicode MS", pointsize = 18.5)
par(xpd = T, 
    mar = mars[[1]])
layout(rbind(
  c(1,1,1,1,4),
  c(2,2,3,3,5),
  c(6,6,6,6,6)
), heights = c(1,1,0.6))

#initialize blank plot
plot(1, xaxt="n",yaxt="n",bty="n",pch="",ylab="",xlab="", main="", sub="", 
     xlim= c(-5,ncol(table_to_use)), ylim = c(-5,nrow(table_to_use)))

if(group_by_tissue_type){
  tissue_cats_bars_xlocs <- sapply(tissue_cats, function(tc) max(strwidth(tc, units = "user"))) + 0.2
  tissue_cats_bars_xlocs <- rep(max(tissue_cats_bars_xlocs), length(tissue_cats_bars_xlocs))
}

for(ri in 1:nrow(table_to_use)){
  # text(x = 0.5, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0), pos = 2, labels = rownames(table_to_use)[ri], col = MotrpacRatTraining6moData::TISSUE_COLORS[rownames(table_to_use)[ri]])
  text(x = 0.5, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0), pos = 2, labels = rownames(table_to_use)[ri], col = 1)
  for(ci in 1:ncol(table_to_use)){
    if(ri == 1){
      #trait names
      text(x = ci+0.5, y = -0.9, pos = 2, srt = 45, cex = 0.85,
           labels = gsub("craneal", "cranial", gsub("_Scatter", "", trait_categories$new_Phenotype[match(colnames(table_to_use)[ci], 
                                                                              trait_categories$Tag)])))
      trait_category <- trait_categories$Category[match(colnames(table_to_use)[ci], trait_categories$Tag)]
      rect(xleft = ci + 1/2,
           xright = ci - 1/2,
           ybottom = ri - 1/2 - 1,
           ytop =  ri + 1/2 - 1,
           col = category_colors[trait_category])
      text(labels = substr(trait_category, 1, 2), y = ri - 1, x = ci,
           col = "white", cex = 0.8)
    }
    
    #vertical total # options
    if(ri == nrow(table_to_use)){
      text(x = ci-0.5, y = ri+1 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), pos = 4, srt = 45,
           labels = sig_twas_by_trait_genes_label[colnames(table_to_use)[ci]], cex = 0.9)
    }
    
    #horiz total # of options
    if(ci == ncol(table_to_use)){
      text(x = ci+0.45, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0), pos = 4,
           labels = n_genes_in_nodes_label[rownames(table_to_use)[ri]])
    }
    
    #the actual cells
    if(use_tissue_cols_for_cols){
      rect(xleft = ci + 1/2,
           xright = ci - 1/2,
           ybottom = ri - 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0),
           ytop =  ri + 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0),
           col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[rownames(table_to_use)[ri]], (table_to_use[ri, ci] / max(table_to_use, na.rm = T))^opacity_power_scaler ))  
    } else {
      rect(xleft = ci + 1/2,
           xright = ci - 1/2,
           ybottom = ri - 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0),
           ytop =  ri + 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0),
           col = heatmap_cols[table_to_use[ri, ci]])
    }
    
    
    
    # rect(xleft = ci + 1/2,
    #      xright = ci - 1/2,
    #      ybottom = ri - 0.475,
    #      ytop =  ri + 0.475,
    #      col = heatmap_cols[table_to_use[ri, ci]],
    #      border = MotrpacRatTraining6moData::TISSUE_COLORS[rownames(table_to_use)[ri]])
    
    #text inside of cells
    if(table_to_use[ri, ci] != 0){
      text_in_cell <- table_to_use[ri, ci]
      ndigs <- nchar(text_in_cell)
      if(incl_significance){
        signif_dir <- c("₋", "","⁺")[match(signif_matrix_to_use[ri, ci], -1:1)]
        # text_in_cell <- paste0(text_in_cell, signif_dir)
        
        #try using corner triangles
        ctr <- corner_triangle_ratio <- 1/3
        cxr = ci + 1/2
        cxl = ci - 1/2
        cyb = ri - 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0)
        cyt =  ri + 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0)
        if(signif_dir == "⁺"){
          polygon(x = c(cxl * (ctr) + cxr * (1-ctr), cxr, cxr),
                  y = c(cyt, cyt, cyt * (1-ctr) + cyb * ctr),
                  col = "red")  
        }
        if(signif_dir == "₋"){
          polygon(x = c(cxl * (ctr) + cxr * (1-ctr), cxr, cxr),
                  y = c(cyb, cyb, cyb * (1-ctr) + cyt * ctr),
                  col = "lightblue")  
        }
      }
      text(text_in_cell, x = ci, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0), 
           col = ifelse((table_to_use[ri, ci] / max(table_to_use, na.rm = T))^opacity_power_scaler > opacity_white_threshold, "black", "white"), 
           cex = 0.85 + ifelse(subset_to_traits, 0.1, 0) - (ndigs-1)/6)
      if(incl_cell_totals){
        text(sig_twas_by_trait_genes_matrix[,colnames(table_to_use)][ri, ci], 
             x = ci - 0.7375, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0) - 0.325, 
             col = "black", cex = 0.2, pos = 4, srt = 90)
        text(n_genes_in_nodes_matrix[,colnames(table_to_use)][ri, ci], 
             x = ci - 0.7, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0) - 0.4, 
             col = "black", cex = 0.2, pos = 4, srt = 0)
      }
    }
    
  }
}

#tissue category bars
if(group_by_tissue_type){
  for(bi in 1:length(tissue_cats_bars_xlocs)){
    segments(x0 = -tissue_cats_bars_xlocs[bi],
             x1 = -tissue_cats_bars_xlocs[bi],
             y0 = tissue_cats_bars_ylocs[bi,1],
             y1 = tissue_cats_bars_ylocs[bi,2],
             lwd = 3)
    bracket_length <- 0.2
    segments(x0 = -tissue_cats_bars_xlocs[bi],
             x1 = -tissue_cats_bars_xlocs[bi]+bracket_length,
             y0 = tissue_cats_bars_ylocs[bi,1],
             y1 = tissue_cats_bars_ylocs[bi,1],
             lwd = 3)
    segments(x0 = -tissue_cats_bars_xlocs[bi],
             x1 = -tissue_cats_bars_xlocs[bi]+bracket_length,
             y0 = tissue_cats_bars_ylocs[bi,2],
             y1 = tissue_cats_bars_ylocs[bi,2],
             lwd = 3)
    
    text(x = -tissue_cats_bars_xlocs[bi], y = mean(tissue_cats_bars_ylocs[bi,]) + ifelse(grepl("_", names(tissue_cats)[bi]), -0.25, 0), pos = 2, 
         labels = gsub("Gi", "GI", stringr::str_to_title(gsub("_", "\n", names(tissue_cats)[bi]))), cex = 1.25)
  }  
}

#legend for heatmap
x_adj <- 2.25
y_adj <- 0
yb_adjust <- ifelse(incl_significance, 3, 0)
n_legend_rects_to_use <- 30
n_legend_labels_to_use <- 10
legend_yvals <- round(seq(0, max(table_to_use), length.out = n_legend_labels_to_use))
legend_ylocs <- seq(yb_adjust, 1 + nrow(table_to_use) + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), length.out = n_legend_rects_to_use)
legend_ycols <- round(seq(1, max(table_to_use), length.out = n_legend_rects_to_use))
for(i in 1:(n_legend_rects_to_use-1)){
  yb = legend_ylocs[i]
  yt = legend_ylocs[i+1]
  # print(paste(c(yb,yt)))
  rect(xleft = ncol(table_to_use) + x_adj + 2 + 1/2,
       xright = ncol(table_to_use) + x_adj + 2 - 1/2,
       ybottom = yb,
       ytop =  yt,
       col = heatmap_cols[legend_ycols[i]], border = NA)
}
for(i in 1:n_legend_labels_to_use){
  text(labels = legend_yvals[i], x = ncol(table_to_use) + x_adj + 2.4, pos = 4, cex = 0.75,
       y = -0.25 + yb_adjust + legend_yvals[i] / max(table_to_use) * (1+nrow(table_to_use) - yb_adjust + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0)))
}

#overall rect for 0
rect(xleft = ncol(table_to_use) + x_adj + 2 + 1/2,
     xright = ncol(table_to_use) + x_adj + 2 - 1/2,
     ybottom = min(legend_ylocs) - diff(range(legend_ylocs))/50,
     ytop =  max(legend_ylocs))

#legend for significance
if(incl_significance){
  # text(labels = paste0("Pr(Enr.) > ", (1 - signif_threshold), " : X⁺\nPr(Dep.) > ", (1 - signif_threshold), " : X₋"),
  #      x = ncol(table_to_use) + x_adj - 1.75,
  #      y = yb_adjust-3.25, pos = 4, cex = 0.75)
  
  signif_label <- paste0("Pr(Dep.) > ", (1 - signif_threshold), " : ")
  ctr <- 0.5
  cxr = ncol(table_to_use) + x_adj - 2.25 + strwidth(signif_label) + 0.45
  cxl = ncol(table_to_use) + x_adj - 2.25 + strwidth(signif_label) - 0.45
  cyb = yb_adjust - 3.75 - 0.45
  cyt = yb_adjust - 3.75 + 0.45
  
  text(labels = signif_label,
       x = ncol(table_to_use) + x_adj - 1.75,
       y = cyb * 0.325 + cyt * 0.675, pos = 4, cex = 0.75)
  rect(xleft = cxl,
       xright = cxr,
       ybottom = cyb,
       ytop =  cyt)
  polygon(x = c(cxl * (ctr) + cxr * (1-ctr), cxr, cxr),
          y = c(cyb, cyb, cyb * (1-ctr) + cyt * ctr),
          col = "lightblue")  
  
  cyb <- cyb + 1.25
  cyt <- cyt + 1.25
  # cxl <- cxl - 0.125
  # cxr <- cxr - 0.125
  
  signif_label2 <- paste0("Pr(Enr.) > ", (1 - signif_threshold), " : ")
  text(labels = signif_label2,
       x = ncol(table_to_use) + x_adj - 1.75,
       y = cyb * 0.675 + cyt * 0.325, pos = 4, cex = 0.75 * strwidth(signif_label) / strwidth(signif_label2))
  rect(xleft = cxl,
       xright = cxr,
       ybottom = cyb,
       ytop =  cyt)
  polygon(x = c(cxl * (ctr) + cxr * (1-ctr), cxr, cxr),
          y = c(cyt, cyt, cyt * (1-ctr) + cyb * ctr),
          col = "red")  
  
} else {
  text(labels = paste0("* indicate \nIHW-adj. \np-vals < 0.05 \nfrom Fisher's Exact Test"),
       x = ncol(table_to_use) + x_adj - 1.75,
       y = yb_adjust-1.25, pos = 4, cex = 0.75, )
}



#legend for trait categories
if(trait_category_legend_below){
  x_adj2 <- 0
  y_adj2 <- -0.5
  for(i in 1:length(categories_represented)){
    rect(xleft = -1/2 + i + cumsum(c(0,strwidth(categories_represented, units = "user") + 1))[i] + x_adj2,
         xright = 1/2 + i + cumsum(c(0,strwidth(categories_represented, units = "user") + 1))[i] + x_adj2,
         ybottom = -11 + y_adj2,
         ytop = -10 + y_adj2,
         col = category_colors[categories_represented[i]], border = 1)
    text(labels = categories_represented[i], pos = 4, y = -10.5 + y_adj2, 
         x = x_adj2 + 0.35 + i + cumsum(c(0,strwidth(categories_represented, units = "user") + 1))[i])
    text(labels = substr(categories_represented, 1, 2)[i], y = -10.5 + y_adj2, 
         x = i + cumsum(c(0,strwidth(categories_represented, units = "user") + 1))[i] + x_adj2,
         col = "white", cex = 0.8)
    #draw circle?
    arc(t1 = 3*pi/2, t2 = pi/2, r1 = (10.5-y_adj2) / 2, r2 = (10.5-y_adj2) / 2, center = c(0,(-10.5 + y_adj2)/2), 
        lwd = 2, res = 100, adjx = ifelse(order_by_counts, 1, 1.1))
    points(0, 0, pch = -9658, cex = 2)
  }
} else {
  x_adj2 <- x_adj - 2.5
  y_adj2 <- y_adj - 2.5
  for(i in 1:length(category_colors)){
    rect(xleft = 0 + i,
         xright = 1 + i,
         ybottom = -10,
         ytop = -9,
         col = category_colors[i], border = 1)
    text(labels = names(category_colors)[i], pos = 4, y = i + y_adj2 - 7, x = ncol(table_to_use) - 1/2 + x_adj2 + 1/3)
  }
}

#labels
#horiz label for total
segments(x0 = -2, x1 = 0.5, y0 = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0),
         y1 = nrow(table_to_use) + 1 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), lwd = 2)
segments(x0 = -2, x1 = 0.5, y0 = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), 
         y1 = nrow(table_to_use) + 2.5 + ifelse(use_range_for_maginal_labels, 0, -0.5) + 
           ifelse(order_by_counts, 0, -0.8) + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0) + 
           ifelse(subset_to_traits, 1, 1.5), lwd = 2)
text(x = -2, y = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), pos = 2, 
     labels = paste0(ifelse(use_range_for_maginal_labels, "total", "max"), " # of\nPrediXcan hits"))

#vertical label for total
segments(x0 = ncol(table_to_use) + 2.5, x1 = ncol(table_to_use) + 1, 
         y0 = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), 
         y1 = nrow(table_to_use) + 0.75 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), lwd = 2)
segments(x0 = ncol(table_to_use) + 2.5, x1 = ncol(table_to_use) + 2 + 
           ifelse(subset_to_traits, 0, 1.25), 
         y0 = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), 
         y1 = nrow(table_to_use) + 0.75 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), lwd = 2)
text(x = ncol(table_to_use) + 2 + ifelse(order_by_counts, 0, 2), y = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), pos = 3, 
     labels = paste0(ifelse(use_range_for_maginal_labels, "total", "max"), " # of DEGs"))


#legend and title labels
text(labels = ifelse(use_counts, latex2exp::TeX("$n_{intersect}$"), latex2exp::TeX("‰_{ intersect}")), pos = 3, font = 2, cex = 1.25,
     x = ncol(table_to_use) + x_adj + 2, y = nrow(table_to_use) + y_adj + 0.875 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0))
if(use_counts){
  text(latex2exp::TeX(paste0("number of genes in 8w - F↓M↓+ 8w - F↑M↑ with IHW-significant PrediXcan at $\\alpha$ = 0.05")), 
       x = 2 + ifelse(subset_to_traits, 0, 20), 
       y = nrow(table_to_use) + 3.25 + 
         ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0) + 
         ifelse(subset_to_traits, 0, 1), pos = 4, cex = 1.5, font = 2)
} else {
  text(latex2exp::TeX(paste0("proportion (‰) of IHW significant TWAS hits at $\\alpha$ = 0.05 in 8w - F↓M↓ or 8w - F↑M↑")), 
       x = 0, y = nrow(table_to_use) + 4.25 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), pos = 4, cex = 2.35, font = 2)
}

fig_label("a)", shrinkX = 0.85, cex = 2, shrinkY = 0.96, xpd = NA)


# ~~~~~~~~~~~~ #
# scatterplots #
# ~~~~~~~~~~~~ #

use_focal_vs_compl <- T

data1 <- data.frame(count = as.integer(unlist(c(n_deg_sigtwas_intersect))))
data1$tissue <- rep(rownames(n_deg_sigtwas_intersect), ncol(n_deg_sigtwas_intersect))
data1$trait <- unlist(lapply(colnames(n_deg_sigtwas_intersect), function(tri) rep(tri, nrow(n_deg_sigtwas_intersect))))
data1$total <- n_genes_in_nodes[data1$tissue]
data1$TWAS_Hit <- "YES"

traits <- unique(data1$trait)
tissues <- unique(data1$tissue)
d <- list(cell_count = data1$count,
          total = sapply(1:nrow(data1), function(i) total_number_of_possible_hits_matrix[data1$tissue[i], data1$trait[i]]),
          row_count = sapply(1:nrow(data1), function(i) n_genes_in_nodes_matrix[data1$tissue[i], data1$trait[i]]),
          col_count = sapply(1:nrow(data1), function(i) sig_twas_by_trait_genes_matrix[data1$tissue[i], data1$trait[i]]),
          row_index = match(data1$tissue, tissues),
          col_index = match(data1$trait, traits),
          row_n = length(unique(data1$tissue)),
          col_n = length(unique(data1$trait)))




# category_shapes <- setNames(15:19, salient.categories)
category_shapes <- setNames(15:19, c("Cardiometabolic", "Aging", "Anthropometric", 
                                     "Immune", "Psychiatric-neurologic")
)
category_shapes <- setNames(c(24,19,15,43,23,25)[1:length(focal_traitcats)], focal_traitcats)



par(mar = mars[[2]])
# layout(t(c(rep(1,10), rep(2,10), 3)))
# layout(t(c(1,2,3)), widths = c(1,1,0.1))
#first plot


# sapply(1:nrow(signif_df), function(i) adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index[i]]], as.numeric(signif_df$signif[i] != 0) * 0.625 + 0.125))
if(use_focal_vs_compl){
  xvals <- (d$col_count - d$cell_count) / (d$total - d$row_count)
  yvals <- d$cell_count / d$row_count
} else {
  xvals <- d$row_count / d$total * d$col_count / d$total  
  yvals <- d$cell_count / d$total
}

# plot(logit(xvals), logit(yvals), xlim = ifelse2(use_focal_vs_compl, c(-14,-1), c(-15,-4.25)), 
#      ylim = ifelse2(use_focal_vs_compl, c(-6, -1), c(-9.5,-4.25)),
#      pch = category_shapes[traitwise_partitions$Category[match(traits[d$col_index], traitwise_partitions$Tag)]], 
#      col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5), bg = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5),
#      cex = 1.5,
#      xlab = "", ylab = "", yaxt = 'n', xaxt = "n")
plot(logit(xvals), logit(yvals), xlim = ifelse2(use_focal_vs_compl, c(-14,-1), c(-15,-4.25)), 
     ylim = ifelse2(use_focal_vs_compl, c(-6, -1), c(-9.5,-4.25)),
     xlab = "", ylab = "", yaxt = 'n', xaxt = "n", type = "n")
#can try rasterizing points?
# raster_points(x = logit(xvals), y =  logit(yvals),
#               pch = category_shapes[traitwise_partitions$Category[match(traits[d$col_index], traitwise_partitions$Tag)]], 
#               col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5), bg = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5),
#               cex = 1.5, raster_res = 200)
#smaller file with regular points, not that much data fundamentally
points(x = logit(xvals), y =  logit(yvals),
       pch = category_shapes[traitwise_partitions$Category[match(traits[d$col_index], traitwise_partitions$Tag)]], 
       col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5), bg = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5),
       cex = 1.5)
xw <- diff(par("usr")[1:2])
yh <- diff(par("usr")[3:4])
box(which = "plot", lwd = 2)
#axes
text(labels = "focal logit(frequency)", x = par("usr")[2] + diff(par("usr")[1:2])/8, srt = 270, pos = 1, 
     y = mean(par("usr")[3:4]), xpd = NA, cex = 1.25)
text(labels = "complement logit(frequency)", x = mean(par("usr")[1:2]), srt = 0, pos = 1, 
     y = par("usr")[3] - diff(par("usr")[3:4])/12, xpd = NA, cex = 1.25)
yaxlocs <- ifelse2(use_focal_vs_compl, c(-6:-1), c(-5:-9))
segments(x0 = par("usr")[2], x1 = par("usr")[2] + diff(par("usr")[1:2])/100, y0 = yaxlocs, y1 = yaxlocs, lwd = 2, xpd = NA)
text(labels = yaxlocs, x = par("usr")[2] + diff(par("usr")[1:2])/200, srt = 0, pos = 4, 
     y = yaxlocs, xpd = NA, cex = 1.25)
xaxlocs <- ifelse2(use_focal_vs_compl, -7:-1*2, -7:-3*2)
segments(x0 = xaxlocs, x1 = xaxlocs, y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4])/80, lwd = 2, xpd = NA)
text(labels = xaxlocs, x = xaxlocs, srt = 0, pos = 1, 
     y = par("usr")[3] - diff(par("usr")[3:4])/80, xpd = NA, cex = 1.25)
#1 to 1 line
abline(0,1, lty = 2, lwd = 4, col = adjustcolor(1,0.75), xpd = F)

#add arrows for enriched vs depleted
arrow_locs_enrich <- c(-3, -3, -2.25, -1.5)
grad_arrow_curve(arrow_locs_enrich, direc = "v", cols = c("white", "red"), 
                 prop_head_width = 3, prop_shaft_length = 0.2, taper_ratio = 0.7, outline_col = "white",
                 w = 0.6)
text(x = mean(arrow_locs_enrich[1:2]), y = sum(arrow_locs_enrich[3:4] * c(0.35,0.65)), 
     labels = "enriched", col = "white", srt = 90, cex = 0.5)

arrow_locs_deplete <- c(-2, -2, -3.5, -4.25)
grad_arrow_curve(arrow_locs_deplete, direc = "v", cols = c("white", "blue"), 
                 prop_head_width = 3, prop_shaft_length = 0.2, taper_ratio = 0.7, outline_col = "white",
                 w = 0.6)

text(x = mean(arrow_locs_deplete[1:2]), y = sum(arrow_locs_deplete[3:4] * c(0.35,0.65)), 
     labels = "depleted", col = "white", srt = 270, cex = 0.5)



#secondplot
xl <- par("usr")[1]
xr <- par("usr")[1] + diff(par("usr")[1:2]) / 2
yt <- par("usr")[4]
yb <- par("usr")[4] - diff(par("usr")[3:4]) / 2

buffer <- 0.075
bx <- (xr-xl) * buffer
by <- (yt-yb) * buffer
rect(xleft = xl, ybottom = yb, xright = xr, ytop = yt, lwd = 2, xpd = NA)

points(x = ((xvals - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
       y = ((yvals - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, 
       pch = category_shapes[traitwise_partitions$Category[match(traits[d$col_index], traitwise_partitions$Tag)]], 
       col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5), bg = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5), 
       cex = 1)
one2one_coords <- c(x0 = 0 * (xr - xl - 2*bx) + xl + bx,
                    y0 = 0 * (yt - yb - 2*by) + yb + by,
                    x1 = ((max(yvals, na.rm = T) - 
                             min(xvals, na.rm = T)) / 
                            max(xvals, na.rm = T)) * 
                      (xr - xl - 2*bx) + xl + bx,
                    y1 = ((max(yvals, na.rm = T) - 
                             min(yvals, na.rm = T)) / 
                            max(yvals, na.rm = T)) * 
                      (yt - yb - 2*by) + yb + by)
one2one_slope <- (one2one_coords["y1"] - one2one_coords["y0"]) / 
  (one2one_coords["x1"] - one2one_coords["x0"])
one2one_coords["y1"] <- one2one_coords["y1"] - (one2one_coords["x1"] - xr) * one2one_slope
one2one_coords["x1"] <- xr

segments(x0 = one2one_coords["x0"],
         y0 = one2one_coords["y0"],
         x1 = one2one_coords["x1"],
         y1 = one2one_coords["y1"],
         lty = 2, lwd = 2, col = adjustcolor(1,0.75))

#axes
text(labels = "focal frequency", x = par("usr")[1] - diff(par("usr")[1:2])/ifelse2(use_focal_vs_compl, 8, 7), srt = 90, pos = 1, 
     y = par("usr")[4] - diff(par("usr")[3:4])/5, xpd = NA, cex = 1)
text(labels = "complement frequency", x = par("usr")[1] + diff(par("usr")[1:2])/4, srt = 0, pos = 1, 
     y = par("usr")[4] + diff(par("usr")[3:4])/8, xpd = NA, cex = 1)
xaxlocs <- ifelse2(use_focal_vs_compl, c(0:4/20), c(0:6/500))
yaxlocs <- ifelse2(use_focal_vs_compl, c(0:6/20), c(0:6/500))

segments(x0 = par("usr")[1], x1 = par("usr")[1] - diff(par("usr")[1:2])/100, 
         y0 = ((yaxlocs - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, 
         y1 = ((yaxlocs - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, 
         lwd = 2, xpd = NA)
segments(y0 = par("usr")[4], y1 = par("usr")[4] + diff(par("usr")[3:4])/100, 
         x0 = ((xaxlocs - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
         x1 = ((xaxlocs - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
         lwd = 2, xpd = NA)
text(labels = xaxlocs, x = ((xaxlocs - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
     srt = 0, pos = 3, 
     y = par("usr")[4], xpd = NA, cex = 0.75)
text(labels = yaxlocs, x = par("usr")[1], 
     srt = 0, pos = 2, 
     y = ((yaxlocs - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, xpd = NA, cex = 0.75)
# text(labels = 0:6/500, x = -7:-3*2, srt = 0, pos = 1, 
#      y = par("usr")[3] - diff(par("usr")[3:4])/80, xpd = NA, cex = 1.25)
# abline(0,1, lty = 2, lwd = 6, col = adjustcolor(1,0.75))

#legends
legend(x = xl, y = yb, legend = names(category_shapes), pch = category_shapes, col = adjustcolor(1, 0.5), pt.bg = adjustcolor(1, 0.5), cex = 0.75, pt.cex = 1.4)
legend(x = xr, y = yt, legend = tissues, pch = 19, col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues], 0.5), cex = 0.55, pt.cex = 1.2)
segments(x0 = par("usr")[1] + diff(par("usr")[1:2]) / 1.55, 
         y0 = par("usr")[4] - diff(par("usr")[3:4])/50,
         x1 = par("usr")[1] + diff(par("usr")[1:2]) / 1.45, 
         y1 = par("usr")[4] - diff(par("usr")[3:4])/50, 
         lty = 2, lwd = 2, col = adjustcolor(1,0.75))
text(x = par("usr")[1] + diff(par("usr")[1:2]) / 1.45, 
     y = par("usr")[4] - diff(par("usr")[3:4])/45, 
     pos = 4, cex = 0.75, labels = "1-to-1 line")


#make second plot with proportion disease-like effects data
par(mar = mars[[3]])
data <- fread("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/prop_positive_effects_in_DEG-TWAS_intersection.txt")
xvals <- data$count_all / data$total_all
yvals <- data$count_inters / data$total_inters


cex_pow <- 1/2
plot(xvals, yvals, xlim = c(0.48, 0.52), 
     ylim = c(0,1),
     pch = category_shapes[traitwise_partitions$Category[match(data$trait, traitwise_partitions$Tag)]], 
     col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[data$tissue], 0.5), bg = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[data$tissue], 0.5), 
     cex = (data$total_inters / max(data$total_inters))^cex_pow * 5,
     xlab = "", ylab = "", yaxt = 'n', xaxt = "n")
box(which = "plot", lwd = 2)

#axes
text(labels = "focal frequency (+ effects)", x = par("usr")[1] - diff(par("usr")[1:2])/7, srt = 90, pos = 1, 
     y = mean(par("usr")[3:4]), xpd = NA, cex = 1.25)
text(labels = "complement frequency (+ effects)", x = mean(par("usr")[1:2]), srt = 0, pos = 1, 
     y = par("usr")[3] - diff(par("usr")[3:4])/12, xpd = NA, cex = 1.25)
yaxlocs <- seq(0,1,by=0.1)
segments(x0 = par("usr")[1], x1 = par("usr")[1] - diff(par("usr")[1:2])/100, y0 = yaxlocs, y1 = yaxlocs, lwd = 2, xpd = NA)
text(labels = yaxlocs, x = par("usr")[1] - diff(par("usr")[1:2])/200, srt = 0, pos = 2, 
     y = yaxlocs, xpd = NA, cex = 1.25)
xaxlocs <- round(seq(par("usr")[1], par("usr")[2], by = 0.01), 2)
segments(x0 = xaxlocs, x1 = xaxlocs, y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4])/80, lwd = 2, xpd = NA)
text(labels = xaxlocs, x = xaxlocs, srt = 0, pos = 1, 
     y = par("usr")[3] - diff(par("usr")[3:4])/80, xpd = NA, cex = 1.25)
# abline(0,1, lty = 2, lwd = 6, col = adjustcolor(1,0.75))
abline(h=0.5,lwd=2,lty=2,col=adjustcolor(1,0.5),xpd=F)
abline(v=0.5,lwd=2,lty=2,col=adjustcolor(1,0.5),xpd=F)

#legend
xl <- par("usr")[2]
xr <- par("usr")[2] + diff(par("usr")[1:2]) / 2
yt <- par("usr")[4]
yb <- par("usr")[4] - diff(par("usr")[3:4]) / 2.5
legend(x = xl, y = yb + (yt-yb)/20, legend = names(category_shapes), bty = "n",
       pch = category_shapes, col = adjustcolor(1, 0.5), pt.bg = adjustcolor(1, 0.5), cex = 0.75, pt.cex = 1.4, xpd = NA)
legend(x = xl + (xr-xl)/100, y = yt, legend = tissues, pch = 19, bty = "n", xpd = NA,
       col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues], 0.5), cex = 0.55, pt.cex = 1.2)

legcexs <- c(0.5, 1:5)
legvals <- (round((legcexs / 5) ^ (1/cex_pow) * max(data$total_inters)))
legcexs <- (legvals / max(data$total_inters))^cex_pow * 5
# points(x = rep(xl + (xr-xl)/11, 6) + (rev(legcexs)) / 4000, y = yb - (yt-yb) * 1.4 + 
#          cumsum(rep((yt-yb)/60, length(legcexs))) + cumsum(legcexs + c(0,legcexs[-length(legcexs)])) / 120, 
#        cex = legcexs, xpd = NA, pch = 19, col = adjustcolor(1, 0.5))
text(labels = latex2exp::TeX("$n_{intersect}$"), x = xl - (xr-xl)/50, y = yb - (yt-yb) * 0.6, pos = 4, xpd = NA)
fixed_incr <- (yt-yb)/25
y_disp <- (yt-yb) * 1.5
points(x = rep(xl + (xr-xl)/10, 6), y = yb - y_disp + 
         cumsum(rep(fixed_incr, length(legcexs))) + cumsum(legcexs + c(0,legcexs[-length(legcexs)])) / 120, 
       cex = legcexs, xpd = NA, pch = 19, col = adjustcolor(1, 0.5))
text(x = rep(xl + (xr-xl)/10, 6), y = yb - y_disp + 
       cumsum(rep(fixed_incr, length(legcexs))) + cumsum(legcexs + c(0,legcexs[-length(legcexs)])) / 120, 
     cex = legcexs / 5, xpd = NA, labels = legvals, col = "white")
text(labels = legvals[1:2], x = xl + (xr-xl)/15 + cumsum(legcexs)[1:2]/4000, y = yb - y_disp + 
       cumsum(rep(fixed_incr, length(legcexs)))[1:2] + cumsum(legcexs + c(0,legcexs[-length(legcexs)]))[1:2] / 120, 
     pos = 4, xpd = NA, cex = 0.5)

#add arrows for enriched vs depleted
arrow_locs_enrich <- c(0.519, 0.519, 0.525, 0.525 + 0.75 * diff(par("usr")[3:4]) / yh)
grad_arrow_curve(arrow_locs_enrich, direc = "v", cols = c("white", "red"), 
                 prop_head_width = 3, prop_shaft_length = 0.2, taper_ratio = 0.7, outline_col = "white",
                 w = 0.6 * diff(par("usr")[1:2]) / xw)
text(x = mean(arrow_locs_enrich[1:2]), y = sum(arrow_locs_enrich[3:4] * c(0.35,0.65)), 
     labels = "enriched", col = "white", srt = 90, cex = 0.5)

arrow_locs_deplete <- c(0.519, 0.519, 0.475, 0.475 - 0.75 * diff(par("usr")[3:4]) / yh)
grad_arrow_curve(arrow_locs_deplete, direc = "v", cols = c("white", "blue"), 
                 prop_head_width = 3, prop_shaft_length = 0.2, taper_ratio = 0.7, outline_col = "white",
                 w = 0.6 * diff(par("usr")[1:2]) / xw)

text(x = mean(arrow_locs_deplete[1:2]), y = sum(arrow_locs_deplete[3:4] * c(0.35,0.65)), 
     labels = "depleted", col = "white", srt = 270, cex = 0.5)

fig_label("b)", shrinkX = 0.865, cex = 2, shrinkY = 0.95, xpd = NA)
fig_label("c)", shrinkX = 0.99, shrinkY = 0.95, cex = 2)


# ~~~~~~~~~~~~ #
# violin plots #
# ~~~~~~~~~~~~ #

par(mar = mars[[4]])
focal_samps <- subset_samps("row_bias", c("raw", "sd"), samps = samps) + samps$overall_bias
n_to_trim <- round(nrow(focal_samps) * 0.002)
focal_samps <- apply(focal_samps, 2, trim_n, n_to_trim)

tissues <- tissues_intersect.model
colnames(focal_samps) <- tissues
tord <- order(apply(focal_samps, 2, mean))
qi_95 <- apply(focal_samps, 2, quantile, probs = c(0.05, 0.95))[,tord]
qi_100 <- apply(focal_samps, 2, quantile, probs = c(0.0, 1.0))[,tord]
posterior_means <- apply(focal_samps, 2, mean)[tord]
inside_cols <- rep("black", length(tord))
inside_cols[overlaps_with_zero(qi_95)] <- "white"
tmp <- vioplot::vioplot(x = focal_samps[,tord], T, 
                        col = MotrpacRatTraining6moData::TISSUE_COLORS[colnames(focal_samps)][tord], outline=FALSE, yaxt = "n",
                        names = colnames(focal_samps)[tord], range = 0, ylab = "", lineCol = MotrpacRatTraining6moData::TISSUE_COLORS[colnames(focal_samps)][tord],
                        xlab = "", cex.lab = 2, plotCentre = "point", 
                        colMed = MotrpacRatTraining6moData::TISSUE_COLORS[colnames(focal_samps)][tord],
                        horizontal = T) 


#axes
xtickvals <- seq3(range(qi_100), 5, 0)
segments(x0 = xtickvals, x1 = xtickvals, y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4])/100, xpd = NA)
text(x = xtickvals, y = par("usr")[3] - diff(par("usr")[3:4])/30, labels = xtickvals)
text(x = mean(par("usr")[1:2]), y = par("usr")[3]-diff(par("usr")[3:4])/10, 
     labels = "multilevel tissue\nenrichment (logodds-scale)", cex = 1.25, xpd = NA)
tick <- seq_along(colnames(focal_samps))
axis(2, at = tick, labels = F)
for(i in 1:length(colnames(focal_samps))){
  segments(y0 = i, y1 = i, x0 = qi_95[1,i], x1 = qi_95[2,i], lwd = 4, col = inside_cols[i])
  points(y = i, x = posterior_means[i], col = inside_cols[i], pch = 19, cex = 1.5)
  
  #guiding lines
  segments(y0 = i, y1 = i, x0 = min(qi_100), x1 = qi_100[1,i], lwd = 1, col = "black", xpd = T, lty = 3)
  segments(y0 = i, y1 = i, x0 = qi_100[2,i], x1 = max(qi_100), lwd = 1, col = "black", xpd = T, lty = 3)
  
}
#group labels
text(y = tick, x = rep(par("usr")[1] - 0.01, length(tick)), colnames(focal_samps)[tord], srt = 0, xpd = T, pos = 2,
     col = MotrpacRatTraining6moData::TISSUE_COLORS[colnames(focal_samps)[tord]], xpd = NA, cex = 1.25)

abline(v=0,lwd=3,lty=2, col = adjustcolor(1,0.5), xpd = F)

#figure label
if(T){
  fig_label("d)", cex = 2, shrinkX = 0.9, shrinkY = 0.97)
}

#legend for violin plots
xval <- -1.25
yval <- 0.5 #2
yh <- 4
yt <- yval + yh/2
yb <- yval - yh/2
xr <- -30:30/10
xvals_disp <- (dnorm(xr) / 6) * 2
xvals <- cbind(xval + xvals_disp, xval - xvals_disp)
yvals <- seq(yb,yt,length.out = nrow(xvals))
# lines(xvals[,1], yvals, xpd = NA)
# lines(xvals[,2], yvals, xpd = NA)
polygon(x = c(xvals[,1], xvals[,2]), y = c(yvals, rev(yvals)), xpd = NA, col = adjustcolor(1,0.15))
ybloc <- yb + (qnorm(p = c(0.05,0.95)) - min(xr)) / diff(range(xr)) * yh
segments(x0 = xval, x1 = xval, y0 = ybloc[1], y1 = ybloc[2], xpd = NA, lwd = 2)
points(xval, yval, pch = 19, xpd = NA)
xlab_loc <- max(xvals) + diff(range(xvals)) / 5
segments(x0 = xval, x1 = xlab_loc, 
         y0 = c(max(yvals), min(yvals), yval, ybloc),
         y1 = c(max(yvals), min(yvals), yval, ybloc),
         xpd = NA, lty = 2)
violabs <- paste0("$P_{", c("99.9", "0.01", "\\mu", "5", "95"), "}$")
violabs[[3]] <- "$\\mu$"
text(latex2exp::TeX(violabs), x = xlab_loc - diff(range(xvals)) / 10, y = c(max(yvals), min(yvals), yval, ybloc), pos = 4, xpd = NA, cex = 0.875)

segments(x0 = min(xvals) - diff(range(xvals))/100, x1 = max(xvals) + diff(range(xvals))/100, 
         y0 = min(yvals) - diff(range(yvals)) / 8, y1 = min(yvals) - diff(range(yvals)) / 8, xpd = NA, lwd = 3.25)
points(xval, min(yvals) - diff(range(yvals)) / 8, pch = 19, xpd = NA, cex = 1.2)
segments(x0 = min(xvals), x1 = max(xvals), y0 = min(yvals) - diff(range(yvals)) / 8, y1 = min(yvals) - diff(range(yvals)) / 8, xpd = NA, lwd = 2, col = "white")
points(xval, min(yvals) - diff(range(yvals)) / 8, pch = 19, xpd = NA, cex = 1, col = "white")
text(latex2exp::TeX("$CI_{95}$ includes 0$"), x = max(xvals) + diff(range(xvals))/100, 
     y = min(yvals) - diff(range(yvals)) / 7, pos = 4, xpd = NA, cex = 0.875)

segments(x0 = min(xvals) - diff(range(xvals))/100, x1 = max(xvals) + diff(range(xvals))/100, 
         y0 = min(yvals) - diff(range(yvals)) / 8 * 2, y1 = min(yvals) - diff(range(yvals)) / 8 * 2, xpd = NA, lwd = 3.25)
points(xval, min(yvals) - diff(range(yvals)) / 8 * 2, pch = 19, xpd = NA, cex = 1.2)
text(latex2exp::TeX("$CI_{95}$ excludes 0$"), x = max(xvals) + diff(range(xvals))/100, 
     y = min(yvals) - diff(range(yvals)) / 7 - diff(range(yvals)) / 8, 
     pos = 4, xpd = NA, cex = 0.875)


###
#colcat bias violin plots
###

salient.categories <- trait_cats
if(length(salient.categories) < 9){
  category_colors <- RColorBrewer::brewer.pal(length(salient.categories), "Dark2")
  names(category_colors) <- sort(salient.categories)  
} else {
  category_colors <- RColorBrewer::brewer.pal(5, "Dark2")
  names(category_colors) <- sort(c("Cardiometabolic", "Aging", "Anthropometric", 
                                   "Immune", "Psychiatric-neurologic"))
  cats_leftover <- setdiff(salient.categories, names(category_colors))
  cols_leftover <- setdiff(RColorBrewer::brewer.pal(8, "Dark2"), category_colors)
  n_cols_needed <- length(cats_leftover) - length(cols_leftover)
  if(n_cols_needed > 0){
    extra_cols <- c(cols_leftover, RColorBrewer::brewer.pal(n_cols_needed, "Set2"))
  } else {
    extra_cols <- cols_leftover
  }
  names(extra_cols) <- cats_leftover
  category_colors <- c(category_colors, extra_cols)
}
cols$category <- category_colors

par(mar = mars[[5]])
focal_samps <- subset_samps("colcat_bias", c("raw", "sd"), samps = samps) + samps$overall_bias
focal_samps <- apply(focal_samps, 2, trim_n, n_to_trim)
colnames(focal_samps) <- trait_cats
tord <- order(apply(focal_samps, 2, mean))
qi_95 <- apply(focal_samps, 2, quantile, probs = c(0.05, 0.95))[,tord]
qi_100 <- apply(focal_samps, 2, quantile, probs = c(0.0, 1.0))[,tord]
posterior_means <- apply(focal_samps, 2, mean)[tord]
inside_cols <- rep("black", length(tord))
inside_cols[overlaps_with_zero(qi_95)] <- "white"

par(xpd=F)
tmp <- vioplot::vioplot(x = focal_samps[,tord], T, 
                        col = cols$category[colnames(focal_samps)][tord], outline=FALSE, yaxt = "n",
                        names = colnames(focal_samps)[tord], range = 0, ylab = "", 
                        lineCol = cols$category[colnames(focal_samps)][tord],
                        xlab = "", cex.lab = 2, plotCentre = "point", 
                        colMed = cols$category[colnames(focal_samps)][tord],
                        horizontal = T)

#axes
xtickvals <- seq3(range(qi_100), 5, 0)
segments(x0 = xtickvals, x1 = xtickvals, y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4])/100, xpd = NA)
text(x = xtickvals, y = par("usr")[3] - diff(par("usr")[3:4])/30, labels = xtickvals)
text(x = mean(par("usr")[1:2]), y = par("usr")[3]-diff(par("usr")[3:4])/10, 
     labels = "multilevel category\nenrichment (logodds-scale)", cex = 1.25, xpd = NA)
tick <- seq_along(colnames(focal_samps))
axis(2, at = tick, labels = F)
for(i in 1:length(colnames(focal_samps))){
  segments(y0 = i, y1 = i, x0 = qi_95[1,i], x1 = qi_95[2,i], lwd = 4, col = inside_cols[i])
  points(y = i, x = posterior_means[i], col = inside_cols[i], pch = 19, cex = 1.5)
  
  #guiding lines
  segments(y0 = i, y1 = i, x0 = min(qi_100), x1 = qi_100[1,i], lwd = 1, col = "black", xpd = T, lty = 3)
  segments(y0 = i, y1 = i, x0 = qi_100[2,i], x1 = max(qi_100), lwd = 1, col = "black", xpd = T, lty = 3)
  
}
catlabs <- colnames(focal_samps)[tord]
# catlabs <- gsub("system ", "system\n", gsub("-", "-\n", catlabs))
catlabs <- gsub("system ", "", gsub("-neurologic", "", catlabs))
text(y = tick, x = rep(par("usr")[1] - 0.01, length(tick)), 
     labels = catlabs, 
     srt = 0, xpd = T, pos = 2, cex = 1.25,
     col = cols$category[colnames(focal_samps)[tord]], xpd = NA)


gsub("_", "\n", colnames(focal_samps))
gsub("Gi", "GI", stringr::str_to_title(gsub("_", "\n", colnames(focal_samps))))

abline(v=0,lwd=3,lty=2, col = adjustcolor(1,0.5), xpd = F)

#figure label
if(T){
  fig_label("e)", cex = 2, shrinkX = 0.905, shrinkY = 0.985)
}

##
# traits
##

par(mar = mars[[6]])
trait_subset_bool <- traitwise_partitions$Category[match(traits, traitwise_partitions$Tag)] %in% salient.categories
focal_samps <- subset_samps("col_bias", c("raw", "sd"), samps = samps) + samps$overall_bias +
  subset_samps("colcat_bias", c("raw", "sd"), samps = samps)[,match(trait_categories$Category[match(traits, trait_categories$Tag)], trait_cats)]
focal_samps <- apply(focal_samps, 2, trim_n, n_to_trim)
colnames(focal_samps) <- traits
tord <- order(apply(focal_samps, 2, mean))
qi_95 <- apply(focal_samps, 2, quantile, probs = c(0.05, 0.95))[,tord]
qi_100 <- apply(focal_samps, 2, quantile, probs = c(0.0, 1.0))[,tord]
posterior_means <- apply(focal_samps, 2, mean)[tord]
inside_cols <- rep("black", length(tord))
inside_cols[overlaps_with_zero(qi_95)] <- "white"

tmp <- vioplot::vioplot(x = focal_samps[,tord], T,
                        col = cols$category[trait_categories$Category[match(colnames(focal_samps)[trait_subset_bool][tord], trait_categories$Tag)]],
                        outline=FALSE, yaxt = "n",
                        names = colnames(focal_samps)[tord], range = 0, ylab = "",
                        lineCol = cols$category[trait_categories$Category[match(colnames(focal_samps)[trait_subset_bool][tord], trait_categories$Tag)]],
                        xlab = "", cex.lab = 2, plotCentre = "point",
                        colMed = cols$category[colnames(focal_samps)][tord],
                        horizontal = F,
                        xlim = c(4,ncol(focal_samps)-3))

ytickvals <- seq3(range(qi_100), 5, 0)
segments(y0 = ytickvals, y1 = ytickvals, x0 = par("usr")[1], x1 = par("usr")[1] - diff(par("usr")[1:2])/200, xpd = NA)
text(y = ytickvals, x = par("usr")[1] - diff(par("usr")[1:2])/75, labels = ytickvals)
text(y = mean(par("usr")[3:4]), x = par("usr")[1]-diff(par("usr")[1:2])/25,
     labels = "multilevel trait enrichment\n(logodds-scale)", cex = 1.25, xpd = NA, srt = 90)
tick <- seq_along(colnames(focal_samps))


tick <- seq_along(traits[trait_subset_bool])
axis(1, at = tick, labels = F)
for(i in 1:length(traits[trait_subset_bool])){
  segments(x0 = i, x1 = i, y0 = qi_95[1,i], y1 = qi_95[2,i], lwd = 2, col = inside_cols[i])
  points(x = i, y = posterior_means[i], col = inside_cols[i], pch = 19)
  #guiding lines
  segments(x0 = i, x1 = i, y0 = min(qi_100) - diff(range(qi_100)) / 30, y1 = qi_100[1,i], lwd = 1, col = "black", xpd = T, lty = 3)
  segments(x0 = i, x1 = i, y0 = qi_100[2,i], y1 = max(qi_100) * 0.75 + 1 * 0.25, lwd = 1, col = "black", xpd = T, lty = 3)
}
trait_labs <- trait_categories$new_Phenotype[match(colnames(focal_samps)[trait_subset_bool][tord], trait_categories$Tag)]
trait_labs[which.max(nchar(trait_labs))]
trait_labs[trait_labs == "Attention_Deficit_Hyperactivity_Disorder"] <- "ADHD"
trait_labs <- gsub("craneal", "cranial", trait_labs)
text(tick + 0.55, par("usr")[3] - diff(par("usr")[3:4])/12, cex = 0.8,
     labels = trait_labs, 
     srt = 45, xpd = T, pos = 2,
     col = cols$category[trait_categories$Category[match(colnames(focal_samps)[trait_subset_bool][tord], trait_categories$Tag)]], xpd = NA)
fig_label("f)", cex = 2, shrinkX = 0.4)
abline(h=0,lwd=3,lty=2, col = adjustcolor(1,0.5), xpd = F)

dev.off()

}

#### diff frequentist test for reviewer ####

#load in some meadata
trait_categories <- read.csv("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/gwas_metadata.csv", header = T)
trait_categories$Category[trait_categories$new_Phenotype == "Multiple_Sclerosis"] <- "Psychiatric-neurologic" #correct original coding from "Cardiometabolic"

library(fgsea)
trait_name_key <- setNames(trait_categories$new_Phenotype, trait_categories$Tag)

#make directory
gsea_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/gsea_output/"
if(!dir.exists(gsea_dir)){dir.create(gsea_dir)}

if(file.exists(paste0(gsea_dir, "gsea_trait_x_tissue.txt"))){
  gsea_output <- fread(paste0(gsea_dir, "gsea_trait_x_tissue.txt"))
} else {
  
  gsea_output <- do.call(rbind, lapply(tissues, function(tiss){
    
    
    gsea_results <- parallel::mclapply(split(twas_by_tissue[[tiss]], twas_by_tissue[[tiss]]$trait), 
                                       function(trait_twas){
                                         
                                         tryCatch({
                                           
                                           trait <- trait_twas$trait[1]
                                           mcprint(paste0(which(tissues == tiss), ".", tiss, ".", trait))
                                           trait_twas$log10_pvalue <- -log10(trait_twas$pvalue)
                                           ranked_genes <- setNames(object = trait_twas$log10_pvalue, 
                                                                    nm = gsub("\\..*", "", trait_twas$gene))
                                           DE_genes <- list(intersect(DE_genes_x_tissue[[tiss]], names(ranked_genes)))
                                           names(DE_genes) <- paste0(tiss, ".", trait)
                                           gsea_out <- fgsea(pathways = DE_genes, stats = ranked_genes, scoreType = "pos")
                                           
                                           if("data.frame" %in% class(gsea_out)){
                                             gsea_out$tissue <- tiss
                                             gsea_out$trait <- trait  
                                           }
                                           
                                           gsea_out
                                           
                                         }, error = function(e){
                                           print(paste("Error msg: ", e$message))
                                           return(NULL)
                                         })
                                       }, mc.cores = 8)
    
    # gsea_results <- gsea_results[lapply(gsea_results, class) != "try-error"]
    gsea_results <- do.call(rbind, gsea_results)
    gsea_results
  }))
  
  gsea_output$trait_name <- trait_name_key[gsea_output$trait]
  gsea_output$pathway_name <- paste0(gsea_output$tissue, ".", gsea_output$trait_name)
  hist(gsea_output$pval, breaks = 0:100/100)
  # gsea_output$adj.pval <- p.adjust(gsea_output$pval, method = "bh")
  gsea_output$adj.pval <- p.adjust(gsea_output$pval, method = "bonf")
  
  fwrite(gsea_output, paste0(gsea_dir, "gsea_trait_x_tissue.txt"))
  
}

hist(gsea_output$adj.pval, breaks = 0:100/100)
head(gsea_output[order(gsea_output$pval),c("trait_name", "tissue", "pval", "adj.pval")], 20)

#try meta-analysis (using HMP?)
library(harmonicmeanp)
L <- nrow(gsea_output)
tissue_ps <- sapply(split(gsea_output$pval, gsea_output$tissue), function(tissps){
  (length(tissps)/L)/sum((1/L)/tissps)
})
tissue_thresh <- qharmonicmeanp(0.05, L) * table(gsea_output$tissue)/L
tissue_thresh <- tissue_thresh[names(tissue_ps)]
any(tissue_ps < tissue_thresh)
head(sort(tissue_ps), n = 10)

trait_ps <- sapply(split(gsea_output$pval, gsea_output$trait_name), function(traitps){
  (length(traitps)/L)/sum((1/L)/traitps)
})
trait_thresh <- qharmonicmeanp(0.05, L) * table(gsea_output$trait_name)/L
trait_thresh <- trait_thresh[names(trait_ps)]
any(trait_ps < trait_thresh)
head(sort(trait_ps), n = 10)
sort(trait_ps / (trait_thresh / 0.05))


trait_cat_key <- setNames(trait_categories$Category, 
                       trait_categories$new_Phenotype)[names(trait_ps)]
L_traits <- length(trait_ps)
trait_category_ps <- sapply(split(trait_ps, trait_cat_key), function(traitcatps){
  (length(traitcatps)/L_traits)/sum((1/L_traits)/traitcatps)
})
trait_category_thresh <- qharmonicmeanp(0.05, L_traits) * 
  table(trait_cat_key) / L_traits
trait_category_thresh <- trait_category_thresh[names(trait_category_ps)]
any(trait_category_ps < trait_category_thresh)
sort(trait_category_ps / (trait_category_thresh / 0.05))

sort(trait_category_ps)

# #try just set unioning the tissue DE genes and combining the predixcan pvals? with HMP also?
tissue_x_twas <- lapply(setNames(unique(some.twas$trait),
                                 unique(some.twas$trait)),
                        function(traiti) {
                          print(traiti)
                          some.twas[some.twas$trait == traiti & compatible_twas_genes,]
                        })
# all_DE_genes <- unlist(lapply(seq_along(DE_genes_x_tissue), function(i) paste0(names(DE_genes_x_tissue)[i], ".", DE_genes_x_tissue[[i]])))
# gsea_output_pooled <- parallel::mclapply(tissue_x_twas, 
#                                          function(trait_twas){
#                                            
#                                            tryCatch({
#                                              
#                                              trait <- trait_twas$trait[1]
#                                              mcprint(trait)
#                                              
#                                              trait_twas$log10_pvalue <- -log10(trait_twas$pvalue)
#                                              ranked_genes <- setNames(object = trait_twas$log10_pvalue, 
#                                                                       nm = paste0(tissue_abbr[trait_twas$tissue], ".", gsub("\\..*", "", trait_twas$gene)))
#                                              # DE_genes <- list(intersect(unique(unlist(DE_genes_x_tissue)), names(ranked_genes)))
#                                              # DE_genes <- list(intersect(names(table(unlist(DE_genes_x_tissue))[table(unlist(DE_genes_x_tissue)) >= 3]), 
#                                              #                            names(ranked_genes)))
#                                              
#                                              #only use 1-to-1 genes here?
#                                              DE_genes <- list(intersect(all_DE_genes, names(ranked_genes)))
#                                              names(DE_genes) <- trait
#                                              gsea_out <- fgseaMultilevel(pathways = DE_genes, 
#                                                                          stats = ranked_genes, 
#                                                                          scoreType = "pos",
#                                                                          eps = 1E-10)
#                                              
#                                              if("data.frame" %in% class(gsea_out)){
#                                                gsea_out$trait <- trait  
#                                              }
#                                              
#                                              return(gsea_out)
#                                              
#                                            }, error = function(e){
#                                              print(paste("Error msg: ", e$message))
#                                              return(NULL)
#                                            })
#                                            
#                                          }, mc.cores = 8)
# 
# gsea_output_pooled <- gsea_output_pooled[lapply(gsea_output_pooled, class) != "NULL"]
# gsea_output_pooled <- do.call(rbind, gsea_output_pooled)
# gsea_output_pooled$trait_name <- setNames(trait_categories$new_Phenotype, trait_categories$Tag)[
#   gsea_output_pooled$trait]
# gsea_output_pooled$trait_fullname <- setNames(trait_categories$Phenotype, trait_categories$Tag)[
#   gsea_output_pooled$trait]
# gsea_output_pooled$adj.pval <- p.adjust(gsea_output_pooled$pval, method = "BH")
# hist(gsea_output_pooled$adj.pval)
# head(gsea_output_pooled[order(gsea_output_pooled$pval),
#                         c("trait_name", "adj.pval")], 20)

#or try combining predixcan probabilities?
if(file.exists(paste0(gsea_dir, "gsea_trait_pooled.txt"))){
  gsea_output_pooled_alt <- fread(paste0(gsea_dir, "gsea_trait_pooled.txt"))
} else {
  gsea_output_pooled_alt <- parallel::mclapply(tissue_x_twas, 
                                               function(trait_twas){
                                                 
                                                 tryCatch({
                                                   
                                                   trait <- trait_twas$trait[1]
                                                   mcprint(trait)
                                                   
                                                   trait_twas$log10_pvalue <- -log10(trait_twas$pvalue)
                                                   trait_twas$ensembl <- gsub("\\..*", "", trait_twas$gene)
                                                   tissue_x_genes <- split(trait_twas$pvalue, trait_twas$ensembl)
                                                   L_txg <- length(unlist(tissue_x_genes))
                                                   ranked_genes <- -log10(sapply(tissue_x_genes, function(x){
                                                     (length(x)/L_txg)/sum((1/L_txg)/x)
                                                   }))
                                                   DE_genes <- list(intersect(names(table(unlist(DE_genes_x_tissue))[table(unlist(DE_genes_x_tissue)) >= 3]),
                                                                              names(ranked_genes)))
                                                   names(DE_genes) <- trait
                                                   gsea_out <- fgseaMultilevel(pathways = DE_genes, 
                                                                               stats = ranked_genes, 
                                                                               scoreType = "pos",
                                                                               eps = 1E-10)
                                                   
                                                   if("data.frame" %in% class(gsea_out)){
                                                     gsea_out$trait <- trait  
                                                   }
                                                   
                                                   return(gsea_out)
                                                   
                                                 }, error = function(e){
                                                   print(paste("Error msg: ", e$message))
                                                   return(NULL)
                                                 })
                                                 
                                               }, mc.cores = 8)
  
  gsea_output_pooled_alt <- gsea_output_pooled_alt[lapply(gsea_output_pooled_alt, class) != "NULL"]
  gsea_output_pooled_alt <- do.call(rbind, gsea_output_pooled_alt)
  gsea_output_pooled_alt$trait_name <- setNames(trait_categories$new_Phenotype, trait_categories$Tag)[
    gsea_output_pooled_alt$trait]
  gsea_output_pooled_alt$adj.pval <- p.adjust(gsea_output_pooled_alt$pval, method = "BH")

  fwrite(gsea_output_pooled_alt, paste0(gsea_dir, "gsea_trait_pooled.txt"))
  
  
}
hist(gsea_output_pooled_alt$adj.pval)
head(gsea_output_pooled_alt[order(gsea_output_pooled_alt$pval),
                            c("trait_name", "pval", "adj.pval")], 20)


#ANOTHER IDEA: create a supplemental figure that compares signed negative log10(pvalues) to log-odds of positive effect

#### supp fig 3 preprocessing ####
#first load in the MCMC output

base = paste0("deviation_from_expected_logodds_split_the_difference_informative-MVN-matrix-normal_priors") 
load(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/mcmc_output/", 
            base, 
            ifelse(T, "_allCats", ""), 
            ifelse(F, "_randomgenes", ""),
            ifelse(T, "", "_intersecting-genes-marginal-corrmats"),
            ifelse(T, "_pairwise-interactions", ""),
            ".cmdStanR.fit"))
samps <- data.frame(as_draws_df(out$draws()))
prop_greater_than_0 <- function(x) mean(x>0)

celltotalbias <- apply(subset_samps("cell_total_prob_bias", c("raw", "sd", "logodds", "L"), samps = samps), 2, prop_greater_than_0)
celltotalbias <- setNames(celltotalbias, paste0(data_subset$tissue, ".", 
                                                trait_name_key[data_subset$trait]))
celltotalbias <- sort(celltotalbias, T)
#adjust 0s and 1s to min and max non-0 and non-1 values
celltotalbias <- celltotalbias + (celltotalbias == 0) * 1/nrow(samps) + (celltotalbias == 1) * -1/nrow(samps)
head(celltotalbias, 20)

rowbias <- apply(subset_samps("row_bias", c("raw", "sd", "L"), samps = samps) + samps$overall_bias, 2, prop_greater_than_0)
rowbias <- setNames(sort(rowbias), tissues[order(rowbias)])
rowbias

colbias <- apply(subset_samps("col_bias", c("raw", "sd", "L"), samps = samps) + 
                   samps$overall_bias +
                   subset_samps("colcat_bias", c("raw", "sd", "L"), 
                                samps = samps)[,match(traitwise_partitions$Category[
                                  match(traits, traitwise_partitions$Tag)], trait_cats)], 
                 2, prop_greater_than_0)
trait_key <- setNames(trait_categories$new_Phenotype, trait_categories$Tag)
colbias <- setNames(sort(colbias), trait_key[traits[order(colbias)]])
tail(colbias, 20)

# colcatbias <- apply(subset_samps("colcat_bias", c("raw", "sd", "L"), samps = samps), 2, prop_greater_than_0)
colcatbias <- apply(subset_samps("colcat_bias", c("raw", "sd", "L"), samps = samps) + samps$overall_bias, 
                    2, prop_greater_than_0)
colcatbias <- setNames(sort(colcatbias), trait_cats[order(colcatbias)])
colcatbias

#### plotting for supp fig 3 ####

# plot(logit(xvals), logit(yvals), xlim = ifelse2(use_focal_vs_compl, c(-14,-1), c(-15,-4.25)), 
#      ylim = ifelse2(use_focal_vs_compl, c(-6, -1), c(-9.5,-4.25)),
#      pch = category_shapes[traitwise_partitions$Category[match(traits[d$col_index], traitwise_partitions$Tag)]], 
#      col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5), bg = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5),
#      cex = 1.5,
#      xlab = "", ylab = "", yaxt = 'n', xaxt = "n")

#xlab = "", ylab = "", yaxt = 'n', xaxt = "n"
# box(which = "plot", lwd = 2)

#start the plotting!

plot.scatter <- function(x, y, xlabel = "", ylabel = "", main.title = "",
                         axis.lab.dev = 0.1, lab.cex = 1.25, pt.labels = NULL,
                         pts.to.label = NULL, pch = 19, pt.col = adjustcolor(1,0.5),
                         pt.pos = NULL, pt.cex = 1, prop_label_push = NULL, 
                         prop_rect_expand = 0.05, ...){
  
  #clean data
  valid_inds <- which(apply(is.na(cbind(x,y)), 1, sum) == 0)
  invalid_inds <- setdiff(1:length(x), valid_inds)
  x <- x[valid_inds]
  y <- y[valid_inds]
  
  #find useful values
  xr <- range(x)
  yr <- range(y)
  
  # initialize blank plot
  plot(NULL, xaxt="n",yaxt="n",bty="n",ylab="",xlab="", main="", sub="",
       xlim = xr + c(-1,1) * diff(xr) * 0.2, 
       ylim = yr + c(-1,1) * diff(yr) * 0.2)
  
  #draw line @ 0.05?
  # abline(h = -log10(0.05), col = adjustcolor(1, 0.25), lwd = 2, lty = 3, xpd = F)
  
  #preprocess pt labels
  if(!is.null(pt.labels)){
    pt.labels <- pt.labels[valid_inds]
    if(is.null(pt.pos)){
      pt.pos <- 3
    }
    if(is.null(pts.to.label)){
      pts.to.label <- 1:length(pt.labels)
    } else {
      pts.to.label <- unlist(lapply(pts.to.label, function(z){
        if(z %in% invalid_inds){
          return(numeric(0))
        } else {
          z - sum(z > invalid_inds)  
        }
      }))
    }
    labx <- x[pts.to.label]
    laby <- y[pts.to.label]
    pt.lab <- pt.labels[pts.to.label]
    
    #find and expand bounding rectangles for text
    word_rects <- text.rects(x = labx, y = laby, pos = pt.pos, labels = pt.lab)
    word_rects <- lapply(word_rects, expand.rect, prop = prop_rect_expand)
    
    #get boxes a bit away from their corresponding pts
    if(is.null(prop_label_push)){
      prop_label_push <- rep(0.04, length(word_rects))
    }
    word_rects <- lapply(seq_along(word_rects), function(i){
      push_rects(word_rects[[i]], c(labx[i], laby[i]), prop_push = prop_label_push[i])
    })
    
    #get those bounding boxes inside the plot
    word_rects <- lapply(word_rects, respect_boundaries, extra = 0.03)
    
    #find closest side to connect to point
    closest_side_pts <- lapply(seq_along(word_rects), function(i){
      closest_side(word_rects[[i]], c(labx[i], laby[i]))
    })
    
  }
  
  #create frame
  box("plot", lwd = 2)
  
  #axes
  tickl <- diff(par("usr")[1:2])/100
  yaxlocs <- nice_ticks(yr, 5)
  yaxlocs <- c(yaxlocs, seq(tail(yaxlocs,1), par("usr")[4], by = diff(yaxlocs[1:2]))[-1])
  segments(x0 = par("usr")[1], 
           x1 = par("usr")[1] - tickl, 
           y0 = yaxlocs, 
           y1 = yaxlocs, 
           lwd = 2, 
           xpd = NA)
  text(labels = yaxlocs, 
       x = par("usr")[1] - tickl, 
       srt = 0, pos = 2, y = yaxlocs, xpd = NA, cex = lab.cex)
  
  xaxlocs <- nice_ticks(xr, 5)
  xaxlocs <- c(seq(head(xaxlocs,1), par("usr")[1], by = -diff(xaxlocs[1:2]))[-1],
               xaxlocs, 
               seq(tail(xaxlocs,1), par("usr")[2], by = diff(xaxlocs[1:2]))[-1])
  segments(y0 = par("usr")[3], 
           y1 = par("usr")[3] - tickl / xyrat(), 
           x0 = xaxlocs, 
           x1 = xaxlocs, 
           lwd = 2, 
           xpd = NA)
  text(labels = xaxlocs, 
       y = par("usr")[3] - tickl / xyrat(), 
       srt = 0, pos = 1, x = xaxlocs, xpd = NA, cex = lab.cex)
  
  #axis labels
  axis.lab.dist <- diff(par("usr")[1:2]) * axis.lab.dev
  maxsh <- max(strheight(xaxlocs, cex = lab.cex))
  maxsw <- max(strwidth(yaxlocs, cex = lab.cex)) / xyrat()
  text(labels = ylabel, 
       x = par("usr")[1] - axis.lab.dist * 1.05, 
       y = mean(par("usr")[3:4]),
       srt = 90, pos = 3, xpd = NA, cex = lab.cex)
  text(labels = xlabel, 
       y = par("usr")[3] - axis.lab.dist / xyrat() + (maxsw - maxsh), 
       x = mean(par("usr")[1:2]),
       srt = 0, pos = 1, xpd = NA, cex = lab.cex)
  
  #title
  text(labels = main.title, pos = 3, cex = lab.cex * 1.25,
       x = mean(par("usr")[1:2]), y = par("usr")[4], xpd = NA)
  
  if(!is.null(pt.labels)){
    
    #draw line segments and rectangles for labels
    for(i in 1:length(word_rects)){
      segments(x0 = closest_side_pts[[i]][1], 
               y0 = closest_side_pts[[i]][2], 
               x1 = labx[i], 
               y1 = laby[i],
               lwd = 2) 
    }
    
    for(i in 1:length(word_rects)){
      rect(xleft = word_rects[[i]][1], 
           xright = word_rects[[i]][2], 
           ybottom = word_rects[[i]][3], 
           ytop = word_rects[[i]][4], 
           col = adjustcolor("white", 1), border = "black",
           lwd = 2)
      
    }
    
    #compute rect midpoints and write labels
    labx <- sapply(word_rects, function(z) mean(z[1:2]))
    laby <- sapply(word_rects, function(z) mean(z[3:4]))
    
    text(labx, laby, labels = pt.lab)
  }
  
  #plot points
  pt_rand_order <- sample(1:length(x))
  points(x, y, pch = pch, col = "black", cex = pt.cex*1.2) #create a nice outer border
  points(x, y, pch = pch, col = "white", cex = pt.cex) #handle the line overlaps
  points(x[pt_rand_order], y[pt_rand_order], 
         pch = pch, col = pt.col[valid_inds][pt_rand_order], 
         cex = pt.cex)
  
  #upper margin Spearman corr
  rho <- cor.test(x = x, y = y, method = 'spearman')
  rho_pval <- strsplit(formatC(rho$p.value, format = "e", digits = 1), "e")[[1]]
  print(summary(rho))
  rholab <- c(paste0("Spearman's $\\rho$ = ", 
                   round(rho$estimate, 2), ","), 
                   paste0("p-value = ", 
                   rho_pval[1], " × 10$^{", as.integer(rho_pval[2]), "}$"))
  text3(x = par("usr")[1] + diff(par("usr")[1:2])/100, 
        y = par("usr")[4] - diff(par("usr")[3:4])/100, pos = c(1,4),
        labels = rholab[1], cex = 1.25)
  text3(x = par("usr")[1] + diff(par("usr")[1:2])/100, 
        y = par("usr")[4] - diff(par("usr")[3:4])/100 -
          strheight(latex2exp::TeX(rholab[2]), cex = 1.25), pos = c(1,4),
        labels = rholab[2], cex = 1.25)
  
}

#snag category colors
{
  category_colors <- RColorBrewer::brewer.pal(5, "Dark2")
  names(category_colors) <- sort(c("Cardiometabolic", "Aging", "Anthropometric", 
                                   "Immune", "Psychiatric-neurologic"))
  cats_leftover <- setdiff(salient.categories, names(category_colors))
  cols_leftover <- setdiff(RColorBrewer::brewer.pal(8, "Dark2"), category_colors)
  n_cols_needed <- length(cats_leftover) - length(cols_leftover)
  if(n_cols_needed > 0){
    extra_cols <- c(cols_leftover, RColorBrewer::brewer.pal(n_cols_needed, "Set2"))
  } else {
    extra_cols <- cols_leftover
  }
  names(extra_cols) <- cats_leftover
  category_colors <- c(category_colors, extra_cols)
}

cairo_pdf(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/figures/figS3_intersect-enrichment_Bayes_vs_freq.pdf"), 
          width = 2200 / 72, 
          height = 2000 / 72, 
          family="Arial Unicode MS", pointsize = 28.5)
par(xpd = NA, mar = c(5,2,4,1), pty="s")
layout(rbind(
  c(1,1,2,2),
  c(1,1,2,2),
  c(3,3,4,6),
  c(3,3,5,6)
), heights = c(1,1,1))

#pairs
trait_x_tissue_ps <- setNames(gsea_output$pval, gsea_output$pathway_name)
# trait_x_tissue_adj_ps <- -log10(setNames(gsea_output$adj.pval, gsea_output$pathway_name))
trait_x_tissue_adj_ps <- -log10(setNames(gsea_output$pval * length(gsea_output$pval), gsea_output$pathway_name))

#from locator()... tried a programmatic solution but not enough time
# temp <- locator(); cat(paste0("piks <- cbind(\n\tx = c(", paste0(temp$x, collapse = ", "), "),\n", "\ty = c(", paste0(temp$y, collapse = ", "), ")\n)"))
#for BH-adjusted values
piks <- cbind(
  x = c(-7.04596421988288, -6.0123057215065, -2.19572049673214, -2.43425707328054, -1.40059857490415, 1.26305986321962, 1.81964520849922, 9.25403517759093, 9.21427908149953),
  y = c(-3.18977144364532, -3.24332748037729, -1.70805442739391, -1.40457021924603, -0.0121132642145909, -0.137077349922541, -0.529821619290384, -2.86843522325345)
)
# pt.pos = c(3,1,3,1,1,3,1), pt.cex = 1.75,
# prop_label_push = c(4,4,4,4,4,9,4)/100

xy <- cbind(logit(celltotalbias), trait_x_tissue_adj_ps[names(celltotalbias)])
dists <- as.matrix(dist(rbind(piks, xy)))[1:nrow(piks), nrow(piks)+1:nrow(xy)]
pts.to.label <- c(apply(dists, 1, which.min), 16, 14)

#plot trait x tissue
plot.scatter(x = logit(celltotalbias), 
             y = trait_x_tissue_adj_ps[names(celltotalbias)],
             cex = 1.5, pch = 19, 
             xlabel = latex2exp::TeX("logit(difference term posterior mass > 0)"),
             ylabel = latex2exp::TeX("-log$_{10}$(adjusted p-value)"),
             axis.lab.dev = 0.0875, lab.cex = 1.5,
             main.title = latex2exp::TeX("pairwise (\\textit{\\textbf{trait}} x \\textit{\\textbf{tissue}}) enrichment comparison"),
             pt.labels = names(celltotalbias),
             # pt.labels = NULL,
             pts.to.label = pts.to.label,
             pt.pos = c(3,1,2,3,2,
                        3,2,1,3,3,
                        4,1), 
             pt.cex = 1.75,
             prop_label_push = c(6,4,4,6,4,
                                 6,4,8,4,13,
                                 6,8,10)/100,
             pt.col = adjustcolor(TISSUE_COLORS[gsub("\\..*", "", names(celltotalbias))], 0.5))
# text(xy[pts.to.label,1], xy[pts.to.label,2], labels = 1:length(pts.to.label), col = 4, cex = 2)




#trait-wise aggregation, HMP
trait_adj_ps <- -log10(trait_ps / (trait_thresh / 0.05))
# temp <- locator(); cat(paste0("piks <- cbind(\n\tx = c(", paste0(temp$x, collapse = ", "), "),\n", "\ty = c(", paste0(temp$y, collapse = ", "), ")\n)"))
piks <- cbind(
  x = c(-7.35450884554205, -7.32513250990435, -6.26758442694711, -3.47683254136552, -2.8893058286115, 2.10467122979766, 2.42781092181237, 4.36664907390063, 4.83667044410385, 2.95658496329099, 0.0777040707962936),
  y = c(-1.97611823905796, -1.69248666124853, -1.45818231436247, -0.188006118085423, 0.86019753903641, -0.311324195393874, -0.668946619588382, -0.594955773203312, -1.17455073655303, -1.72948208444106, -2.27208162459825)
)
xy <- cbind(logit(colbias), trait_adj_ps[names(colbias)])
dists <- as.matrix(dist(rbind(piks, xy)))[1:nrow(piks), nrow(piks)+1:nrow(xy)]
pts.to.label <- apply(dists, 1, which.min)

plot.scatter(x = logit(colbias), 
             y = trait_adj_ps[names(colbias)],
             cex = 1.5, pch = 19, 
             xlabel = latex2exp::TeX("logit(difference term posterior mass > 0)"),
             ylabel = latex2exp::TeX("-log$_{10}$(adjusted p-value)"),
             axis.lab.dev = 0.0875, lab.cex = 1.5,
             main.title = latex2exp::TeX("\\textit{\\textbf{trait}} enrichment comparison (HMP)"),
             pt.labels = names(colbias),
             # pt.labels = NULL,
             pts.to.label = pts.to.label,
             pt.pos = c(1,3,3,2,1,2,3,3,1,1,1), pt.cex = 1.75,
             prop_label_push = c(4,8,10,4,4,3,2,8,4,9,4)/100,
             pt.col = adjustcolor(category_colors[trait_cat_key[names(colbias)]], 0.5))
# text(xy[pts.to.label,1], xy[pts.to.label,2], labels = 1:length(pts.to.label), col = 4, cex = 3)

#genes in 3 tissues, HMP for TWAS
# trait_adj_ps <- -log10(setNames(gsea_output_pooled_alt$adj.pval, gsea_output_pooled_alt$trait_name))
trait_adj_ps <- -log10(setNames(gsea_output_pooled_alt$pval * length(gsea_output_pooled_alt$pval), gsea_output_pooled_alt$trait_name))
# temp <- locator(); cat(paste0("piks <- cbind(\n\tx = c(", paste0(temp$x, collapse = ", "), "),\n", "\ty = c(", paste0(temp$y, collapse = ", "), ")\n)"))
piks <- cbind(
  x = c(-7.3179032742577, -7.20766574694757, -6.25962301208041, -5.64229285914366, -2.53359458899785, -1.80602690875097, 2.5814266781924, 3.00032928197091, 4.23498958784442, 4.69798720254698, 3.2428518420532, 2.80190173281267),
  y = c(-2.05367997178105, -1.82772236457882, -1.61158900116798, -1.43475261292275, -0.776528278898839, -0.452328233782585, 1.30621140487831, 0.431853707443563, 0.0192354682046946, -0.206722138997544, -0.648813109610618, -1.52317080704536)
)
xy <- cbind(logit(colbias), trait_adj_ps[names(colbias)])
dists <- as.matrix(dist(rbind(piks, xy)))[1:nrow(piks), nrow(piks)+1:nrow(xy)]
pts.to.label <- apply(dists, 1, which.min)

plot.scatter(x = logit(colbias), 
             y = trait_adj_ps[names(colbias)],
             cex = 1.5, pch = 19, 
             xlabel = latex2exp::TeX("logit(difference term posterior mass > 0)"),
             ylabel = latex2exp::TeX("-log$_{10}$(adjusted p-value)"),
             axis.lab.dev = 0.0875, lab.cex = 1.5,
             main.title = latex2exp::TeX("\\textit{\\textbf{trait}} enrichment comparison (multi-tissue)"),
             pt.labels = names(colbias),
             # pt.labels = NULL,
             pts.to.label = pts.to.label,
             pt.pos = c(1,3,3,3,3,
                        3,2,2,3,3,
                        2,4), pt.cex = 1.75,
             prop_label_push = c(4,11,11,13,8,
                                 8,4,4,10,6,
                                 2,4)/100,
             pt.col = adjustcolor(category_colors[trait_cat_key[names(colbias)]], 0.5))
# text(xy[pts.to.label,1], xy[pts.to.label,2], labels = 1:length(pts.to.label), col = 4, cex = 3)

#change size for smaller plots
par(xpd = NA, mar = c(4,1,3,0), pty="s")

#tissue aggregation
tissue_adj_ps <- -log10(tissue_ps / (tissue_thresh / 0.05))
pts.to.label <- 1:length(tissue_adj_ps)
xy <- cbind(x = logit(rowbias), 
        y = tissue_adj_ps[names(rowbias)])
plot.scatter(x = logit(rowbias), 
             y = tissue_adj_ps[names(rowbias)],
             cex = 1.5, pch = 19, 
             xlabel = latex2exp::TeX("logit(difference term posterior mass > 0)"),
             ylabel = latex2exp::TeX("-log$_{10}$(adjusted p-value)"),
             axis.lab.dev = 0.225, lab.cex = 1.5,
             main.title = latex2exp::TeX("\\textit{\\textbf{tissue}} enrichment comparison (HMP)"),
             pt.labels = tolower(names(rowbias)),
             # pt.labels = NULL,
             pts.to.label = pts.to.label,
             pt.pos = c(3,2,1,1,2,
                        1,2,3,1,4,
                        3,1,3,1,1), 
             pt.cex = 1.75,
             prop_label_push = c(2.5,12,4,14,4,
                                 4,4,25,8,4,
                                 4,17,4,8,4)/100,
             pt.col = adjustcolor(TISSUE_COLORS[names(rowbias)], 0.5),
             prop_rect_expand = 0.1)
# text(xy[pts.to.label,1], xy[pts.to.label,2], labels = 1:length(pts.to.label), col = 4, cex = 1)

#trait category aggregation
trait_category_adj_ps <- -log10(trait_category_ps / (trait_category_thresh / 0.05))
pts.to.label <- 1:length(trait_category_adj_ps)
xy <- cbind(x = logit(colcatbias),
            y = trait_category_adj_ps[names(colcatbias)])
plot.scatter(x = logit(colcatbias),
             y = trait_category_adj_ps[names(colcatbias)],
             cex = 1.5, pch = 19,
             xlabel = latex2exp::TeX("logit(difference term posterior mass > 0)"),
             ylabel = latex2exp::TeX("-log$_{10}$(adjusted p-value)"),
             axis.lab.dev = 0.225, lab.cex = 1.5,
             main.title = latex2exp::TeX("\\textit{\\textbf{trait category}} enrichment comparison (HMP)"),
             pt.labels = gsub("( |\\-).*", "", names(colcatbias)),
             # pt.labels = NULL,
             pts.to.label = pts.to.label,
             pt.pos = c(1,2,2,4,3,
                        2,2,1,1,4,
                        1,3), 
             pt.cex = 1.75,
             prop_label_push = c(6,6,15,6,6,
                                 6,6,3,6,6,
                                 3,6)/100,
             pt.col = adjustcolor(category_colors[names(colcatbias)], 0.5),
             prop_rect_expand = 0.1)
# text(xy[pts.to.label,1], xy[pts.to.label,2], labels = 1:length(pts.to.label), col = 4, cex = 1)

#legends
legend(x = par("usr")[2] + diff(par("usr")[1:2]) / 2, 
       y = par("usr")[3] + diff(par("usr")[3:4]) * 1.1, 
       legend = names(category_colors), pch = 21, 
       col = 1, pt.bg = adjustcolor(category_colors, 0.5), 
       cex = 1.5, pt.cex = 1.75, bty = "n", title = latex2exp::TeX("\\textbf{Trait Categories}"))
legend(x = par("usr")[2] + diff(par("usr")[1:2]) / 2, 
       y = par("usr")[3] + diff(par("usr")[3:4]) * 2.55, 
       legend = tissues, pch = 21, ncol = 2,
       col = 1, pt.bg = adjustcolor(TISSUE_COLORS[tissues], 0.5), 
       cex = 1.5, pt.cex = 1.75, bty = "n", title = latex2exp::TeX("\\textbf{Tissues}"))

#figure subpanel labels
fig_label("a)", cex = 2.5, shrinkX = 30.5, shrinkY = 8.5)
fig_label("b)", cex = 2.5, shrinkX = 30.5, shrinkY = 3.5)
fig_label("c)", cex = 2.5, shrinkX = 0.5, shrinkY = 8.5)
fig_label("d)", cex = 2.5, shrinkX = 3.5, shrinkY = 3.75)
fig_label("e)", cex = 2.5, shrinkX = 3.5, shrinkY = 1.25)

# legend(x = xr, y = yt, legend = tissues, pch = 19, col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues], 0.5), cex = 0.55, pt.cex = 1.2)
# segments(x0 = par("usr")[1] + diff(par("usr")[1:2]) / 1.55, 
#          y0 = par("usr")[4] - diff(par("usr")[3:4])/50,
#          x1 = par("usr")[1] + diff(par("usr")[1:2]) / 1.45, 
#          y1 = par("usr")[4] - diff(par("usr")[3:4])/50, 
#          lty = 2, lwd = 2, col = adjustcolor(1,0.75))
# text(x = par("usr")[1] + diff(par("usr")[1:2]) / 1.45, 
#      y = par("usr")[4] - diff(par("usr")[3:4])/45, 
#      pos = 4, cex = 0.75, labels = "1-to-1 line")



dev.off()


#### prop pos hits bayesian model ####
twas_with_hits <- colnames(prop_degs_are_twas)
tissue_code <- MotrpacBicQC::bic_animal_tissue_code
tissue_code <- tissue_code[,c("tissue_name_release", "abbreviation")]
tissue_code <- tissue_code[tissue_code$tissue_name_release != "" & !is.na(tissue_code$abbreviation)]
tissue_code_vec <- setNames(tissue_code$tissue_name_release, tissue_code$abbreviation)
tissue_code_vec <- MotrpacRatTraining6moData::TISSUE_ABBREV_TO_CODE
possible_genes <- intersect(all_orthologs_tested, all_twas_genes_tested)
compatible_twas_genes <- some.twas$gene_name %in% possible_genes
twas_by_tissue <- lapply(setNames(unique(some.twas$tissue), MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[unique(some.twas$tissue)]), function(tiss) {
  print(tiss)
  some.twas[some.twas$tissue == tiss & compatible_twas_genes,]
})

#get "null hypothesis" genes for overall probabilities
rna_dea <- MotrpacRatTraining6mo::combine_da_results(assays = "TRNSCRPT")
rna_dea <- rna_dea[rna_dea$comparison_group == "8w",]
rna_dea$tissue_abbreviation <- rna_dea$tissue
rna_dea_null <- lapply(setNames(unique(rna_dea$tissue_abbreviation), unique(rna_dea$tissue_abbreviation)), function(tiss) {
  print(tiss)
  subset <- rna_dea[rna_dea$tissue_abbreviation == tiss,]
  subset$gene <- gene_map$HUMAN_ORTHOLOG_SYMBOL[match(subset$feature_ID, gene_map$RAT_ENSEMBL_ID)]
  subset <-subset[subset$gene %in% possible_genes,]
  genes_in_sub <- table(subset$feature_ID)
  genes_in_sub <- names(genes_in_sub)[genes_in_sub == 2]
  if(length(genes_in_sub) == 0){return(NULL)}
  subset <- subset[subset$feature_ID %in% genes_in_sub,c("feature_ID", "sex", "logFC")]
  subset <- tidyr::spread(subset, key = "sex", value = "logFC")
  subset <- subset[(sign(subset$female) * sign(subset$male)) == 1,]
  x <- sign(subset$female)
  names(x) <- gene_map$HUMAN_ORTHOLOG_SYMBOL[match(subset$feature_ID, gene_map$RAT_ENSEMBL_ID)]
  return(x)
})
do.call(rbind, lapply(rna_dea_null, function(x) table(x)))
rm(rna_dea)

tissues <- tissue_code$abbreviation[match(names(motrpac_gtex_map), tissue_code$tissue_name_release)]
tissues <- setdiff(tissues, c("HYPOTH", "TESTES", "OVARY"))
adj_pvalue_alpha <- 0.05
data <- lapply(setNames(tissues, tissues), function(tiss){
  print(tiss)
  tws <- twas_by_tissue[[tiss]]
  
  #get focal genesets
  sig_tws <- tws[tws$adj_pvalue < adj_pvalue_alpha,]
  DE_genes_in_Nodes <- node_metadata_list[["8w"]]$human_gene_symbol[node_metadata_list[["8w"]]$tissue == tiss]
  DE_genes_in_Nodes <- DE_genes_in_Nodes[!is.na(DE_genes_in_Nodes)]
  DE_genes_in_Nodes <- unique(DE_genes_in_Nodes[DE_genes_in_Nodes %in% possible_genes])
  
  motr <- rna_dea_null[[tiss]]
  if(is.null(motr) | is.null(tws)){
    return(NULL)
  }
  compatible_genes <- intersect(names(motr), tws$gene_name)
  motr <- motr[compatible_genes]
  tws <- tws[tws$gene_name %in% compatible_genes,]
  
  tws$motrpac_sign <- motr[tws$gene_name]
  tws$twas_sign <- sign(tws$zscore)
  tws <- tws[tws$twas_sign != 0,]
  tws$effect_dir <- tws$twas_sign * tws$motrpac_sign
  output <- as.data.frame(t(sapply(setNames(twas_with_hits, twas_with_hits), function(trait_i) {
    trait_subset <- tws[tws$trait == trait_i, c("gene_name", "effect_dir")]
    trait_subset <- setNames(trait_subset$effect_dir, trait_subset$gene_name)
    twas_genes <- sig_tws$gene_name[sig_tws$trait == trait_i & sig_tws$gene_name %in% all_orthologs_tested]
    twas_genes <- twas_genes[twas_genes %in% names(trait_subset)]
    motr_genes <- DE_genes_in_Nodes[DE_genes_in_Nodes %in% names(trait_subset)]
    c(count_twas = sum(trait_subset[twas_genes] > 0), 
      count_motr = sum(trait_subset[motr_genes] > 0),
      count_inters = sum(trait_subset[intersect(twas_genes, motr_genes)] > 0),
      count_compl = sum(trait_subset[setdiff(motr_genes, twas_genes)] > 0),
      count_compl_all = sum(trait_subset[setdiff(names(trait_subset), intersect(twas_genes, motr_genes))] > 0),
      count_all = sum(trait_subset > 0),
      total_twas = length(twas_genes),
      total_motr = length(motr_genes),
      total_inters = length(intersect(twas_genes, motr_genes)),
      total_compl = length(setdiff(motr_genes, twas_genes)),
      total_compl_all = length(setdiff(names(trait_subset), intersect(twas_genes, motr_genes))),
      total_all = length(trait_subset))
  })))
  output$trait <- rownames(output)
  output$tissue <- tiss
  rownames(output) <- NULL
  output
})
data <- do.call(rbind, data)
fwrite(data, "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/prop_positive_effects_in_DEG-TWAS_intersection.txt")

# plot((data$count_twas / data$total_twas + data$count_motr / data$total_motr)/2, data$count_inters / data$total_inters,
#      main = paste0("prop above line = ", round(mean((data$count_twas / data$total_twas + data$count_motr / data$total_motr)/2 < 
#                                                       (data$count_inters / data$total_inters), na.rm = T), 2)),
#      cex = data$total_inters / max(data$total_inters) * 5, pch = 19, col = adjustcolor(1,0.5))
plot(data$count_compl / data$total_compl, data$count_inters / data$total_inters,
     main = paste0("prop above line = ", round(mean((data$count_inters / data$total_inters > data$count_compl / data$total_compl), na.rm = T), 2)),
     cex = data$total_inters / max(data$total_inters) * 5, pch = 19, col = adjustcolor(1,0.5))
abline(0,1,lwd=2,lty=2,col=adjustcolor(2,0.5))
abline(h=0.5,lwd=2,lty=2,col=adjustcolor(1,0.5))
plot(data$count_all / data$total_all, data$count_inters / data$total_inters,
     main = paste0("prop above line = ", round(mean((data$count_inters / data$total_inters > data$count_all / data$total_all), na.rm = T), 2)),
     cex = data$total_inters / max(data$total_inters) * 5, pch = 19, col = adjustcolor(1,0.5))
abline(0,1,lwd=2,lty=2,col=adjustcolor(2,0.5))
abline(h=0.5,lwd=2,lty=2,col=adjustcolor(1,0.5))

#make slightly fancier plot
pchs <- category_shapes[traitwise_partitions$Category[match(data$trait, traitwise_partitions$Tag)]]
pchs[is.na(pchs)] <- 1
plot(data$count_all / data$total_all, data$count_inters / data$total_inters,
     main = "", 
     pch = pchs,
     cex = data$total_inters / max(data$total_inters) * 5, col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[data$tissue],0.5))
abline(h=0.5,lwd=2,lty=2,col=adjustcolor(1,0.5))
abline(v=0.5,lwd=2,lty=2,col=adjustcolor(1,0.5))


basic_posteriors_masses <- 1 - pbeta(q = 0.5, shape1 = 1 + data$count_inters, shape2 = 1 + data$total_inters - data$count_inters)
sum(basic_posteriors_masses > 0.95)
sum(basic_posteriors_masses < 0.05)
head(cbind(data, basic_posteriors_masses)[order(abs(basic_posteriors_masses - 0.5), decreasing = T),
                                          c("trait", "tissue", "basic_posteriors_masses", "count_inters", "total_inters")], 30)
# table(data_not_null[which(basic_posteriors_masses > 0.95),"trait"])
# table(data_not_null[which(basic_posteriors_masses < 0.05),"trait"])
# 
# hist(sapply(tissues, function(tiss) diff(range(data_null$count_twas[data_null$tissue == tiss] / data_null$total[data_null$tissue == tiss]))))
# hist(sapply(traits, function(trait) diff(range(data_null$count_twas[data_null$trait == trait] / data_null$total[data_null$trait == trait]))))
# hist(sapply(tissues, function(tiss) diff(range(data_null$count_motr[data_null$tissue == tiss] / data_null$total[data_null$tissue == tiss]))))
# hist(sapply(traits, function(trait) diff(range(data_null$count_motr[data_null$trait == trait] / data_null$total[data_null$trait == trait]))))

#alternative prob model
# tissues <- tissues
# traits <- twas_with_hits
# d <- list(count = data$count_inters,
#           total = data$total_inters,
#           count_compl = data$count_compl_all,
#           total_compl = data$total_compl_all,
#           row_index = match(data$tissue, tissues),
#           col_index = match(data$trait, traits),
#           row_n = length(tissues),
#           col_n = length(traits)
# )

# stan_program <- '
# data {
#     int<lower=1> row_n;
#     int<lower=1> col_n;
#     int<lower=1,upper=row_n> row_index[row_n * col_n];
#     int<lower=1,upper=col_n> col_index[row_n * col_n];
#     int<lower=0> count[row_n * col_n];
#     int<lower=0> total[row_n * col_n];
#     int<lower=0> count_compl[row_n * col_n];
#     int<lower=0> total_compl[row_n * col_n];
# }
# transformed data {
#     int<lower=1> n = row_n * col_n;
# }
# parameters {
#     //main parameters
#     real col_mean;
#     real<lower=0> col_sd;
#     vector[col_n] raw_col_logodds;
#     
#     vector[n] raw_mean_logodds;
#     real<lower=0> sd_mean_logodds;
#     
#     vector[n] raw_logodds;
#     vector[n] raw_logodds_compl;
#     real<lower=0> logodds_sd;
#     
#     //biases in deviations terms
#     vector[row_n] raw_row_bias;
#     vector[col_n] raw_col_bias;
# 
#     //multilevel deviation term params
#     real<lower=0> row_bias_sd;
#     real<lower=0> col_bias_sd;
# }
# transformed parameters {
#     //recenter params
#     vector[col_n] col_logodds = raw_col_logodds * col_sd + col_mean;
#     vector[n] mean_logodds = raw_mean_logodds * sd_mean_logodds + col_logodds[col_index];
#     vector[row_n] row_bias = raw_row_bias * row_bias_sd;
#     vector[col_n] col_bias = raw_col_bias * col_bias_sd;
#     vector[n] logodds = raw_logodds * logodds_sd + mean_logodds + 
#                         row_bias[row_index] / 2 + col_bias[col_index] / 2;
#     vector[n] logodds_compl = raw_logodds_compl * logodds_sd + mean_logodds - 
#                               row_bias[row_index] / 2 - col_bias[col_index] / 2;
#     
# }
# model {
#     //high-level cell params
#     col_mean ~ normal(0,2);
#     col_sd ~ std_normal();
#     raw_col_logodds ~ std_normal();
#     raw_mean_logodds ~ std_normal();
#     sd_mean_logodds ~ std_normal();
#     
#     //bias params
#     raw_row_bias ~ std_normal();
#     raw_col_bias ~ std_normal();
#     row_bias_sd ~ std_normal();
#     col_bias_sd ~ std_normal();
#     
#     //cell params
#     raw_logodds ~ std_normal();
#     raw_logodds_compl ~ std_normal();
#     logodds_sd ~ std_normal();
#     
#     //likelihood
#     count ~ binomial_logit(total, logodds);
#     count_compl ~ binomial_logit(total_compl, logodds_compl);
# }
# generated quantities {
#     vector[n] cell_total_prob_bias = inv_logit(logodds) - inv_logit(logodds_compl);
# }
# '
# 
# 
# stan_program <- '
# data {
#     int<lower=1> row_n;
#     int<lower=1> col_n;
#     int<lower=1,upper=row_n> row_index[row_n * col_n];
#     int<lower=1,upper=col_n> col_index[row_n * col_n];
#     int<lower=0> count[row_n * col_n];
#     int<lower=0> total[row_n * col_n];
#     int<lower=0> count_compl[row_n * col_n];
#     int<lower=0> total_compl[row_n * col_n];
# }
# transformed data {
#     int<lower=1> n = row_n * col_n;
# }
# parameters {
#     //main parameters
#     real col_mean;
#     real<lower=0> col_sd;
#     vector[col_n] raw_col_logodds;
#     
#     vector[n] raw_mean_logodds;
#     real<lower=0> sd_mean_logodds;
#     
#     vector[n] raw_logodds;
#     vector[n] raw_logodds_compl;
#     real<lower=0> logodds_sd;
#     
#     //biases in deviations terms
#     vector[col_n] raw_col_bias;
# 
#     //multilevel deviation term params
#     real<lower=0> col_bias_sd;
# }
# transformed parameters {
#     //recenter params
#     vector[col_n] col_logodds = raw_col_logodds * col_sd + col_mean;
#     vector[n] mean_logodds = raw_mean_logodds * sd_mean_logodds + col_logodds[col_index];
#     vector[col_n] col_bias = raw_col_bias * col_bias_sd;
#     vector[n] logodds = raw_logodds * logodds_sd + mean_logodds + col_bias[col_index] / 2;
#     vector[n] logodds_compl = raw_logodds_compl * logodds_sd + mean_logodds - col_bias[col_index] / 2;
#     
# }
# model {
#     //high-level cell params
#     col_mean ~ normal(0,2);
#     col_sd ~ std_normal();
#     raw_col_logodds ~ std_normal();
#     raw_mean_logodds ~ std_normal();
#     sd_mean_logodds ~ std_normal();
#     
#     //bias params
#     raw_col_bias ~ std_normal();
#     col_bias_sd ~ std_normal();
#     
#     //cell params
#     raw_logodds ~ std_normal();
#     raw_logodds_compl ~ std_normal();
#     logodds_sd ~ std_normal();
#     
#     //likelihood
#     count ~ binomial_logit(total, logodds);
#     count_compl ~ binomial_logit(total_compl, logodds_compl);
# }
# generated quantities {
#     vector[n] cell_total_prob_bias = inv_logit(logodds) - inv_logit(logodds_compl);
# }
# '
# 
# stan_program <- '
# data {
#     int<lower=1> row_n;
#     int<lower=1> col_n;
#     int<lower=1,upper=row_n> row_index[row_n * col_n];
#     int<lower=1,upper=col_n> col_index[row_n * col_n];
#     int<lower=0> count[row_n * col_n];
#     int<lower=0> total[row_n * col_n];
#     int<lower=0> count_compl[row_n * col_n];
#     int<lower=0> total_compl[row_n * col_n];
# }
# transformed data {
#     int<lower=1> n = row_n * col_n;
# }
# parameters {
#     //main parameters
#     real col_mean;
#     real<lower=0> col_sd;
#     vector[col_n] raw_col_logodds;
#     
#     vector[n] raw_mean_logodds;
#     real<lower=0> sd_mean_logodds;
#     
#     vector[n] raw_logodds;
#     vector[n] raw_logodds_compl;
#     real<lower=0> logodds_sd;
#     
#     //biases in deviations terms
#     vector[col_n] raw_col_bias;
# 
#     //multilevel deviation term mean params
#     real<lower=0> col_bias_sd;
#     
#     //multilevel deviation term dispersion params
#     vector[col_n] raw_col_logodds_logsd;
#     real<lower=0> logodds_logsd_sd;
#     
# }
# transformed parameters {
#     //recenter params
#     vector[n] col_dispersion_bias = logodds_sd * exp(raw_col_logodds_logsd[col_index] * logodds_logsd_sd);
#     vector[col_n] col_logodds = raw_col_logodds * col_sd + col_mean;
#     vector[n] mean_logodds = raw_mean_logodds * sd_mean_logodds + col_logodds[col_index];
#     vector[col_n] col_bias = raw_col_bias * col_bias_sd;
#     vector[n] logodds = raw_logodds .* col_dispersion_bias + mean_logodds + col_bias[col_index] / 2;
#     vector[n] logodds_compl = raw_logodds_compl .* col_dispersion_bias + mean_logodds - col_bias[col_index] / 2;
#     
# }
# model {
#     //high-level cell params
#     col_mean ~ normal(0,2);
#     col_sd ~ std_normal();
#     raw_col_logodds ~ std_normal();
#     raw_mean_logodds ~ std_normal();
#     sd_mean_logodds ~ std_normal();
#     
#     //bias params
#     raw_col_bias ~ std_normal();
#     col_bias_sd ~ std_normal();
#     
#     //cell params
#     raw_logodds ~ std_normal();
#     raw_logodds_compl ~ std_normal();
#     logodds_sd ~ std_normal();
#     logodds_logsd_sd ~ std_normal();
#     raw_col_logodds_logsd ~ std_normal();
#     
#     //likelihood
#     count ~ binomial_logit(total, logodds);
#     count_compl ~ binomial_logit(total_compl, logodds_compl);
# }
# generated quantities {
#     vector[n] cell_total_prob_bias = inv_logit(logodds) - inv_logit(logodds_compl);
# }
# '
# 
# 
# stan_program <- '
# data {
#     int<lower=1> row_n;
#     int<lower=1> col_n;
#     int<lower=1,upper=row_n> row_index[row_n * col_n];
#     int<lower=1,upper=col_n> col_index[row_n * col_n];
#     int<lower=0> count[row_n * col_n];
#     int<lower=0> total[row_n * col_n];
# }
# transformed data {
#     int<lower=1> n = row_n * col_n;
# }
# parameters {
#     real<lower=0> col_sd_sd;
#     vector<lower=0>[col_n] col_sd;
#     vector[n] raw_logodds;
# }
# transformed parameters {
#     vector[n] logodds = raw_logodds .* col_sd[col_index] * col_sd_sd + 0.5;
# }
# model {
#     col_sd_sd ~ std_normal();
#     col_sd ~ std_normal();
#     raw_logodds ~ std_normal();
#     
#     count ~ binomial_logit(total, logodds);
# }
# generated quantities {
#     vector[n] cell_total_prob_bias = inv_logit(logodds) - 0.5;
# }
# '
# 

# stan_program <- '
# data {
#     int<lower=1> row_n;
#     int<lower=1> col_n;
#     int<lower=1,upper=row_n> row_index[row_n * col_n];
#     int<lower=1,upper=col_n> col_index[row_n * col_n];
#     int<lower=0> count[row_n * col_n];
#     int<lower=0> total[row_n * col_n];
#     int<lower=0> count_compl[row_n * col_n];
#     int<lower=0> total_compl[row_n * col_n];
# }
# transformed data {
#     int<lower=1> n = row_n * col_n;
# }
# parameters {
#     real<lower=0> col_sd_sd;
#     vector<lower=0>[col_n] col_sd;
#     vector[n] raw_logodds;
#     vector[col_n] raw_col_bias;
#     real<lower=0> col_bias_sd;
#     
#     real col_mean_compl;
#     real<lower=0> col_sd_compl;
#     vector[col_n] raw_col_means_compl;
#     vector<lower=0>[col_n] col_sds_compl;
#     real<lower=0> col_sds_compl_sd;
#     vector[n] raw_logodds_compl;
# }
# transformed parameters {
#     vector[col_n] col_means_compl = raw_col_means_compl * col_sd_compl + col_mean_compl;
#     vector[n] logodds_compl = raw_logodds_compl .* col_sds_compl[col_index] * col_sds_compl_sd + col_means_compl[col_index];
#     vector[col_n] col_bias = raw_col_bias * col_bias_sd;
#     vector[n] logodds = raw_logodds .* col_sd[col_index] * col_sd_sd + logodds_compl + col_bias[col_index];
# }
# model {
#     //compl params
#     col_mean_compl ~ normal(0,2);
#     col_sd_compl ~ std_normal();
#     raw_col_means_compl ~ std_normal();
#     col_sds_compl ~ std_normal();
#     col_sds_compl_sd ~ std_normal();
#     raw_logodds_compl ~ std_normal();
#  
#     col_sd_sd ~ std_normal();
#     col_sd ~ std_normal();
#     raw_logodds ~ std_normal();
#     raw_col_bias ~ std_normal();
#     col_bias_sd ~ std_normal();
#     
#     count ~ binomial_logit(total, logodds);
#     count_compl ~ binomial_logit(total_compl, logodds_compl);
# }
# generated quantities {
#     vector[n] cell_total_prob_bias = inv_logit(logodds) - inv_logit(logodds_compl);
# }
# '

#incorporate a multivariate-normal bias term?
# tissues <- tissues
# traits <- twas_with_hits
# trait_cats <- salient.categories
# d <- list(count = data$count_inters,
#           total = data$total_inters,
#           count_compl = data$count_compl_all,
#           total_compl = data$total_compl_all,
#           row_index = match(data$tissue, tissues),
#           col_index = match(data$trait, traits),
#           colcat_index = match(traitwise_partitions$Category[match(traits, traitwise_partitions$Tag)] , trait_cats),
#           row_n = length(tissues),
#           col_n = length(traits),
#           colcat_n = length(trait_cats)
# )


# load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/est_gcor_mat.RData")
# tissues <- tissues
# traits <- twas_with_hits
# # traits <- twas_with_hits[twas_with_hits %in% trait_categories$Tag[trait_categories$Category == "Cardiometabolic"]]
# gcor <- gcor_mat[traits, traits]
# gcor[is.na(gcor)] <- 0
# gcor <- as.matrix(Matrix::nearPD(gcor, corr = T)$mat)
# L <- t(chol(gcor))
# subdata <- data[data$trait %in% traits,]
# d <- list(count = subdata$count_inters,
#           total = subdata$total_inters,
#           count_compl = subdata$count_compl_all,
#           total_compl = subdata$total_compl_all,
#           row_index = match(subdata$tissue, tissues),
#           col_index = match(subdata$trait, traits),
#           row_n = length(tissues),
#           col_n = length(traits),
#           L = L
# )
# 
# base = "deviation_from_expected_propposhits"
# stan_program <- '
# data {
#     int<lower=1> row_n;
#     int<lower=1> col_n;
#     int<lower=1,upper=row_n> row_index[row_n * col_n];
#     int<lower=1,upper=col_n> col_index[row_n * col_n];
#     int<lower=0> count[row_n * col_n];
#     int<lower=0> total[row_n * col_n];
#     int<lower=0> count_compl[row_n * col_n];
#     int<lower=0> total_compl[row_n * col_n];
#     cholesky_factor_corr[col_n] L;
# }
# transformed data {
#     int<lower=1> n = row_n * col_n;
# }
# parameters {
#     real<lower=0> col_sd_sd;
#     vector<lower=0>[col_n] raw_col_sd;
#     vector[n] raw_logodds;
#     vector[col_n] raw_col_bias;
#     real<lower=0> col_bias_sd;
# 
#     real col_mean_compl;
#     real<lower=0> col_sd_compl;
#     vector[col_n] raw_col_means_compl;
#     vector<lower=0>[col_n] col_sds_compl;
#     real<lower=0> col_sds_compl_sd;
#     vector[n] raw_logodds_compl;
# }
# transformed parameters {
#     vector[col_n] col_means_compl = raw_col_means_compl * col_sd_compl + col_mean_compl;
#     vector[n] logodds_compl = raw_logodds_compl .* col_sds_compl[col_index] * col_sds_compl_sd + col_means_compl[col_index];
#     vector[col_n] col_bias = L * raw_col_bias * col_bias_sd;
#     vector[n] logodds = raw_logodds .* raw_col_sd[col_index] * col_sd_sd + logodds_compl + col_bias[col_index];
# }
# model {
#     //compl params
#     col_mean_compl ~ normal(0,2);
#     col_sd_compl ~ std_normal();
#     raw_col_means_compl ~ std_normal();
#     col_sds_compl ~ std_normal();
#     col_sds_compl_sd ~ std_normal();
#     raw_logodds_compl ~ std_normal();
# 
#     col_sd_sd ~ std_normal();
#     raw_col_sd ~ std_normal();
#     raw_logodds ~ std_normal();
#     raw_col_bias ~ std_normal();
#     col_bias_sd ~ std_normal();
# 
#     count ~ binomial_logit(total, logodds);
#     count_compl ~ binomial_logit(total_compl, logodds_compl);
# }
# generated quantities {
#     vector[n] cell_total_prob_bias = inv_logit(logodds) - inv_logit(logodds_compl);
# }
# '

# # now let's try the split the difference approach again?!
# load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/est_gcor_mat.RData")
# tissues <- tissues
# traits <- twas_with_hits
# # traits <- twas_with_hits[twas_with_hits %in% trait_categories$Tag[trait_categories$Category == "Cardiometabolic"]]
# gcor <- gcor_mat[traits, traits]
# gcor[is.na(gcor)] <- 0
# gcor <- as.matrix(Matrix::nearPD(gcor, corr = T)$mat)
# L <- t(chol(gcor))
# subdata <- data[data$trait %in% traits,]
# d <- list(count_focal = subdata$count_inters,
#           total_focal = subdata$total_inters,
#           count_compl = subdata$count_compl_all,
#           total_compl = subdata$total_compl_all,
#           row_index = match(subdata$tissue, tissues),
#           col_index = match(subdata$trait, traits),
#           row_n = length(tissues),
#           col_n = length(traits),
#           L = L
# )
# 
# base = "deviation_from_expected_propposhits_splitthediff"
# stan_program <- '
# data {
#     int<lower=1> row_n;
#     int<lower=1> col_n;
#     int<lower=1,upper=row_n> row_index[row_n * col_n];
#     int<lower=1,upper=col_n> col_index[row_n * col_n];
#     int<lower=0> count_focal[row_n * col_n];
#     int<lower=0> total_focal[row_n * col_n];
#     int<lower=0> count_compl[row_n * col_n];
#     int<lower=0> total_compl[row_n * col_n];
#     cholesky_factor_corr[col_n] L;
# }
# transformed data {
#     int<lower=1> n = row_n * col_n;
# }
# parameters {
#     real col_mean;
#     real<lower=0> col_sd;
#     vector[col_n] raw_col_means;
#     vector<lower=0>[col_n] raw_col_sds;
#     real<lower=0> col_sds_sd;
#     vector[n] raw_logodds;
# 
#     real col_bias_mean;
#     real<lower=0> col_bias_sd;
#     vector[col_n] raw_col_bias_means;
#     vector<lower=0>[col_n] raw_col_bias_sds;
#     real<lower=0> col_bias_sds_sd;
#     vector[n] raw_logodds_bias;
# 
# }
# transformed parameters {
#     vector[col_n] col_means = raw_col_means * col_sd + col_mean;
#     vector[n] logodds = raw_logodds .* raw_col_sds[col_index] * col_sds_sd + col_means[col_index];
# 
#     vector[col_n] col_bias_means = raw_col_bias_means * col_bias_sd + col_bias_mean;
#     vector[n] logodds_bias = raw_logodds_bias .* raw_col_bias_sds[col_index] * col_bias_sds_sd + col_bias_means[col_index];
# 
#     vector[n] logodds_focal = logodds + logodds_bias / 2;
#     vector[n] logodds_compl = logodds - logodds_bias / 2;
# }
# model {
#     //priors
#     col_mean ~ normal(0,2);
#     col_sd ~ std_normal();
#     raw_col_means ~ std_normal();
#     raw_col_sds ~ std_normal();
#     col_sds_sd ~ std_normal();
#     raw_logodds ~ std_normal();
# 
#     col_bias_mean ~ std_normal();
#     col_bias_sd ~ std_normal();
#     raw_col_bias_means ~ std_normal();
#     raw_col_bias_sds ~ std_normal();
#     col_bias_sds_sd ~ std_normal();
#     raw_logodds_bias ~ std_normal();
# 
#     count_focal ~ binomial_logit(total_focal, logodds_focal);
#     count_compl ~ binomial_logit(total_compl, logodds_compl);
# }
# generated quantities {
#     vector[n] cell_total_prob_bias = inv_logit(logodds_focal) - inv_logit(logodds_compl);
# }
# '

# #"split the difference" approach to comparing complement to focal
# stan_program <- '
# data {
#     int<lower=1> row_n;
#     int<lower=1> col_n;
#     int<lower=1,upper=row_n> row_index[row_n * col_n];
#     int<lower=1,upper=col_n> col_index[row_n * col_n];
#     int<lower=0> count[row_n * col_n];
#     int<lower=0> total[row_n * col_n];
#     int<lower=0> count_compl[row_n * col_n];
#     int<lower=0> total_compl[row_n * col_n];
#     cholesky_factor_corr[col_n] L;
# }
# transformed data {
#     int<lower=1> n = row_n * col_n;
# }
# parameters {
#     //overall and column means
#     real col_mean;
#     real<lower=0> col_sd;
#     vector[col_n] raw_col_means;
#     real<lower=0> logodds_sd;
#     vector[n] raw_logodds;
# 
#     //col bias terms
#     real col_bias_mean;
#     real<lower=0> col_bias_sd;
#     vector[col_n] raw_col_bias_means;
#     vector<lower=0>[col_n] raw_logodds_bias_sd;
#     real<lower=0> logodds_bias_sd_sd;
#     vector[n] raw_bias_logodds;
# 
# }
# transformed parameters {
#     vector[col_n] col_means = raw_col_means * col_sd + col_mean;
#     vector[n] logodds = raw_logodds * logodds_sd + col_means[col_index];
# 
#     vector[col_n] col_bias_means = L * raw_col_bias_means * col_bias_sd + col_bias_mean;
#     vector[n] bias_logodds = raw_bias_logodds .* raw_logodds_bias_sd[col_index] * logodds_bias_sd_sd + col_bias_means[col_index];
# 
#     vector[n] logodds_focal = logodds + bias_logodds / 2;
#     vector[n] logodds_compl = logodds - bias_logodds / 2;
# }
# model {
#     //priors
#     col_mean ~ normal(0,2);
#     col_sd ~ std_normal();
#     raw_col_means ~ std_normal();
#     logodds_sd ~ std_normal();
#     raw_logodds ~ std_normal();
# 
#     col_bias_mean ~ std_normal();
#     col_bias_sd ~ std_normal();
#     raw_col_bias_means ~ std_normal();
#     raw_bias_logodds ~ std_normal();
#     raw_logodds_bias_sd ~ std_normal();
#     logodds_bias_sd_sd ~ std_normal();
# 
#     //likelihood
#     count ~ binomial_logit(total, logodds_focal);
#     count_compl ~ binomial_logit(total_compl, logodds_compl);
# }
# generated quantities {
#     vector[n] cell_total_prob_bias = inv_logit(logodds_focal) - inv_logit(logodds_compl);
# }
# '
# 
# #a model where we try to partially pool across trait categories & incorporate goodness
# 
# load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/est_gcor_mat.RData")
# tissues <- tissues
# traits <- twas_with_hits[twas_with_hits %in% trait_categories$Tag[trait_categories$new_Phenotype %in% names(trait_goodness)[trait_goodness != 0]]]
# gcor <- gcor_mat[traits, traits]
# gcor[is.na(gcor)] <- 0
# gcor <- as.matrix(Matrix::nearPD(gcor, corr = T)$mat)
# L <- t(chol(gcor))
# 
# #generate subdata matrix
# subdata <- data[data$trait %in% traits,]
# harmful_traits <- trait_categories$Tag[trait_categories$new_Phenotype %in% names(trait_goodness)[trait_goodness == -1]]
# subdata_harmful <- subdata[subdata$trait %in% harmful_traits,]
# subdata_helpful <- subdata[!(subdata$trait %in% harmful_traits),]
# subdata_harmful$count_inters <- subdata_harmful$total_inters - subdata_harmful$count_inters
# subdata_harmful$count_compl_all <- subdata_harmful$total_compl_all - subdata_harmful$count_compl_all
# subdata <- rbind(subdata_harmful, subdata_helpful)
# 
# trait_cats <- salient.categories
# 
# d <- list(count = subdata$count_inters,
#           total = subdata$total_inters,
#           count_compl = subdata$count_compl_all,
#           total_compl = subdata$total_compl_all,
#           row_index = match(subdata$tissue, tissues),
#           col_index = match(subdata$trait, traits),
#           row_n = length(tissues),
#           col_n = length(traits),
#           L = L,
#           colcat_index = match(traitwise_partitions$Category[match(traits, traitwise_partitions$Tag)] , trait_cats),
#           colcat_n = length(trait_cats)
# )
# 
# 
# stan_program <- '
# data {
#     int<lower=1> row_n;
#     int<lower=1> col_n;
#     int<lower=1> colcat_n;
#     int<lower=1,upper=row_n> row_index[row_n * col_n];
#     int<lower=1,upper=col_n> col_index[row_n * col_n];
#     int<lower=0> count[row_n * col_n];
#     int<lower=0> total[row_n * col_n];
#     int<lower=0> count_compl[row_n * col_n];
#     int<lower=0> total_compl[row_n * col_n];
#     cholesky_factor_corr[col_n] L;
#     int<lower=1,upper=colcat_n> colcat_index[col_n];
# }
# transformed data {
#     int<lower=1> n = row_n * col_n;
# }
# parameters {
#     //overall and column means
#     real col_mean;
#     real<lower=0> col_sd;
#     vector[col_n] raw_col_means;
#     real<lower=0> logodds_sd;
#     vector[n] raw_logodds;
# 
#     //col bias terms
#     real col_bias_mean;
#     real<lower=0> col_bias_sd;
#     vector[colcat_n] raw_colcat_bias_means;
#     real<lower=0> colcat_bias_sd;
#     vector[col_n] raw_col_bias_means;
#     vector<lower=0>[col_n] raw_logodds_bias_sd;
#     real<lower=0> logodds_bias_sd_sd;
#     vector[n] raw_bias_logodds;
# 
# }
# transformed parameters {
#     vector[col_n] col_means = raw_col_means * col_sd + col_mean;
#     vector[n] logodds = raw_logodds * logodds_sd + col_means[col_index];
# 
#     vector[colcat_n] colcat_bias_means = raw_colcat_bias_means * colcat_bias_sd + col_bias_mean;
#     vector[col_n] col_bias_means = L * raw_col_bias_means * col_bias_sd + colcat_bias_means[colcat_index];
#     vector[n] bias_logodds = raw_bias_logodds .* raw_logodds_bias_sd[col_index] * logodds_bias_sd_sd + col_bias_means[col_index];
# 
#     vector[n] logodds_focal = logodds + bias_logodds / 2;
#     vector[n] logodds_compl = logodds - bias_logodds / 2;
# }
# model {
#     //priors
#     col_mean ~ normal(0,2);
#     col_sd ~ std_normal();
#     raw_col_means ~ std_normal();
#     logodds_sd ~ std_normal();
#     raw_logodds ~ std_normal();
#     
#     //bias priors
#     col_bias_mean ~ std_normal();
#     col_bias_sd ~ std_normal();
#     raw_colcat_bias_means ~ std_normal();
#     colcat_bias_sd ~ std_normal();
#     raw_col_bias_means ~ std_normal();
#     raw_bias_logodds ~ std_normal();
#     raw_logodds_bias_sd ~ std_normal();
#     logodds_bias_sd_sd ~ std_normal();
# 
#     //likelihood
#     count ~ binomial_logit(total, logodds_focal);
#     count_compl ~ binomial_logit(total_compl, logodds_compl);
# }
# generated quantities {
#     vector[n] cell_total_prob_bias = inv_logit(logodds) - inv_logit(logodds_compl);
# }
# '
# 
# #also incorporate a tissue bias?
# stan_program <- '
# data {
#     int<lower=1> row_n;
#     int<lower=1> col_n;
#     int<lower=1> colcat_n;
#     int<lower=1,upper=row_n> row_index[row_n * col_n];
#     int<lower=1,upper=col_n> col_index[row_n * col_n];
#     int<lower=0> count[row_n * col_n];
#     int<lower=0> total[row_n * col_n];
#     int<lower=0> count_compl[row_n * col_n];
#     int<lower=0> total_compl[row_n * col_n];
#     cholesky_factor_corr[col_n] L;
#     int<lower=1,upper=colcat_n> colcat_index[col_n];
# }
# transformed data {
#     int<lower=1> n = row_n * col_n;
# }
# parameters {
#     //overall and column means
#     real col_mean;
#     real<lower=0> col_sd;
#     vector[col_n] raw_col_means;
#     real<lower=0> logodds_sd;
#     vector[n] raw_logodds;
# 
#     //col bias terms
#     real col_bias_mean;
#     real<lower=0> col_bias_sd;
#     vector[colcat_n] raw_colcat_bias_means;
#     real<lower=0> colcat_bias_sd;
#     vector[col_n] raw_col_bias_means;
#     vector<lower=0>[col_n] raw_logodds_bias_sd;
#     real<lower=0> logodds_bias_sd_sd;
#     vector[n] raw_bias_logodds;
#     
#     //tissue bias terms
#     real<lower=0> row_bias_sd;
#     vector[row_n] raw_row_bias_means;
# 
# }
# transformed parameters {
#     vector[col_n] col_means = raw_col_means * col_sd + col_mean;
#     vector[n] logodds = raw_logodds * logodds_sd + col_means[col_index];
# 
#     vector[colcat_n] colcat_bias_means = raw_colcat_bias_means * colcat_bias_sd + col_bias_mean;
#     vector[col_n] col_bias_means = L * raw_col_bias_means * col_bias_sd + colcat_bias_means[colcat_index];
#     vector[row_n] row_bias_means = raw_row_bias_means * row_bias_sd;
#     vector[n] bias_logodds = raw_bias_logodds .* raw_logodds_bias_sd[col_index] * logodds_bias_sd_sd + col_bias_means[col_index] + row_bias_means[row_index];
# 
#     vector[n] logodds_focal = logodds + bias_logodds / 2;
#     vector[n] logodds_compl = logodds - bias_logodds / 2;
# }
# model {
#     //priors
#     col_mean ~ normal(0,2);
#     col_sd ~ std_normal();
#     raw_col_means ~ std_normal();
#     logodds_sd ~ std_normal();
#     raw_logodds ~ std_normal();
#     
#     //bias priors
#     col_bias_mean ~ std_normal();
#     col_bias_sd ~ std_normal();
#     raw_colcat_bias_means ~ std_normal();
#     colcat_bias_sd ~ std_normal();
#     raw_col_bias_means ~ std_normal();
#     raw_bias_logodds ~ std_normal();
#     raw_logodds_bias_sd ~ std_normal();
#     logodds_bias_sd_sd ~ std_normal();
#     
#     row_bias_sd ~ std_normal();
#     raw_row_bias_means ~ std_normal();
# 
#     //likelihood
#     count ~ binomial_logit(total, logodds_focal);
#     count_compl ~ binomial_logit(total_compl, logodds_compl);
# }
# generated quantities {
#     vector[n] cell_total_prob_bias = inv_logit(logodds) - inv_logit(logodds_compl);
# }
# '


#let's try a more basic flavor of model
# load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/est_gcor_mat.RData")
# tissues <- tissues
# traits <- twas_with_hits
# trait_cats <- salient.categories
# # traits <- twas_with_hits[twas_with_hits %in% trait_categories$Tag[trait_categories$Category == "Cardiometabolic"]]
# gcor <- gcor_mat[traits, traits]
# gcor[is.na(gcor)] <- 0
# gcor <- as.matrix(Matrix::nearPD(gcor, corr = T)$mat)
# L <- t(chol(gcor))
# subdata <- data[data$trait %in% traits,]
# d <- list(count = subdata$count_inters,
#           total = subdata$total_inters,
#           row_index = match(subdata$tissue, tissues),
#           col_index = match(subdata$trait, traits),
#           colcat_index = match(traitwise_partitions$Category[match(traits, traitwise_partitions$Tag)], trait_cats),
#           row_n = length(tissues),
#           col_n = length(traits),
#           colcat_n = length(trait_cats),
#           L = L
# )
# 
# 
# base = "deviation_from_half"
# stan_program <- '
# data {
#     int<lower=1> row_n;
#     int<lower=1> col_n;
#     int<lower=1> colcat_n;
#     int<lower=1,upper=row_n> row_index[row_n * col_n];
#     int<lower=1,upper=col_n> col_index[row_n * col_n];
#     int<lower=1,upper=colcat_n> colcat_index[col_n];
#     int<lower=0> count[row_n * col_n];
#     int<lower=0> total[row_n * col_n];
#     cholesky_factor_corr[col_n] L;
# }
# transformed data {
#     int<lower=1> n = row_n * col_n;
# }
# parameters {
#     vector[col_n] raw_col_bias;
#     real<lower=0> col_bias_sd;
#     vector[colcat_n] raw_colcat_logsd;
#     real<lower=0> colcat_logsd_sd;
#     
#     vector[n] raw_logodds;
#     real<lower=0> col_sd_sd;
#     vector[col_n] raw_col_logsd;
#     real<lower=0> col_logsd_sd;
# }
# transformed parameters {
#     vector[col_n] col_bias = raw_col_bias * col_bias_sd.* exp(raw_colcat_logsd[colcat_index] * colcat_logsd_sd);
#     vector[n] logodds = raw_logodds * col_sd_sd .* exp(raw_col_logsd[col_index] * col_logsd_sd) + 0 + col_bias[col_index];
# }
# model {
#     raw_col_bias ~ std_normal();
#     col_bias_sd ~ std_normal();
#     raw_colcat_logsd ~ std_normal();
#     colcat_logsd_sd ~ std_normal();
#     
#     raw_logodds ~ std_normal();
#     col_sd_sd ~ std_normal();
#     raw_col_logsd ~ std_normal();
#     col_logsd_sd ~ std_normal();
# 
#     count ~ binomial_logit(total, logodds);
# }
# generated quantities {
#     vector[n] cell_total_prob_bias = inv_logit(logodds) - 0.5;
# }
# '
# 
# #incorporate non-independence via cholesky factor
# stan_program <- '
# data {
#     int<lower=1> row_n;
#     int<lower=1> col_n;
#     int<lower=1> colcat_n;
#     int<lower=1,upper=row_n> row_index[row_n * col_n];
#     int<lower=1,upper=col_n> col_index[row_n * col_n];
#     int<lower=1,upper=colcat_n> colcat_index[col_n];
#     int<lower=0> count[row_n * col_n];
#     int<lower=0> total[row_n * col_n];
#     cholesky_factor_corr[col_n] L;
# }
# transformed data {
#     int<lower=1> n = row_n * col_n;
# }
# parameters {
#     vector[col_n] raw_col_bias;
#     real<lower=0> col_bias_sd;
#     vector[colcat_n] raw_colcat_logsd;
#     real<lower=0> colcat_logsd_sd;
#     
#     vector[n] raw_logodds;
#     real<lower=0> col_sd_sd;
#     vector[col_n] raw_col_logsd;
#     real<lower=0> col_logsd_sd;
# }
# transformed parameters {
#     vector[col_n] col_bias = L * raw_col_bias * col_bias_sd.* exp(raw_colcat_logsd[colcat_index] * colcat_logsd_sd);
#     vector[n] logodds = raw_logodds * col_sd_sd .* exp(raw_col_logsd[col_index] * col_logsd_sd) + 0 + col_bias[col_index];
# }
# model {
#     raw_col_bias ~ std_normal();
#     col_bias_sd ~ std_normal();
#     raw_colcat_logsd ~ std_normal();
#     colcat_logsd_sd ~ std_normal();
#     
#     raw_logodds ~ std_normal();
#     col_sd_sd ~ std_normal();
#     raw_col_logsd ~ std_normal();
#     col_logsd_sd ~ std_normal();
# 
#     count ~ binomial_logit(total, logodds);
# }
# generated quantities {
#     vector[n] cell_total_prob_bias = inv_logit(logodds) - 0.5;
# }
# '

#incorporate non-independence, but allow for independence too
load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/est_gcor_mat.RData")
tissues <- tissues
traits <- twas_with_hits
trait_cats <- salient.categories
# traits <- twas_with_hits[twas_with_hits %in% trait_categories$Tag[trait_categories$Category == "Cardiometabolic"]]
gcor <- gcor_mat[traits, traits]
gcor[is.na(gcor)] <- 0
gcor <- as.matrix(Matrix::nearPD(gcor, corr = T)$mat)
L <- t(chol(gcor))
subdata <- data[data$trait %in% traits,]
d <- list(count = subdata$count_inters,
          total = subdata$total_inters,
          row_index = match(subdata$tissue, tissues),
          col_index = match(subdata$trait, traits),
          colcat_index = match(traitwise_partitions$Category[match(traits, traitwise_partitions$Tag)], trait_cats),
          row_n = length(tissues),
          col_n = length(traits),
          colcat_n = length(trait_cats),
          gcor = gcor
)

base = "deviation_from_half"
stan_program <- '
data {
    int<lower=1> row_n;
    int<lower=1> col_n;
    int<lower=1> colcat_n;
    int<lower=1,upper=row_n> row_index[row_n * col_n];
    int<lower=1,upper=col_n> col_index[row_n * col_n];
    int<lower=1,upper=colcat_n> colcat_index[col_n];
    int<lower=0> count[row_n * col_n];
    int<lower=0> total[row_n * col_n];
    corr_matrix[col_n] gcor;
}
transformed data {
    int<lower=1> n = row_n * col_n;
}
parameters {
    real<lower=0, upper=1> prop_gcor;

    vector[col_n] raw_col_bias;
    real<lower=0> col_bias_sd;
    vector[colcat_n] raw_colcat_logsd;
    real<lower=0> colcat_logsd_sd;
    
    vector[n] raw_logodds;
    real<lower=0> col_sd_sd;
    vector[col_n] raw_col_logsd;
    real<lower=0> col_logsd_sd;
}
transformed parameters {
    corr_matrix[col_n] averaged_corrmat = gcor * prop_gcor + diag_matrix(rep_vector(1.0, col_n)) * (1-prop_gcor);
    vector[col_n] col_bias = cholesky_decompose(averaged_corrmat) * raw_col_bias * col_bias_sd.* exp(raw_colcat_logsd[colcat_index] * colcat_logsd_sd);
    vector[n] logodds = raw_logodds * col_sd_sd .* exp(raw_col_logsd[col_index] * col_logsd_sd) + 0 + col_bias[col_index];
}
model {
    raw_col_bias ~ std_normal();
    col_bias_sd ~ normal(0,2);
    raw_colcat_logsd ~ std_normal();
    colcat_logsd_sd ~ std_normal();
    
    raw_logodds ~ std_normal();
    col_sd_sd ~ normal(0,2);
    raw_col_logsd ~ std_normal();
    col_logsd_sd ~ std_normal();

    count ~ binomial_logit(total, logodds);
}
generated quantities {
    vector[n] cell_total_prob_bias = inv_logit(logodds) - 0.5;
}
'

# stan_program <- '
# data {
#     int<lower=1> row_n;
#     int<lower=1> col_n;
#     int<lower=1> colcat_n;
#     int<lower=1,upper=row_n> row_index[row_n * col_n];
#     int<lower=1,upper=col_n> col_index[row_n * col_n];
#     int<lower=1,upper=colcat_n> colcat_index[col_n];
#     int<lower=0> count[row_n * col_n];
#     int<lower=0> total[row_n * col_n];
#     corr_matrix[col_n] gcor;
# }
# transformed data {
#     int<lower=1> n = row_n * col_n;
# }
# parameters {
#     real<lower=0, upper=1> prop_gcor;
# 
#     vector[col_n] raw_col_bias;
#     vector<lower=0>[colcat_n] raw_colcat_bias_sd;
#     real<lower=0> colcat_bias_sd2;
#     
#     vector[n] raw_logodds;
#     vector<lower=0>[col_n] raw_within_col_sd;
#     real<lower=0> within_col_sd2;
# }
# transformed parameters {
#     corr_matrix[col_n] averaged_corrmat = gcor * prop_gcor + diag_matrix(rep_vector(1.0, col_n)) * (1-prop_gcor);
#     vector[col_n] col_bias = cholesky_decompose(averaged_corrmat) * raw_col_bias .* raw_colcat_bias_sd[colcat_index] * colcat_bias_sd2;
#     vector[n] logodds = raw_logodds .* raw_within_col_sd[col_index] * within_col_sd2;
# }
# model {
# 
#     raw_col_bias ~ std_normal();
#     raw_colcat_bias_sd ~ std_normal();
#     colcat_bias_sd2 ~ normal(0,2);
#     
#     raw_logodds ~ std_normal();
#     raw_within_col_sd ~ std_normal();
#     within_col_sd2 ~ normal(0,2);
# 
#     count ~ binomial_logit(total, logodds);
# }
# generated quantities {
#     vector[n] cell_total_prob_bias = inv_logit(logodds) - 0.5;
# }
# '

fit_model <- F
if(fit_model){
  
  #compile model
  if(!exists("curr_stan_program") || stan_program != curr_stan_program){
    curr_stan_program <- stan_program
    f <- write_stan_file(stan_program)
  }
  mod <- cmdstan_model(f)
  
  #write model
  write_stan_file(stan_program, dir = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/mcmc_output/", basename = paste0(base, ifelse(use_all_cats, "_allCats", "")))
  write_stan_json(d, paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/mcmc_output/", paste0(base, ifelse(use_all_cats, "_allCats", "")),".json"))
  
  #fit model
  out <- mod$sample(chains = 4, iter_sampling = 2.5E3, iter_warmup = 2.5E3, data = d, parallel_chains = 4, 
                    adapt_delta = 0.99, refresh = 50, init = 0.1, max_treedepth = 20, thin = 1)
  summ <- out$summary()
  summ[order(summ$ess_bulk),]
  summ[order(summ$rhat, decreasing = T),]
  save(out, file = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/mcmc_output/", paste0(base, ifelse(use_all_cats, "_allCats", "")),".cmdStanR.fit"))
  save(summ, file = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/mcmc_output/", paste0(base, ifelse(use_all_cats, "_allCats", "")),".cmdStanR.summ"))
} else {
  load(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/mcmc_output/", paste0(base, ifelse(use_all_cats, "_allCats", "")),".cmdStanR.fit"))
}

samps <- data.frame(as_draws_df(out$draws()))


# samps[,grep("logit_prop_21\\.", colnames(samps))]
# sum(apply(samps[,grep("logit_prop_21\\.", colnames(samps))], 2, prop_greater_than_0) > 0.95)
# sum(apply(samps[,grep("logit_prop_21\\.", colnames(samps))], 2, prop_greater_than_0) < 0.05)

dev.off()
hist(apply(subset_samps("cell_total_prob_bias", c("raw", "sd"), samps = samps), 2, mean))
hist(apply(subset_samps("col_bias", c("raw", "sd"), samps = samps), 2, mean), breaks = 10)
# hist(samps$col_sd_comp)
# hist(invlogit(apply(subset_samps("col_means_compl", c("raw", "sd"), samps = samps), 2, mean)))
# hist(invlogit(apply(subset_samps("logodds_compl", c("raw", "sd"), samps = samps), 2, mean)))
# hist(invlogit(apply(subset_samps("logodds_focal", c("raw", "sd"), samps = samps), 2, mean)))
hist(invlogit(apply(subset_samps("logodds", c("raw", "sd"), samps = samps), 2, mean)))

hist(apply(subset_samps("logodds", c("raw", "compl", "sd", "mean", "bias"), samps = samps), 2, mean))
prop_greater_than_0 <- function(x) mean(x>0)
sum((apply(subset_samps("cell_total_prob_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)) > 0.95)
sum((apply(subset_samps("cell_total_prob_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)) < 0.05)

subdata[which((apply(subset_samps("cell_total_prob_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)) > 0.9),
        c("count_inters", "total_inters", "tissue", "trait")]
subdata[which((apply(subset_samps("cell_total_prob_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)) <0.1),
        c("count_inters", "total_inters", "tissue", "trait")]


cbind(subdata[order(((apply(subset_samps("cell_total_prob_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)))), 
              c("count_inters", "total_inters", "tissue", "trait")], 
      sort(apply(subset_samps("cell_total_prob_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)))
plot(d$count / d$total, apply(subset_samps("cell_total_prob_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0))
# plot(d$count / d$total, invlogit(apply(subset_samps("logodds_focal", c("raw", "sd", "compl"), samps = samps), 2, mean)), cex = d$total / max(d$total) * 5); abline(0,1,col=2,lty=2)

sum((apply(subset_samps("col_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)) > 0.9)
sum((apply(subset_samps("col_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)) < 0.1)
trait_key <- setNames(trait_categories$new_Phenotype, trait_categories$Tag)

alpha = 0.05
#cell enrichments
subdata[which((apply(subset_samps("cell_total_prob_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)) > 1-alpha), 
        c("count_inters", "total_inters", "tissue", "trait")]
#cell depletions
subdata[which((apply(subset_samps("cell_total_prob_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)) < alpha),
        c("count_inters", "total_inters", "tissue", "trait")]
#multilevel trait enrichments
trait_key[traits[which((apply(subset_samps("col_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)) > 1-alpha)]]
#multilevel trait depletions
trait_key[traits[which((apply(subset_samps("col_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)) < alpha)]]
#category enrichments
cbind(trait_cats[which((apply(subset_samps("colcat_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)) < alpha)],
      apply(subset_samps("colcat_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0))
#gcor proportion
hist(samps$prop_gcor, probability = T)
mean(samps$prop_gcor > 0.9)


#process this into a format for plotting
posterior_mean_probs <- data.frame(trait = data$trait, 
                                   tissue = data$tissue, 
                                   prob = invlogit(apply(subset_samps("logodds", c("raw", "sd", "compl"), samps = samps), 2, mean)))
deg_sigtwas_proportion_posterior_mean <- reshape(posterior_mean_probs, idvar = "tissue", timevar = "trait", direction = "wide")
rownames(deg_sigtwas_proportion_posterior_mean) <- deg_sigtwas_proportion_posterior_mean$tissue
deg_sigtwas_proportion_posterior_mean <- deg_sigtwas_proportion_posterior_mean[,-1]
colnames(deg_sigtwas_proportion_posterior_mean) <- gsub("prob\\.", "", colnames(deg_sigtwas_proportion_posterior_mean))

posterior_mean_prob_bias <- data.frame(trait = data$trait, 
                                       tissue = data$tissue, 
                                       prob = apply(subset_samps("cell_total_prob_bias", c("raw", "sd", "compl"), samps = samps), 2, prop_greater_than_0))
deg_sigtwas_proportion_posterior_prob_bias <- reshape(posterior_mean_prob_bias, idvar = "tissue", timevar = "trait", direction = "wide")
rownames(deg_sigtwas_proportion_posterior_prob_bias) <- deg_sigtwas_proportion_posterior_prob_bias$tissue
deg_sigtwas_proportion_posterior_prob_bias <- deg_sigtwas_proportion_posterior_prob_bias[,-1]
colnames(deg_sigtwas_proportion_posterior_prob_bias) <- gsub("prob\\.", "", colnames(deg_sigtwas_proportion_posterior_prob_bias))

#check to make sure this is consistent with model
# raw_logodds * col_sd_sd .* exp(raw_col_logsd[col_index] * col_logsd_sd) + 0 + col_bias[col_index]
apply(subset_samps("col_bias", c("raw", "sd", "compl"), samps = samps), 2, prop_greater_than_0)
prop_greater_than_0(subset_samps("raw_logodds\\.1\\.", c("sd", "compl"), samps = samps) * subset_samps("col_sd_sd", c("compl"), samps = samps) * 
                      exp(subset_samps("raw_col_logsd\\.1\\.", c("compl"), samps = samps) * subset_samps("col_logsd_sd", c("compl"), samps = samps)) + 
                      subset_samps("col_bias\\.1\\.", c("raw", "sd", "compl"), samps = samps))
posterior_mean_prob_bias[1,]
#YEP!

# trait_means <- invlogit(apply(subset_samps("col_bias", c("raw", "sd", "compl"), samps = samps) + 
#                                 subset_samps("col_means_compl", c("raw", "sd"), samps = samps), 
#                               2, mean))
trait_means <- invlogit(apply(subset_samps("col_bias", c("raw", "sd", "compl"), samps = samps), 
                              2, mean))
names(trait_means) <- traits
trait_bias_probs <- apply(subset_samps("col_bias", c("raw", "sd", "compl"), samps = samps), 2, prop_greater_than_0)
names(trait_bias_probs) <- traits


#### proportion of DEGs & TWAS hits in + & - directions ####
twas_with_hits <- colnames(prop_degs_are_twas)
deg_sigtwas_proportion <- array(NA, dim = c(length(motrpac_gtex_map), length(twas_with_hits), 4, 3, 2), 
                                dimnames = list(names(motrpac_gtex_map), twas_with_hits, paste0(2^(0:3), "w"), c("p", "n", "genes"), c("male", "female")))
nuniq <- function(x) length(unique(x))
tissue_abbr_rev <- MotrpacRatTraining6moData::TISSUE_ABBREV_TO_CODE

load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/relative_effect_sizes_deg_eqtl_list.RData")
trace_8w_backwards <- T
# paths <- data.table::fread("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/feature_repfdr_states_20220117.tsv")
paths <- MotrpacRatTraining6moData::GRAPH_STATES
paths <- paths[paths$ome == "TRNSCRPT",]
paths <- lapply(setNames(unique(paths$tissue),unique(paths$tissue)), function(tiss) paths[paths$tissue == tiss,])

for(tissue_i in names(motrpac_gtex_map)){
  
  print(tissue_i)
  
  sig_twas_by_trait <- lapply(setNames(twas_with_hits, twas_with_hits), 
                              function(trait_i) sig_twas_by_tissue[[tissue_i]][sig_twas_by_tissue[[tissue_i]]$trait == trait_i,])
  sig_twas_by_trait_genes <- lapply(setNames(twas_with_hits, twas_with_hits), 
                                    function(trait_i) sig_twas_by_trait[[trait_i]]$gene_name)
  sig_twas_by_trait_signs <- lapply(setNames(twas_with_hits, twas_with_hits), 
                                    function(trait_i){
                                      signs <- sign(sig_twas_by_trait[[trait_i]]$zscore)
                                      names(signs) <- sig_twas_by_trait_genes[[trait_i]]
                                      signs
                                    })
  
  paths_tiss <- as.data.frame(paths[[MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[tissue_i]]])
  
  for(timepoint in paste0(2^(0:3), "w")){
    
    for(sex in c("male", "female")){
      
      if(trace_8w_backwards){
        timepoint_to_use <- "8w"
      } else {
        timepoint_to_use <- timepoint
      }
      
      #get human gene symbols
      DE_genes_in_Nodes <- node_metadata_list[[timepoint_to_use]]$human_gene_symbol[
        node_metadata_list[[timepoint_to_use]]$tissue == tissue_code$abbreviation[tissue_code$tissue_name_release == tissue_i]]
      #get rat ensembl genes
      DE_genes_in_Nodes_ensembl <- node_metadata_list[[timepoint_to_use]]$ensembl_gene[
        node_metadata_list[[timepoint_to_use]]$tissue == tissue_code$abbreviation[tissue_code$tissue_name_release == tissue_i]]
      if(length(DE_genes_in_Nodes) == 0){next()}
      
      #snag signs
      DE_genes_in_Nodes_sign <- cbind(as.data.frame(do.call(rbind, strsplit(node_metadata_list[[timepoint_to_use]]$node[
        node_metadata_list[[timepoint_to_use]]$tissue == tissue_code$abbreviation[tissue_code$tissue_name_release == tissue_i]], split = "_"))),
        gene = DE_genes_in_Nodes)
      colnames(DE_genes_in_Nodes_sign) <- c("time", "female", "male", "gene")
      DE_genes_in_Nodes_sign$female[grep(pattern = "F1", DE_genes_in_Nodes_sign$female)] <- 1
      DE_genes_in_Nodes_sign$female[grep(pattern = "F-1", DE_genes_in_Nodes_sign$female)] <- -1
      DE_genes_in_Nodes_sign$male[grep(pattern = "M1", DE_genes_in_Nodes_sign$male)] <- 1
      DE_genes_in_Nodes_sign$male[grep(pattern = "M-1", DE_genes_in_Nodes_sign$male)] <- -1
      DE_genes_in_Nodes_sign$ensembl_gene <- DE_genes_in_Nodes_ensembl
      
      #remove genes with no known human gene symbol
      DE_genes_in_Nodes_sign <- DE_genes_in_Nodes_sign[!is.na(DE_genes_in_Nodes_sign$gene),]
      
      #filter down to sex
      DE_genes_signs <- as.integer(DE_genes_in_Nodes_sign[,sex])
      names(DE_genes_signs) <- DE_genes_in_Nodes_sign$gene
      
      #get relative expression scores
      sex_i <- sex
      relative_effect_sizes <- deg_eqtl_list[[tissue_i]][deg_eqtl_list[[tissue_i]]$sex == sex_i & 
                                                           deg_eqtl_list[[tissue_i]]$comparison_group == timepoint,]
      
      
      if(trace_8w_backwards){
        new_nodes <- paths_tiss[match(DE_genes_in_Nodes_sign$ensembl_gene, paths_tiss$feature_ID), paste0("state_", timepoint)]
        
        new_DE_genes_in_Nodes_sign <- cbind(as.data.frame(do.call(rbind, strsplit(new_nodes, split = "_"))), gene = DE_genes_in_Nodes_sign$gene)
        colnames(new_DE_genes_in_Nodes_sign) <- c("female", "male", "gene")
        new_DE_genes_in_Nodes_sign$female[grep(pattern = "F1", new_DE_genes_in_Nodes_sign$female)] <- 1
        new_DE_genes_in_Nodes_sign$female[grep(pattern = "F-1", new_DE_genes_in_Nodes_sign$female)] <- -1
        new_DE_genes_in_Nodes_sign$male[grep(pattern = "M1", new_DE_genes_in_Nodes_sign$male)] <- 1
        new_DE_genes_in_Nodes_sign$male[grep(pattern = "M-1", new_DE_genes_in_Nodes_sign$male)] <- -1
        new_DE_genes_in_Nodes_sign$female[grep(pattern = "F0", new_DE_genes_in_Nodes_sign$female)] <- 0
        new_DE_genes_in_Nodes_sign$male[grep(pattern = "M0", new_DE_genes_in_Nodes_sign$male)] <- 0
        
        #NA nodes get a 0 too
        new_DE_genes_in_Nodes_sign[is.na(new_DE_genes_in_Nodes_sign)] <- 0
        
        DE_genes_signs <- setNames(as.integer(new_DE_genes_in_Nodes_sign[,sex]), DE_genes_in_Nodes_sign$gene)
        
      }
      
      
      results <- t(sapply(setNames(twas_with_hits, twas_with_hits), function(trait_i) {
        twas_signs <- sig_twas_by_trait_signs[[trait_i]]
        shared_genes <- intersect(names(DE_genes_signs), names(twas_signs))
        sign_of_effect <- twas_signs[shared_genes] * DE_genes_signs[shared_genes]
        
        #get proportion of positive hits
        npos <- sum(sign_of_effect > 0.5)
        nneg <- sum(sign_of_effect < -0.5)
        nzero <- length(sign_of_effect) - npos - nneg
        prop_positive <- (npos + nzero * 0.5) / length(sign_of_effect)
        #prop_positive <- mean(sign_of_effect > 0)
        
        sign_of_effect_symbol <- rep("+", length(shared_genes))
        sign_of_effect_symbol[sign_of_effect == -1] <- "-"
        eff_size <- round(abs(relative_effect_sizes$phenotypic_expression_Z[match(shared_genes, relative_effect_sizes$gene_name)]), 2)
        shared_gene_names <- paste0(shared_genes, " (",  sign_of_effect_symbol, ", ", eff_size, ")")
        
        return(list(prop_positive = prop_positive, n = length(shared_genes), genes = paste0(shared_gene_names, collapse = " ~ ")))
      }))
      
      # DELT <- as.data.frame(deg_eqtl_list_TWAS_cluster_subset[[tissue_i]])
      # 
      # #subset to time, sex, and unique gene entries
      # DELT <- DELT[DELT$comparison_group == time,]
      # DELT <- DELT[DELT$sex == sex,]
      # # DELT <- DELT[-which(is.na(DELT$gene_name)),] #get rid of na genes
      # if(nrow(DELT) == 0){next()}
      # DELT <- DELT[match(unique(DELT$gene_name), DELT$gene_name),] #pull out only first gene entry
      # 
      # 
      # DE_inds <- which(DELT$adj_p_value <= 1.1)
      # # TWAS_inds <- apply(DELT[,grep(colnames(DELT), pattern = "BH_PValue")] < 0.05, 2, which)
      # TWAS_inds <- apply(log(DELT[,grep(colnames(DELT), pattern = ".pvalue")]) <= (log(0.05) - log(n_twas_comparisons_clusters)), 2, which)
      # if(length(TWAS_inds) == 0){next()}
      # 
      # intersect_inds <- lapply(TWAS_inds, function(twi) intersect(DE_inds, twi))
      # # intersect_inds <- lapply(intersect_inds, function(ii) ii[DELT$gene_name[ii] %in% cluster_genes]) #subset to just monotonic sex-homogenous clusters
      # intersect_genes <- lapply(intersect_inds, function(ii) unique(DELT$gene_name[ii]))
      # 
      # 
      # intersect_sign_logFC <- lapply(intersect_inds, function(ii) sign(DELT$logFC[ii]))
      # names(intersect_sign_logFC) <- gsub(names(intersect_sign_logFC), pattern = ".pvalue", replacement = "")
      # intersect_sign_TWAS <- lapply(gsub(names(intersect_inds), pattern = ".pvalue", replacement = ""), function(ii) 
      #                        sign(DELT[intersect_inds[paste0(ii, ".pvalue")][[1]], paste0(ii, ".zscore")]))
      # intersect_sign_match <- lapply(1:length(intersect_sign_logFC), function(i) intersect_sign_logFC[[i]]*intersect_sign_TWAS[[i]])
      # 
      # intersect_genes_strings <- lapply(1:length(intersect_sign_logFC), function(i) 
      #   paste0(intersect_genes[[i]], " (", c("-", "", "+")[intersect_sign_match[[i]] + 2], ")", collapse = " ~ "))
      # intersect_genes_strings <- unlist(intersect_genes_strings)
      # intersect_genes_strings[intersect_genes_strings == " ()"] <- ""
      # 
      # intersect_sign_match <- sapply(intersect_sign_match, function(x) mean(x == 1))
      # intersect_sign_match_n <- sapply(intersect_inds, function(x) length(x))
      # deg_sigtwas_proportion[tissue_i,names(n_intersect),time, "p", sex] <- intersect_sign_match
      # deg_sigtwas_proportion[tissue_i,names(n_intersect),time, "n", sex] <- intersect_sign_match_n
      # deg_sigtwas_proportion[tissue_i,names(n_intersect),time, "genes", sex] <- intersect_genes_strings
      
      deg_sigtwas_proportion[tissue_i,rownames(results),timepoint, "p", sex] <- unlist(results[,"prop_positive"])
      deg_sigtwas_proportion[tissue_i,rownames(results),timepoint, "n", sex] <- unlist(results[,"n"])
      deg_sigtwas_proportion[tissue_i,rownames(results),timepoint, "genes", sex] <- unlist(results[,"genes"])
      
    }
    
  }
  
}
save(deg_sigtwas_proportion, file = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/deg_sigtwas_proportion.txt")


#### figure 5 preprocessing ####

#load relevant data
base = "deviation_from_half"
load(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/mcmc_output/", paste0(base, ifelse(use_all_cats, "_allCats", "")),".cmdStanR.fit"))
samps <- data.frame(as_draws_df(out$draws()))

# samps[,grep("logit_prop_21\\.", colnames(samps))]
# sum(apply(samps[,grep("logit_prop_21\\.", colnames(samps))], 2, prop_greater_than_0) > 0.95)
# sum(apply(samps[,grep("logit_prop_21\\.", colnames(samps))], 2, prop_greater_than_0) < 0.05)

alpha = 0.05
#process this into a format for plotting
posterior_mean_probs <- data.frame(trait = data$trait, 
                                   tissue = data$tissue, 
                                   prob = invlogit(apply(subset_samps("logodds", c("raw", "sd", "compl"), samps = samps), 2, mean)))
deg_sigtwas_proportion_posterior_mean <- reshape(posterior_mean_probs, idvar = "tissue", timevar = "trait", direction = "wide")
rownames(deg_sigtwas_proportion_posterior_mean) <- deg_sigtwas_proportion_posterior_mean$tissue
deg_sigtwas_proportion_posterior_mean <- deg_sigtwas_proportion_posterior_mean[,-1]
colnames(deg_sigtwas_proportion_posterior_mean) <- gsub("prob\\.", "", colnames(deg_sigtwas_proportion_posterior_mean))

posterior_mean_prob_bias <- data.frame(trait = data$trait, 
                                       tissue = data$tissue, 
                                       prob = apply(subset_samps("cell_total_prob_bias", c("raw", "sd", "compl"), samps = samps), 2, prop_greater_than_0))
deg_sigtwas_proportion_posterior_prob_bias <- reshape(posterior_mean_prob_bias, idvar = "tissue", timevar = "trait", direction = "wide")
rownames(deg_sigtwas_proportion_posterior_prob_bias) <- deg_sigtwas_proportion_posterior_prob_bias$tissue
deg_sigtwas_proportion_posterior_prob_bias <- deg_sigtwas_proportion_posterior_prob_bias[,-1]
colnames(deg_sigtwas_proportion_posterior_prob_bias) <- gsub("prob\\.", "", colnames(deg_sigtwas_proportion_posterior_prob_bias))

#check to make sure this is consistent with model
# raw_logodds * col_sd_sd .* exp(raw_col_logsd[col_index] * col_logsd_sd) + 0 + col_bias[col_index]
# apply(subset_samps("col_bias", c("raw", "sd", "compl"), samps = samps), 2, prop_greater_than_0)
prop_greater_than_0(subset_samps("raw_logodds\\.1\\.", c("sd", "compl"), samps = samps) * subset_samps("col_sd_sd", c("compl"), samps = samps) * 
                      exp(subset_samps("raw_col_logsd\\.1\\.", c("compl"), samps = samps) * subset_samps("col_logsd_sd", c("compl"), samps = samps)) + 
                      subset_samps("col_bias\\.1\\.", c("raw", "sd", "compl"), samps = samps))
posterior_mean_prob_bias[1,]
#YEP!

# trait_means <- invlogit(apply(subset_samps("col_bias", c("raw", "sd", "compl"), samps = samps) + 
#                                 subset_samps("col_means_compl", c("raw", "sd"), samps = samps), 
#                               2, mean))
trait_means <- invlogit(apply(subset_samps("col_bias", c("raw", "sd", "compl"), samps = samps), 
                              2, mean))
names(trait_means) <- traits
trait_bias_probs <- apply(subset_samps("col_bias", c("raw", "sd", "compl"), samps = samps), 2, prop_greater_than_0)
names(trait_bias_probs) <- traits

#set a few graphical parameters
if(length(salient.categories) < 9){
  category_colors <- RColorBrewer::brewer.pal(length(salient.categories), "Dark2")
  names(category_colors) <- sort(salient.categories)  
} else {
  category_colors <- RColorBrewer::brewer.pal(5, "Dark2")
  names(category_colors) <- sort(c("Cardiometabolic", "Aging", "Anthropometric", 
                                   "Immune", "Psychiatric-neurologic"))
  cats_leftover <- setdiff(salient.categories, names(category_colors))
  cols_leftover <- setdiff(RColorBrewer::brewer.pal(8, "Dark2"), category_colors)
  n_cols_needed <- length(cats_leftover) - length(cols_leftover)
  if(n_cols_needed > 0){
    extra_cols <- c(cols_leftover, RColorBrewer::brewer.pal(n_cols_needed, "Set2"))
  } else {
    extra_cols <- cols_leftover
  }
  names(extra_cols) <- cats_leftover
  category_colors <- c(category_colors, extra_cols)
}

#process significance files
# deg_sigtwas_proportion_posterior_mean
alpha_post = 0.1
cell_sigs <- deg_sigtwas_proportion_posterior_prob_bias < alpha_post | deg_sigtwas_proportion_posterior_prob_bias > (1-alpha_post)
# trait_means 
trait_sigs <- trait_bias_probs < alpha_post | trait_bias_probs > (1-alpha_post)

#### figure 5 plotting ####
load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/deg_sigtwas_proportion.txt")
subset_to_traits <- T
cairo_pdf(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/figures/prop_pos_hits_bayesian_scatterplot", 
                 ifelse(subset_to_traits, "_sub.pdf", "")), 
          width = 2000 / 72 * length(twas_with_hits) / 2 / 73, height = 700 / 72 * 2, family="Arial Unicode MS")
par(mfrow = c(2,1), xpd = NA, 
    mar = c(19,
            9,
            5,
            7 + ifelse(subset_to_traits, 2, 0)))
lwd <- 3
tissue_name_cex <- 1.2


for(subset_i in 1:2){
  if(subset_to_traits){
    trait_subset <- colnames(n_deg_sigtwas_intersect)[
      traitwise_partitions$Category[match(colnames(n_deg_sigtwas_intersect), traitwise_partitions$Tag)] %in% 
        ifelse2(subset_to_traits, c("Cardiometabolic", "Psychiatric-neurologic"), salient.categories)]
    
    #alternatively
    if(subset_i == 1){
      trait_subset <- names(trait_means)[trait_means > 0.5]  
    } else {
      trait_subset <- names(trait_means)[trait_means < 0.5]
    }
    
    traits_to_plot <- trait_subset
    
  } else {
    traits_to_plot <- twas_with_hits
  }
  categories_represented <- unique(traitwise_partitions$Category[match(traits_to_plot, traitwise_partitions$Tag)])
  categories_represented <- gsub("\\ .*", "", gsub("\\-.*", "", categories_represented))
  category_colors <- c(category_colors, setNames(category_colors[sapply(categories_represented, function(cati) grep(cati, names(category_colors))[1])], 
                                                 categories_represented)[!(categories_represented %in% names(category_colors))])
  
  cols = list(Tissue=MotrpacRatTraining6moData::TISSUE_COLORS[names(motrpac_gtex_map)], 
              Time=MotrpacRatTraining6moData::GROUP_COLORS[paste0(c(1,2,4,8), "w")],
              Sex=MotrpacRatTraining6moData::SEX_COLORS[c('male','female')])
  cols$Tissue<- cols$Tissue[order(match(MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[names(cols$Tissue)], MotrpacRatTraining6moData::TISSUE_ORDER))]
  
  tissue_names <- sapply(strsplit(names(cols$Tissue), "-"), function(x) paste0(x[ifelse(length(x) == 2, c(2), list(2:3))[[1]]], collapse = " "))
  
  #other plotting params
  trait_category_legend_below <- T
  
  #do the plotting
  plot(100,100, xlim = c(3,length(traits_to_plot)), ylim = range(deg_sigtwas_proportion_posterior_mean[,traits_to_plot]), xpd=NA, 
       ylab = "", xlab = "", xaxt = "n", 
       yaxt = "n", bty="n", cex.lab = 2, cex.axis = 1.25)
  ry <- diff(par("usr")[3:4])
  
  text(label = "Posterior Mean Difference (logit-scale)\nin Positive Effects on GWAS Trait", x = par("usr")[1] - diff(par("usr")[1:2])/15, 
       y = mean(par("usr")[3:4]), srt = 90, pos = 3,
       cex = 1.75)
  
  for(sex in c("male", "female")[1]){
    
    # mean_freq <- apply(sapply(traits_to_plot, function(trait) as.numeric(deg_sigtwas_proportion[, trait,"8w","p",sex])), 2, weighted.mean, na.rm = T)
    
    mean_freq <- trait_means[traits_to_plot]
    
    order_traits_to_plot <- order(mean_freq, decreasing = T)
    
    #plot faded positive and negative regions
    yseq <- round(seq(min(deg_sigtwas_proportion_posterior_mean[,traits_to_plot]) - 0.0025, max(deg_sigtwas_proportion_posterior_mean[,traits_to_plot]) + 0.0025, length.out = 9), 2)
    rect(xl = 1, xr = length(traits_to_plot) + 2, yb = 0.5, ytop = max(yseq),
         col = grDevices::adjustcolor("red", 0.1), border = NA)
    rect(xl = 1, xr = length(traits_to_plot) + 2, yb = min(yseq), ytop = 0.5,
         col = grDevices::adjustcolor("blue", 0.1), border = NA)
    
    #vertical axis
    segments(x0 = 1 - length(traits_to_plot) * 0.005, x1 = 1 - 3 * 0.04, y0 = yseq, y1 = yseq, lwd = 2)
    segments(x0 = 1 - 3 * 0.02, x1 = 1 - 3 * 0.02, y0 = min(yseq), y1 = max(yseq), lwd = 2)
    text(x = 1 - 3 * 0.05, y = yseq, labels = round(logit(yseq), 2), pos = 2, cex = 1.5)
    
    #horizontal axis
    text(1:length(traits_to_plot) + 1.5, min(yseq) - ry / 10, srt = 45, pos = 2,
         labels = gsub("raneal", "ranial", paste0(ifelse(trait_sigs[traits_to_plot][order_traits_to_plot], "♦ ", ""), 
                         trait_categories$new_Phenotype[match(traits_to_plot[order_traits_to_plot], trait_categories$Tag)])), 
         cex = 1.125, col = 1)
    segments(1:length(traits_to_plot) + 1, min(yseq) - ry / 15, 1:length(traits_to_plot) + 1, max(yseq), col = adjustcolor(1, 0.2), lty = 2)
    
    #plot points
    points(1:length(traits_to_plot)+1, y = sort(mean_freq[traits_to_plot], decreasing = T), pch = "*", cex = 3)
    
    for(trait_i in (1:length(traits_to_plot))){
      
      #plot category blocks
      rect(xleft = which(order_traits_to_plot == trait_i) + 1/2, 
           xright = which(order_traits_to_plot == trait_i) + 3/2,
           ybottom = min(yseq)- ry / 30, 
           ytop = min(yseq)-ry / 14,
           col = category_colors[trait_categories$Category[match(traits_to_plot[trait_i], trait_categories$Tag)]])
      text(labels = substr(trait_categories$Category[match(traits_to_plot[trait_i], trait_categories$Tag)], 1, 2), 
           x = which(order_traits_to_plot == trait_i) + 1,
           y = min(yseq)-ry / 19.5, 
           col = "white", cex = 1)
      
      
      d <- unlist(deg_sigtwas_proportion_posterior_mean[,traits_to_plot[trait_i]])
      d_sig <- unlist(cell_sigs[,traits_to_plot[trait_i]])
      tissue_abbr_rev <- MotrpacRatTraining6moData::TISSUE_ABBREV_TO_CODE
      names(d) <- tissue_abbr_rev[rownames(deg_sigtwas_proportion_posterior_mean)]
      
      dn <- as.numeric(deg_sigtwas_proportion[, traits_to_plot[trait_i],"8w","n",sex])
      names(dn) <- rownames(deg_sigtwas_proportion)
      dn <- dn[names(d)]
      
      #plot white points first
      points(x = rep(which(order_traits_to_plot == trait_i) + 1, length(d)), y = d, pch = 19 - d_sig, col = "white", cex = dn^0.25)
      #and then the actual points
      points(x = rep(which(order_traits_to_plot == trait_i) + 1, length(d)), y = d, pch = 19 - d_sig, 
             col = adjustcolor(cols$Tissue[names(d)], 0.5), cex = dn^0.25)
      #and then the significant points
      points(x = rep(which(order_traits_to_plot == trait_i) + 1, length(d)), y = d, pch = 5, 
             col = 1, cex = dn^0.25*d_sig*0.8)
      
    }
    
    
    #legend for tissues
    tisscols <- cols$Tissue[!(names(cols$Tissue) %in% c("t64-ovaries", "t63-testes"))]
    points(x = rep(length(traits_to_plot) + 3, length(tisscols)), 
           y = seq(mean(yseq)-ry / 15, max(yseq), length.out = length(tisscols)), 
           col = adjustcolor(tisscols, 0.5),
           pch = 19, cex = 2)
    text(x = rep(length(traits_to_plot) + 3.25, length(tisscols)), 
         y = seq(mean(yseq)-ry / 15, max(yseq), length.out = length(tisscols)), 
         pos = 4, pch = 19, cex = 1,
         labels = MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[names(tisscols)])

    #legend for tissue size
    n_pts <- 5
    pt_size_logs <- seq(1, log(max(as.numeric(deg_sigtwas_proportion[,,"8w","n",sex]), 
                                   na.rm = T)) / log(2), length.out = n_pts)
    pt_size_legend <- round(2^pt_size_logs)
    pty <-  par("usr")[3] + cumsum(pt_size_legend^0.25/800) * ry * 15 + 
      cumsum(rep(ry / 200, n_pts)) - ry / 6
    text(x = length(traits_to_plot) + 2.25, 
         y = max(pty) + diff(pty)[4], 
         pos = 4, pch = 19, cex = 1.1,
         labels = "Sample Size")
    points(x = rep(length(traits_to_plot) + 3.25, n_pts), 
           y = pty, 
           col = adjustcolor(1, 0.5),
           pch = 19, cex = pt_size_legend^0.25)
    text(x = rep(length(traits_to_plot) + 3.25, n_pts) + pt_size_legend^0.25/10, 
         y = pty, 
         pos = 4, pch = 19, cex = 1,
         labels = pt_size_legend)
    
    #legend for stars
    points(x = length(traits_to_plot) + 3, 
           y = par("usr")[4]-ry / 1.525, 
           col = "black", pch = "*", cex = 3)
    text(x = length(traits_to_plot) + 3.25, cex = 1.1,
         y = par("usr")[4]-ry / 1.5, 
         pos = 4, labels = "Posterior\nMean")
    
    
    #legend for diamonds
    text(x = length(traits_to_plot) * 1.055, cex = 1,
         y = par("usr")[3] + ry / 4.9, 
         pos = 4, labels = paste0("♦ indicates\nPr(enr. or dep.)\n> ", 1-alpha_post))
    
    
    #legend for categories
    if(trait_category_legend_below){
      yadj <- 0
      yloc_factor <- ifelse(subset_i == 1, 1.8, 2)
      for(i in 1:length(categories_represented)){
        rect(xleft = 0 + i + cumsum(c(0,strwidth(categories_represented, cex = 1.25, units = "user") + 1))[i],
             xright = 1 + i + cumsum(c(0,strwidth(categories_represented, cex = 1.25, units = "user") + 1))[i],
             ybottom = min(yseq) - ry / yloc_factor - ry / 14,
             ytop = min(yseq) - ry / yloc_factor - ry / 30,
             col = category_colors[categories_represented[i]], border = 1)
        text(labels = categories_represented[i], pos = 4, 
             y = min(yseq) - ry / yloc_factor - ry / 18, 
             cex = 1.25, 
             x = 0.85 + i + cumsum(c(0,strwidth(categories_represented, cex = 1.25, units = "user") + 1))[i])
        text(labels = substr(categories_represented[i], 1, 2), 
             x = 0.5 + i + cumsum(c(0,strwidth(categories_represented, cex = 1.25, units = "user") + 1))[i],
             y = min(yseq) - ry / yloc_factor - ry / 19.5, 
             col = "white", cex = 1)
        #draw circle?
        ends_y <- c(min(yseq) - ry / yloc_factor - ry / 18, min(yseq) - ry / 18)
        arc(t1 = 3*pi/2, t2 = pi/2+1E-6, r1 = diff(range(ends_y))/2, r2 = diff(range(ends_y))/2, 
            center = c(0.5, mean(ends_y)), lwd = 3, res = 50, adjx = 210 + ifelse(subset_i == 1, 50, 70))
        points(0.5, ends_y[2], pch = -9658, cex = 3)
      }
    } else {
      x_adj2 <- x_adj - 2.5
      y_adj2 <- y_adj - 2.5
      for(i in 1:length(categories_represented)){
        rect(xleft = 0 + i,
             xright = 1 + i,
             ybottom = -10,
             ytop = -9,
             col = category_colors[categories_represented[i]], border = 1)
        text(labels = categories_represented[i], pos = 4, y = i + y_adj2 - 7, x = ncol(table_to_use) - 1/2 + x_adj2 + 1/3)
      }
    }
    
  }
  
  title_text <- paste0("Traits \\textbf{", ifelse(subset_i == 1, "More", "Less"), 
                       "} Similar to Transcriptomic Adaptation to Exercise")
  text(latex2exp::TeX(title_text), 
       x = (par("usr")[1]), y = par("usr")[4] + ry / 15, 
       cex = 2.25, pos = 4)
  rect(xleft = mean(par("usr")[1]) + strwidth(latex2exp::TeX("Traits  "), cex = 2.25),
       xright = mean(par("usr")[1]) + strwidth(latex2exp::TeX(paste0("Traits  \\textbf{", 
                                                                     ifelse(subset_i == 1, "More", "Less"), "}")), 
                                               cex = 2.25),
       ybottom = par("usr")[4] + ry / 20,
       ytop = par("usr")[4] + ry / 20 + 
         strheight(latex2exp::TeX(paste0("\\textbf{", ifelse(subset_i == 1, "More", "Less"), "}")), cex = 2.8),
       col = "white", border = "white")
  text(latex2exp::TeX(paste0("\\textbf{", ifelse(subset_i == 1, "More", "Less"), "}")), 
       x = (par("usr")[1]) + strwidth("Traits ", cex = 2.25), 
       y = par("usr")[4] + ry / 15 + (
         strheight(latex2exp::TeX(title_text), cex = 2.25) - 
           strheight(latex2exp::TeX(paste0("\\textbf{", ifelse(subset_i == 1, "More", "Less"), "}")), cex = 2.25)
       ) / 2, 
       cex = 2.25, pos = 4, col = adjustcolor(c("red", "blue")[subset_i], "0.75"))
  
  
}
dev.off()


#### END Script ####

#### try doing just the frequentist test ####
FET_pvals <- sapply(1:nrow(data1), function(i){
  two_by_two <- matrix(0,2,2, dimnames = list(c("DEG", "nDEG"), c("TWAS", "nTWAS")))
  two_by_two["DEG", "TWAS"] <-  data1$count[i]
  two_by_two["nDEG", "TWAS"] <-  sig_twas_by_trait_genes_matrix[data1$tissue[i], data1$trait[i]] - data1$count[i]
  two_by_two["DEG", "nTWAS"] <-  n_genes_in_nodes_matrix[data1$tissue[i], data1$trait[i]] - data1$count[i]
  two_by_two["nDEG", "nTWAS"] <-  total_number_of_possible_hits_matrix[data1$tissue[i], data1$trait[i]] - sum(two_by_two)
  fisher.test(two_by_two)$p.value
})
hist(FET_pvals, breaks = 20, main = "nominal pvals")

data1$FET_pval <- FET_pvals
data_ihw <- data1
# data_ihw <- data_ihw[data1$count != 0,]
# data_ihw <- data_ihw[!(data_ihw$trait %in% names(which(table(data_ihw$trait) < 3))),]
data_ihw$trait <- as.factor(as.character(data_ihw$trait))
ihw_results_FET <- IHW::ihw(FET_pval ~ trait, data = data_ihw, alpha = 0.05)

par(mfrow = c(2,1))
hist(data_ihw$FET_pval, breaks = 20, main = "nominal pvals")
hist(ihw_results_FET@df$adj_pvalue)
sum(ihw_results_FET@df$adj_pvalue < 0.05)
data_ihw$TWAS_total <- as.numeric(sapply(1:nrow(data_ihw),function(i){sig_twas_by_trait_genes_matrix[data_ihw$tissue[i], data_ihw$trait[i]]}))
ihw_results_FET <- IHW::ihw(FET_pval ~ TWAS_total, data = data_ihw, alpha = 0.05, nbins = 2)  
sum(ihw_results_FET@df$adj_pvalue < 0.05)

sapply(1:20, function(n_bins) sum(IHW::ihw(FET_pval ~ TWAS_total, data = data_ihw, alpha = 0.05, nbins = n_bins)@df$adj_pvalue < 0.05))

#### quick plot of exp vs obs freqs ####
use_focal_vs_compl <- F

traits <- unique(data1$trait)
tissues <- unique(data1$tissue)
d <- list(cell_count = data1$count,
          total = sapply(1:nrow(data1), function(i) total_number_of_possible_hits_matrix[data1$tissue[i], data1$trait[i]]),
          row_count = sapply(1:nrow(data1), function(i) n_genes_in_nodes_matrix[data1$tissue[i], data1$trait[i]]),
          col_count = sapply(1:nrow(data1), function(i) sig_twas_by_trait_genes_matrix[data1$tissue[i], data1$trait[i]]),
          row_index = match(data1$tissue, tissues),
          col_index = match(data1$trait, traits),
          row_n = length(unique(data1$tissue)),
          col_n = length(unique(data1$trait)))


cairo_pdf(paste0("~/Documents/Documents - nikolai/motrpac_companion/figures/pass1b_fig8_DEG-TWAS_logodds.pdf"), 
          width = 1000 / 72, height = 950 / 72, family="Arial Unicode MS", pointsize = 25)

category_shapes <- setNames(15:19, salient.categories)
category_shapes <- setNames(15:19, c("Cardiometabolic", "Aging", "Anthropometric", 
                                     "Immune", "Psychiatric-neurologic")
)

par(mar = c(5,4,3,4))
#first plot


sapply(1:nrow(signif_df), function(i) adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index[i]]], as.numeric(signif_df$signif[i] != 0) * 0.625 + 0.125))
if(use_focal_vs_compl){
  xvals <- (d$col_count - d$cell_count) / (d$total - d$row_count)
  yvals <- d$cell_count / d$row_count
} else {
  xvals <- d$row_count / d$total * d$col_count / d$total  
  yvals <- d$cell_count / d$total
}




plot(logit(xvals), logit(yvals), xlim = ifelse2(use_focal_vs_compl, c(-14,-1), c(-15,-4.25)), 
     ylim = ifelse2(use_focal_vs_compl, c(-6, -1), c(-9.5,-4.25)),
     pch = category_shapes[traitwise_partitions$Category[match(traits[d$col_index], traitwise_partitions$Tag)]], 
     col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5), 
     cex = 1.5,
     xlab = "", ylab = "", yaxt = 'n', xaxt = "n")
box(which = "plot", lwd = 2)
#axes
text(labels = "observed logit(frequency)", x = par("usr")[2] + diff(par("usr")[1:2])/8, srt = 270, pos = 1, 
     y = mean(par("usr")[3:4]), xpd = NA, cex = 1.25)
text(labels = "expected logit(frequency)", x = mean(par("usr")[1:2]), srt = 0, pos = 1, 
     y = par("usr")[3] - diff(par("usr")[3:4])/10, xpd = NA, cex = 1.25)
yaxlocs <- ifelse2(use_focal_vs_compl, c(-6:-1), c(-5:-9))
segments(x0 = par("usr")[2], x1 = par("usr")[2] + diff(par("usr")[1:2])/100, y0 = yaxlocs, y1 = yaxlocs, lwd = 4, xpd = NA)
text(labels = yaxlocs, x = par("usr")[2] + diff(par("usr")[1:2])/200, srt = 0, pos = 4, 
     y = yaxlocs, xpd = NA, cex = 1.25)
xaxlocs <- ifelse2(use_focal_vs_compl, -7:-1*2, -7:-3*2)
segments(x0 = xaxlocs, x1 = xaxlocs, y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4])/80, lwd = 4, xpd = NA)
text(labels = xaxlocs, x = xaxlocs, srt = 0, pos = 1, 
     y = par("usr")[3] - diff(par("usr")[3:4])/80, xpd = NA, cex = 1.25)
abline(0,1, lty = 2, lwd = 6, col = adjustcolor(1,0.75))

#secondplot
xl <- par("usr")[1]
xr <- par("usr")[1] + diff(par("usr")[1:2]) / 2
yt <- par("usr")[4]
yb <- par("usr")[4] - diff(par("usr")[3:4]) / 2

buffer <- 0.075
bx <- (xr-xl) * buffer
by <- (yt-yb) * buffer
rect(xleft = xl, ybottom = yb, xright = xr, ytop = yt, lwd = 2, xpd = NA)

points(x = ((xvals - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
       y = ((yvals - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, 
       pch = category_shapes[traitwise_partitions$Category[match(traits[d$col_index], traitwise_partitions$Tag)]], 
       col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5), 
       cex = 1)
segments(x0 = 0 * (xr - xl - 2*bx) + xl, 
         y0 = 0 * (yt - yb - 2*by) + yb,
         x1 = ((max(yvals, na.rm = T) - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
         y1 = ((max(yvals, na.rm = T) - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, 
         lty = 2, lwd = 4, col = adjustcolor(1,0.75))

#axes
text(labels = "observed frequency", x = par("usr")[1] - diff(par("usr")[1:2])/ifelse2(use_focal_vs_compl, 8, 7), srt = 90, pos = 1, 
     y = par("usr")[4] - diff(par("usr")[3:4])/5, xpd = NA, cex = 1)
text(labels = "expected frequency", x = par("usr")[1] + diff(par("usr")[1:2])/4, srt = 0, pos = 1, 
     y = par("usr")[4] + diff(par("usr")[3:4])/8, xpd = NA, cex = 1)
xaxlocs <- ifelse2(use_focal_vs_compl, c(0:4/20), c(0:6/500))
yaxlocs <- ifelse2(use_focal_vs_compl, c(0:6/20), c(0:6/500))

segments(x0 = par("usr")[1], x1 = par("usr")[1] - diff(par("usr")[1:2])/100, 
         y0 = ((yaxlocs - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, 
         y1 = ((yaxlocs - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, 
         lwd = 4, xpd = NA)
segments(y0 = par("usr")[4], y1 = par("usr")[4] + diff(par("usr")[3:4])/100, 
         x0 = ((xaxlocs - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
         x1 = ((xaxlocs - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
         lwd = 4, xpd = NA)
text(labels = xaxlocs, x = ((xaxlocs - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
     srt = 0, pos = 3, 
     y = par("usr")[4], xpd = NA, cex = 0.75)
text(labels = yaxlocs, x = par("usr")[1], 
     srt = 0, pos = 2, 
     y = ((yaxlocs - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, xpd = NA, cex = 0.75)
# text(labels = 0:6/500, x = -7:-3*2, srt = 0, pos = 1, 
#      y = par("usr")[3] - diff(par("usr")[3:4])/80, xpd = NA, cex = 1.25)
abline(0,1, lty = 2, lwd = 6, col = adjustcolor(1,0.75))

#legends
legend(x = xl, y = yb, legend = names(category_shapes), pch = category_shapes, col = adjustcolor(1, 0.5), cex = 0.75, pt.cex = 1.4)
legend(x = xr, y = yt, legend = tissues, pch = 19, col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues], 0.5), cex = 0.55, pt.cex = 1.2)
segments(x0 = par("usr")[1] + diff(par("usr")[1:2]) / 1.55, 
         y0 = par("usr")[4] - diff(par("usr")[3:4])/50,
         x1 = par("usr")[1] + diff(par("usr")[1:2]) / 1.45, 
         y1 = par("usr")[4] - diff(par("usr")[3:4])/50, 
         lty = 2, lwd = 4, col = adjustcolor(1,0.75))
text(x = par("usr")[1] + diff(par("usr")[1:2]) / 1.45, 
     y = par("usr")[4] - diff(par("usr")[3:4])/45, 
     pos = 4, cex = 0.75, labels = "1-to-1 line")

# plot(d$row_count / d$total * d$col_count / d$total, d$cell_count / d$total)

dev.off()


#### plot rowbias as a boxplot ####
if(length(salient.categories) < 9){
  category_colors <- RColorBrewer::brewer.pal(length(salient.categories), "Dark2")
  names(category_colors) <- sort(salient.categories)  
} else {
  category_colors <- RColorBrewer::brewer.pal(5, "Dark2")
  names(category_colors) <- sort(c("Cardiometabolic", "Aging", "Anthropometric", 
                                   "Immune", "Psychiatric-neurologic"))
  cats_leftover <- setdiff(salient.categories, names(category_colors))
  cols_leftover <- setdiff(RColorBrewer::brewer.pal(8, "Dark2"), category_colors)
  n_cols_needed <- length(cats_leftover) - length(cols_leftover)
  if(n_cols_needed > 0){
    extra_cols <- c(cols_leftover, RColorBrewer::brewer.pal(n_cols_needed, "Set2"))
  } else {
    extra_cols <- cols_leftover
  }
  names(extra_cols) <- cats_leftover
  category_colors <- c(category_colors, extra_cols)
}
cols$category <- category_colors
CI_width <- 0.8
CI_range <- c((1-CI_width)/2, CI_width+(1-CI_width)/2)


grDevices::cairo_pdf(filename = paste0("~/Documents/Documents - nikolai/motrpac_companion/figures/figure_intersect_enrichment_posterior_summary.pdf"), 
                     width = 1400 / 72, height = 800 / 72, family="Arial Unicode MS", pointsize = 20)
layout(rbind(c(1,1,1,1,2,2,2), c(3,3,3,3,3,3,3)))

#plot rowbias as a violin plot

par(mar = c(5,4.5,2,2))
tord <- order(apply(subset_samps("row_bias", c("raw", "sd"), samps = samps), 2, mean))
CI_locs <- apply(subset_samps("row_bias", c("raw", "sd"), samps = samps) + samps$overall_bias, 2, quantile, probs = CI_range)[,tord]
qi_100 <- apply(subset_samps("row_bias", c("raw", "sd"), samps = samps) + samps$overall_bias, 2, quantile, probs = c(0.0, 1.0))[,tord]
posterior_means <- apply(subset_samps("row_bias", c("raw", "sd"), samps = samps) + samps$overall_bias, 2, mean)[tord]
inside_cols <- rep("black", length(tord))
inside_cols[overlaps_with_zero(CI_locs)] <- "white"
tmp <- vioplot::vioplot(x = as.matrix(subset_samps("row_bias", c("raw", "sd"), samps = samps) + samps$overall_bias)[,tord], T, 
                        col = MotrpacRatTraining6moData::TISSUE_COLORS[tissues][tord], outline=FALSE, xaxt = "n",
                        names = tissues[tord], srt = 90, range = 0, xlab = "", lineCol = MotrpacRatTraining6moData::TISSUE_COLORS[tissues][tord],
                        ylab = "multilevel tissue enrichment (logodds-scale)", cex.lab = 1, plotCentre = "point", 
                        colMed = MotrpacRatTraining6moData::TISSUE_COLORS[tissues][tord])
tick <- seq_along(tissues)
axis(1, at = tick, labels = F)
for(i in 1:length(tissues)){
  segments(x0 = i, x1 = i, y0 = CI_locs[1,i], y1 = CI_locs[2,i], lwd = 4, col = inside_cols[i])
  points(x = i, y = posterior_means[i], col = inside_cols[i], pch = 19, cex = 1.5)
  #guiding lines
  segments(x0 = i, x1 = i, y0 = min(qi_100) - diff(range(qi_100)) / 30, y1 = qi_100[1,i], lwd = 1, col = "black", xpd = T, lty = 3)
  segments(x0 = i, x1 = i, y0 = qi_100[2,i], y1 = max(qi_100) * 0.93 + 1 * 0.07, lwd = 1, col = "black", xpd = T, lty = 3)
}
text(tick + 0.2, par("usr")[3] - 0.05, tissues[tord], srt = 45, xpd = T, pos = 2)
fig_label("a)", cex = 2, shrinkX = 0.4)
abline(h=0,lwd=3,lty=2, col = adjustcolor(1,0.5))


### plot colcatbias as a boxplot ###

par(mar = c(7,5.5,2,2))
tord <- order(apply(subset_samps("colcat_bias", c("raw", "sd"), samps = samps), 2, mean))
CI_locs <- apply(subset_samps("colcat_bias", c("raw", "sd"), samps = samps) + samps$overall_bias, 2, quantile, probs = CI_range)[,tord]
qi_100 <- apply(subset_samps("colcat_bias", c("raw", "sd"), samps = samps) + samps$overall_bias, 2, quantile, probs = c(0.0, 1.0))[,tord]
posterior_means <- apply(subset_samps("colcat_bias", c("raw", "sd"), samps = samps) + samps$overall_bias, 2, mean)[tord]
inside_cols <- rep("black", length(tord))
inside_cols[overlaps_with_zero(CI_locs)] <- "white"
vioplot::vioplot(x = as.matrix(subset_samps("colcat_bias", c("raw", "sd"), samps = samps) + samps$overall_bias)[,tord], T, 
                 col = category_colors[trait_cats][tord], outline=FALSE, xaxt = "n",
                 names = trait_cats[tord], srt = 90, range = 0, xlab = "", lineCol = category_colors[trait_cats][tord],
                 ylab = "multilevel category enrichment\n(logodds-scale)", cex.lab = 1, plotCentre = "point", 
                 colMed = category_colors[trait_cats][tord])
tick <- seq_along(trait_cats)
axis(1, at = tick, labels = F)
for(i in 1:length(trait_cats)){
  segments(x0 = i, x1 = i, y0 = CI_locs[1,i], y1 = CI_locs[2,i], lwd = 4, col = inside_cols[i])
  points(x = i, y = posterior_means[i], col = inside_cols[i], pch = 19, cex = 1.5)
  #guiding lines
  segments(x0 = i, x1 = i, y0 = min(qi_100) - diff(range(qi_100)) / 30, y1 = qi_100[1,i], lwd = 1, col = "black", xpd = T, lty = 3)
  segments(x0 = i, x1 = i, y0 = qi_100[2,i], y1 = max(qi_100) * 0.8 + 1 * 0.2, lwd = 1, col = "black", xpd = T, lty = 3)
}
text(tick + 0.35, par("usr")[3] - 0.1, trait_cats[tord], srt = 45, xpd = T, pos = 2, xpd = NA)
abline(h=0,lwd=3,lty=2, col = adjustcolor(1,0.5))


#### plot colbias as a boxplot ###
par(mar = c(12,5.5,2,2))
trait_subset_bool <- traitwise_partitions$Category[match(traits, traitwise_partitions$Tag)] %in% c("Cardiometabolic")
trait_subset_bool <- traitwise_partitions$Category[match(traits, traitwise_partitions$Tag)] %in% salient.categories
tord <- order(apply(subset_samps("col_bias", c("raw", "sd"), samps = samps), 2, mean)[trait_subset_bool])
CI_locs <- apply((subset_samps("col_bias", c("raw", "sd"), samps = samps) + samps$overall_bias)[,trait_subset_bool], 2, quantile, probs = CI_range)[,tord]
qi_100 <- apply((subset_samps("col_bias", c("raw", "sd"), samps = samps) + samps$overall_bias)[,trait_subset_bool], 2, quantile, probs = c(0.00, 1.00))[,tord]
posterior_means <- apply(subset_samps("col_bias", c("raw", "sd"), samps = samps)[,trait_subset_bool] + samps$overall_bias, 2, mean)[tord]
inside_cols <- rep("black", length(tord))
inside_cols[overlaps_with_zero(CI_locs)] <- "white"
tmp <- vioplot::vioplot(x = as.matrix((subset_samps("col_bias", c("raw", "sd"), samps = samps) + samps$overall_bias)[,trait_subset_bool])[,tord], T, 
                        col = category_colors[traitwise_partitions$Category[match(traits[trait_subset_bool], traitwise_partitions$Tag)]][tord], 
                        outline=FALSE, xaxt = "n", xlim = c(4,sum(trait_subset_bool)-3),
                        names = tissues[tord], srt = 90, range = 0, xlab = "", 
                        lineCol = category_colors[traitwise_partitions$Category[match(traits[trait_subset_bool], traitwise_partitions$Tag)]][tord],
                        ylab = "multilevel trait enrichment\n(logodds-scale)", cex.lab = 1, plotCentre = "point", 
                        colMed = category_colors[traitwise_partitions$Category[match(traits[trait_subset_bool], traitwise_partitions$Tag)]][tord])
tick <- seq_along(traits[trait_subset_bool])
axis(1, at = tick, labels = F)
for(i in 1:length(traits[trait_subset_bool])){
  segments(x0 = i, x1 = i, y0 = CI_locs[1,i], y1 = CI_locs[2,i], lwd = 2, col = inside_cols[i])
  points(x = i, y = posterior_means[i], col = inside_cols[i], pch = 19)
  #guiding lines
  segments(x0 = i, x1 = i, y0 = min(qi_100) - diff(range(qi_100)) / 30, y1 = qi_100[1,i], lwd = 1, col = "black", xpd = T, lty = 3)
  segments(x0 = i, x1 = i, y0 = qi_100[2,i], y1 = max(qi_100) * 0.75 + 1 * 0.25, lwd = 1, col = "black", xpd = T, lty = 3)
}
text(tick + 0.55, par("usr")[3] - 0.125, cex = 0.7,
     trait_categories$new_Phenotype[match(traits[trait_subset_bool][tord], trait_categories$Tag)], srt = 45, xpd = T, pos = 2)
fig_label("c)", cex = 2, shrinkX = 0.3)
abline(h=0,lwd=3,lty=2, col = adjustcolor(1,0.5))

dev.off()

#### plot rowcatbias as a boxplot ####

par(mar = c(8,5.5,2,2))
tissue_category_colors <- setNames(RColorBrewer::brewer.pal(length(tissue_cats), "Set1"), tissue_cats)
tord <- order(apply(subset_samps("rowcat_bias", c("raw", "sd"), samps = samps), 2, mean))
CI_locs <- apply(subset_samps("rowcat_bias", c("raw", "sd"), samps = samps) + samps$overall_bias, 2, quantile, probs = c(0.05, 0.95))[,tord]
posterior_means <- apply(subset_samps("row_bias", c("raw", "sd"), samps = samps) + samps$overall_bias, 2, mean)[tord]
inside_cols <- rep("black", length(tord))
inside_cols[overlaps_with_zero(CI_locs)] <- "white"
vioplot::vioplot(x = as.matrix(subset_samps("rowcat_bias", c("raw", "sd"), samps = samps) + samps$overall_bias)[,tord], T,
                 col = tissue_category_colors[tissue_cats][tord], outline=FALSE, xaxt = "n",
                 names = tissue_cats[tord], srt = 90, range = 0, xlab = "", lineCol = tissue_category_colors[tissue_cats][tord],
                 ylab = "multilevel category enrichment\n(logodds-scale)", cex.lab = 1, plotCentre = "point", colMed = 1)
tick <- seq_along(tissue_cats)
axis(1, at = tick, labels = F)
for(i in 1:length(tissue_cats)){
  segments(x0 = i, x1 = i, y0 = CI_locs[1,i], y1 = CI_locs[2,i], lwd = 4)
}
text(tick + 0.2, par("usr")[3] - 0.075, tissue_cats[tord], srt = 90, xpd = T, pos = 2)

abline(h=0,lwd=3,lty=2, col = adjustcolor(1,0.5))

#### now plot the intersect of DEGs & TWAS hits ####

#counts or props?
incl_significance <- T
incl_cell_totals <- F
trait_category_legend_below = T
use_tissue_cols_for_cols <- T
opacity_power_scaler <- 0.25
opacity_white_threshold <- 1
use_counts <- T
prop_TWAS <- T
order_by_counts <- F
order_by_posterior_enrichment <- T
group_by_tissue_type <- T
use_range_for_maginal_labels <- F

subset_to_traits <- T
if(subset_to_traits){
  trait_subset <- colnames(n_deg_sigtwas_intersect)[
    traitwise_partitions$Category[match(colnames(n_deg_sigtwas_intersect), traitwise_partitions$Tag)] %in% 
      c("Cardiometabolic", "Psychiatric-neurologic")]
  categories_represented <- unique(traitwise_partitions$Category[match(trait_subset, traitwise_partitions$Tag)])
} else {
  trait_subset <- colnames(n_deg_sigtwas_intersect)
}
categories_represented <- unique(traitwise_partitions$Category[match(trait_subset, traitwise_partitions$Tag)])

if(use_counts){
  table_to_use <- n_deg_sigtwas_intersect  
} else {
  if(prop_TWAS){
    table_to_use <- round(prop_twas_are_degs * 1000)  
  } else{
    table_to_use <- round(prop_degs_are_twas * 1000)  
  }
}

if(subset_to_traits){
  table_to_use <- table_to_use[,trait_subset]
  signif_matrix_to_use <- signif_matrix_to_use[,trait_subset]
}


if(order_by_counts){
  table_to_use <- table_to_use[,colnames(n_deg_sigtwas_intersect)]
  signif_matrix_to_use <- signif_matrix[,colnames(n_deg_sigtwas_intersect)]
} else if(order_by_posterior_enrichment) {
  trait_order <- traits[order(colbias, decreasing = T)]
  trait_order <- trait_order[trait_order %in% colnames(table_to_use)]
  table_to_use <- table_to_use[,trait_order]
  signif_matrix_to_use <- signif_matrix[,match(trait_order, colnames(signif_matrix))]
} else {
  table_to_use <- table_to_use[,colnames(prop_twas_are_degs)]
  signif_matrix_to_use <- signif_matrix[,colnames(prop_twas_are_degs)]
}

if(group_by_tissue_type){
  tissue_cats <- list(circulation = c("BLOOD", "HEART", "SPLEEN"),
                      skeletal_muscle = c("SKM-GN", "SKM-VL"),
                      adipose = c("WAT−SC"),
                      other = rev(c("ADRNL", "KIDNEY", "LUNG", "LIVER")),
                      brain = c("CORTEX", "HYPOTH", "HIPPOC"),
                      GI = c("SMLINT", "COLON"))
  tissue_cats <- rev(tissue_cats)
  disp_amount <- 0.5
  tissue_disps <- unlist(lapply(1:length(tissue_cats), function(tci) rep(disp_amount * (tci), length(tissue_cats[[tci]]))))
  tissue_cats_bars_ylocs <- cbind(start = (c(0, cumsum(unlist(lapply(tissue_cats, function(tc) length(tc))) + disp_amount)) + 2 * disp_amount)[-(length(tissue_cats)+1)], 
                                  end = cumsum(unlist(lapply(tissue_cats, function(tc) length(tc))) + disp_amount) + disp_amount)
  
}

if(use_range_for_maginal_labels){
  n_genes_in_nodes_label <- apply(apply(n_genes_in_nodes_matrix, 1, range), 2, paste0, collapse = " - ")
  sig_twas_by_trait_genes_label <- apply(apply(sig_twas_by_trait_genes_matrix, 2, range), 2, paste0, collapse = " - ")
} else {
  n_genes_in_nodes_label <- apply(n_genes_in_nodes_matrix, 1, max)
  sig_twas_by_trait_genes_label <- apply(sig_twas_by_trait_genes_matrix, 2, max)
}



if(length(salient.categories) < 9){
  category_colors <- RColorBrewer::brewer.pal(length(salient.categories), "Dark2")
  names(category_colors) <- sort(salient.categories)  
} else {
  category_colors <- RColorBrewer::brewer.pal(5, "Dark2")
  names(category_colors) <- sort(c("Cardiometabolic", "Aging", "Anthropometric", 
                                   "Immune", "Psychiatric-neurologic"))
  cats_leftover <- setdiff(salient.categories, names(category_colors))
  cols_leftover <- setdiff(RColorBrewer::brewer.pal(8, "Dark2"), category_colors)
  n_cols_needed <- length(cats_leftover) - length(cols_leftover)
  if(n_cols_needed > 0){
    extra_cols <- c(cols_leftover, RColorBrewer::brewer.pal(n_cols_needed, "Set2"))
  } else {
    extra_cols <- cols_leftover
  }
  names(extra_cols) <- cats_leftover
  category_colors <- c(category_colors, extra_cols)
}
cols$category <- category_colors

cairo_pdf(paste0("~/Documents/Documents - nikolai/pass1b_fig8_DEG-TWAS_Intersect", 
                 ifelse(use_counts, "_counts", "_permille"), 
                 ifelse(subset_to_traits, "_subset-traits", "_all-traits"),".pdf"), 
          width = 2100 / 72 * ncol(table_to_use) / 80, 
          height = 500 / 72 + ifelse(group_by_tissue_type, disp_amount * 0.75, 0), 
          family="Arial Unicode MS")
par(xpd = T, 
    mar = c(6,
            0 + ifelse(subset_to_traits, 4, 0),
            6 + ifelse(group_by_tissue_type, disp_amount * 4.5, 0),
            5.5 + ifelse(subset_to_traits, 1, 0)))
plot(1, xaxt="n",yaxt="n",bty="n",pch="",ylab="",xlab="", main="", sub="", xlim= c(-5,ncol(table_to_use)), ylim = c(-5,nrow(table_to_use)))

if(group_by_tissue_type){
  tissue_cats_bars_xlocs <- sapply(tissue_cats, function(tc) max(strwidth(tc, units = "user"))) + 0.2
  tissue_cats_bars_xlocs <- rep(max(tissue_cats_bars_xlocs), length(tissue_cats_bars_xlocs))
}


if(use_tissue_cols_for_cols){
  heatmap_cols <- sapply((1:max(table_to_use) / max(table_to_use))^opacity_power_scaler, function(opcty) 
    adjustcolor("black", opcty))
} else {
  heatmap_cols <- viridisLite::viridis(n = max(table_to_use, na.rm = T)*100+1)
  heatmap_cols <- heatmap_cols[round(log(1:max(table_to_use, na.rm = T)) / log(max(table_to_use, na.rm = T)) * max(table_to_use, na.rm = T) * 100 + 1)]
}
for(ri in 1:nrow(table_to_use)){
  text(x = 0.5, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0), pos = 2, labels = rownames(table_to_use)[ri])
  for(ci in 1:ncol(table_to_use)){
    if(ri == 1){
      text(x = ci+0.5, y = -0.9, pos = 2, srt = 45,
           labels = gsub("_Scatter", "", trait_categories$new_Phenotype[match(colnames(table_to_use)[ci], trait_categories$Tag)]))
      rect(xleft = ci + 1/2,
           xright = ci - 1/2,
           ybottom = ri - 1/2 - 1,
           ytop =  ri + 1/2 - 1,
           col = category_colors[trait_categories$Category[match(colnames(table_to_use)[ci], trait_categories$Tag)]])
    }
    
    #vertical total # options
    if(ri == nrow(table_to_use)){
      text(x = ci-0.5, y = ri+1 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), pos = 4, srt = 45,
           labels = sig_twas_by_trait_genes_label[colnames(table_to_use)[ci]])
    }
    
    #horiz total # of options
    if(ci == ncol(table_to_use)){
      text(x = ci+0.45, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0), pos = 4,
           labels = n_genes_in_nodes_label[rownames(table_to_use)[ri]])
    }
    
    #the actual cells
    if(use_tissue_cols_for_cols){
      rect(xleft = ci + 1/2,
           xright = ci - 1/2,
           ybottom = ri - 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0),
           ytop =  ri + 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0),
           col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[rownames(table_to_use)[ri]], (table_to_use[ri, ci] / max(table_to_use, na.rm = T))^opacity_power_scaler ))  
    } else {
      rect(xleft = ci + 1/2,
           xright = ci - 1/2,
           ybottom = ri - 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0),
           ytop =  ri + 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0),
           col = heatmap_cols[table_to_use[ri, ci]])
    }
    
    
    
    # rect(xleft = ci + 1/2,
    #      xright = ci - 1/2,
    #      ybottom = ri - 0.475,
    #      ytop =  ri + 0.475,
    #      col = heatmap_cols[table_to_use[ri, ci]],
    #      border = MotrpacRatTraining6moData::TISSUE_COLORS[rownames(table_to_use)[ri]])
    
    #text inside of cells
    if(table_to_use[ri, ci] != 0){
      text_in_cell <- table_to_use[ri, ci]
      if(incl_significance){
        signif_dir <- c("₋", "","⁺")[match(signif_matrix_to_use[ri, ci], -1:1)]
        text_in_cell <- paste0(text_in_cell, signif_dir)
      }
      if(use_tissue_cols_for_cols){
        text(text_in_cell, x = ci, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0), 
             col = ifelse((table_to_use[ri, ci] / max(table_to_use, na.rm = T))^opacity_power_scaler > opacity_white_threshold, "white", "black"), 
             cex = 0.85 + ifelse(subset_to_traits, 0.1, 0))  
      } else {
        text(text_in_cell, x = ci, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0), col = "white", cex = 0.85 + ifelse(subset_to_traits, 0.5, 0))
      }
      if(incl_cell_totals){
        text(sig_twas_by_trait_genes_matrix[,colnames(table_to_use)][ri, ci], 
             x = ci - 0.7375, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0) - 0.325, 
             col = "black", cex = 0.2, pos = 4, srt = 90)
        text(n_genes_in_nodes_matrix[,colnames(table_to_use)][ri, ci], 
             x = ci - 0.7, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0) - 0.4, 
             col = "black", cex = 0.2, pos = 4, srt = 0)
      }
    }
    
  }
}

#tissue category bars
if(group_by_tissue_type){
  for(bi in 1:length(tissue_cats_bars_xlocs)){
    segments(x0 = -tissue_cats_bars_xlocs[bi],
             x1 = -tissue_cats_bars_xlocs[bi],
             y0 = tissue_cats_bars_ylocs[bi,1],
             y1 = tissue_cats_bars_ylocs[bi,2],
             lwd = 3)
    bracket_length <- 0.2
    segments(x0 = -tissue_cats_bars_xlocs[bi],
             x1 = -tissue_cats_bars_xlocs[bi]+bracket_length,
             y0 = tissue_cats_bars_ylocs[bi,1],
             y1 = tissue_cats_bars_ylocs[bi,1],
             lwd = 3)
    segments(x0 = -tissue_cats_bars_xlocs[bi],
             x1 = -tissue_cats_bars_xlocs[bi]+bracket_length,
             y0 = tissue_cats_bars_ylocs[bi,2],
             y1 = tissue_cats_bars_ylocs[bi,2],
             lwd = 3)
    
    text(x = -tissue_cats_bars_xlocs[bi], y = mean(tissue_cats_bars_ylocs[bi,]), pos = 2, 
         labels = gsub("Gi", "GI", stringr::str_to_title(gsub("_", " ", names(tissue_cats)[bi]))), cex = 1.25)
  }  
}

#legend for heatmap
x_adj <- 2.25
y_adj <- 0
yb_adjust <- ifelse(incl_significance, 3, 0)
n_legend_rects_to_use <- 30
n_legend_labels_to_use <- 10
legend_yvals <- round(seq(0, max(table_to_use), length.out = n_legend_labels_to_use))
legend_ylocs <- seq(yb_adjust, 1 + nrow(table_to_use) + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), length.out = n_legend_rects_to_use)
legend_ycols <- round(seq(1, max(table_to_use), length.out = n_legend_rects_to_use))
for(i in 1:(n_legend_rects_to_use-1)){
  yb = legend_ylocs[i]
  yt = legend_ylocs[i+1]
  print(paste(c(yb,yt)))
  rect(xleft = ncol(table_to_use) + x_adj + 2 + 1/2,
       xright = ncol(table_to_use) + x_adj + 2 - 1/2,
       ybottom = yb,
       ytop =  yt,
       col = heatmap_cols[legend_ycols[i]], border = NA)
}
for(i in 1:n_legend_labels_to_use){
  text(labels = legend_yvals[i], x = ncol(table_to_use) + x_adj + 2.4, pos = 4, cex = 0.75,
       y = -0.25 + yb_adjust + legend_yvals[i] / max(table_to_use) * (1+nrow(table_to_use) - yb_adjust + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0)))
}

#overall rect for 0
rect(xleft = ncol(table_to_use) + x_adj + 2 + 1/2,
     xright = ncol(table_to_use) + x_adj + 2 - 1/2,
     ybottom = min(legend_ylocs) - diff(range(legend_ylocs))/50,
     ytop =  max(legend_ylocs))

#legend for significance
if(incl_significance){
  text(labels = paste0("Pr(Enr.) > ", (1 - signif_threshold), " : X⁺\nPr(Dep.) > ", (1 - signif_threshold), " : X₋"),
       x = ncol(table_to_use) + x_adj - 1.75,
       y = yb_adjust-3.25, pos = 4, cex = 0.75)
} else {
  text(labels = paste0("* indicate \nIHW-adj. \np-vals < 0.05 \nfrom Fisher's Exact Test"),
       x = ncol(table_to_use) + x_adj - 1.75,
       y = yb_adjust-1.25, pos = 4, cex = 0.75, )
}



#legend for trait categories
if(trait_category_legend_below){
  x_adj2 <- 0
  y_adj2 <- 2.5
  for(i in 1:length(categories_represented)){
    rect(xleft = -1/2 + i + cumsum(c(0,strwidth(categories_represented, units = "user") + 1))[i] + x_adj2,
         xright = 1/2 + i + cumsum(c(0,strwidth(categories_represented, units = "user") + 1))[i] + x_adj2,
         ybottom = -11,
         ytop = -10,
         col = category_colors[categories_represented[i]], border = 1)
    text(labels = categories_represented[i], pos = 4, y = -10.5, x = x_adj2 + 0.35 + i + cumsum(c(0,strwidth(categories_represented, units = "user") + 1))[i])
    #draw circle?
    arc(t1 = 3*pi/2, t2 = pi/2, r1 = 10.5 / 2, r2 = 10.5 / 2, center = c(0,-10.5/2), lwd = 2, res = 100, adjx = ifelse(order_by_counts, 1, 1))
    points(0, 0, pch = -9658, cex = 2)
  }
} else {
  x_adj2 <- x_adj - 2.5
  y_adj2 <- y_adj - 2.5
  for(i in 1:length(category_colors)){
    rect(xleft = 0 + i,
         xright = 1 + i,
         ybottom = -10,
         ytop = -9,
         col = category_colors[i], border = 1)
    text(labels = names(category_colors)[i], pos = 4, y = i + y_adj2 - 7, x = ncol(table_to_use) - 1/2 + x_adj2 + 1/3)
  }
}

#labels
#horiz label for total
segments(x0 = -2, x1 = 0.5, y0 = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0),
         y1 = nrow(table_to_use) + 1 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), lwd = 2)
segments(x0 = -2, x1 = 0.5, y0 = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), 
         y1 = nrow(table_to_use) + 2.5 + ifelse(use_range_for_maginal_labels, 0, -0.5) + 
           ifelse(order_by_counts, 0, -0.8) + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0) + 
           ifelse(subset_to_traits, 1, 1.5), lwd = 2)
text(x = -2, y = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), pos = 2, 
     labels = paste0(ifelse(use_range_for_maginal_labels, "total", "max"), " # of\nTWAS hits"))

#vertical label for total
segments(x0 = ncol(table_to_use) + 2, x1 = ncol(table_to_use) + 0.75, 
         y0 = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), 
         y1 = nrow(table_to_use) + 0.75 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), lwd = 2)
segments(x0 = ncol(table_to_use) + 2, x1 = ncol(table_to_use) + 1.75 + 
           ifelse(subset_to_traits, 0, 1.25), 
         y0 = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), 
         y1 = nrow(table_to_use) + 0.75 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), lwd = 2)
text(x = ncol(table_to_use) + 2 + ifelse(order_by_counts, 0, 2), y = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), pos = 3, 
     labels = paste0(ifelse(use_range_for_maginal_labels, "total", "max"), " # of DEGs"))


#legend and title labels
text(labels = ifelse(use_counts, latex2exp::TeX("$n_{intersect}$"), latex2exp::TeX("‰_{ intersect}")), pos = 3, font = 2, cex = 1.25,
     x = ncol(table_to_use) + x_adj + 2, y = nrow(table_to_use) + y_adj + 0.875 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0))
if(use_counts){
  text(latex2exp::TeX(paste0("number of genes in 8w - F↓M↓+ 8w - F↑M↑ with IHW significant TWAS at $\\alpha$ = 0.05")), 
       x = 2 + ifelse(subset_to_traits, 0, 20), 
       y = nrow(table_to_use) + 3.25 + 
         ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0) + 
         ifelse(subset_to_traits, 0, 1), pos = 4, cex = 1.5, font = 2)
} else {
  text(latex2exp::TeX(paste0("proportion (‰) of IHW significant TWAS hits at $\\alpha$ = 0.05 in 8w - F↓M↓ or 8w - F↑M↑")), 
       x = 0, y = nrow(table_to_use) + 4.25 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), pos = 4, cex = 2.35, font = 2)
}
dev.off()





#### plot a few Beta distributions ####

xr <- 0:1000/1000
dens <- dbeta(xr, 5,15)
par(mfrow = c(1,1))
plot(xr, dens, type = "l", lwd = 2, frame = F, yaxt = "n", xlab = "", ylab = "", cex.axis = 1.25)
polygon(xr, dens, lwd = 2, col = adjustcolor(1,0.2))

par(mfrow = c(3,5), mar = c(2,1,0,1))
for(i in 1:15){
  dens <- dbeta(xr, sample(5:15, 1), sample(10:30, 1))
  plot(xr, dens, type = "l", lwd = 2, frame = F, yaxt = "n", xlab = "", ylab = "", cex.axis = 1.25)
  polygon(xr, dens, lwd = 2, col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[rownames(table_to_use)[i]],0.2))
  
}

par(mfrow = c(10,20), mar = c(0.1,0.1,0.1,0.1))
for(i in 1:200){
  dens <- dbeta(xr, sample(30:40, 1), sample(80:120, 1))
  plot(xr, dens, type = "l", lwd = 0.5, frame = F, yaxt = "n", xaxt = "n", xlab = "", ylab = "")
  polygon(xr, dens, lwd = 0.5, col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[sample(rownames(table_to_use),1)],1))
  
}

#### compute binomial model counts to test for enrichment in + hits ####
twas_with_hits <- colnames(prop_degs_are_twas)
# tissue_code <- MotrpacBicQC::bic_animal_tissue_code
# tissue_code <- tissue_code[,c("tissue_name_release", "abbreviation")]
# tissue_code <- tissue_code[tissue_code$tissue_name_release != "" & !is.na(tissue_code$abbreviation)]
tissue_code <- data.frame(tissue_name_release = names(MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV),
                          abbreviation = MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV)

#get "null hypothesis" genes for comparison group
# load('/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/transcript_rna_seq_20211008.RData')
# rna_dea <- transcript_rna_seq$timewise_dea
rna_dea <- MotrpacRatTraining6mo::combine_da_results(assays = "TRNSCRPT")
rna_dea <- rna_dea[rna_dea$comparison_group == "8w" & rna_dea$selection_fdr > 0.05,]
rna_dea$tissue_abbreviation <- rna_dea$tissue
rna_dea_null <- lapply(setNames(unique(rna_dea$tissue_abbreviation), unique(rna_dea$tissue_abbreviation)), function(tiss) {
  print(tiss)
  subset <- rna_dea[rna_dea$tissue_abbreviation == tiss,]
  genes_in_sub <- table(subset$feature_ID)
  genes_in_sub <- names(genes_in_sub)[genes_in_sub == 2]
  if(length(genes_in_sub) == 0){return(NULL)}
  subset <- subset[subset$feature_ID %in% genes_in_sub,c("feature_ID", "sex", "logFC")]
  subset <- tidyr::spread(subset, key = "sex", value = "logFC")
  subset <- subset[(sign(subset$female) * sign(subset$male)) == 1,]
  x <- sign(subset$female)
  names(x) <- subset$feature_ID
  return(x)
})
do.call(rbind, lapply(rna_dea_null, function(x) table(x)))
# rm(transcript_rna_seq)
rm(rna_dea)

data <- list()
for(tissue_abbrev in names(twas_intersect)){
  
  tissue_i <- tissue_code$tissue_name_release[tissue_code$abbreviation == tissue_abbrev]
  
  print(tissue_i)
  
  sig_twas_by_trait <- lapply(setNames(twas_with_hits, twas_with_hits), 
                              function(trait_i) sig_twas_by_tissue[[tissue_i]][sig_twas_by_tissue[[tissue_i]]$trait == trait_i,])
  sig_twas_by_trait_genes <- lapply(setNames(twas_with_hits, twas_with_hits), 
                                    function(trait_i) sig_twas_by_trait[[trait_i]]$gene_name)
  sig_twas_by_trait_signs <- lapply(setNames(twas_with_hits, twas_with_hits), 
                                    function(trait_i){
                                      signs <- sign(sig_twas_by_trait[[trait_i]]$zscore)
                                      names(signs) <- sig_twas_by_trait_genes[[trait_i]]
                                      signs
                                    })
  
  timepoint = "8w"
  sex = "male"
  
  DE_genes_in_Nodes <- node_metadata_list[[timepoint]]$human_gene_symbol[
    node_metadata_list[[timepoint]]$tissue == tissue_code$abbreviation[tissue_code$tissue_name_release == tissue_i]]
  
  if(length(DE_genes_in_Nodes) == 0){next()}
  
  DE_genes_in_Nodes_sign <- cbind(as.data.frame(do.call(rbind, strsplit(node_metadata_list[[timepoint]]$node[
    node_metadata_list[[timepoint]]$tissue == tissue_code$abbreviation[tissue_code$tissue_name_release == tissue_i]], split = "_"))),
    gene = DE_genes_in_Nodes)
  
  colnames(DE_genes_in_Nodes_sign) <- c("time", "female", "male", "gene")
  DE_genes_in_Nodes_sign$female[grep(pattern = "F1", DE_genes_in_Nodes_sign$female)] <- 1
  DE_genes_in_Nodes_sign$female[grep(pattern = "F-1", DE_genes_in_Nodes_sign$female)] <- -1
  DE_genes_in_Nodes_sign$male[grep(pattern = "M1", DE_genes_in_Nodes_sign$male)] <- 1
  DE_genes_in_Nodes_sign$male[grep(pattern = "M-1", DE_genes_in_Nodes_sign$male)] <- -1
  DE_genes_in_Nodes_sign <- DE_genes_in_Nodes_sign[!is.na(DE_genes_in_Nodes_sign$gene),]
  
  DE_genes_signs <- as.integer(DE_genes_in_Nodes_sign[,sex])
  names(DE_genes_signs) <- DE_genes_in_Nodes_sign$gene
  
  nDE_genes_signs <- rna_dea_null[[tissue_abbrev]]
  names(nDE_genes_signs) <- gene_map$HUMAN_ORTHOLOG_SYMBOL[match(names(nDE_genes_signs), gene_map$RAT_ENSEMBL_ID)]
  nDE_genes_signs <- nDE_genes_signs[!is.na(names(nDE_genes_signs))]
  nDE_genes_signs <- nDE_genes_signs[!is.na(nDE_genes_signs)]
  
  results_1 <- as.data.frame(t(sapply(setNames(twas_with_hits, twas_with_hits), function(trait_i) {
    twas_signs <- sig_twas_by_trait_signs[[trait_i]]
    shared_genes <- intersect(names(DE_genes_signs), names(twas_signs))
    if(length(shared_genes) != 0){
      shared_genes_signs <- twas_signs[shared_genes] * DE_genes_signs[shared_genes]
      output <- c(count = sum(shared_genes_signs > 0), total = length(shared_genes_signs))
      return(output)
    } else {return(c(0,0))}
  })))
  
  results_2 <- as.data.frame(t(sapply(setNames(twas_with_hits, twas_with_hits), function(trait_i) {
    twas_signs <- sig_twas_by_trait_signs[[trait_i]]
    shared_genes <- intersect(names(nDE_genes_signs), names(twas_signs))
    if(length(shared_genes) != 0){
      shared_genes_signs <- twas_signs[shared_genes] * nDE_genes_signs[shared_genes]
      output <- c(count = sum(shared_genes_signs > 0), total = length(shared_genes_signs))
      return(output)
    } else {return(c(0,0))}
  })))
  
  results_1$tissue <- tissue_abbrev
  results_1$trait <- rownames(results_1)
  results_1$group <- "focal"
  
  results_2$tissue <- tissue_abbrev
  results_2$trait <- rownames(results_2)
  results_2$group <- "compl"
  
  data[[tissue_abbrev]] <- rbind(results_1, results_2)
  
}

hist(log2(unlist(lapply(data, function(x) (x$count / x$total)[x$group == "focal"] / (x$count / x$total)[x$group == "compl"]))), breaks = c(-3,0, 3))
x <- log2(unlist(lapply(data, function(x) (x$count / x$total)[x$group == "focal"] / (x$count / x$total)[x$group == "compl"])))
x <- x[x != Inf & x != -Inf & !is.na(x)]
hist(x, breaks = c(-max(abs(x)),0,max(abs(x))))
mean(x > 0)

#### compute binomial model counts to test for enrichment in + hits, alternate version ####
twas_with_hits <- colnames(prop_degs_are_twas)
# tissue_code <- MotrpacBicQC::bic_animal_tissue_code
# tissue_code <- tissue_code[,c("tissue_name_release", "abbreviation")]
# tissue_code <- tissue_code[tissue_code$tissue_name_release != "" & !is.na(tissue_code$abbreviation)]
tissue_code <- data.frame(tissue_name_release=names(MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV),
                          abbreviation=MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV)
possible_genes <- intersect(all_orthologs_tested, all_twas_genes_tested)
compatible_twas_genes <- some.twas$gene_name %in% possible_genes
twas_by_tissue <- lapply(setNames(unique(some.twas$tissue), MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[unique(some.twas$tissue)]), function(tiss) {
  print(tiss)
  some.twas[some.twas$tissue == tiss & compatible_twas_genes,]
})

#get "null hypothesis" genes for overall probabilities
# load('/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/transcript_rna_seq_20211008.RData')
# rna_dea <- transcript_rna_seq$timewise_dea
rna_dea <- MotrpacRatTraining6mo::combine_da_results(assays = "TRNSCRPT")
rna_dea <- rna_dea[rna_dea$comparison_group == "8w",]
rna_dea$tissue_abbreviation <- rna_dea$tissue
rna_dea_null <- lapply(setNames(unique(rna_dea$tissue_abbreviation), unique(rna_dea$tissue_abbreviation)), function(tiss) {
  print(tiss)
  subset <- rna_dea[rna_dea$tissue_abbreviation == tiss,]
  subset$gene <- gene_map$HUMAN_ORTHOLOG_SYMBOL[match(subset$feature_ID, gene_map$RAT_ENSEMBL_ID)]
  subset <-subset[subset$gene %in% possible_genes,]
  genes_in_sub <- table(subset$feature_ID)
  genes_in_sub <- names(genes_in_sub)[genes_in_sub == 2]
  if(length(genes_in_sub) == 0){return(NULL)}
  subset <- subset[subset$feature_ID %in% genes_in_sub,c("feature_ID", "sex", "logFC")]
  subset <- tidyr::spread(subset, key = "sex", value = "logFC")
  subset <- subset[(sign(subset$female) * sign(subset$male)) == 1,]
  x <- sign(subset$female)
  names(x) <- gene_map$HUMAN_ORTHOLOG_SYMBOL[match(subset$feature_ID, gene_map$RAT_ENSEMBL_ID)]
  return(x)
})
do.call(rbind, lapply(rna_dea_null, function(x) table(x)))
rm(rna_dea)

tissues <- tissue_code$abbreviation[match(names(motrpac_gtex_map), tissue_code$tissue_name_release)]
tissues <- setdiff(tissues, c("HYPOTH", "TESTES", "OVARY"))
data_null <- lapply(setNames(tissues, tissues), function(tiss){
  print(tiss)
  tws <- twas_by_tissue[[tiss]]
  motr <- rna_dea_null[[tiss]]
  if(is.null(motr) | is.null(tws)){
    return(NULL)
  }
  compatible_genes <- intersect(names(motr), tws$gene_name)
  motr <- motr[compatible_genes]
  tws <- tws[tws$gene_name %in% compatible_genes,]
  tws$motrpac_sign <- motr[tws$gene_name]
  tws$twas_sign <- sign(tws$zscore)
  tws <- tws[tws$twas_sign != 0,]
  output <- as.data.frame(t(sapply(setNames(twas_with_hits, twas_with_hits), function(trait_i) {
    trait_subset <- tws[tws$trait == trait_i, c("motrpac_sign", "twas_sign")]
    c(count_twas = sum(trait_subset$twas_sign > 0), count_motr = sum(trait_subset$motrpac_sign > 0), total = nrow(trait_subset))
  })))
  output$trait <- rownames(output)
  output$tissue <- tiss
  rownames(output) <- NULL
  output
})
data_null <- do.call(rbind, data_null)
hist(data_null$count_twas / data_null$total)
hist(data_null$count_motr / data_null$total)

adj_pvalue_alpha <- 0.05
sig_twas_by_tissue <- lapply(setNames(names(motrpac_gtex_map), names(motrpac_gtex_map)), function(tissue){
  inds <- which(some.twas$tissue == tissue & some.twas$adj_pvalue < adj_pvalue_alpha)
  some.twas[inds,]
})
data_not_null <- lapply(setNames(tissues, tissues), function(tissue_abbrev){
  tissue_i <- tissue_code$tissue_name_release[tissue_code$abbreviation == tissue_abbrev]
  
  print(tissue_i)
  
  sig_twas_by_trait <- lapply(setNames(twas_with_hits, twas_with_hits), 
                              function(trait_i) sig_twas_by_tissue[[tissue_i]][sig_twas_by_tissue[[tissue_i]]$trait == trait_i,])
  sig_twas_by_trait_genes <- lapply(setNames(twas_with_hits, twas_with_hits), 
                                    function(trait_i) sig_twas_by_trait[[trait_i]]$gene_name)
  sig_twas_by_trait_signs <- lapply(setNames(twas_with_hits, twas_with_hits), 
                                    function(trait_i){
                                      signs <- sign(sig_twas_by_trait[[trait_i]]$zscore)
                                      names(signs) <- sig_twas_by_trait_genes[[trait_i]]
                                      signs
                                    })
  
  timepoint = "8w"
  sex = "male"
  
  DE_genes_in_Nodes <- node_metadata_list[[timepoint]]$human_gene_symbol[
    node_metadata_list[[timepoint]]$tissue == tissue_code$abbreviation[tissue_code$tissue_name_release == tissue_i]]
  
  if(length(DE_genes_in_Nodes) == 0){return(NULL)}
  
  DE_genes_in_Nodes_sign <- cbind(as.data.frame(do.call(rbind, strsplit(node_metadata_list[[timepoint]]$node[
    node_metadata_list[[timepoint]]$tissue == tissue_code$abbreviation[tissue_code$tissue_name_release == tissue_i]], split = "_"))),
    gene = DE_genes_in_Nodes)
  
  colnames(DE_genes_in_Nodes_sign) <- c("time", "female", "male", "gene")
  DE_genes_in_Nodes_sign$female[grep(pattern = "F1", DE_genes_in_Nodes_sign$female)] <- 1
  DE_genes_in_Nodes_sign$female[grep(pattern = "F-1", DE_genes_in_Nodes_sign$female)] <- -1
  DE_genes_in_Nodes_sign$male[grep(pattern = "M1", DE_genes_in_Nodes_sign$male)] <- 1
  DE_genes_in_Nodes_sign$male[grep(pattern = "M-1", DE_genes_in_Nodes_sign$male)] <- -1
  DE_genes_in_Nodes_sign <- DE_genes_in_Nodes_sign[!is.na(DE_genes_in_Nodes_sign$gene),]
  
  DE_genes_signs <- as.integer(DE_genes_in_Nodes_sign[,sex])
  names(DE_genes_signs) <- DE_genes_in_Nodes_sign$gene
  
  output <- as.data.frame(t(sapply(setNames(twas_with_hits, twas_with_hits), function(trait_i) {
    twas_signs <- sig_twas_by_trait_signs[[trait_i]]
    shared_genes <- intersect(names(DE_genes_signs), names(twas_signs))
    if(length(shared_genes) != 0){
      shared_genes_signs <- twas_signs[shared_genes] * DE_genes_signs[shared_genes]
      output <- c(count = sum(shared_genes_signs > 0), total = length(shared_genes_signs))
      return(output)
    } else {return(c(0,0))}
  })))
  
  output$trait <- rownames(output)
  output$tissue <- tissue_abbrev
  rownames(output) <- NULL
  output
})
data_not_null <- do.call(rbind, data_not_null)

#confirm matchup
all(paste0(data_null$trait, "~", data_null$tissue) == paste0(data_not_null$trait, "~", data_not_null$tissue))

#some quick eda
plot(data_null$count_twas / data_null$total * data_null$count_motr / data_null$total + 
       (1-data_null$count_twas / data_null$total) * (1-data_null$count_motr / data_null$total),
     (data_not_null$count+1) / (data_not_null$total+2), pch = 19, col = adjustcolor(1, 0.5),
     cex = data_not_null$total / max(data_not_null$total) * 5,
     xlab = "expected proportion positives under null",
     ylab = "posterior mean proportion updated from Beta(1,1)")
# points((data_null$count_twas / data_null$total * data_null$count_motr / data_null$total + 
#        (1-data_null$count_twas / data_null$total) * (1-data_null$count_motr / data_null$total))[c(3, 403, 563, 803)],
#      (data_not_null$count / data_not_null$total)[c(3, 403, 563, 803)], pch = 19, col = adjustcolor(3, 0.5),
#      cex = (data_not_null$total / max(data_not_null$total) * 5)[c(3, 403, 563, 803)])
abline(h = 0.5, lty = 2, col = adjustcolor(2, 0.75), lwd = 2)
abline(v = 0.5, lty = 2, col = adjustcolor(2, 0.75), lwd = 2)

#check iid flat beta model
basic_posteriors_masses <- 1 - pbeta(q = data_null$count_twas / data_null$total * data_null$count_motr / data_null$total + 
                                       (1-data_null$count_twas / data_null$total) * (1-data_null$count_motr / data_null$total),
                                     shape1 = 1 + data_not_null$count, shape2 = 1 + data_not_null$total - data_not_null$count)
sum(basic_posteriors_masses > 0.95)
sum(basic_posteriors_masses < 0.05)
head(data_not_null[order(abs(basic_posteriors_masses - 0.5), decreasing = T),], 20)
table(data_not_null[which(basic_posteriors_masses > 0.95),"trait"])
table(data_not_null[which(basic_posteriors_masses < 0.05),"trait"])

hist(sapply(tissues, function(tiss) diff(range(data_null$count_twas[data_null$tissue == tiss] / data_null$total[data_null$tissue == tiss]))))
hist(sapply(traits, function(trait) diff(range(data_null$count_twas[data_null$trait == trait] / data_null$total[data_null$trait == trait]))))
hist(sapply(tissues, function(tiss) diff(range(data_null$count_motr[data_null$tissue == tiss] / data_null$total[data_null$tissue == tiss]))))
hist(sapply(traits, function(trait) diff(range(data_null$count_motr[data_null$trait == trait] / data_null$total[data_null$trait == trait]))))
trait
data_null$count_motr[data_null$trait == trait] / data_null$total[data_null$trait == trait]
data_null[data_null$trait == trait,]$count_motr / data_null[data_null$trait == trait,]$total
#potential source of heterogeneity here! 


#bayesian model
traits <- unique(data_null$trait)
tissues <- unique(data_null$tissue)
d <- list(intersect_count = data_not_null$count,
          intersect_total = data_not_null$total,
          twas_count = data_null$count_twas,
          motr_count = data_null$count_motr,
          twas_motr_total = data_null$total,
          trait = match(data_null$trait, traits),
          tissue = match(data_null$tissue, tissues),
          n_trait = length(traits),
          n_tissue = length(tissues)
)

# stan_program <- '
# data {
#     int<lower=1> n_trait;
#     int<lower=1> n_tissue;
#     int<lower=1,upper=n_trait> trait[n_trait * n_tissue];
#     int<lower=1,upper=n_tissue> tissue[n_trait * n_tissue];
#     int<lower=0> intersect_count[n_trait * n_tissue];
#     int<lower=0> intersect_total[n_trait * n_tissue];
#     int<lower=0> twas_count[n_trait * n_tissue];
#     int<lower=0> motr_count[n_trait * n_tissue];
#     int<lower=0> twas_motr_total[n_trait * n_tissue];
# }
# transformed data {
#     int<lower=1> n = n_trait * n_tissue;
# }
# parameters {
#     //main parameters
#     real trait_mean;
#     real<lower=0> trait_sd;
#     vector[n_trait] raw_trait_logodds;
# 
#     real tissue_mean;
#     real<lower=0> tissue_sd;
#     vector[n_tissue] raw_tissue_logodds;
#     
#     vector[n] raw_intersect_logodds;
#     real<lower=0> intersect_sd;
#     
#     //biases in deviations terms
#     vector[n_tissue] tissue_bias_mean;
#     vector<lower=0>[n_tissue] tissue_bias_sd;
#     vector[n_trait] trait_bias_mean;
#     vector<lower=0>[n_trait] trait_bias_sd;
#     
#     //multilevel deviation term params
#     real<lower=0> tissue_bias_mean_sd;
#     real<lower=0> trait_bias_mean_sd;
# }
# transformed parameters {
#     //recenter params
#     vector[n_trait] trait_logodds = raw_trait_logodds * trait_sd + trait_mean;
#     vector[n_tissue] tissue_logodds = raw_tissue_logodds * tissue_sd + tissue_mean;
#     
#     vector[n_trait] neg_trait_logodds = logit(1 - inv_logit(trait_logodds));
#     vector[n_tissue] neg_tissue_logodds = logit(1 - inv_logit(tissue_logodds));
#     
#     vector[n] intersect_mean = logit(inv_logit(trait_logodds[trait]) .* inv_logit(tissue_logodds[tissue]) +
#                                      inv_logit(neg_trait_logodds[trait]) .* inv_logit(neg_tissue_logodds[tissue]));
#     vector[n] intersect_logodds = raw_intersect_logodds * intersect_sd .* tissue_bias_sd[tissue] .* trait_bias_sd[trait] +
#                                   intersect_mean + tissue_bias_mean[tissue] * tissue_bias_mean_sd + 
#                                   trait_bias_mean[trait] * trait_bias_mean_sd;
# }
# model {
#     //priors / hyperpriors
#     raw_trait_logodds ~ std_normal();
#     trait_mean ~ normal(0,2);
#     trait_sd ~ normal(0,2);
#     
#     raw_tissue_logodds ~ std_normal();
#     tissue_mean ~ normal(0,2);
#     tissue_sd ~ normal(0,2);
#     
#     raw_intersect_logodds ~ std_normal();
#     intersect_sd ~ normal(0,2);
#     
#     tissue_bias_mean ~ std_normal();
#     tissue_bias_sd ~ std_normal();
#     trait_bias_mean ~ std_normal();
#     trait_bias_sd ~ std_normal();
#     tissue_bias_mean_sd ~ std_normal();
#     trait_bias_mean_sd ~ std_normal();
#     
#     //likelihood
#     twas_count ~ binomial_logit(twas_motr_total, trait_logodds[trait]);
#     motr_count ~ binomial_logit(twas_motr_total, tissue_logodds[tissue]);
#     intersect_count ~ binomial_logit(intersect_total, intersect_logodds);
# }
# generated quantities {
#     vector[n] difference_in_props = inv_logit(intersect_logodds) - inv_logit(intersect_mean);
# }
# '

# stan_program <- '
# data {
#     int<lower=1> n_trait;
#     int<lower=1> n_tissue;
#     int<lower=1,upper=n_trait> trait[n_trait * n_tissue];
#     int<lower=1,upper=n_tissue> tissue[n_trait * n_tissue];
#     int<lower=0> intersect_count[n_trait * n_tissue];
#     int<lower=0> intersect_total[n_trait * n_tissue];
#     int<lower=0> twas_count[n_trait * n_tissue];
#     int<lower=0> motr_count[n_trait * n_tissue];
#     int<lower=0> twas_motr_total[n_trait * n_tissue];
# }
# transformed data {
#     int<lower=1> n = n_trait * n_tissue;
# }
# parameters {
#     vector[n] raw_intersect_logodds;
#     real intersect_logodds_log_sd;
# }
# transformed parameters {
#     vector[n] intersect_logodds = raw_intersect_logodds * exp(intersect_logodds_log_sd);
# 
# }
# model {
#     raw_intersect_logodds ~ std_normal();
#     intersect_logodds_log_sd ~ normal(0,2);
#     intersect_count ~ binomial_logit(intersect_total, intersect_logodds);
# }
# generated quantities {
#     vector[n] difference_in_props = inv_logit(intersect_logodds) - 0.5;
# }
# '

stan_program <- '
data {
    int<lower=1> n_trait;
    int<lower=1> n_tissue;
    int<lower=1,upper=n_trait> trait[n_trait * n_tissue];
    int<lower=1,upper=n_tissue> tissue[n_trait * n_tissue];
    int<lower=0> intersect_count[n_trait * n_tissue];
    int<lower=0> intersect_total[n_trait * n_tissue];
    int<lower=0> twas_count[n_trait * n_tissue];
    int<lower=0> motr_count[n_trait * n_tissue];
    int<lower=0> twas_motr_total[n_trait * n_tissue];
}
transformed data {
    int<lower=1> n = n_trait * n_tissue;
}
parameters {
    //main parameters
    real trait_mean;
    real<lower=0> trait_sd;
    vector[n_trait] raw_trait_logodds;
    real<lower=0> within_trait_sd;
    vector[n] raw_within_trait_logodds;

    real tissue_mean;
    real<lower=0> tissue_sd;
    vector[n_tissue] raw_tissue_logodds;
    real<lower=0> within_tissue_sd;
    vector[n] raw_within_tissue_logodds;

    vector[n] raw_intersect_logodds;
    real<lower=0> intersect_sd;
    
    //biases in deviations terms
    vector[n_tissue] tissue_bias_mean;
    vector[n_tissue] tissue_bias_log_sd;
    vector[n_trait] trait_bias_mean;
    vector[n_trait] trait_bias_log_sd;
    
    //multilevel deviation term params
    real<lower=0> tissue_bias_log_sd_sd;
    real<lower=0> trait_bias_log_sd_sd;
    real<lower=0> tissue_bias_mean_sd;
    real<lower=0> trait_bias_mean_sd;
}
transformed parameters {
    //recenter params
    vector[n_trait] trait_logodds = raw_trait_logodds * trait_sd + trait_mean;
    vector[n_tissue] tissue_logodds = raw_tissue_logodds * tissue_sd + tissue_mean;
    
    vector[n] within_trait_logodds = raw_within_trait_logodds * within_trait_sd + trait_logodds[trait];
    vector[n] within_tissue_logodds = raw_within_tissue_logodds * within_tissue_sd + tissue_logodds[tissue];
    
    vector<lower=0, upper=1>[n] within_trait_prob = inv_logit(within_trait_logodds);
    vector<lower=0, upper=1>[n] within_tissue_prob = inv_logit(within_tissue_logodds);
    
    vector[n] intersect_mean = logit(within_trait_prob .* within_tissue_prob +
                                     (1-within_trait_prob) .* (1-within_tissue_prob));
    vector[n] intersect_logodds = raw_intersect_logodds * intersect_sd .* exp(tissue_bias_log_sd[tissue] * tissue_bias_log_sd_sd) .* 
                                  exp(trait_bias_log_sd[trait] * trait_bias_log_sd_sd) +
                                  intersect_mean + tissue_bias_mean[tissue] * tissue_bias_mean_sd + 
                                  trait_bias_mean[trait] * trait_bias_mean_sd;
}
model {
    //priors and hyperpriors
    raw_trait_logodds ~ std_normal();
    trait_mean ~ normal(0,2);
    trait_sd ~ normal(0,2);
    raw_within_trait_logodds ~ std_normal();
    within_trait_sd ~ normal(0,2);
    
    raw_tissue_logodds ~ std_normal();
    tissue_mean ~ normal(0,2);
    tissue_sd ~ normal(0,2);
    raw_within_tissue_logodds ~ std_normal();
    within_tissue_sd ~ normal(0,2);
    
    raw_intersect_logodds ~ std_normal();
    intersect_sd ~ normal(0,2);
    
    tissue_bias_mean ~ std_normal();
    tissue_bias_log_sd ~ std_normal();
    tissue_bias_log_sd_sd ~ std_normal();
    trait_bias_mean ~ std_normal();
    trait_bias_log_sd ~ std_normal();
    trait_bias_log_sd_sd ~ std_normal();
    tissue_bias_mean_sd ~ std_normal();
    trait_bias_mean_sd ~ std_normal();
    
    //likelihood
    twas_count ~ binomial_logit(twas_motr_total, within_trait_logodds);
    motr_count ~ binomial_logit(twas_motr_total, within_tissue_logodds);
    intersect_count ~ binomial_logit(intersect_total, intersect_logodds);
}
generated quantities {
    vector[n] difference_in_props = inv_logit(intersect_logodds) - inv_logit(intersect_mean);
}
'

fit_model <- F
if(fit_model){
  
  #compile model
  if(!exists("curr_stan_program") || stan_program != curr_stan_program){
    curr_stan_program <- stan_program
    f <- write_stan_file(stan_program)
  }
  mod <- cmdstan_model(f)
  
  #write model
  write_stan_file(stan_program, dir = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/mcmc_output/", basename = paste0(base, ifelse(use_all_cats, "_allCats", "")))
  write_stan_json(d, paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/mcmc_output/", paste0(base, ifelse(use_all_cats, "_allCats", "")),".json"))
  
  #fit model
  out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = d, parallel_chains = 4, 
                    adapt_delta = 0.9, refresh = 50, init = 0.1, max_treedepth = 20, thin = 5)
  # out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = d, parallel_chains = 4, adapt_delta = 0.9, refresh = 10, init = 0.1, max_treedepth = 15)
  summ <- out$summary()
  summ[order(summ$ess_bulk),]
  summ[order(summ$rhat, decreasing = T),]
  save(out, file = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/mcmc_output/", paste0(base, ifelse(use_all_cats, "_allCats", "")),".cmdStanR.fit"))
} else {
  load(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/mcmc_output/", paste0(base, ifelse(use_all_cats, "_allCats", "")),".cmdStanR.fit"))
}

samps <- data.frame(as_draws_df(out$draws()))

#take a look at output
hist(apply(samps[,grep("difference_in_props", colnames(samps))], 2, mean))
prop_greater_than_0 <- function(x) mean(x>0)
sum((apply(samps[,grep("difference_in_props", colnames(samps))], 2, prop_greater_than_0)) > 0.95)
sum((apply(samps[,grep("difference_in_props", colnames(samps))], 2, prop_greater_than_0)) < 0.05)
hist(invlogit(samps$intersect_logodds.2.) - invlogit(samps$intersect_mean.2.))
data_not_null[which((apply(samps[,grep("difference_in_props", colnames(samps))], 2, prop_greater_than_0)) > 0.95),]
hist(invlogit(samps$intersect_mean.1.))
hist(invlogit(samps$intersect_logodds.1.))

#compare estimated mean probs t

sum((apply(samps[,grep("tissue_bias_mean\\.", colnames(samps))], 2, prop_greater_than_0)) > 0.95)
traits[which((apply(samps[,grep("trait_bias_mean\\.", colnames(samps))], 2, prop_greater_than_0)) > 0.95)]

plot(data_null$count_twas / data_null$total,
     apply(invlogit(samps[,setdiff(grep("within_trait_logodds", colnames(samps)), grep("raw", colnames(samps)))]), 2, mean),
     pch = 19, col = adjustcolor(match(data_null$trait, traits), 0.75))
abline(0,1,col=2,lty=2)

plot(data_null$count_motr / data_null$total,
     apply(invlogit(samps[,setdiff(grep("within_tissue_logodds", colnames(samps)), grep("raw", colnames(samps)))]), 2, mean))
abline(0,1,col=2,lty=2)

plot((data_not_null$count + 1) / (data_not_null$total + 2), 
     apply(invlogit(samps[,setdiff(grep("intersect_logodds", colnames(samps)), grep("raw", colnames(samps)))]), 2, mean),
     col = adjustcolor(1, 0.25), pch = 19)
points((data_not_null$count / data_not_null$total)[which((apply(samps[,grep("difference_in_props", colnames(samps))], 2, prop_greater_than_0)) > 0.95)], 
       (apply(invlogit(samps[,setdiff(grep("intersect_logodds", colnames(samps)), grep("raw", colnames(samps)))]), 2, mean))[which((apply(samps[,grep("difference_in_props", colnames(samps))], 2, prop_greater_than_0)) > 0.95)],
       col = adjustcolor(3, 1), pch = 19)
abline(0,1,col=2,lty=2)

hist(apply(invlogit(samps[,setdiff(grep("intersect_logodds", colnames(samps)), grep("raw", colnames(samps)))]), 2, mean))
hist(apply(invlogit(samps[,setdiff(grep("intersect_logodds", colnames(samps)), grep("raw", colnames(samps)))]), 2, mean))
data_not_null[which(abs(apply(invlogit(samps[,setdiff(grep("intersect_logodds", colnames(samps)), grep("raw", colnames(samps)))]), 2, mean)-0.5) > 0.1),]

hist(apply(invlogit(samps[,setdiff(grep("intersect_mean", colnames(samps)), grep("raw", colnames(samps)))]), 2, mean))

#### make a summary table for Nicole ####
twas_intersect <- lapply(setNames(rownames(n_deg_sigtwas_intersect), rownames(n_deg_sigtwas_intersect)), function(x) NULL)
twas_with_hits <- colnames(prop_degs_are_twas)
prop_pos <- lapply(setNames(rownames(n_deg_sigtwas_intersect), rownames(n_deg_sigtwas_intersect)), function(i) NULL)
# tissue_code <- MotrpacBicQC::bic_animal_tissue_code
# tissue_code <- tissue_code[,c("tissue_name_release", "abbreviation")]
# tissue_code <- tissue_code[tissue_code$tissue_name_release != "" & !is.na(tissue_code$abbreviation)]
tissue_code <- data.frame(tissue_name_release = names(MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV),
                          abbreviation = MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV)
for(tissue_abbrev in names(twas_intersect)){
  
  tissue_i <- tissue_code$tissue_name_release[tissue_code$abbreviation == tissue_abbrev]
  
  print(tissue_i)
  
  sig_twas_by_trait <- lapply(setNames(twas_with_hits, twas_with_hits), 
                              function(trait_i) sig_twas_by_tissue[[tissue_i]][sig_twas_by_tissue[[tissue_i]]$trait == trait_i,])
  sig_twas_by_trait_genes <- lapply(setNames(twas_with_hits, twas_with_hits), 
                                    function(trait_i) sig_twas_by_trait[[trait_i]]$gene_name)
  sig_twas_by_trait_signs <- lapply(setNames(twas_with_hits, twas_with_hits), 
                                    function(trait_i){
                                      signs <- sign(sig_twas_by_trait[[trait_i]]$zscore)
                                      names(signs) <- sig_twas_by_trait_genes[[trait_i]]
                                      signs
                                    })
  
  timepoint = "8w"
  sex = "male"
  
  DE_genes_in_Nodes <- node_metadata_list[[timepoint]]$human_gene_symbol[
    node_metadata_list[[timepoint]]$tissue == tissue_code$abbreviation[tissue_code$tissue_name_release == tissue_i]]
  
  if(length(DE_genes_in_Nodes) == 0){next()}
  
  DE_genes_in_Nodes_sign <- cbind(as.data.frame(do.call(rbind, strsplit(node_metadata_list[[timepoint]]$node[
    node_metadata_list[[timepoint]]$tissue == tissue_code$abbreviation[tissue_code$tissue_name_release == tissue_i]], split = "_"))),
    gene = DE_genes_in_Nodes)
  
  colnames(DE_genes_in_Nodes_sign) <- c("time", "female", "male", "gene")
  DE_genes_in_Nodes_sign$female[grep(pattern = "F1", DE_genes_in_Nodes_sign$female)] <- 1
  DE_genes_in_Nodes_sign$female[grep(pattern = "F-1", DE_genes_in_Nodes_sign$female)] <- -1
  DE_genes_in_Nodes_sign$male[grep(pattern = "M1", DE_genes_in_Nodes_sign$male)] <- 1
  DE_genes_in_Nodes_sign$male[grep(pattern = "M-1", DE_genes_in_Nodes_sign$male)] <- -1
  DE_genes_in_Nodes_sign <- DE_genes_in_Nodes_sign[!is.na(DE_genes_in_Nodes_sign$gene),]
  
  DE_genes_signs <- as.integer(DE_genes_in_Nodes_sign[,sex])
  names(DE_genes_signs) <- DE_genes_in_Nodes_sign$gene
  
  results <- lapply(setNames(twas_with_hits, twas_with_hits), function(trait_i) {
    twas_signs <- sig_twas_by_trait_signs[[trait_i]]
    shared_genes <- intersect(names(DE_genes_signs), names(twas_signs))
    if(length(shared_genes) != 0){
      shared_genes_signs <- twas_signs[shared_genes] * DE_genes_signs[shared_genes]
      shared_genes_rat <- gene_map$RAT_SYMBOL[match(shared_genes, gene_map$HUMAN_ORTHOLOG_SYMBOL)]
      shared_genes_rat_ensembl_id <- gene_map$RAT_ENSEMBL_ID[match(shared_genes, gene_map$HUMAN_ORTHOLOG_SYMBOL)]
      
      enrichment_test <- signif_df[which(signif_df$tissue == tissue_abbrev & signif_df$trait == trait_i),]
      
      output <- data.frame(tissue = tissue_abbrev,
                           trait = trait_i,
                           human_ortholog_symbol = shared_genes,
                           rat_symbol = shared_genes_rat,
                           rat_ensembl_ID = shared_genes_rat_ensembl_id,
                           direction_of_DE_on_trait = shared_genes_signs,
                           Pr_Geneset_DiffProp_is_Pos = enrichment_test$prob_diff_is_positive,
                           geneset_enrichment = enrichment_test$signif)
      return(output)
    } else {return(NULL)}
  })
  
  results <- do.call(rbind, results)
  
  twas_intersect[[tissue_abbrev]] <- results
  
}

twas_intersect <- do.call(rbind, twas_intersect)
rownames(twas_intersect) <- NULL
twas_intersect <- as.data.table(twas_intersect)
fwrite(x = twas_intersect, file = "~/data/twas_intersect_table_uncalibrated_Pr.txt")
fread(file = "~/data/twas_intersect_table_uncalibrated_Pr.txt")

# #### add in "all tissues" group ####
# 
# str(deg_sigtwas_proportion)
# 
# all_tissues_array <- deg_sigtwas_proportion[1,,,,]
# 
# list(names(motrpac_gtex_map), twas_with_hits, paste0(2^(0:3), "w"), c("p", "n", "genes"), c("male", "female"))
# 
# for(tissue_i in names(motrpac_gtex_map)[-1]){
#   
#   for(timepoint in paste0(2^(0:3), "w")){
#     
#     for(sex in c("male", "female")){
#     
#       all_tissues_array[tissue_i,,,"n",] <- all_tissues_array[,,"n",] + deg_sigtwas_proportion[tissue_i,,,"n",]
#             
#     }
#   }
# }
# 
# deg_sigtwas_proportion <- abind::abind(deg_sigtwas_proportion, all_tissues_array, along = 1)
# deg_sigtwas_proportion <- deg_sigtwas_proportion[-18,,,,]




#### the actual trajectory plotting ####
deg_sigtwas_proportion[,twas_with_hits[grep(twas_with_hits, pattern = "type_1_diabetes", ignore.case = T)],,"p",sex]
deg_sigtwas_proportion[,twas_with_hits[grep(twas_with_hits, pattern = "reported_hypertension", ignore.case = T)],,"p",sex]
deg_sigtwas_proportion[,twas_with_hits[grep(twas_with_hits, pattern = "heart_problems", ignore.case = T)],,"n","male"]
deg_sigtwas_proportion[,twas_with_hits[grep(twas_with_hits, pattern = "Body_mass_index_BMI", ignore.case = T)],,"n",sex]
deg_sigtwas_proportion[,twas_with_hits[grep(twas_with_hits, pattern = "Body_fat_percentage", ignore.case = T)],,"p",sex]
deg_sigtwas_proportion[,twas_with_hits[grep(twas_with_hits, pattern = "Body_fat_percentage", ignore.case = T)],,"p",sex]
deg_sigtwas_proportion[,twas_with_hits[grep(twas_with_hits, pattern = "UKB_20002_1462_self_reported_crohns_disease", ignore.case = T)],,"p",sex]
deg_sigtwas_proportion[,twas_with_hits[grep(twas_with_hits, pattern = "tag.evrsmk.tbl", ignore.case = T)],,"p",sex]

deg_sigtwas_proportion[,twas_with_hits[grep(twas_with_hits, pattern = "UKB_50_Standing_height", ignore.case = T)],,"n",sex]

#plot lines for proportions
trait <- twas_with_hits[grep(twas_with_hits, pattern = "reported_hypertension", ignore.case = T)][1]
trait <- twas_with_hits[grep(twas_with_hits, pattern = "heart_problems", ignore.case = T)][1]
trait <- twas_with_hits[grep(twas_with_hits, pattern = "Body_mass_index_BMI", ignore.case = T)][1]

trait <- twas_with_hits[grep(twas_with_hits, pattern = "ldl", ignore.case = T)][1]
trait <- twas_with_hits[grep(twas_with_hits, pattern = "CAD", ignore.case = T)][1]
trait <- twas_with_hits[grep(twas_with_hits, pattern = "alzhei", ignore.case = T)][1]

trait <- twas_with_hits[grep(twas_with_hits, pattern = "Body_fat_percentage", ignore.case = T)][1]
trait <- twas_with_hits[grep(twas_with_hits, pattern = "cholesterol", ignore.case = T)][1]
trait <- twas_with_hits[grep(twas_with_hits, pattern = "asthma", ignore.case = T)][1]

trait_patt <- trait_patts[4]
trait <- intersect(twas_with_hits, trait_categories$Tag[grep(trait_categories$new_Phenotype, pattern = trait_patt, ignore.case = T)])[1]

cols = list(Tissue=MotrpacRatTraining6moData::TISSUE_COLORS[names(motrpac_gtex_map)], 
            Time=MotrpacRatTraining6moData::GROUP_COLORS[paste0(c(1,2,4,8), "w")],
            Sex=MotrpacRatTraining6moData::SEX_COLORS[c('male','female')])
tissue_names <- sapply(strsplit(names(cols$Tissue), "-"), function(x) paste0(x[ifelse(length(x) == 2, c(2), list(2:3))[[1]]], collapse = " "))

#plotting params
plot_gene_names <- T
trait_patts <- c("triglycerides", "rheumatoid_arthritis", "high_cholesterol", "asthma_ukbs")
# trait_patts <- c("triglycerides")


cairo_pdf(paste0("~/Documents/Documents - nikolai/pass1b_fig8_DE_protective_effect.pdf"), 
          width = 1400 / 72, height = 500 / 72 * length(trait_patts), family="Arial Unicode MS", pointsize = 18)
par(mfrow = c(length(trait_patts),2), xpd = NA, mar = c(5,5,5,17))
lwd <- 3
tissue_name_cex <- 1.2

for(trait_patt in trait_patts){
  trait <- intersect(twas_with_hits, trait_categories$Tag[grep(trait_categories$new_Phenotype, pattern = trait_patt, ignore.case = T)])[1]
  
  
  for(sex in c("male", "female")){
    
    #retrieve data
    d <- apply(deg_sigtwas_proportion[, trait,,"p",sex], 2, as.numeric)
    dn <-  apply(deg_sigtwas_proportion[, trait,,"n",sex], 2, as.numeric)
    colnames(d) <- colnames(dn) <- colnames(deg_sigtwas_proportion[, trait,,"p",sex])
    rownames(d) <- rownames(dn) <- rownames(deg_sigtwas_proportion[, trait,,"p",sex])
    dn <- dn[!apply(apply(d, 1, is.na), 2, all),]
    d <- d[!apply(apply(d, 1, is.na), 2, all),]
    d <- d[!is.nan(d[,"8w"]),]
    d[is.nan(d)] <- 0.5
    dg <- deg_sigtwas_proportion[, trait,,"genes",sex]
    dg <- lapply(dg[,"8w"], function(x) strsplit(x, " ~ ")[[1]])
    
    #iterate through d to make identical lines parallel
    line_thickness <- lwd / 96 / par("pin")[2] * (par("usr")[4] - par("usr")[3])
    need_to_increment <- matrix(T, nrow = nrow(d), ncol = ncol(d), dimnames = list(rownames(d), colnames(d)))
    while(any(need_to_increment)){
      need_to_increment <- matrix(F, nrow = nrow(d), ncol = ncol(d), dimnames = list(rownames(d), colnames(d)))
      for(coi in 1:(ncol(d)-1)){
        unchanging_tissues <- rownames(unique(d[,c(coi,coi+1)]))
        need_to_increment[setdiff(rownames(d), unchanging_tissues),c(coi,coi+1)] <- T
      }
      d[need_to_increment] <- d[need_to_increment] + line_thickness * 1.05
    }
    
    #find coordinates to plot tissue names
    if(sex == "male"){
      ylocs_scale <- 100
      xylocs_tissue_names <- FField::FFieldPtRep(coords = cbind(rep(4 + 0.25, nrow(d)),
                                                                c((d[,"8w"] + 1:length(d[,"8w"])/1000 - length(d[,"8w"])/2000) * ylocs_scale)),
                                                 rep.fact = 20, adj.max = 15)
      xylocs_tissue_names$y <- xylocs_tissue_names$y / ylocs_scale
      # xylocs_tissue_names <- as.data.frame(cbind(x = rep(4 + 0.25, nrow(d)), y = c(d[,"8w"])))
    }
    
    #find coordinates to plot gene names
    if(plot_gene_names){
      xylocs_tissues_genes <- data.frame(gene = unlist(dg[rownames(xylocs_tissue_names)[order(xylocs_tissue_names$y, decreasing = T)]]))
      xylocs_tissues_genes$tissue <- rownames(xylocs_tissue_names)[sapply(rownames(xylocs_tissues_genes), function(x) 
        grep(strsplit(x, "-")[[1]][1], rownames(xylocs_tissue_names)))]
      xylocs_tissues_genes$x <- 5.25
      xylocs_tissues_genes$y <- xylocs_tissue_names$y[match(xylocs_tissues_genes$tissue, rownames(xylocs_tissue_names))]
      xylocs_tissues_genes$y <- xylocs_tissues_genes$y + 1:length(xylocs_tissues_genes$y)/1E4
      xylocs_tissues_genes$y <- seq(1.15, -0.15, length.out = nrow(xylocs_tissues_genes))
      xylocs_tissues_genes_locs <- FField::FFieldPtRep(coords = cbind(xylocs_tissues_genes$x,
                                                                      xylocs_tissues_genes$y * 500),
                                                       rep.fact = 30, adj.max = 5, iter.max = 5E3)
      xylocs_tissues_genes$y <- xylocs_tissues_genes_locs$y / 500
    }
    
    
    plot(100,100,xlim = c(1,4), ylim = c(0,1), xpd=NA, 
         ylab = "Proportion Positive Effects on GWAS Trait", xlab = "Timepoint", xaxt = "n", 
         yaxt = "n", bty="n", cex.lab = 1.5, cex.axis = 1.25)
    
    #plot faded positive and negative regions
    rect(xl = 1, xr = 4, yb = 0.5, ytop = 1,
         col = grDevices::adjustcolor("red", 0.1), border = NA)
    rect(xl = 1, xr = 4, yb = 0, ytop = 0.5,
         col = grDevices::adjustcolor("blue", 0.1), border = NA)
    
    text(paste0(trait_categories$new_Phenotype[trait_categories$Tag == trait], ","), col = 1, cex = 2, font = 2, pos = 3,
         x = par("usr")[2] * 0.45 + par("usr")[1] * 0.55, y = par("usr")[3] * 0 + par("usr")[4] * 1)
    
    text(c("\u2642", "\u2640")[c("male", "female") == sex], col = cols$Sex[sex], cex = 3.5, font = 3, pos = 3,
         x = par("usr")[2] * 0.45 + par("usr")[1] * 0.55 + strwidth(paste0(trait_categories$new_Phenotype[trait_categories$Tag == trait], ",     "), units = "use")*2/2, 
         y = par("usr")[3] * 0 + par("usr")[4] * 1)
    
    
    #horizontal axis
    segments(x0 = 1:4, x1 = 1:4, y0 = - 0.02, y1 = - 0.04, lwd = 2)
    segments(x0 = 1, x1 = 4, y0 = - 0.02, y1 = - 0.02, lwd = 2)
    text(x = 1:4, y = - 0.07, labels = paste0(2^(0:3), "w"), pos = , cex = 1.251)
    
    #vertical axis
    segments(x0 = 1 - 3 * 0.02, x1 = 1 - 3 * 0.04, y0 = 0:5/5, y1 = 0:5/5, lwd = 2)
    segments(x0 = 1 - 3 * 0.02, x1 = 1 - 3 * 0.02, y0 = 0, y1 = 1, lwd = 2)
    text(x = 1 - 3 * 0.035, y = 0:5/5, labels = 0:5/5, pos = 2, cex = 1.25)
    
    
    for(tissue in rownames(d)){
      
      if(all(is.na(d[tissue,]))){
        next()
      }
      
      lines(1:4, d[tissue,], lwd = lwd, col = cols$Tissue[tissue])
      
      #plot gene names
      if(plot_gene_names){
        tissue_genes <- xylocs_tissues_genes[xylocs_tissues_genes$tissue == tissue,]
        
        #break up gene names to allow for different colors
        genes_df <- do.call(rbind, strsplit(tissue_genes$gene, "(?=[\\(,)])", perl = TRUE))
        tissue_genes$gene <- tissue_genes$gene[order(genes_df[,3], decreasing = T)]
        genes_df <- genes_df[order(genes_df[,3], decreasing = T),]
        if(nrow(tissue_genes) == 1){
          genes_df <- t(as.matrix(genes_df))
        }
        max_ef <- 5
        colors_mat <- cbind(cols$Tissue[tissue], 
                            1, 
                            ifelse3(genes_df[,3] == "+", "red", "blue"), 
                            1, 
                            rev(viridisLite::cividis(n = 100))[min2(ceiling(as.numeric(genes_df[,5]) / max_ef * 100), 100)], 
                            1)
        text_indivcolor(labels_mat = genes_df, xloc = tissue_genes$x, y = tissue_genes$y, 
                        colors_mat = colors_mat, pos = 4, cex = 0.75)
        
        # #all same color text
        # text(labels = tissue_genes$gene,
        #      x = tissue_genes$x,
        #      y = tissue_genes$y,
        #      cex = 0.75, col = cols$Tissue[tissue], pos = 4)
        
        #plot connecting lines
        for(gene_i in 1:nrow(tissue_genes)){
          segments(x0 = xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue] + 
                     strwidth(paste0(MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[tissue], " (", dn[tissue,"8w"], ")   "), cex = tissue_name_cex, units = "user"), 
                   y0 = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue],
                   x1 = tissue_genes$x[gene_i] + 0.75*strwidth("  ", cex = tissue_name_cex, units = "user"),
                   y1 = tissue_genes$y[gene_i],
                   col = cols$Tissue[tissue], lty = 3)
        }
      }
      
    }
    
    #plot tissue names
    for(tissue in rownames(d)){
      text(x = xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue], cex = tissue_name_cex,
           y = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue], 
           labels = paste0(MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[tissue], " (", dn[tissue,"8w"], ")"), col = cols$Tissue[tissue], pos = 4)
      segments(x0 = 4, y0 = d[tissue,"8w"], 
               x1 = xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue] + 0.075, 
               y1 = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue],
               col = cols$Tissue[tissue], lty = 3)
    }
    
    # legend(x = 1, y = 1.5, legend = tissue_names,
    #        col = cols$Tissue, lwd = 3, ncol = 4, cex = 1, border = NA, seg.len = 1, bg = NA, bty = "n", x.intersp = 0.25, text.width = 0.65)
    # segments(x0 = 0.9, y0 = 0.5, x1 = 4.1, y1 = 0.5, lty = 3, lwd = 3, col = "lightgrey")
    
  }
  
}

dev.off()

pdftools::pdf_combine(c("~/Documents/Documents - nikolai/pass1b_fig8_DE_protective_effect_in_UKB_20002_1111_self_reported_asthma.pdf",
                        "~/Documents/Documents - nikolai/pass1b_fig8_DE_protective_effect_in_UKB_20002_1473_self_reported_high_cholesterol.pdf",
                        "~/Documents/Documents - nikolai/pass1b_fig8_DE_protective_effect_in_GLGC_Mc_TG.pdf",
                        "~/Documents/Documents - nikolai/pass1b_fig8_DE_protective_effect_in_RA_OKADA_TRANS_ETHNIC.pdf"),
                      output = paste0("~/Documents/Documents - nikolai/motrpac_companion/figures/multi-trait-trajectories.pdf"))

#### trajectories, but sexes opposing each other ####

#plotting params
plot_gene_names <- T
add_in_alltiss <- T
trait_patts <- c("triglycerides", "rheumatoid_arthritis", "high_cholesterol", "asthma_ukbs")[c(3,4)]
n_deg_sigtwas_intersect
# trait_patts <- c("body", "standing")
# trait_patts <- c("triglycerides")


cairo_pdf(paste0("~/Documents/Documents - nikolai/motrpac_companion/figures/pass1b_fig8_DE_protective_effect_opposing.pdf"), 
          width = 1400 / 72, height = 475 / 72 * length(trait_patts), family="Arial Unicode MS", pointsize = 12.5)
par(mfrow = c(length(trait_patts),1), xpd = NA, mar = c(4,3,4,17.75))
lwd <- 3
tissue_name_cex <- 1.2
gene_cex <- 0.75
gene_location <- 5.5
f_offset <- 7.15

my_tissue_abbr <- MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV
my_tissue_abbr["all_tissues"] <- "ALL"
for(trait_patt in trait_patts){
  trait <- intersect(twas_with_hits, 
                     trait_categories$Tag[grep(trait_categories$new_Phenotype, pattern = trait_patt, ignore.case = T)])[1]
  
  for(sex in c("male", "female")){
    
    #retrieve data
    d <- apply(deg_sigtwas_proportion[, trait,,"p",sex], 2, as.numeric)
    dn <-  apply(deg_sigtwas_proportion[, trait,,"n",sex], 2, as.numeric)
    colnames(d) <- colnames(dn) <- colnames(deg_sigtwas_proportion[, trait,,"p",sex])
    rownames(d) <- rownames(dn) <- rownames(deg_sigtwas_proportion[, trait,,"p",sex])
    
    #add in summary of all tissues
    if(add_in_alltiss){
      at_prop <- apply(d * dn, 2, sum, na.rm = T)
      dn <- rbind(dn, all_tissues = apply(dn, 2, sum, na.rm = T))
      d <- rbind(d, all_tissues = at_prop / dn["all_tissues",])
      # d <- rbind(d, all_tissues = apply(d, 2, mean, na.rm = T))
    }
    
    dn <- dn[!apply(apply(d, 1, is.na), 2, all),]
    d <- d[!apply(apply(d, 1, is.na), 2, all),]
    d <- d[!is.nan(d[,"8w"]),]
    d[is.nan(d)] <- 0.5
    dg <- deg_sigtwas_proportion[, trait,,"genes",sex]
    dg <- lapply(dg[,"8w"], function(x) strsplit(x, " ~ ")[[1]])
    
    #munge gene names into composite strings
    dgm <- deg_sigtwas_proportion[, trait,,"genes","male"]
    dgm <- lapply(dgm[,"8w"], function(x) strsplit(x, " ~ ")[[1]])
    dgf <- deg_sigtwas_proportion[, trait,,"genes","female"]
    dgf <- lapply(dgf[,"8w"], function(x) strsplit(x, " ~ ")[[1]])
    dg <- lapply(setNames(names(dgm), names(dgm)), function(tiss){
      if(is.na(dgm[[tiss]][1])){return(NA)}
      mvs <- do.call(rbind, strsplit(dgm[[tiss]], " "))
      fvs <- do.call(rbind, strsplit(dgf[[tiss]], " "))
      if(nrow(mvs) > 1){
        apply(cbind(mvs[,2:3], mvs[,1], fvs[match(mvs[,1], fvs[,1]), 2:3]), 1, paste0, collapse = " ")  
      } else {
        mvs <- c(mvs)
        fvs <- c(fvs)
        paste0(c(mvs[2:3], mvs[1], fvs[2:3]), collapse = " ")
      }
    })
    
    
    #iterate through d to make identical lines parallel
    line_thickness <- lwd / 96 / par("pin")[2] * (par("usr")[4] - par("usr")[3])
    need_to_increment <- matrix(T, nrow = nrow(d), ncol = ncol(d), dimnames = list(rownames(d), colnames(d)))
    while(any(need_to_increment)){
      need_to_increment <- matrix(F, nrow = nrow(d), ncol = ncol(d), dimnames = list(rownames(d), colnames(d)))
      for(coi in 1:(ncol(d)-1)){
        unchanging_tissues <- rownames(unique(d[,c(coi,coi+1)]))
        need_to_increment[setdiff(rownames(d), unchanging_tissues),c(coi,coi+1)] <- T
      }
      d[need_to_increment] <- d[need_to_increment] + line_thickness * 1.05
    }
    
    #find coordinates to plot tissue names
    if(sex == "male"){
      # ylocs_scale <- 100
      xylocs_tissue_names <- FField::FFieldPtRep(coords = cbind(rep(4 + 0.25, nrow(d)),
                                                                c((d[,"8w"] + 1:length(d[,"8w"])/1000 - length(d[,"8w"])/2000) * ylocs_scale)),
                                                 rep.fact = 20, adj.max = 15)
      xylocs_tissue_names$y <- xylocs_tissue_names$y / ylocs_scale
      # xylocs_tissue_names <- as.data.frame(cbind(x = rep(4 + 0.25, nrow(d)), y = c(d[,"8w"])))
    }
    
    #find coordinates to plot gene names
    if(plot_gene_names){
      xylocs_tissues_genes <- data.frame(gene = unlist(dg[rownames(xylocs_tissue_names)[order(xylocs_tissue_names$y, decreasing = T)]]))
      xylocs_tissues_genes$tissue <- rownames(xylocs_tissue_names)[sapply(rownames(xylocs_tissues_genes), function(x) 
        grep(strsplit(x, "-")[[1]][1], rownames(xylocs_tissue_names)))]
      xylocs_tissues_genes$x <- gene_location
      xylocs_tissues_genes$y <- xylocs_tissue_names$y[match(xylocs_tissues_genes$tissue, rownames(xylocs_tissue_names))]
      xylocs_tissues_genes$y <- xylocs_tissues_genes$y + 1:length(xylocs_tissues_genes$y)/1E4
      xylocs_tissues_genes$y <- seq(1.15, -0.15, length.out = nrow(xylocs_tissues_genes))
      # xylocs_tissues_genes_locs <- FField::FFieldPtRep(coords = cbind(xylocs_tissues_genes$x,
      #                                                                 xylocs_tissues_genes$y * 500),
      #                                                  rep.fact = 30, adj.max = 5, iter.max = 5E3)
      # xylocs_tissues_genes$y <- xylocs_tissues_genes_locs$y / 500
    }
    
    if(sex == "male"){
      
      plot(100,100,xlim = c(1,9.5), ylim = c(0,1), xpd=NA, 
           ylab = "", xlab = "", xaxt = "n", 
           yaxt = "n", bty="n", cex.lab = 1.5, cex.axis = 1.25)
      
      
      text("Timepoint", x = 2.5 + c(0,f_offset), y = -0.125, pos = 1, cex = 1.5)
      if(trait_goodness[trait_categories$new_Phenotype[trait_categories$Tag == trait]] == -1){
        vlab <- "Proportion Risk Enhancing Effects"
      } else if(trait_goodness[trait_categories$new_Phenotype[trait_categories$Tag == trait]] == 0){
        vlab <- "Proportion Positive Effects on GWAS Trait"
      } else if(trait_goodness[trait_categories$new_Phenotype[trait_categories$Tag == trait]] == 1){
        vlab <- "Proportion Risk Reducing Effects"
      }
      
      text(vlab, x = 0.5, y = 0.5, pos = 3, cex = 1.5, srt = 90, xpd = NA)
      text(vlab, x = 4 + f_offset + 0.6, y = 0.5, pos = 1, cex = 1.5, srt = 270, xpd = NA)
      
      
      #plot faded positive and negative regions
      rect(xl = 1, xr = 4, yb = 0.5, ytop = 1,
           col = grDevices::adjustcolor("red", 0.1), border = NA)
      rect(xl = 1, xr = 4, yb = 0, ytop = 0.5,
           col = grDevices::adjustcolor("blue", 0.1), border = NA)
      
      #trait and sex name
      text(paste0(trait_categories$new_Phenotype[trait_categories$Tag == trait], ","), col = 1, cex = 2, font = 2, pos = 3,
           x = par("usr")[2] * 0.185 + par("usr")[1] * 0.825, y = par("usr")[3] * 0 + par("usr")[4] * 1)
      text(c("\u2642", "\u2640")[c("male", "female") == sex], col = cols$Sex[sex], cex = 3.5, font = 3, pos = 3,
           x = par("usr")[2] * 0.185 + par("usr")[1] * 0.825 + strwidth(paste0(trait_categories$new_Phenotype[trait_categories$Tag == trait], ",     "), units = "use")*2/2, 
           y = par("usr")[3] * 0 + par("usr")[4] * 1)
      
      
      #horizontal axis
      segments(x0 = 1:4, x1 = 1:4, y0 = - 0.02, y1 = - 0.04, lwd = 2)
      segments(x0 = 1, x1 = 4, y0 = - 0.02, y1 = - 0.02, lwd = 2)
      text(x = 1:4, y = - 0.07, labels = paste0(2^(0:3), "w"), pos = , cex = 1.251)
      
      #vertical axis
      segments(x0 = 1 - 3 * 0.02, x1 = 1 - 3 * 0.04, y0 = 0:5/5, y1 = 0:5/5, lwd = 2)
      segments(x0 = 1 - 3 * 0.02, x1 = 1 - 3 * 0.02, y0 = 0, y1 = 1, lwd = 2)
      text(x = 1 - 3 * 0.035, y = 0:5/5, labels = 0:5/5, pos = 2, cex = 1.25)
      
      
      for(tissue in rownames(d)){
        
        if(all(is.na(d[tissue,]))){
          next()
        }
        
        lines(1:4, d[tissue,], lwd = lwd, col = ifelse(is.na(cols$Tissue[tissue]), 1, cols$Tissue[tissue]))
        
        #plot gene names
        maxwidth_genename <- max(strwidth(xylocs_tissues_genes$gene, units = "user", cex = gene_cex))
        if(plot_gene_names & tissue != "all_tissues"){
          tissue_genes <- xylocs_tissues_genes[xylocs_tissues_genes$tissue == tissue,]
          
          #break up gene names to allow for different colors
          genes_df <- do.call(rbind, strsplit(tissue_genes$gene, "(?=[\\(,)])", perl = TRUE))
          tissue_genes$gene <- tissue_genes$gene[order(genes_df[,3], decreasing = T)]
          genes_df <- genes_df[order(genes_df[,3], decreasing = T),]
          if(nrow(tissue_genes) == 1){
            genes_df <- t(as.matrix(genes_df))
          }
          max_ef <- 5
          colors_mat <- cbind(1, 
                              ifelse3(genes_df[,2] == "+", "red", "blue"), 
                              1, 
                              rev(viridisLite::cividis(n = 100))[min2(ceiling(as.numeric(genes_df[,4]) / max_ef * 100), 100)], 
                              1,
                              cols$Tissue[tissue], 
                              1, 
                              ifelse3(genes_df[,8] == "+", "red", "blue"), 
                              1, 
                              rev(viridisLite::cividis(n = 100))[min2(ceiling(as.numeric(genes_df[,10]) / max_ef * 100), 100)], 
                              1)
          text_indivcolor(labels_mat = genes_df, xloc = tissue_genes$x + (maxwidth_genename - strwidth(tissue_genes$gene, units = "user", cex = gene_cex)) / 2, y = tissue_genes$y, 
                          colors_mat = colors_mat, pos = 4, cex = gene_cex)
          
          
          # #all same color text
          # text(labels = tissue_genes$gene,
          #      x = tissue_genes$x,
          #      y = tissue_genes$y,
          #      cex = gene_cex, col = cols$Tissue[tissue], pos = 4)
          
          #plot connecting lines
          for(gene_i in 1:nrow(tissue_genes)){
            segments(x0 = xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue] + 
                       strwidth(paste0(my_tissue_abbr[tissue], " (", dn[tissue,"8w"], ")   "), cex = tissue_name_cex, units = "user"), 
                     y0 = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue],
                     x1 = tissue_genes$x[gene_i] + gene_cex*strwidth("  ", cex = tissue_name_cex, units = "user"),
                     y1 = tissue_genes$y[gene_i],
                     col = cols$Tissue[tissue], lty = 3)
            segments(x0 = tissue_genes$x[gene_i] + gene_cex*strwidth("  ", cex = tissue_name_cex, units = "user"), 
                     y0 = tissue_genes$y[gene_i],
                     x1 = tissue_genes$x[gene_i] + gene_cex*strwidth("  ", cex = tissue_name_cex, units = "user") + 
                       (maxwidth_genename - strwidth(tissue_genes$gene[gene_i], units = "user", cex = gene_cex)) / 2,
                     y1 = tissue_genes$y[gene_i],
                     col = cols$Tissue[tissue], lty = 3)
          }
        }
        
      }
      
      #plot tissue names
      for(tissue in rownames(d)){
        text(x = xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue], cex = tissue_name_cex,
             y = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue], 
             labels = paste0(my_tissue_abbr[tissue], " (", dn[tissue,"8w"], ")"), col = cols$Tissue[tissue], pos = 4)
        segments(x0 = 4, y0 = d[tissue,"8w"], 
                 x1 = xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue] + 0.075, 
                 y1 = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue],
                 col = ifelse(is.na(cols$Tissue[tissue]), 1, cols$Tissue[tissue]), lty = 3)
      }
      
      # legend(x = 1, y = 1.5, legend = tissue_names,
      #        col = cols$Tissue, lwd = 3, ncol = 4, cex = 1, border = NA, seg.len = 1, bg = NA, bty = "n", x.intersp = 0.25, text.width = 0.65)
      # segments(x0 = 0.9, y0 = 0.5, x1 = 4.1, y1 = 0.5, lty = 3, lwd = 3, col = "lightgrey")
    } else if(sex == "female") {
      
      #find desired offset for female plot
      # f_offset <- (gene_location - 4) * 2 + maxwidth_genename
      
      #plot faded positive and negative regions
      rect(xl = 1+f_offset, xr = 4+f_offset, yb = 0.5, ytop = 1,
           col = grDevices::adjustcolor("red", 0.1), border = NA)
      rect(xl = 1+f_offset, xr = 4+f_offset, yb = 0, ytop = 0.5,
           col = grDevices::adjustcolor("blue", 0.1), border = NA)
      
      #trait and sex name
      text(paste0(trait_categories$new_Phenotype[trait_categories$Tag == trait], ","), col = 1, cex = 2, font = 2, pos = 3,
           x = par("usr")[2] * 0.185 + par("usr")[1] * 0.825+f_offset, y = par("usr")[3] * 0 + par("usr")[4] * 1)
      text(c("\u2642", "\u2640")[c("male", "female") == sex], col = cols$Sex[sex], cex = 3.5, font = 3, pos = 3,
           x = par("usr")[2] * 0.185 + par("usr")[1] * 0.825 + strwidth(paste0(trait_categories$new_Phenotype[trait_categories$Tag == trait], ",     "), units = "use")*2/2+f_offset, 
           y = par("usr")[3] * 0 + par("usr")[4] * 1)
      
      #horizontal axis
      segments(x0 = 1:4+f_offset, x1 = 1:4+f_offset, y0 = - 0.02, y1 = - 0.04, lwd = 2)
      segments(x0 = 1+f_offset, x1 = 4+f_offset, y0 = - 0.02, y1 = - 0.02, lwd = 2)
      text(x = 1:4+f_offset, y = - 0.07, labels = rev(paste0(2^(0:3), "w")), pos = , cex = 1.251)
      
      #vertical axis
      segments(x0 = 1 + 3 * 0.02+f_offset+3, x1 = 1 + 3 * 0.04+f_offset+3, y0 = 0:5/5, y1 = 0:5/5, lwd = 2)
      segments(x0 = 1 + 3 * 0.02+f_offset+3, x1 = 1 + 3 * 0.02+f_offset+3, y0 = 0, y1 = 1, lwd = 2)
      text(x = 1 + 3 * 0.035+f_offset+3, y = 0:5/5, labels = 0:5/5, pos = 4, cex = 1.25)
      
      
      for(tissue in rownames(d)){
        
        if(all(is.na(d[tissue,]))){
          next()
        }
        
        lines(rev(1:4)+f_offset, d[tissue,], lwd = lwd, col = ifelse(is.na(cols$Tissue[tissue]), 1, cols$Tissue[tissue]))
        
      }
      
      #plot tissue names
      for(tissue in rownames(d)){
        text(x = -(xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue] - 4) + f_offset + 1 - 0.025, cex = tissue_name_cex,
             y = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue], 
             labels = paste0("(", dn[tissue,"8w"], ") ", my_tissue_abbr[tissue]), col = cols$Tissue[tissue], pos = 2)
        segments(x0 = 1+f_offset, y0 = d[tissue,"8w"], 
                 x1 = -(xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue] - 4) + f_offset + 1 - 0.075, 
                 y1 = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue],
                 col = ifelse(is.na(cols$Tissue[tissue]), 1, cols$Tissue[tissue]), lty = 3)
      }
      
      #plot connecting lines
      
      for(tissue in rownames(d)){
        if(tissue == "all_tissues"){next()}
        tissue_genes <- xylocs_tissues_genes[xylocs_tissues_genes$tissue == tissue,]
        
        for(gene_i in 1:nrow(tissue_genes)){
          segments(x0 = -(xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue] - 4) + f_offset + 1 - 0.075 - 
                     strwidth(paste0(my_tissue_abbr[tissue], " (", dn[tissue,"8w"], ")  "), cex = tissue_name_cex, units = "user"), 
                   y0 = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue],
                   x1 = tissue_genes$x[gene_i] + gene_cex*strwidth("  ", cex = tissue_name_cex, units = "user") + 
                     (maxwidth_genename - strwidth(tissue_genes$gene[gene_i], units = "user", cex = gene_cex)) + 
                     strwidth(tissue_genes$gene[gene_i], units = "user", cex = gene_cex),
                   y1 = tissue_genes$y[gene_i],
                   col = cols$Tissue[tissue], lty = 3)
          
          segments(x0 = tissue_genes$x[gene_i] + gene_cex*strwidth("  ", cex = tissue_name_cex, units = "user") + 
                     (maxwidth_genename - strwidth(tissue_genes$gene[gene_i], units = "user", cex = gene_cex)) + 
                     strwidth(tissue_genes$gene[gene_i], units = "user", cex = gene_cex),
                   y0 = tissue_genes$y[gene_i],
                   x1 = tissue_genes$x[gene_i] + gene_cex*strwidth("  ", cex = tissue_name_cex, units = "user") +
                     (maxwidth_genename - strwidth(tissue_genes$gene[gene_i], units = "user", cex = gene_cex)) / 2 +
                     strwidth(tissue_genes$gene[gene_i], units = "user", cex = gene_cex),
                   y1 = tissue_genes$y[gene_i],
                   col = cols$Tissue[tissue], lty = 3)
          
        }
      }
      
    }
  }
  
}

dev.off()

#### plot raw scatterplot ####
category_colors <- RColorBrewer::brewer.pal(length(salient.categories), "Dark2")
names(category_colors) <- sort(salient.categories)

subset_to_traits <- T
trait_subset <- colnames(n_deg_sigtwas_intersect)[
  traitwise_partitions$Category[match(colnames(n_deg_sigtwas_intersect), traitwise_partitions$Tag)] %in% 
    c("Cardiometabolic", "Psychiatric-neurologic")]
categories_represented <- unique(traitwise_partitions$Category[match(trait_subset, traitwise_partitions$Tag)])
if(subset_to_traits){
  traits_to_plot <- trait_subset
} else {
  traits_to_plot <- twas_with_hits
}

cols = list(Tissue=MotrpacRatTraining6moData::TISSUE_COLORS[names(motrpac_gtex_map)], 
            Time=MotrpacRatTraining6moData::GROUP_COLORS[paste0(c(1,2,4,8), "w")],
            Sex=MotrpacRatTraining6moData::SEX_COLORS[c('male','female')])
cols$Tissue<- cols$Tissue[order(match(MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[names(cols$Tissue)], MotrpacRatTraining6moData::TISSUE_ORDER))]

tissue_names <- sapply(strsplit(names(cols$Tissue), "-"), function(x) paste0(x[ifelse(length(x) == 2, c(2), list(2:3))[[1]]], collapse = " "))

#plotting params
trait_category_legend_below <- T
deg_sigtwas_proportion[,order(),,,]

cairo_pdf(paste0("~/Documents/Documents - nikolai/pass1b_fig8_DE_protective_effect_scatterplot.pdf"), 
          width = 2000 / 72 * length(traits_to_plot) / 73, height = 700 / 72, family="Arial Unicode MS")
par(mfrow = c(1,1), xpd = NA, 
    mar = c(14,
            9,
            5,
            7 + ifelse(subset_to_traits, 2, 0)))
lwd <- 3
tissue_name_cex <- 1.2

plot(100,100, xlim = c(3,length(traits_to_plot)), ylim = c(0,1), xpd=NA, 
     ylab = "", xlab = "", xaxt = "n", 
     yaxt = "n", bty="n", cex.lab = 2, cex.axis = 1.25)
text(label = "Proportion Positive Effects on GWAS Trait", x = par("usr")[1] - diff(par("usr")[1:2])/12, y = mean(par("usr")[3:4]), srt = 90, pos = 3,
     cex = 1.75)

for(sex in c("male", "female")[1]){
  
  # mean_freq <- apply(sapply(traits_to_plot, function(trait) as.numeric(deg_sigtwas_proportion[, trait,"8w","p",sex])), 2, weighted.mean, na.rm = T)
  
  mean_freq <- sapply(setNames(1:length(traits_to_plot), traits_to_plot), function(trait_i){
    weighted.mean(as.numeric(deg_sigtwas_proportion[, traits_to_plot[trait_i],"8w","p",sex]), 
                  w = as.numeric(deg_sigtwas_proportion[, traits_to_plot[trait_i],"8w","n",sex]),
                  na.rm = T)})
  
  order_traits_to_plot <- order(mean_freq, decreasing = T)
  
  #plot faded positive and negative regions
  rect(xl = 1, xr = length(traits_to_plot) + 2, yb = 0.5, ytop = 1,
       col = grDevices::adjustcolor("red", 0.1), border = NA)
  rect(xl = 1, xr = length(traits_to_plot) + 2, yb = 0, ytop = 0.5,
       col = grDevices::adjustcolor("blue", 0.1), border = NA)
  
  #vertical axis
  segments(x0 = 1 - length(traits_to_plot) * 0.005, x1 = 1 - 3 * 0.04, y0 = 0:5/5, y1 = 0:5/5, lwd = 2)
  segments(x0 = 1 - 3 * 0.02, x1 = 1 - 3 * 0.02, y0 = 0, y1 = 1, lwd = 2)
  text(x = 1 - 3 * 0.05, y = 0:5/5, labels = 0:5/5, pos = 2, cex = 1.5)
  
  #horizontal axis
  text(1:length(traits_to_plot) + 1.5, -0.07, srt = 45, pos = 2,
       trait_categories$new_Phenotype[match(traits_to_plot[order_traits_to_plot], trait_categories$Tag)])
  segments(1:length(traits_to_plot) + 1, 0, 1:length(traits_to_plot) + 1, 1, col = adjustcolor(1, 0.2), lty = 2)
  
  #plot points
  points(1:length(traits_to_plot)+1, y = sort(mean_freq, decreasing = T), pch = "*", cex = 3)
  
  for(trait_i in (1:length(traits_to_plot))){
    
    #plot category blocks
    rect(xleft = which(order_traits_to_plot == trait_i) + 1/2, xright = which(order_traits_to_plot == trait_i) + 3/2,
         ybottom = -0.06, ytop = -0.02,
         col = category_colors[trait_categories$Category[match(traits_to_plot[trait_i], trait_categories$Tag)]])
    
    
    d <- as.numeric(deg_sigtwas_proportion[, traits_to_plot[trait_i],"8w","p",sex])
    dn <- as.numeric(deg_sigtwas_proportion[, traits_to_plot[trait_i],"8w","n",sex])
    names(d) <- names(dn) <- rownames(deg_sigtwas_proportion)
    dn <- dn[!is.na(d)]
    d <- d[!is.na(d)]
    #plot white points first
    points(x = rep(which(order_traits_to_plot == trait_i) + 1, length(d)), y = d, pch = 19, col = "white", cex = dn^0.25)
    #and then the actual points
    points(x = rep(which(order_traits_to_plot == trait_i) + 1, length(d)), y = d, pch = 19, col = adjustcolor(cols$Tissue[names(d)], 0.5), cex = dn^0.25)
    
  }
  
  
  #legend for tissues
  points(x = rep(length(traits_to_plot) + 3, length(cols$Tissue)), 
         y = seq(0.5, 0.99, length.out = length(cols$Tissue)), 
         col = adjustcolor(cols$Tissue, 0.5),
         pch = 19, cex = 2)
  text(x = rep(length(traits_to_plot) + 3.25, length(cols$Tissue)), 
       y = seq(0.5, 0.99, length.out = length(cols$Tissue)), 
       pos = 4, pch = 19, cex = 1,
       labels = MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[names(cols$Tissue)])
  
  points(x = length(traits_to_plot) + 3, 
         y = min(seq(0.5, 0.99, length.out = length(cols$Tissue))) - 0.0675, 
         col = "black", pch = "*", cex = 3)
  text(x = length(traits_to_plot) + 3.25, cex = 1.1,
       y = min(seq(0.5, 0.99, length.out = length(cols$Tissue))) - 0.075, 
       pos = 4, labels = "Weighted\nMean")
  
  #legend for tissue size
  n_pts <- 5
  pt_size_logs <- seq(1, log(max(as.numeric(deg_sigtwas_proportion[,,"8w","n",sex]), na.rm = T)) / log(2), length.out = n_pts)
  pt_size_legend <- round(2^pt_size_logs)
  text(x = length(traits_to_plot) + 2.25, 
       y = 0.1875 + 0.01 * n_pts + sum(pt_size_legend^0.25/100), 
       pos = 4, pch = 19, cex = 1.1,
       labels = "Sample Size")
  points(x = rep(length(traits_to_plot) + 3.25, n_pts), 
         y = 0.15 + cumsum(pt_size_legend^0.25/100) + cumsum(rep(0.01, n_pts)), 
         col = adjustcolor(1, 0.5),
         pch = 19, cex = pt_size_legend^0.25)
  text(x = rep(length(traits_to_plot) + 3.25, n_pts) + pt_size_legend^0.25/10, 
       y = 0.15 + cumsum(pt_size_legend^0.25/100) + cumsum(rep(0.01, n_pts)), 
       pos = 4, pch = 19, cex = 1,
       labels = pt_size_legend)
  
  #legend for categories
  if(trait_category_legend_below){
    yadj <- -0.15
    for(i in 1:length(categories_represented)){
      rect(xleft = 0 + i + cumsum(c(0,strwidth(categories_represented, cex = 1.25, units = "user") + 1))[i],
           xright = 1 + i + cumsum(c(0,strwidth(categories_represented, cex = 1.25, units = "user") + 1))[i],
           ybottom = -0.375 + yadj,
           ytop = -0.335 + yadj,
           col = category_colors[categories_represented[i]], border = 1)
      text(labels = categories_represented[i], pos = 4, 
           y = -0.355 + yadj, cex = 1.25, x = 0.85 + i + cumsum(c(0,strwidth(categories_represented, cex = 1.25, units = "user") + 1))[i])
      #draw circle?
      arc(t1 = 3*pi/2, t2 = pi/2, r1 = (0.355-0.04-yadj) / 2, r2 = (0.355-0.04-yadj) / 2, center = c(0.5,(-0.355 + yadj -0.04)/2), lwd = 3, res = 50, adjx = 20)
      points(0.5, -0.04, pch = -9658, cex = 3)
    }
  } else {
    x_adj2 <- x_adj - 2.5
    y_adj2 <- y_adj - 2.5
    for(i in 1:length(categories_represented)){
      rect(xleft = 0 + i,
           xright = 1 + i,
           ybottom = -10,
           ytop = -9,
           col = category_colors[categories_represented[i]], border = 1)
      text(labels = categories_represented[i], pos = 4, y = i + y_adj2 - 7, x = ncol(table_to_use) - 1/2 + x_adj2 + 1/3)
    }
  }
  
}

text("Proportion of Positive Effects Across Tissues and Traits at 8W Timepoint", x = mean(par("usr")[1:2]), y = 1.075, cex = 2.25, font = 2)

dev.off()

#### stitch pdfs together ####

pdftools::pdf_combine(paste0("~/Documents/Documents - nikolai/", 
                             c("pass1b_fig8_DEG-TWAS_Intersect_counts",
                               "pass1b_fig8_DEG-TWAS_Intersect_permille", 
                               "pass1b_fig8_DE_protective_effect_scatterplot"),".pdf"),
                      output = paste0("~/Documents/Documents - nikolai/", ifelse(use_tissue_cols_for_cols, "option_1", "option_2"), ".pdf"))



#### plot bayesian posterior scatterplot ####
if(length(salient.categories) < 9){
  category_colors <- RColorBrewer::brewer.pal(length(salient.categories), "Dark2")
  names(category_colors) <- sort(salient.categories)  
} else {
  category_colors <- RColorBrewer::brewer.pal(5, "Dark2")
  names(category_colors) <- sort(c("Cardiometabolic", "Aging", "Anthropometric", 
                                   "Immune", "Psychiatric-neurologic"))
  cats_leftover <- setdiff(salient.categories, names(category_colors))
  cols_leftover <- setdiff(RColorBrewer::brewer.pal(8, "Dark2"), category_colors)
  n_cols_needed <- length(cats_leftover) - length(cols_leftover)
  if(n_cols_needed > 0){
    extra_cols <- c(cols_leftover, RColorBrewer::brewer.pal(n_cols_needed, "Set2"))
  } else {
    extra_cols <- cols_leftover
  }
  names(extra_cols) <- cats_leftover
  category_colors <- c(category_colors, extra_cols)
}

#process significance files
# deg_sigtwas_proportion_posterior_mean
alpha_post = 0.1
cell_sigs <- deg_sigtwas_proportion_posterior_prob_bias < alpha_post | deg_sigtwas_proportion_posterior_prob_bias > (1-alpha_post)
# trait_means 
trait_sigs <- trait_bias_probs < alpha_post | trait_bias_probs > (1-alpha_post)

subset_to_traits <- F
if(subset_to_traits){
  trait_subset <- colnames(n_deg_sigtwas_intersect)[
    traitwise_partitions$Category[match(colnames(n_deg_sigtwas_intersect), traitwise_partitions$Tag)] %in% 
      ifelse2(subset_to_traits, c("Cardiometabolic", "Psychiatric-neurologic"), salient.categories)]
  
  #alternatively
  trait_subset <- names(trait_means)[trait_means > 0.5]
  traits_to_plot <- trait_subset
  
} else {
  traits_to_plot <- twas_with_hits
}
categories_represented <- unique(traitwise_partitions$Category[match(traits_to_plot, traitwise_partitions$Tag)])


cols = list(Tissue=MotrpacRatTraining6moData::TISSUE_COLORS[names(motrpac_gtex_map)], 
            Time=MotrpacRatTraining6moData::GROUP_COLORS[paste0(c(1,2,4,8), "w")],
            Sex=MotrpacRatTraining6moData::SEX_COLORS[c('male','female')])
cols$Tissue<- cols$Tissue[order(match(MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[names(cols$Tissue)], MotrpacRatTraining6moData::TISSUE_ORDER))]

tissue_names <- sapply(strsplit(names(cols$Tissue), "-"), function(x) paste0(x[ifelse(length(x) == 2, c(2), list(2:3))[[1]]], collapse = " "))

#other plotting params
trait_category_legend_below <- T

#do the plotting
cairo_pdf(paste0("~/Documents/Documents - nikolai/motrpac_companion/figures/prop_pos_hits_bayesian_scatterplot", ifelse(subset_to_traits, "_sub.pdf", "")), 
          width = 2000 / 72 * length(traits_to_plot) / 73, height = 700 / 72, family="Arial Unicode MS")
par(mfrow = c(1,1), xpd = NA, 
    mar = c(20,
            12,
            9,
            7 + ifelse(subset_to_traits, 2, 0)))
lwd <- 3
tissue_name_cex <- 1.2

plot(100,100, xlim = c(3,length(traits_to_plot)), ylim = range(deg_sigtwas_proportion_posterior_mean), xpd=NA, 
     ylab = "", xlab = "", xaxt = "n", 
     yaxt = "n", bty="n", cex.lab = 2, cex.axis = 1.25)
text(label = "Proportion Positive Effects on GWAS Trait", x = par("usr")[1] - diff(par("usr")[1:2])/12, y = mean(par("usr")[3:4]), srt = 90, pos = 3,
     cex = 1.75)

for(sex in c("male", "female")[1]){
  
  # mean_freq <- apply(sapply(traits_to_plot, function(trait) as.numeric(deg_sigtwas_proportion[, trait,"8w","p",sex])), 2, weighted.mean, na.rm = T)
  
  mean_freq <- trait_means[traits_to_plot]
  
  order_traits_to_plot <- order(mean_freq, decreasing = T)
  
  #plot faded positive and negative regions
  yseq <- round(seq(min(deg_sigtwas_proportion_posterior_mean) - 0.025, max(deg_sigtwas_proportion_posterior_mean) + 0.025, length.out = 6), 2)
  rect(xl = 1, xr = length(traits_to_plot) + 2, yb = 0.5, ytop = max(yseq),
       col = grDevices::adjustcolor("red", 0.1), border = NA)
  rect(xl = 1, xr = length(traits_to_plot) + 2, yb = min(yseq), ytop = 0.5,
       col = grDevices::adjustcolor("blue", 0.1), border = NA)
  
  #vertical axis
  segments(x0 = 1 - length(traits_to_plot) * 0.005, x1 = 1 - 3 * 0.04, y0 = yseq, y1 = yseq, lwd = 2)
  segments(x0 = 1 - 3 * 0.02, x1 = 1 - 3 * 0.02, y0 = min(yseq), y1 = max(yseq), lwd = 2)
  text(x = 1 - 3 * 0.05, y = yseq, labels = yseq, pos = 2, cex = 1.5)
  
  #horizontal axis
  text(1:length(traits_to_plot) + 1.5, min(yseq)-0.0275, srt = 45, pos = 2,
       labels = paste0(ifelse(trait_sigs[traits_to_plot][order_traits_to_plot], "♦ ", ""), 
                       trait_categories$new_Phenotype[match(traits_to_plot[order_traits_to_plot], trait_categories$Tag)]), 
       cex = 1.125, col = 1)
  segments(1:length(traits_to_plot) + 1, min(yseq)-0.01, 1:length(traits_to_plot) + 1, max(yseq), col = adjustcolor(1, 0.2), lty = 2)
  
  #plot points
  points(1:length(traits_to_plot)+1, y = sort(mean_freq[traits_to_plot], decreasing = T), pch = "*", cex = 3)
  
  for(trait_i in (1:length(traits_to_plot))){
    
    #plot category blocks
    rect(xleft = which(order_traits_to_plot == trait_i) + 1/2, xright = which(order_traits_to_plot == trait_i) + 3/2,
         ybottom = min(yseq)-0.02, ytop = min(yseq)-0.01,
         col = category_colors[trait_categories$Category[match(traits_to_plot[trait_i], trait_categories$Tag)]])
    
    
    d <- unlist(deg_sigtwas_proportion_posterior_mean[,traits_to_plot[trait_i]])
    d_sig <- unlist(cell_sigs[,traits_to_plot[trait_i]])
    tissue_abbr_rev <- MotrpacRatTraining6moData::TISSUE_ABBREV_TO_CODE
    names(d) <- tissue_abbr_rev[rownames(deg_sigtwas_proportion_posterior_mean)]
    
    dn <- as.numeric(deg_sigtwas_proportion[, traits_to_plot[trait_i],"8w","n",sex])
    names(dn) <- rownames(deg_sigtwas_proportion)
    dn <- dn[names(d)]
    
    #plot white points first
    points(x = rep(which(order_traits_to_plot == trait_i) + 1, length(d)), y = d, pch = 19 - d_sig, col = "white", cex = dn^0.25)
    #and then the actual points
    points(x = rep(which(order_traits_to_plot == trait_i) + 1, length(d)), y = d, pch = 19 - d_sig, 
           col = adjustcolor(cols$Tissue[names(d)], 0.5), cex = dn^0.25)
    #and then the significant points
    points(x = rep(which(order_traits_to_plot == trait_i) + 1, length(d)), y = d, pch = 5, 
           col = 1, cex = dn^0.25*d_sig*0.8)
    
  }
  
  
  #legend for tissues
  points(x = rep(length(traits_to_plot) + 3, length(cols$Tissue)), 
         y = seq(mean(yseq), max(yseq), length.out = length(cols$Tissue)), 
         col = adjustcolor(cols$Tissue, 0.5),
         pch = 19, cex = 2)
  text(x = rep(length(traits_to_plot) + 3.25, length(cols$Tissue)), 
       y = seq(mean(yseq), max(yseq), length.out = length(cols$Tissue)), 
       pos = 4, pch = 19, cex = 1,
       labels = MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[names(cols$Tissue)])
  
  points(x = length(traits_to_plot) + 3, 
         y = min(seq(mean(yseq), max(yseq), length.out = length(cols$Tissue))) - 0.025, 
         col = "black", pch = "*", cex = 3)
  text(x = length(traits_to_plot) + 3.25, cex = 1.1,
       y = min(seq(mean(yseq), max(yseq), length.out = length(cols$Tissue))) - 0.0275, 
       pos = 4, labels = "Posterior\nMean")
  
  #legend for tissue size
  n_pts <- 5
  pt_size_logs <- seq(1, log(max(as.numeric(deg_sigtwas_proportion[,,"8w","n",sex]), na.rm = T)) / log(2), length.out = n_pts)
  pt_size_legend <- round(2^pt_size_logs)
  per_pt_incr <- 0.01
  text(x = length(traits_to_plot) + 2.25, 
       y = 0.2875 + per_pt_incr * n_pts + sum(pt_size_legend^0.25/100), 
       pos = 4, pch = 19, cex = 1.1,
       labels = "Sample Size")
  points(x = rep(length(traits_to_plot) + 3.25, n_pts), 
         y = 0.362 + cumsum(pt_size_legend^0.25/800) + cumsum(rep(per_pt_incr, n_pts)), 
         col = adjustcolor(1, 0.5),
         pch = 19, cex = pt_size_legend^0.25)
  text(x = rep(length(traits_to_plot) + 3.25, n_pts) + pt_size_legend^0.25/10, 
       y = 0.362 + cumsum(pt_size_legend^0.25/800) + cumsum(rep(per_pt_incr, n_pts)), 
       pos = 4, pch = 19, cex = 1,
       labels = pt_size_legend)
  
  #legend for diamonds
  text(x = length(traits_to_plot) + 2.25, cex = 1,
       y = min(yseq) + 0.01, 
       pos = 4, labels = paste0("♦ indicates\nPr(enr. or dep.)\n> ", 1-alpha_post))
  
  
  #legend for categories
  if(trait_category_legend_below){
    yadj <- 0.4875
    for(i in 1:length(categories_represented)){
      rect(xleft = 0 + i + cumsum(c(0,strwidth(categories_represented, cex = 1.25, units = "user") + 1))[i],
           xright = 1 + i + cumsum(c(0,strwidth(categories_represented, cex = 1.25, units = "user") + 1))[i],
           ybottom = -0.345 + yadj,
           ytop = -0.335 + yadj,
           col = category_colors[categories_represented[i]], border = 1)
      text(labels = categories_represented[i], pos = 4, 
           y = -0.34 + yadj, cex = 1.25, x = 0.85 + i + cumsum(c(0,strwidth(categories_represented, cex = 1.25, units = "user") + 1))[i])
      #draw circle?
      arc(t1 = 3*pi/2, t2 = pi/2-1E-6, r1 = (0.34-0.015-yadj) / 2, r2 = (0.34-0.025-yadj) / 2, 
          center = c(0.5, (-0.34 + yadj + min(yseq) - 0.01)/2), lwd = 3, res = 50, adjx = 40)
      points(0.5, min(yseq)-0.015, pch = -9658, cex = 3)
    }
  } else {
    x_adj2 <- x_adj - 2.5
    y_adj2 <- y_adj - 2.5
    for(i in 1:length(categories_represented)){
      rect(xleft = 0 + i,
           xright = 1 + i,
           ybottom = -10,
           ytop = -9,
           col = category_colors[categories_represented[i]], border = 1)
      text(labels = categories_represented[i], pos = 4, y = i + y_adj2 - 7, x = ncol(table_to_use) - 1/2 + x_adj2 + 1/3)
    }
  }
  
}

text("Proportion of Positive Effects Across Tissues and Traits at 8W Timepoint", x = mean(par("usr")[1:2]), y = max(yseq) + 0.02, cex = 2.25, font = 2)

dev.off()

#### identify genes in intersect & x-reference against relative effect size ####
twas_intersect <- lapply(setNames(rownames(n_deg_sigtwas_intersect), rownames(n_deg_sigtwas_intersect)), function(x) NULL)
twas_with_hits <- colnames(prop_degs_are_twas)
prop_pos <- lapply(setNames(rownames(n_deg_sigtwas_intersect), rownames(n_deg_sigtwas_intersect)), function(i) NULL)
# tissue_code <- MotrpacBicQC::bic_animal_tissue_code
# tissue_code <- tissue_code[,c("tissue_name_release", "abbreviation")]
# tissue_code <- tissue_code[tissue_code$tissue_name_release != "" & !is.na(tissue_code$abbreviation)]
tissue_code <- data.frame(tissue_name_release = names(MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV),
                          abbreviation = MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV)
for(tissue_abbrev in names(twas_intersect)){
  
  tissue_i <- tissue_code$tissue_name_release[tissue_code$abbreviation == tissue_abbrev]
  
  print(tissue_i)
  
  sig_twas_by_trait <- lapply(setNames(twas_with_hits, twas_with_hits), 
                              function(trait_i) sig_twas_by_tissue[[tissue_i]][sig_twas_by_tissue[[tissue_i]]$trait == trait_i,])
  sig_twas_by_trait_genes <- lapply(setNames(twas_with_hits, twas_with_hits), 
                                    function(trait_i) sig_twas_by_trait[[trait_i]]$gene_name)
  sig_twas_by_trait_signs <- lapply(setNames(twas_with_hits, twas_with_hits), 
                                    function(trait_i){
                                      signs <- sign(sig_twas_by_trait[[trait_i]]$zscore)
                                      names(signs) <- sig_twas_by_trait_genes[[trait_i]]
                                      signs
                                    })
  
  timepoint = "8w"
  sex = "male"
  
  DE_genes_in_Nodes <- node_metadata_list[[timepoint]]$human_gene_symbol[
    node_metadata_list[[timepoint]]$tissue == tissue_code$abbreviation[tissue_code$tissue_name_release == tissue_i]]
  
  if(length(DE_genes_in_Nodes) == 0){next()}
  
  DE_genes_in_Nodes_sign <- cbind(as.data.frame(do.call(rbind, strsplit(node_metadata_list[[timepoint]]$node[
    node_metadata_list[[timepoint]]$tissue == tissue_code$abbreviation[tissue_code$tissue_name_release == tissue_i]], split = "_"))),
    gene = DE_genes_in_Nodes)
  
  colnames(DE_genes_in_Nodes_sign) <- c("time", "female", "male", "gene")
  DE_genes_in_Nodes_sign$female[grep(pattern = "F1", DE_genes_in_Nodes_sign$female)] <- 1
  DE_genes_in_Nodes_sign$female[grep(pattern = "F-1", DE_genes_in_Nodes_sign$female)] <- -1
  DE_genes_in_Nodes_sign$male[grep(pattern = "M1", DE_genes_in_Nodes_sign$male)] <- 1
  DE_genes_in_Nodes_sign$male[grep(pattern = "M-1", DE_genes_in_Nodes_sign$male)] <- -1
  DE_genes_in_Nodes_sign <- DE_genes_in_Nodes_sign[!is.na(DE_genes_in_Nodes_sign$gene),]
  
  DE_genes_signs <- as.integer(DE_genes_in_Nodes_sign[,sex])
  names(DE_genes_signs) <- DE_genes_in_Nodes_sign$gene
  
  results <- lapply(setNames(twas_with_hits, twas_with_hits), function(trait_i) {
    twas_signs <- sig_twas_by_trait_signs[[trait_i]]
    shared_genes <- intersect(names(DE_genes_signs), names(twas_signs))
    if(length(shared_genes) != 0){
      shared_genes_signs <- twas_signs[shared_genes] * DE_genes_signs[shared_genes]
      shared_genes_rat <- gene_map$RAT_SYMBOL[match(shared_genes, gene_map$HUMAN_ORTHOLOG_SYMBOL)]
      shared_genes_rat_ensembl_id <- gene_map$RAT_ENSEMBL_ID[match(shared_genes, gene_map$HUMAN_ORTHOLOG_SYMBOL)]
      
      # enrichment_test <- signif_df[which(signif_df$tissue == tissue_abbrev & signif_df$trait == trait_i),]
      
      # output <- data.frame(tissue = tissue_abbrev,
      #                      trait = trait_i,
      #                      human_ortholog_symbol = shared_genes,
      #                      rat_symbol = shared_genes_rat,
      #                      rat_ensembl_ID = shared_genes_rat_ensembl_id,
      #                      direction_of_DE_on_trait = shared_genes_signs,
      #                      Pr_Geneset_DiffProp_is_Pos = enrichment_test$prob_diff_is_positive,
      #                      geneset_enrichment = enrichment_test$signif)
      
      output <- data.frame(tissue = tissue_abbrev,
                           trait = trait_i,
                           human_ortholog_symbol = shared_genes,
                           rat_symbol = shared_genes_rat,
                           rat_ensembl_ID = shared_genes_rat_ensembl_id,
                           direction_of_DE_on_trait = shared_genes_signs)
      return(output)
    } else {return(NULL)}
  })
  
  results <- do.call(rbind, results)
  
  twas_intersect[[tissue_abbrev]] <- results
  
}



#### quick plot of exp vs obs freqs, incl prop pos hits ####
use_focal_vs_compl <- T

traits <- unique(data1$trait)
tissues <- unique(data1$tissue)
d <- list(cell_count = data1$count,
          total = sapply(1:nrow(data1), function(i) total_number_of_possible_hits_matrix[data1$tissue[i], data1$trait[i]]),
          row_count = sapply(1:nrow(data1), function(i) n_genes_in_nodes_matrix[data1$tissue[i], data1$trait[i]]),
          col_count = sapply(1:nrow(data1), function(i) sig_twas_by_trait_genes_matrix[data1$tissue[i], data1$trait[i]]),
          row_index = match(data1$tissue, tissues),
          col_index = match(data1$trait, traits),
          row_n = length(unique(data1$tissue)),
          col_n = length(unique(data1$trait)))




cairo_pdf(paste0("~/Documents/Documents - nikolai/motrpac_companion/figures/pass1b_fig8_DEG-TWAS_logodds_double.pdf"), 
          width = 1000 / 72 * 2, height = 925 / 72, family="Arial Unicode MS", pointsize = 36)

category_shapes <- setNames(15:19, salient.categories)
category_shapes <- setNames(15:19, c("Cardiometabolic", "Aging", "Anthropometric", 
                                     "Immune", "Psychiatric-neurologic")
)

par(mar = c(5,4,3,4), mfrow = c(1,2))
layout(t(c(rep(1,10), rep(2,10), 3)))
layout(t(c(1,2,3)), widths = c(1,1,0.1))
#first plot


sapply(1:nrow(signif_df), function(i) adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index[i]]], as.numeric(signif_df$signif[i] != 0) * 0.625 + 0.125))
if(use_focal_vs_compl){
  xvals <- (d$col_count - d$cell_count) / (d$total - d$row_count)
  yvals <- d$cell_count / d$row_count
} else {
  xvals <- d$row_count / d$total * d$col_count / d$total  
  yvals <- d$cell_count / d$total
}




plot(logit(xvals), logit(yvals), xlim = ifelse2(use_focal_vs_compl, c(-14,-1), c(-15,-4.25)), 
     ylim = ifelse2(use_focal_vs_compl, c(-6, -1), c(-9.5,-4.25)),
     pch = category_shapes[traitwise_partitions$Category[match(traits[d$col_index], traitwise_partitions$Tag)]], 
     col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5), 
     cex = 1.5,
     xlab = "", ylab = "", yaxt = 'n', xaxt = "n")
box(which = "plot", lwd = 2)
#axes
text(labels = "observed logit(frequency)", x = par("usr")[2] + diff(par("usr")[1:2])/8, srt = 270, pos = 1, 
     y = mean(par("usr")[3:4]), xpd = NA, cex = 1.25)
text(labels = "expected logit(frequency)", x = mean(par("usr")[1:2]), srt = 0, pos = 1, 
     y = par("usr")[3] - diff(par("usr")[3:4])/12, xpd = NA, cex = 1.25)
yaxlocs <- ifelse2(use_focal_vs_compl, c(-6:-1), c(-5:-9))
segments(x0 = par("usr")[2], x1 = par("usr")[2] + diff(par("usr")[1:2])/100, y0 = yaxlocs, y1 = yaxlocs, lwd = 4, xpd = NA)
text(labels = yaxlocs, x = par("usr")[2] + diff(par("usr")[1:2])/200, srt = 0, pos = 4, 
     y = yaxlocs, xpd = NA, cex = 1.25)
xaxlocs <- ifelse2(use_focal_vs_compl, -7:-1*2, -7:-3*2)
segments(x0 = xaxlocs, x1 = xaxlocs, y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4])/80, lwd = 4, xpd = NA)
text(labels = xaxlocs, x = xaxlocs, srt = 0, pos = 1, 
     y = par("usr")[3] - diff(par("usr")[3:4])/80, xpd = NA, cex = 1.25)
abline(0,1, lty = 2, lwd = 6, col = adjustcolor(1,0.75))

#secondplot
xl <- par("usr")[1]
xr <- par("usr")[1] + diff(par("usr")[1:2]) / 2
yt <- par("usr")[4]
yb <- par("usr")[4] - diff(par("usr")[3:4]) / 2

buffer <- 0.075
bx <- (xr-xl) * buffer
by <- (yt-yb) * buffer
rect(xleft = xl, ybottom = yb, xright = xr, ytop = yt, lwd = 2, xpd = NA)

points(x = ((xvals - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
       y = ((yvals - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, 
       pch = category_shapes[traitwise_partitions$Category[match(traits[d$col_index], traitwise_partitions$Tag)]], 
       col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5), 
       cex = 1)
segments(x0 = 0 * (xr - xl - 2*bx) + xl, 
         y0 = 0 * (yt - yb - 2*by) + yb,
         x1 = ((max(yvals, na.rm = T) - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
         y1 = ((max(yvals, na.rm = T) - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, 
         lty = 2, lwd = 4, col = adjustcolor(1,0.75))

#axes
text(labels = "observed frequency", x = par("usr")[1] - diff(par("usr")[1:2])/ifelse2(use_focal_vs_compl, 8, 7), srt = 90, pos = 1, 
     y = par("usr")[4] - diff(par("usr")[3:4])/5, xpd = NA, cex = 1)
text(labels = "expected frequency", x = par("usr")[1] + diff(par("usr")[1:2])/4, srt = 0, pos = 1, 
     y = par("usr")[4] + diff(par("usr")[3:4])/8, xpd = NA, cex = 1)
xaxlocs <- ifelse2(use_focal_vs_compl, c(0:4/20), c(0:6/500))
yaxlocs <- ifelse2(use_focal_vs_compl, c(0:6/20), c(0:6/500))

segments(x0 = par("usr")[1], x1 = par("usr")[1] - diff(par("usr")[1:2])/100, 
         y0 = ((yaxlocs - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, 
         y1 = ((yaxlocs - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, 
         lwd = 4, xpd = NA)
segments(y0 = par("usr")[4], y1 = par("usr")[4] + diff(par("usr")[3:4])/100, 
         x0 = ((xaxlocs - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
         x1 = ((xaxlocs - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
         lwd = 4, xpd = NA)
text(labels = xaxlocs, x = ((xaxlocs - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
     srt = 0, pos = 3, 
     y = par("usr")[4], xpd = NA, cex = 0.75)
text(labels = yaxlocs, x = par("usr")[1], 
     srt = 0, pos = 2, 
     y = ((yaxlocs - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, xpd = NA, cex = 0.75)
# text(labels = 0:6/500, x = -7:-3*2, srt = 0, pos = 1, 
#      y = par("usr")[3] - diff(par("usr")[3:4])/80, xpd = NA, cex = 1.25)
abline(0,1, lty = 2, lwd = 6, col = adjustcolor(1,0.75))

#legends
legend(x = xl, y = yb, legend = names(category_shapes), pch = category_shapes, col = adjustcolor(1, 0.5), cex = 0.75, pt.cex = 1.4)
legend(x = xr, y = yt, legend = tissues, pch = 19, col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues], 0.5), cex = 0.55, pt.cex = 1.2)
segments(x0 = par("usr")[1] + diff(par("usr")[1:2]) / 1.55, 
         y0 = par("usr")[4] - diff(par("usr")[3:4])/50,
         x1 = par("usr")[1] + diff(par("usr")[1:2]) / 1.45, 
         y1 = par("usr")[4] - diff(par("usr")[3:4])/50, 
         lty = 2, lwd = 4, col = adjustcolor(1,0.75))
text(x = par("usr")[1] + diff(par("usr")[1:2]) / 1.45, 
     y = par("usr")[4] - diff(par("usr")[3:4])/45, 
     pos = 4, cex = 0.75, labels = "1-to-1 line")
fig_label("a)", shrinkX = 0.975, cex = 2, shrinkY = 2)

#make second plot with raw proportion data
xvals <- data$count_all / data$total_all
yvals <- data$count_inters / data$total_inters


cex_pow <- 1/2
plot(xvals, yvals, xlim = c(0.48, 0.52), 
     ylim = c(0,1),
     pch = category_shapes[traitwise_partitions$Category[match(data$trait, traitwise_partitions$Tag)]], 
     col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[data$tissue], 0.5), 
     cex = (data$total_inters / max(data$total_inters))^cex_pow * 5,
     xlab = "", ylab = "", yaxt = 'n', xaxt = "n")
box(which = "plot", lwd = 2)
#axes
text(labels = "observed frequency (+ effects)", x = par("usr")[1] - diff(par("usr")[1:2])/7, srt = 90, pos = 1, 
     y = mean(par("usr")[3:4]), xpd = NA, cex = 1.25)
text(labels = "expected frequency (+ effects)", x = mean(par("usr")[1:2]), srt = 0, pos = 1, 
     y = par("usr")[3] - diff(par("usr")[3:4])/12, xpd = NA, cex = 1.25)
yaxlocs <- seq(0,1,by=0.1)
segments(x0 = par("usr")[1], x1 = par("usr")[1] - diff(par("usr")[1:2])/100, y0 = yaxlocs, y1 = yaxlocs, lwd = 4, xpd = NA)
text(labels = yaxlocs, x = par("usr")[1] - diff(par("usr")[1:2])/200, srt = 0, pos = 2, 
     y = yaxlocs, xpd = NA, cex = 1.25)
xaxlocs <- round(seq(par("usr")[1], par("usr")[2], by = 0.01), 2)
segments(x0 = xaxlocs, x1 = xaxlocs, y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4])/80, lwd = 4, xpd = NA)
text(labels = xaxlocs, x = xaxlocs, srt = 0, pos = 1, 
     y = par("usr")[3] - diff(par("usr")[3:4])/80, xpd = NA, cex = 1.25)
# abline(0,1, lty = 2, lwd = 6, col = adjustcolor(1,0.75))
abline(h=0.5,lwd=4,lty=2,col=adjustcolor(1,0.5))
abline(v=0.5,lwd=4,lty=2,col=adjustcolor(1,0.5))

#legened
xl <- par("usr")[2]
xr <- par("usr")[2] + diff(par("usr")[1:2]) / 2
yt <- par("usr")[4]
yb <- par("usr")[4] - diff(par("usr")[3:4]) / 2.75
legend(x = xl, y = yb, legend = names(category_shapes), bty = "n",
       pch = category_shapes, col = adjustcolor(1, 0.5), cex = 0.75, pt.cex = 1.4, xpd = NA)
legend(x = xl + (xr-xl)/100, y = yt, legend = tissues, pch = 19, bty = "n", xpd = NA,
       col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues], 0.5), cex = 0.55, pt.cex = 1.2)


legcexs <- c(0.5, 1:5)
legvals <- (round((legcexs / 5) ^ (1/cex_pow) * max(data$total_inters)))
legcexs <- (legvals / max(data$total_inters))^cex_pow * 5
# points(x = rep(xl + (xr-xl)/11, 6) + (rev(legcexs)) / 4000, y = yb - (yt-yb) * 1.4 + 
#          cumsum(rep((yt-yb)/60, length(legcexs))) + cumsum(legcexs + c(0,legcexs[-length(legcexs)])) / 120, 
#        cex = legcexs, xpd = NA, pch = 19, col = adjustcolor(1, 0.5))
text(labels = latex2exp::TeX("$n_{intersect}$"), x = xl - (xr-xl)/50, y = yb - (yt-yb) * 0.6, pos = 4, xpd = NA)
points(x = rep(xl + (xr-xl)/11, 6), y = yb - (yt-yb) * 1.4 + 
         cumsum(rep((yt-yb)/60, length(legcexs))) + cumsum(legcexs + c(0,legcexs[-length(legcexs)])) / 120, 
       cex = legcexs, xpd = NA, pch = 19, col = adjustcolor(1, 0.5))
text(x = rep(xl + (xr-xl)/11, 6), y = yb - (yt-yb) * 1.4 + 
       cumsum(rep((yt-yb)/60, length(legcexs))) + cumsum(legcexs + c(0,legcexs[-length(legcexs)])) / 120, 
     cex = legcexs / 5, xpd = NA, labels = legvals, col = "white")
text(labels = legvals[1:2], x = xl + (xr-xl)/15 + cumsum(legcexs)[1:2]/4000, y = yb - (yt-yb) * 1.4 + 
       cumsum(rep((yt-yb)/60, length(legcexs)))[1:2] + cumsum(legcexs + c(0,legcexs[-length(legcexs)]))[1:2] / 120, 
     pos = 4, xpd = NA, cex = 0.5)


fig_label("b)", shrinkX = 1, cex = 2)
dev.off()

#### composite figure for the heatmap, posterior output, and scatterplot ####

#load appropriate model fit
if(base != "deviation_from_expected_logodds_split_the_difference"){
  
  base = "deviation_from_expected_logodds_split_the_difference"
  load(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/mcmc_output/", paste0(base, ifelse(use_all_cats, "_allCats", "")),".cmdStanR.fit"))
  
  samps <- data.frame(as_draws_df(out$draws()))
  prop_greater_than_0 <- function(x) mean(x>0)
  cellbias <- apply(subset_samps("cell_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)
  celltotalbias <- apply(subset_samps("cell_total_prob_bias", c("raw", "sd", "logodds"), samps = samps), 2, prop_greater_than_0)
  rowbias <- apply(subset_samps("row_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)
  colbias <- apply(subset_samps("col_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)
  trait_key <- setNames(trait_categories$new_Phenotype, trait_categories$Tag)
  colcatbias <- apply(subset_samps("colcat_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)
  
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
}


cols = list(Tissue=MotrpacRatTraining6moData::TISSUE_COLORS[names(motrpac_gtex_map)], 
            Time=MotrpacRatTraining6moData::GROUP_COLORS[paste0(c(1,2,4,8), "w")],
            Sex=MotrpacRatTraining6moData::SEX_COLORS[c('male','female')])
mars <- list(c(6,
               0 + ifelse(subset_to_traits, 2, 0),
               6 + ifelse(group_by_tissue_type, disp_amount * 4.5, 0),
               6.5 + ifelse(subset_to_traits, 1, 0)),
             c(3.5,7,3,7)+0.5, 
             c(3.5,1,3,14)+0.5,
             c(4.25,2.5,2.5,5.5), 
             c(4.25,2.5,2.5,5.5), 
             c(4.25,7.5,2.5,0.5))

cairo_pdf(paste0("~/Documents/Documents - nikolai/motrpac_companion/figures/fig4_intersect-enrichment.pdf"), 
          width = 2300 / 72 * ncol(table_to_use) / 80, 
          height = 900 / 72 + ifelse(group_by_tissue_type, disp_amount * 0.75, 0), 
          family="Arial Unicode MS", pointsize = 18.5)
par(xpd = T, 
    mar = mars[[1]])

layout(rbind(
  c(1,1,1,1,4,6),
  c(2,2,3,3,5,6)
), heights = c(1,0.9))

#counts or props?
incl_significance <- T
incl_cell_totals <- F
trait_category_legend_below = T
use_tissue_cols_for_cols <- F
opacity_power_scaler <- 0.25
opacity_white_threshold <- 0.85
use_counts <- T
prop_TWAS <- T
order_by_counts <- F
order_by_posterior_enrichment <- T
group_by_tissue_type <- T
use_range_for_maginal_labels <- F

subset_to_traits <- T
focal_traitcats <- c("Cardiometabolic", "Immune", "Endocrine system", "Anthropometric", "Allergy", "Aging")
if(subset_to_traits){
  trait_subset <- colnames(n_deg_sigtwas_intersect)[
    traitwise_partitions$Category[match(colnames(n_deg_sigtwas_intersect), traitwise_partitions$Tag)] %in% 
      focal_traitcats]
  categories_represented <- unique(traitwise_partitions$Category[match(trait_subset, traitwise_partitions$Tag)])
} else {
  trait_subset <- colnames(n_deg_sigtwas_intersect)
}
categories_represented <- unique(traitwise_partitions$Category[match(trait_subset, traitwise_partitions$Tag)])

if(use_counts){
  table_to_use <- n_deg_sigtwas_intersect  
} else {
  if(prop_TWAS){
    table_to_use <- round(prop_twas_are_degs * 1000)  
  } else{
    table_to_use <- round(prop_degs_are_twas * 1000)  
  }
}

if(subset_to_traits){
  table_to_use <- table_to_use[,trait_subset]
  signif_matrix_to_use <- signif_matrix_to_use[,trait_subset]
}


if(order_by_counts){
  table_to_use <- table_to_use[,colnames(n_deg_sigtwas_intersect)]
  signif_matrix_to_use <- signif_matrix[,colnames(n_deg_sigtwas_intersect)]
} else if(order_by_posterior_enrichment) {
  trait_order <- traits[order(colbias, decreasing = T)]
  trait_order <- trait_order[trait_order %in% colnames(table_to_use)]
  table_to_use <- table_to_use[,trait_order]
  signif_matrix_to_use <- signif_matrix[,match(trait_order, colnames(signif_matrix))]
} else {
  table_to_use <- table_to_use[,colnames(prop_twas_are_degs)]
  signif_matrix_to_use <- signif_matrix[,colnames(prop_twas_are_degs)]
}

if(group_by_tissue_type){
  tissue_cats <- list(circulation = c("BLOOD", "HEART", "SPLEEN"),
                      skeletal_muscle = c("SKM-GN", "SKM-VL"),
                      other = rev(c("ADRNL", "KIDNEY", "LUNG", "LIVER")),
                      adipose = c("WAT−SC"),
                      brain = c("CORTEX", "HYPOTH", "HIPPOC"),
                      GI = c("SMLINT", "COLON"))
  tissue_cats <- rev(tissue_cats)
  disp_amount <- 0.5
  tissue_disps <- unlist(lapply(1:length(tissue_cats), function(tci) rep(disp_amount * (tci), length(tissue_cats[[tci]]))))
  tissue_cats_bars_ylocs <- cbind(start = (c(0, cumsum(unlist(lapply(tissue_cats, function(tc) length(tc))) + disp_amount)) + 2 * disp_amount)[-(length(tissue_cats)+1)], 
                                  end = cumsum(unlist(lapply(tissue_cats, function(tc) length(tc))) + disp_amount) + disp_amount)
  
}

if(use_range_for_maginal_labels){
  n_genes_in_nodes_label <- apply(apply(n_genes_in_nodes_matrix, 1, range), 2, paste0, collapse = " - ")
  sig_twas_by_trait_genes_label <- apply(apply(sig_twas_by_trait_genes_matrix, 2, range), 2, paste0, collapse = " - ")
} else {
  n_genes_in_nodes_label <- apply(n_genes_in_nodes_matrix, 1, max)
  sig_twas_by_trait_genes_label <- apply(sig_twas_by_trait_genes_matrix, 2, max)
}



if(length(salient.categories) < 9){
  category_colors <- RColorBrewer::brewer.pal(length(salient.categories), "Dark2")
  names(category_colors) <- sort(salient.categories)  
} else {
  category_colors <- RColorBrewer::brewer.pal(5, "Dark2")
  names(category_colors) <- sort(c("Cardiometabolic", "Aging", "Anthropometric", 
                                   "Immune", "Psychiatric-neurologic"))
  cats_leftover <- setdiff(salient.categories, names(category_colors))
  cols_leftover <- setdiff(RColorBrewer::brewer.pal(8, "Dark2"), category_colors)
  n_cols_needed <- length(cats_leftover) - length(cols_leftover)
  if(n_cols_needed > 0){
    extra_cols <- c(cols_leftover, RColorBrewer::brewer.pal(n_cols_needed, "Set2"))
  } else {
    extra_cols <- cols_leftover
  }
  names(extra_cols) <- cats_leftover
  category_colors <- c(category_colors, extra_cols)
}
cols$category <- category_colors

plot(1, xaxt="n",yaxt="n",bty="n",pch="",ylab="",xlab="", main="", sub="", 
     xlim= c(-5,ncol(table_to_use)), ylim = c(-5,nrow(table_to_use)))

if(group_by_tissue_type){
  tissue_cats_bars_xlocs <- sapply(tissue_cats, function(tc) max(strwidth(tc, units = "user"))) + 0.2
  tissue_cats_bars_xlocs <- rep(max(tissue_cats_bars_xlocs), length(tissue_cats_bars_xlocs))
}


if(use_tissue_cols_for_cols){
  heatmap_cols <- sapply((1:max(table_to_use) / max(table_to_use))^opacity_power_scaler, function(opcty) 
    adjustcolor("black", opcty))
} else {
  heatmap_cols <- viridisLite::viridis(n = max(table_to_use, na.rm = T)*100+1)
  heatmap_cols <- heatmap_cols[round(log(1:max(table_to_use, na.rm = T)) / log(max(table_to_use, na.rm = T)) * max(table_to_use, na.rm = T) * 100 + 1)]
}
for(ri in 1:nrow(table_to_use)){
  # text(x = 0.5, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0), pos = 2, labels = rownames(table_to_use)[ri], col = MotrpacRatTraining6moData::TISSUE_COLORS[rownames(table_to_use)[ri]])
  text(x = 0.5, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0), pos = 2, labels = rownames(table_to_use)[ri], col = 1)
  for(ci in 1:ncol(table_to_use)){
    if(ri == 1){
      #trait names
      text(x = ci+0.5, y = -0.9, pos = 2, srt = 45, cex = 0.85,
           labels = gsub("_Scatter", "", trait_categories$new_Phenotype[match(colnames(table_to_use)[ci], trait_categories$Tag)]))
      rect(xleft = ci + 1/2,
           xright = ci - 1/2,
           ybottom = ri - 1/2 - 1,
           ytop =  ri + 1/2 - 1,
           col = category_colors[trait_categories$Category[match(colnames(table_to_use)[ci], trait_categories$Tag)]])
    }
    
    #vertical total # options
    if(ri == nrow(table_to_use)){
      text(x = ci-0.5, y = ri+1 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), pos = 4, srt = 45,
           labels = sig_twas_by_trait_genes_label[colnames(table_to_use)[ci]], cex = 0.9)
    }
    
    #horiz total # of options
    if(ci == ncol(table_to_use)){
      text(x = ci+0.45, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0), pos = 4,
           labels = n_genes_in_nodes_label[rownames(table_to_use)[ri]])
    }
    
    #the actual cells
    if(use_tissue_cols_for_cols){
      rect(xleft = ci + 1/2,
           xright = ci - 1/2,
           ybottom = ri - 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0),
           ytop =  ri + 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0),
           col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[rownames(table_to_use)[ri]], (table_to_use[ri, ci] / max(table_to_use, na.rm = T))^opacity_power_scaler ))  
    } else {
      rect(xleft = ci + 1/2,
           xright = ci - 1/2,
           ybottom = ri - 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0),
           ytop =  ri + 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0),
           col = heatmap_cols[table_to_use[ri, ci]])
    }
    
    
    
    # rect(xleft = ci + 1/2,
    #      xright = ci - 1/2,
    #      ybottom = ri - 0.475,
    #      ytop =  ri + 0.475,
    #      col = heatmap_cols[table_to_use[ri, ci]],
    #      border = MotrpacRatTraining6moData::TISSUE_COLORS[rownames(table_to_use)[ri]])
    
    #text inside of cells
    if(table_to_use[ri, ci] != 0){
      text_in_cell <- table_to_use[ri, ci]
      ndigs <- nchar(text_in_cell)
      if(incl_significance){
        signif_dir <- c("₋", "","⁺")[match(signif_matrix_to_use[ri, ci], -1:1)]
        # text_in_cell <- paste0(text_in_cell, signif_dir)
        
        #try using corner triangles
        ctr <- corner_triangle_ratio <- 1/3
        cxr = ci + 1/2
        cxl = ci - 1/2
        cyb = ri - 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0)
        cyt =  ri + 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0)
        if(signif_dir == "⁺"){
          polygon(x = c(cxl * (ctr) + cxr * (1-ctr), cxr, cxr),
                  y = c(cyt, cyt, cyt * (1-ctr) + cyb * ctr),
                  col = "red")  
        }
        if(signif_dir == "₋"){
          polygon(x = c(cxl * (ctr) + cxr * (1-ctr), cxr, cxr),
                  y = c(cyb, cyb, cyb * (1-ctr) + cyt * ctr),
                  col = "lightblue")  
        }
      }
      text(text_in_cell, x = ci, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0), 
           col = ifelse((table_to_use[ri, ci] / max(table_to_use, na.rm = T))^opacity_power_scaler > opacity_white_threshold, "black", "white"), 
           cex = 0.85 + ifelse(subset_to_traits, 0.1, 0) - (ndigs-1)/6)
      if(incl_cell_totals){
        text(sig_twas_by_trait_genes_matrix[,colnames(table_to_use)][ri, ci], 
             x = ci - 0.7375, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0) - 0.325, 
             col = "black", cex = 0.2, pos = 4, srt = 90)
        text(n_genes_in_nodes_matrix[,colnames(table_to_use)][ri, ci], 
             x = ci - 0.7, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0) - 0.4, 
             col = "black", cex = 0.2, pos = 4, srt = 0)
      }
    }
    
  }
}

#tissue category bars
if(group_by_tissue_type){
  for(bi in 1:length(tissue_cats_bars_xlocs)){
    segments(x0 = -tissue_cats_bars_xlocs[bi],
             x1 = -tissue_cats_bars_xlocs[bi],
             y0 = tissue_cats_bars_ylocs[bi,1],
             y1 = tissue_cats_bars_ylocs[bi,2],
             lwd = 3)
    bracket_length <- 0.2
    segments(x0 = -tissue_cats_bars_xlocs[bi],
             x1 = -tissue_cats_bars_xlocs[bi]+bracket_length,
             y0 = tissue_cats_bars_ylocs[bi,1],
             y1 = tissue_cats_bars_ylocs[bi,1],
             lwd = 3)
    segments(x0 = -tissue_cats_bars_xlocs[bi],
             x1 = -tissue_cats_bars_xlocs[bi]+bracket_length,
             y0 = tissue_cats_bars_ylocs[bi,2],
             y1 = tissue_cats_bars_ylocs[bi,2],
             lwd = 3)
    
    text(x = -tissue_cats_bars_xlocs[bi], y = mean(tissue_cats_bars_ylocs[bi,]) + ifelse(grepl("_", names(tissue_cats)[bi]), -0.25, 0), pos = 2, 
         labels = gsub("Gi", "GI", stringr::str_to_title(gsub("_", "\n", names(tissue_cats)[bi]))), cex = 1.25)
  }  
}

#legend for heatmap
x_adj <- 2.25
y_adj <- 0
yb_adjust <- ifelse(incl_significance, 3, 0)
n_legend_rects_to_use <- 30
n_legend_labels_to_use <- 10
legend_yvals <- round(seq(0, max(table_to_use), length.out = n_legend_labels_to_use))
legend_ylocs <- seq(yb_adjust, 1 + nrow(table_to_use) + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), length.out = n_legend_rects_to_use)
legend_ycols <- round(seq(1, max(table_to_use), length.out = n_legend_rects_to_use))
for(i in 1:(n_legend_rects_to_use-1)){
  yb = legend_ylocs[i]
  yt = legend_ylocs[i+1]
  print(paste(c(yb,yt)))
  rect(xleft = ncol(table_to_use) + x_adj + 2 + 1/2,
       xright = ncol(table_to_use) + x_adj + 2 - 1/2,
       ybottom = yb,
       ytop =  yt,
       col = heatmap_cols[legend_ycols[i]], border = NA)
}
for(i in 1:n_legend_labels_to_use){
  text(labels = legend_yvals[i], x = ncol(table_to_use) + x_adj + 2.4, pos = 4, cex = 0.75,
       y = -0.25 + yb_adjust + legend_yvals[i] / max(table_to_use) * (1+nrow(table_to_use) - yb_adjust + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0)))
}

#overall rect for 0
rect(xleft = ncol(table_to_use) + x_adj + 2 + 1/2,
     xright = ncol(table_to_use) + x_adj + 2 - 1/2,
     ybottom = min(legend_ylocs) - diff(range(legend_ylocs))/50,
     ytop =  max(legend_ylocs))

#legend for significance
if(incl_significance){
  # text(labels = paste0("Pr(Enr.) > ", (1 - signif_threshold), " : X⁺\nPr(Dep.) > ", (1 - signif_threshold), " : X₋"),
  #      x = ncol(table_to_use) + x_adj - 1.75,
  #      y = yb_adjust-3.25, pos = 4, cex = 0.75)
  
  signif_label <- paste0("Pr(Dep.) > ", (1 - signif_threshold), " : ")
  ctr <- 0.5
  cxr = ncol(table_to_use) + x_adj - 2.25 + strwidth(signif_label) + 0.45
  cxl = ncol(table_to_use) + x_adj - 2.25 + strwidth(signif_label) - 0.45
  cyb = yb_adjust - 3.75 - 0.45
  cyt = yb_adjust - 3.75 + 0.45
  
  text(labels = signif_label,
       x = ncol(table_to_use) + x_adj - 1.75,
       y = cyb * 0.325 + cyt * 0.675, pos = 4, cex = 0.75)
  rect(xleft = cxl,
       xright = cxr,
       ybottom = cyb,
       ytop =  cyt)
  polygon(x = c(cxl * (ctr) + cxr * (1-ctr), cxr, cxr),
          y = c(cyb, cyb, cyb * (1-ctr) + cyt * ctr),
          col = "lightblue")  
  
  cyb <- cyb + 1.25
  cyt <- cyt + 1.25
  # cxl <- cxl - 0.125
  # cxr <- cxr - 0.125
  
  signif_label2 <- paste0("Pr(Enr.) > ", (1 - signif_threshold), " : ")
  text(labels = signif_label2,
       x = ncol(table_to_use) + x_adj - 1.75,
       y = cyb * 0.675 + cyt * 0.325, pos = 4, cex = 0.75 * strwidth(signif_label) / strwidth(signif_label2))
  rect(xleft = cxl,
       xright = cxr,
       ybottom = cyb,
       ytop =  cyt)
  polygon(x = c(cxl * (ctr) + cxr * (1-ctr), cxr, cxr),
          y = c(cyt, cyt, cyt * (1-ctr) + cyb * ctr),
          col = "red")  
  
} else {
  text(labels = paste0("* indicate \nIHW-adj. \np-vals < 0.05 \nfrom Fisher's Exact Test"),
       x = ncol(table_to_use) + x_adj - 1.75,
       y = yb_adjust-1.25, pos = 4, cex = 0.75, )
}



#legend for trait categories
if(trait_category_legend_below){
  x_adj2 <- 0
  y_adj2 <- -0.5
  for(i in 1:length(categories_represented)){
    rect(xleft = -1/2 + i + cumsum(c(0,strwidth(categories_represented, units = "user") + 1))[i] + x_adj2,
         xright = 1/2 + i + cumsum(c(0,strwidth(categories_represented, units = "user") + 1))[i] + x_adj2,
         ybottom = -11 + y_adj2,
         ytop = -10 + y_adj2,
         col = category_colors[categories_represented[i]], border = 1)
    text(labels = categories_represented[i], pos = 4, y = -10.5 + y_adj2, x = x_adj2 + 0.35 + i + cumsum(c(0,strwidth(categories_represented, units = "user") + 1))[i])
    #draw circle?
    arc(t1 = 3*pi/2, t2 = pi/2, r1 = (10.5-y_adj2) / 2, r2 = (10.5-y_adj2) / 2, center = c(0,(-10.5 + y_adj2)/2), 
        lwd = 2, res = 100, adjx = ifelse(order_by_counts, 1, 1.1))
    points(0, 0, pch = -9658, cex = 2)
  }
} else {
  x_adj2 <- x_adj - 2.5
  y_adj2 <- y_adj - 2.5
  for(i in 1:length(category_colors)){
    rect(xleft = 0 + i,
         xright = 1 + i,
         ybottom = -10,
         ytop = -9,
         col = category_colors[i], border = 1)
    text(labels = names(category_colors)[i], pos = 4, y = i + y_adj2 - 7, x = ncol(table_to_use) - 1/2 + x_adj2 + 1/3)
  }
}

#labels
#horiz label for total
segments(x0 = -2, x1 = 0.5, y0 = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0),
         y1 = nrow(table_to_use) + 1 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), lwd = 2)
segments(x0 = -2, x1 = 0.5, y0 = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), 
         y1 = nrow(table_to_use) + 2.5 + ifelse(use_range_for_maginal_labels, 0, -0.5) + 
           ifelse(order_by_counts, 0, -0.8) + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0) + 
           ifelse(subset_to_traits, 1, 1.5), lwd = 2)
text(x = -2, y = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), pos = 2, 
     labels = paste0(ifelse(use_range_for_maginal_labels, "total", "max"), " # of\nPrediXcan hits"))

#vertical label for total
segments(x0 = ncol(table_to_use) + 2.5, x1 = ncol(table_to_use) + 1, 
         y0 = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), 
         y1 = nrow(table_to_use) + 0.75 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), lwd = 2)
segments(x0 = ncol(table_to_use) + 2.5, x1 = ncol(table_to_use) + 2 + 
           ifelse(subset_to_traits, 0, 1.25), 
         y0 = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), 
         y1 = nrow(table_to_use) + 0.75 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), lwd = 2)
text(x = ncol(table_to_use) + 2 + ifelse(order_by_counts, 0, 2), y = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), pos = 3, 
     labels = paste0(ifelse(use_range_for_maginal_labels, "total", "max"), " # of DEGs"))


#legend and title labels
text(labels = ifelse(use_counts, latex2exp::TeX("$n_{intersect}$"), latex2exp::TeX("‰_{ intersect}")), pos = 3, font = 2, cex = 1.25,
     x = ncol(table_to_use) + x_adj + 2, y = nrow(table_to_use) + y_adj + 0.875 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0))
if(use_counts){
  text(latex2exp::TeX(paste0("number of genes in 8w - F↓M↓+ 8w - F↑M↑ with IHW-significant PrediXcan at $\\alpha$ = 0.05")), 
       x = 2 + ifelse(subset_to_traits, 0, 20), 
       y = nrow(table_to_use) + 3.25 + 
         ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0) + 
         ifelse(subset_to_traits, 0, 1), pos = 4, cex = 1.5, font = 2)
} else {
  text(latex2exp::TeX(paste0("proportion (‰) of IHW significant TWAS hits at $\\alpha$ = 0.05 in 8w - F↓M↓ or 8w - F↑M↑")), 
       x = 0, y = nrow(table_to_use) + 4.25 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), pos = 4, cex = 2.35, font = 2)
}

fig_label("a)", shrinkX = 0.85, cex = 2, shrinkY = 0.96, xpd = NA)


# ~~~~~~~~~~~~ #
# scatterplots #
# ~~~~~~~~~~~~ #

use_focal_vs_compl <- T

traits <- unique(data1$trait)
tissues <- unique(data1$tissue)
d <- list(cell_count = data1$count,
          total = sapply(1:nrow(data1), function(i) total_number_of_possible_hits_matrix[data1$tissue[i], data1$trait[i]]),
          row_count = sapply(1:nrow(data1), function(i) n_genes_in_nodes_matrix[data1$tissue[i], data1$trait[i]]),
          col_count = sapply(1:nrow(data1), function(i) sig_twas_by_trait_genes_matrix[data1$tissue[i], data1$trait[i]]),
          row_index = match(data1$tissue, tissues),
          col_index = match(data1$trait, traits),
          row_n = length(unique(data1$tissue)),
          col_n = length(unique(data1$trait)))




# category_shapes <- setNames(15:19, salient.categories)
category_shapes <- setNames(15:19, c("Cardiometabolic", "Aging", "Anthropometric", 
                                     "Immune", "Psychiatric-neurologic")
)
category_shapes <- setNames(c(24,19,15,43,23,25)[1:length(focal_traitcats)], focal_traitcats)



par(mar = mars[[2]])
# layout(t(c(rep(1,10), rep(2,10), 3)))
# layout(t(c(1,2,3)), widths = c(1,1,0.1))
#first plot


sapply(1:nrow(signif_df), function(i) adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index[i]]], as.numeric(signif_df$signif[i] != 0) * 0.625 + 0.125))
if(use_focal_vs_compl){
  xvals <- (d$col_count - d$cell_count) / (d$total - d$row_count)
  yvals <- d$cell_count / d$row_count
} else {
  xvals <- d$row_count / d$total * d$col_count / d$total  
  yvals <- d$cell_count / d$total
}




plot(logit(xvals), logit(yvals), xlim = ifelse2(use_focal_vs_compl, c(-14,-1), c(-15,-4.25)), 
     ylim = ifelse2(use_focal_vs_compl, c(-6, -1), c(-9.5,-4.25)),
     pch = category_shapes[traitwise_partitions$Category[match(traits[d$col_index], traitwise_partitions$Tag)]], 
     col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5), bg = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5),
     cex = 1.5,
     xlab = "", ylab = "", yaxt = 'n', xaxt = "n")
box(which = "plot", lwd = 2)
#axes
text(labels = "focal logit(frequency)", x = par("usr")[2] + diff(par("usr")[1:2])/8, srt = 270, pos = 1, 
     y = mean(par("usr")[3:4]), xpd = NA, cex = 1.25)
text(labels = "complement logit(frequency)", x = mean(par("usr")[1:2]), srt = 0, pos = 1, 
     y = par("usr")[3] - diff(par("usr")[3:4])/12, xpd = NA, cex = 1.25)
yaxlocs <- ifelse2(use_focal_vs_compl, c(-6:-1), c(-5:-9))
segments(x0 = par("usr")[2], x1 = par("usr")[2] + diff(par("usr")[1:2])/100, y0 = yaxlocs, y1 = yaxlocs, lwd = 2, xpd = NA)
text(labels = yaxlocs, x = par("usr")[2] + diff(par("usr")[1:2])/200, srt = 0, pos = 4, 
     y = yaxlocs, xpd = NA, cex = 1.25)
xaxlocs <- ifelse2(use_focal_vs_compl, -7:-1*2, -7:-3*2)
segments(x0 = xaxlocs, x1 = xaxlocs, y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4])/80, lwd = 2, xpd = NA)
text(labels = xaxlocs, x = xaxlocs, srt = 0, pos = 1, 
     y = par("usr")[3] - diff(par("usr")[3:4])/80, xpd = NA, cex = 1.25)
#1 to 1 line
abline(0,1, lty = 2, lwd = 4, col = adjustcolor(1,0.75), xpd = F)

#secondplot
xl <- par("usr")[1]
xr <- par("usr")[1] + diff(par("usr")[1:2]) / 2
yt <- par("usr")[4]
yb <- par("usr")[4] - diff(par("usr")[3:4]) / 2

buffer <- 0.075
bx <- (xr-xl) * buffer
by <- (yt-yb) * buffer
rect(xleft = xl, ybottom = yb, xright = xr, ytop = yt, lwd = 2, xpd = NA)

points(x = ((xvals - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
       y = ((yvals - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, 
       pch = category_shapes[traitwise_partitions$Category[match(traits[d$col_index], traitwise_partitions$Tag)]], 
       col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5), bg = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5), 
       cex = 1)
segments(x0 = 0 * (xr - xl - 2*bx) + xl, 
         y0 = 0 * (yt - yb - 2*by) + yb,
         x1 = ((max(yvals, na.rm = T) - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
         y1 = ((max(yvals, na.rm = T) - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, 
         lty = 2, lwd = 2, col = adjustcolor(1,0.75))

#axes
text(labels = "focal frequency", x = par("usr")[1] - diff(par("usr")[1:2])/ifelse2(use_focal_vs_compl, 8, 7), srt = 90, pos = 1, 
     y = par("usr")[4] - diff(par("usr")[3:4])/5, xpd = NA, cex = 1)
text(labels = "complement frequency", x = par("usr")[1] + diff(par("usr")[1:2])/4, srt = 0, pos = 1, 
     y = par("usr")[4] + diff(par("usr")[3:4])/8, xpd = NA, cex = 1)
xaxlocs <- ifelse2(use_focal_vs_compl, c(0:4/20), c(0:6/500))
yaxlocs <- ifelse2(use_focal_vs_compl, c(0:6/20), c(0:6/500))

segments(x0 = par("usr")[1], x1 = par("usr")[1] - diff(par("usr")[1:2])/100, 
         y0 = ((yaxlocs - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, 
         y1 = ((yaxlocs - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, 
         lwd = 2, xpd = NA)
segments(y0 = par("usr")[4], y1 = par("usr")[4] + diff(par("usr")[3:4])/100, 
         x0 = ((xaxlocs - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
         x1 = ((xaxlocs - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
         lwd = 2, xpd = NA)
text(labels = xaxlocs, x = ((xaxlocs - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
     srt = 0, pos = 3, 
     y = par("usr")[4], xpd = NA, cex = 0.75)
text(labels = yaxlocs, x = par("usr")[1], 
     srt = 0, pos = 2, 
     y = ((yaxlocs - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, xpd = NA, cex = 0.75)
# text(labels = 0:6/500, x = -7:-3*2, srt = 0, pos = 1, 
#      y = par("usr")[3] - diff(par("usr")[3:4])/80, xpd = NA, cex = 1.25)
# abline(0,1, lty = 2, lwd = 6, col = adjustcolor(1,0.75))

#legends
legend(x = xl, y = yb, legend = names(category_shapes), pch = category_shapes, col = adjustcolor(1, 0.5), pt.bg = adjustcolor(1, 0.5), cex = 0.75, pt.cex = 1.4)
legend(x = xr, y = yt, legend = tissues, pch = 19, col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues], 0.5), cex = 0.55, pt.cex = 1.2)
segments(x0 = par("usr")[1] + diff(par("usr")[1:2]) / 1.55, 
         y0 = par("usr")[4] - diff(par("usr")[3:4])/50,
         x1 = par("usr")[1] + diff(par("usr")[1:2]) / 1.45, 
         y1 = par("usr")[4] - diff(par("usr")[3:4])/50, 
         lty = 2, lwd = 2, col = adjustcolor(1,0.75))
text(x = par("usr")[1] + diff(par("usr")[1:2]) / 1.45, 
     y = par("usr")[4] - diff(par("usr")[3:4])/45, 
     pos = 4, cex = 0.75, labels = "1-to-1 line")


#make second plot with raw proportion data
par(mar = mars[[3]])

xvals <- data$count_all / data$total_all
yvals <- data$count_inters / data$total_inters


cex_pow <- 1/2
plot(xvals, yvals, xlim = c(0.48, 0.52), 
     ylim = c(0,1),
     pch = category_shapes[traitwise_partitions$Category[match(data$trait, traitwise_partitions$Tag)]], 
     col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[data$tissue], 0.5), bg = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[data$tissue], 0.5), 
     cex = (data$total_inters / max(data$total_inters))^cex_pow * 5,
     xlab = "", ylab = "", yaxt = 'n', xaxt = "n")
box(which = "plot", lwd = 2)
#axes
text(labels = "focal frequency (+ effects)", x = par("usr")[1] - diff(par("usr")[1:2])/7, srt = 90, pos = 1, 
     y = mean(par("usr")[3:4]), xpd = NA, cex = 1.25)
text(labels = "complement frequency (+ effects)", x = mean(par("usr")[1:2]), srt = 0, pos = 1, 
     y = par("usr")[3] - diff(par("usr")[3:4])/12, xpd = NA, cex = 1.25)
yaxlocs <- seq(0,1,by=0.1)
segments(x0 = par("usr")[1], x1 = par("usr")[1] - diff(par("usr")[1:2])/100, y0 = yaxlocs, y1 = yaxlocs, lwd = 2, xpd = NA)
text(labels = yaxlocs, x = par("usr")[1] - diff(par("usr")[1:2])/200, srt = 0, pos = 2, 
     y = yaxlocs, xpd = NA, cex = 1.25)
xaxlocs <- round(seq(par("usr")[1], par("usr")[2], by = 0.01), 2)
segments(x0 = xaxlocs, x1 = xaxlocs, y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4])/80, lwd = 2, xpd = NA)
text(labels = xaxlocs, x = xaxlocs, srt = 0, pos = 1, 
     y = par("usr")[3] - diff(par("usr")[3:4])/80, xpd = NA, cex = 1.25)
# abline(0,1, lty = 2, lwd = 6, col = adjustcolor(1,0.75))
abline(h=0.5,lwd=2,lty=2,col=adjustcolor(1,0.5),xpd=F)
abline(v=0.5,lwd=2,lty=2,col=adjustcolor(1,0.5),xpd=F)

#legend
xl <- par("usr")[2]
xr <- par("usr")[2] + diff(par("usr")[1:2]) / 2
yt <- par("usr")[4]
yb <- par("usr")[4] - diff(par("usr")[3:4]) / 2.5
legend(x = xl, y = yb + (yt-yb)/20, legend = names(category_shapes), bty = "n",
       pch = category_shapes, col = adjustcolor(1, 0.5), pt.bg = adjustcolor(1, 0.5), cex = 0.75, pt.cex = 1.4, xpd = NA)
legend(x = xl + (xr-xl)/100, y = yt, legend = tissues, pch = 19, bty = "n", xpd = NA,
       col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues], 0.5), cex = 0.55, pt.cex = 1.2)

legcexs <- c(0.5, 1:5)
legvals <- (round((legcexs / 5) ^ (1/cex_pow) * max(data$total_inters)))
legcexs <- (legvals / max(data$total_inters))^cex_pow * 5
# points(x = rep(xl + (xr-xl)/11, 6) + (rev(legcexs)) / 4000, y = yb - (yt-yb) * 1.4 + 
#          cumsum(rep((yt-yb)/60, length(legcexs))) + cumsum(legcexs + c(0,legcexs[-length(legcexs)])) / 120, 
#        cex = legcexs, xpd = NA, pch = 19, col = adjustcolor(1, 0.5))
text(labels = latex2exp::TeX("$n_{intersect}$"), x = xl - (xr-xl)/50, y = yb - (yt-yb) * 0.6, pos = 4, xpd = NA)
fixed_incr <- (yt-yb)/25
y_disp <- (yt-yb) * 1.5
points(x = rep(xl + (xr-xl)/10, 6), y = yb - y_disp + 
         cumsum(rep(fixed_incr, length(legcexs))) + cumsum(legcexs + c(0,legcexs[-length(legcexs)])) / 120, 
       cex = legcexs, xpd = NA, pch = 19, col = adjustcolor(1, 0.5))
text(x = rep(xl + (xr-xl)/10, 6), y = yb - y_disp + 
       cumsum(rep(fixed_incr, length(legcexs))) + cumsum(legcexs + c(0,legcexs[-length(legcexs)])) / 120, 
     cex = legcexs / 5, xpd = NA, labels = legvals, col = "white")
text(labels = legvals[1:2], x = xl + (xr-xl)/15 + cumsum(legcexs)[1:2]/4000, y = yb - y_disp + 
       cumsum(rep(fixed_incr, length(legcexs)))[1:2] + cumsum(legcexs + c(0,legcexs[-length(legcexs)]))[1:2] / 120, 
     pos = 4, xpd = NA, cex = 0.5)


fig_label("b)", shrinkX = 0.865, cex = 2, shrinkY = 0.95, xpd = NA)
fig_label("c)", shrinkX = 0.99, shrinkY = 0.95, cex = 2)


# ~~~~~~~~~~~~ #
# violin plots #
# ~~~~~~~~~~~~ #

par(mar = mars[[4]])
focal_samps <- subset_samps("row_bias", c("raw", "sd"), samps = samps) + samps$overall_bias
n_to_trim <- round(nrow(focal_samps) * 0.002)
focal_samps <- apply(focal_samps, 2, trim_n, n_to_trim)

tissues <- tissues_intersect.model
colnames(focal_samps) <- tissues
tord <- order(apply(focal_samps, 2, mean))
qi_95 <- apply(focal_samps, 2, quantile, probs = c(0.05, 0.95))[,tord]
qi_100 <- apply(focal_samps, 2, quantile, probs = c(0.0, 1.0))[,tord]
posterior_means <- apply(focal_samps, 2, mean)[tord]
inside_cols <- rep("black", length(tord))
inside_cols[overlaps_with_zero(qi_95)] <- "white"
tmp <- vioplot::vioplot(x = focal_samps[,tord], T, 
                        col = MotrpacRatTraining6moData::TISSUE_COLORS[colnames(focal_samps)][tord], outline=FALSE, yaxt = "n",
                        names = colnames(focal_samps)[tord], range = 0, ylab = "", lineCol = MotrpacRatTraining6moData::TISSUE_COLORS[colnames(focal_samps)][tord],
                        xlab = "", cex.lab = 2, plotCentre = "point", 
                        colMed = MotrpacRatTraining6moData::TISSUE_COLORS[colnames(focal_samps)][tord],
                        horizontal = T)

#axes
xtickvals <- -1:2/5
segments(x0 = xtickvals, x1 = xtickvals, y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4])/100, xpd = NA)
text(x = xtickvals, y = par("usr")[3] - diff(par("usr")[3:4])/30, labels = xtickvals)
text(x = mean(par("usr")[1:2]), y = par("usr")[3]-diff(par("usr")[3:4])/10, 
     labels = "multilevel tissue\nenrichment (logodds-scale)", cex = 1.25, xpd = NA)
tick <- seq_along(colnames(focal_samps))
axis(2, at = tick, labels = F)
for(i in 1:length(colnames(focal_samps))){
  segments(y0 = i, y1 = i, x0 = qi_95[1,i], x1 = qi_95[2,i], lwd = 4, col = inside_cols[i])
  points(y = i, x = posterior_means[i], col = inside_cols[i], pch = 19, cex = 1.5)
  
  #guiding lines
  segments(y0 = i, y1 = i, x0 = min(qi_100), x1 = qi_100[1,i], lwd = 1, col = "black", xpd = T, lty = 3)
  segments(y0 = i, y1 = i, x0 = qi_100[2,i], x1 = max(qi_100), lwd = 1, col = "black", xpd = T, lty = 3)
  
}
#group labels
text(y = tick, x = rep(par("usr")[1] - 0.01, length(tick)), colnames(focal_samps)[tord], srt = 0, xpd = T, pos = 2,
     col = MotrpacRatTraining6moData::TISSUE_COLORS[colnames(focal_samps)[tord]], xpd = NA)

abline(v=0,lwd=3,lty=2, col = adjustcolor(1,0.5), xpd = F)

#figure label
if(T){
  fig_label("d)", cex = 2, shrinkX = 0.9, shrinkY = 0.97)
}

#legend for violin plots
xval <- -0.85
yval <- 0.5 #2
yh <- 4
yt <- yval + yh/2
yb <- yval - yh/2
xr <- -30:30/10
xvals_disp <- (dnorm(xr) / 4)
xvals <- cbind(xval + xvals_disp, xval - xvals_disp)
yvals <- seq(yb,yt,length.out = nrow(xvals))
# lines(xvals[,1], yvals, xpd = NA)
# lines(xvals[,2], yvals, xpd = NA)
polygon(x = c(xvals[,1], xvals[,2]), y = c(yvals, rev(yvals)), xpd = NA, col = adjustcolor(1,0.15))
ybloc <- yb + (qnorm(p = c(0.05,0.95)) - min(xr)) / diff(range(xr)) * yh
segments(x0 = xval, x1 = xval, y0 = ybloc[1], y1 = ybloc[2], xpd = NA, lwd = 2)
points(xval, yval, pch = 19, xpd = NA)
xlab_loc <- max(xvals) + diff(range(xvals)) / 5
segments(x0 = xval, x1 = xlab_loc, 
         y0 = c(max(yvals), min(yvals), yval, ybloc),
         y1 = c(max(yvals), min(yvals), yval, ybloc),
         xpd = NA, lty = 2)
violabs <- paste0("$P_{", c("99.9", "0.01", "\\mu", "5", "95"), "}$")
violabs[[3]] <- "$\\mu$"
text(latex2exp::TeX(violabs), x = xlab_loc - diff(range(xvals)) / 10, y = c(max(yvals), min(yvals), yval, ybloc), pos = 4, xpd = NA, cex = 0.875)

segments(x0 = min(xvals) - diff(range(xvals))/100, x1 = max(xvals) + diff(range(xvals))/100, 
         y0 = min(yvals) - diff(range(yvals)) / 8, y1 = min(yvals) - diff(range(yvals)) / 8, xpd = NA, lwd = 3.25)
points(xval, min(yvals) - diff(range(yvals)) / 8, pch = 19, xpd = NA, cex = 1.2)
segments(x0 = min(xvals), x1 = max(xvals), y0 = min(yvals) - diff(range(yvals)) / 8, y1 = min(yvals) - diff(range(yvals)) / 8, xpd = NA, lwd = 2, col = "white")
points(xval, min(yvals) - diff(range(yvals)) / 8, pch = 19, xpd = NA, cex = 1, col = "white")
text(latex2exp::TeX("$CI_{95}$ includes 0$"), x = max(xvals) + diff(range(xvals))/100, 
     y = min(yvals) - diff(range(yvals)) / 7, pos = 4, xpd = NA, cex = 0.875)

segments(x0 = min(xvals) - diff(range(xvals))/100, x1 = max(xvals) + diff(range(xvals))/100, 
         y0 = min(yvals) - diff(range(yvals)) / 8 * 2, y1 = min(yvals) - diff(range(yvals)) / 8 * 2, xpd = NA, lwd = 3.25)
points(xval, min(yvals) - diff(range(yvals)) / 8 * 2, pch = 19, xpd = NA, cex = 1.2)
text(latex2exp::TeX("$CI_{95}$ excludes 0$"), x = max(xvals) + diff(range(xvals))/100, 
     y = min(yvals) - diff(range(yvals)) / 7 - diff(range(yvals)) / 8, 
     pos = 4, xpd = NA, cex = 0.875)


###
#colcat bias violin plots
###

salient.categories <- trait_cats
if(length(salient.categories) < 9){
  category_colors <- RColorBrewer::brewer.pal(length(salient.categories), "Dark2")
  names(category_colors) <- sort(salient.categories)  
} else {
  category_colors <- RColorBrewer::brewer.pal(5, "Dark2")
  names(category_colors) <- sort(c("Cardiometabolic", "Aging", "Anthropometric", 
                                   "Immune", "Psychiatric-neurologic"))
  cats_leftover <- setdiff(salient.categories, names(category_colors))
  cols_leftover <- setdiff(RColorBrewer::brewer.pal(8, "Dark2"), category_colors)
  n_cols_needed <- length(cats_leftover) - length(cols_leftover)
  if(n_cols_needed > 0){
    extra_cols <- c(cols_leftover, RColorBrewer::brewer.pal(n_cols_needed, "Set2"))
  } else {
    extra_cols <- cols_leftover
  }
  names(extra_cols) <- cats_leftover
  category_colors <- c(category_colors, extra_cols)
}
cols$category <- category_colors

par(mar = mars[[5]])
focal_samps <- subset_samps("colcat_bias", c("raw", "sd"), samps = samps) + samps$overall_bias
focal_samps <- apply(focal_samps, 2, trim_n, n_to_trim)
colnames(focal_samps) <- trait_cats
tord <- order(apply(focal_samps, 2, mean))
qi_95 <- apply(focal_samps, 2, quantile, probs = c(0.05, 0.95))[,tord]
qi_100 <- apply(focal_samps, 2, quantile, probs = c(0.0, 1.0))[,tord]
posterior_means <- apply(focal_samps, 2, mean)[tord]
inside_cols <- rep("black", length(tord))
inside_cols[overlaps_with_zero(qi_95)] <- "white"
tmp <- vioplot::vioplot(x = focal_samps[,tord], T, 
                        col = cols$category[colnames(focal_samps)][tord], outline=FALSE, yaxt = "n",
                        names = colnames(focal_samps)[tord], range = 0, ylab = "", lineCol = MotrpacRatTraining6moData::TISSUE_COLORS[colnames(focal_samps)][tord],
                        xlab = "", cex.lab = 2, plotCentre = "point", 
                        colMed = cols$category[colnames(focal_samps)][tord],
                        horizontal = T)

#axes
xtickvals <- -2:2/5
segments(x0 = xtickvals, x1 = xtickvals, y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4])/100, xpd = NA)
text(x = xtickvals, y = par("usr")[3] - diff(par("usr")[3:4])/30, labels = xtickvals)
text(x = mean(par("usr")[1:2]), y = par("usr")[3]-diff(par("usr")[3:4])/10, 
     labels = "multilevel category\nenrichment (logodds-scale)", cex = 1.25, xpd = NA)
tick <- seq_along(colnames(focal_samps))
axis(2, at = tick, labels = F)
for(i in 1:length(colnames(focal_samps))){
  segments(y0 = i, y1 = i, x0 = qi_95[1,i], x1 = qi_95[2,i], lwd = 4, col = inside_cols[i])
  points(y = i, x = posterior_means[i], col = inside_cols[i], pch = 19, cex = 1.5)
  
  #guiding lines
  segments(y0 = i, y1 = i, x0 = min(qi_100), x1 = qi_100[1,i], lwd = 1, col = "black", xpd = T, lty = 3)
  segments(y0 = i, y1 = i, x0 = qi_100[2,i], x1 = max(qi_100), lwd = 1, col = "black", xpd = T, lty = 3)
  
}
text(y = tick, x = rep(par("usr")[1] - 0.01, length(tick)), 
     labels = colnames(focal_samps)[tord], 
     srt = 0, xpd = T, pos = 2,
     col = cols$category[colnames(focal_samps)[tord]], xpd = NA)

gsub("_", "\n", colnames(focal_samps))
gsub("Gi", "GI", stringr::str_to_title(gsub("_", "\n", colnames(focal_samps))))

abline(v=0,lwd=3,lty=2, col = adjustcolor(1,0.5), xpd = F)

#figure label
if(T){
  fig_label("e)", cex = 2, shrinkX = 0.905, shrinkY = 0.985)
}

##
# traits
##

par(mar = mars[[6]])
trait_subset_bool <- traitwise_partitions$Category[match(traits, traitwise_partitions$Tag)] %in% salient.categories
focal_samps <- subset_samps("col_bias", c("raw", "sd"), samps = samps) + samps$overall_bias + 
  subset_samps("colcat_bias", c("raw", "sd"), samps = samps)[,match(trait_categories$Category[match(traits, trait_categories$Tag)], trait_cats)]
focal_samps <- apply(focal_samps, 2, trim_n, n_to_trim)
colnames(focal_samps) <- traits
tord <- order(apply(focal_samps, 2, mean))
qi_95 <- apply(focal_samps, 2, quantile, probs = c(0.05, 0.95))[,tord]
qi_100 <- apply(focal_samps, 2, quantile, probs = c(0.0, 1.0))[,tord]
posterior_means <- apply(focal_samps, 2, mean)[tord]
inside_cols <- rep("black", length(tord))
inside_cols[overlaps_with_zero(qi_95)] <- "white"
tmp <- vioplot::vioplot(x = focal_samps[,tord], T, 
                        col = cols$category[trait_categories$Category[match(colnames(focal_samps)[trait_subset_bool][tord], trait_categories$Tag)]], 
                        outline=FALSE, yaxt = "n",
                        names = colnames(focal_samps)[tord], range = 0, ylab = "", 
                        lineCol = cols$category[trait_categories$Category[match(colnames(focal_samps)[trait_subset_bool][tord], trait_categories$Tag)]],
                        xlab = "", cex.lab = 2, plotCentre = "point", 
                        colMed = cols$category[colnames(focal_samps)][tord],
                        horizontal = T,
                        xlim = c(4,ncol(focal_samps)-3))

#axes
xtickvals <- -2:2/2
segments(x0 = xtickvals, x1 = xtickvals, y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4])/200, xpd = NA)
text(x = xtickvals, y = par("usr")[3] - diff(par("usr")[3:4])/75, labels = xtickvals)
text(x = mean(par("usr")[1:2]), y = par("usr")[3]-diff(par("usr")[3:4])/25, 
     labels = "multilevel trait\nenrichment (logodds-scale)", cex = 1.25, xpd = NA)
tick <- seq_along(colnames(focal_samps))
axis(2, at = tick, labels = F)
for(i in 1:length(colnames(focal_samps))){
  segments(y0 = i, y1 = i, x0 = qi_95[1,i], x1 = qi_95[2,i], lwd = 2, col = inside_cols[i])
  points(y = i, x = posterior_means[i], col = inside_cols[i], pch = 19, cex = 0.75)
  
  #guiding lines
  segments(y0 = i, y1 = i, x0 = min(qi_100), x1 = qi_100[1,i], lwd = 1, col = "black", xpd = T, lty = 3)
  segments(y0 = i, y1 = i, x0 = qi_100[2,i], x1 = max(qi_100), lwd = 1, col = "black", xpd = T, lty = 3)
  
}
text(y = tick, x = rep(par("usr")[1] - 0.01, length(tick)), 
     labels = trait_categories$new_Phenotype[match(colnames(focal_samps)[trait_subset_bool][tord], trait_categories$Tag)], 
     srt = 0, xpd = T, pos = 2, cex = 0.75,
     col = cols$category[trait_categories$Category[match(colnames(focal_samps)[trait_subset_bool][tord], trait_categories$Tag)]], xpd = NA)




abline(v=0,lwd=3,lty=2, col = adjustcolor(1,0.5), xpd = F)

#figure label
if(T){
  fig_label("f)", cex = 2, shrinkX = 0.95, shrinkY = 0.9855)
}
dev.off()



#### end ####
