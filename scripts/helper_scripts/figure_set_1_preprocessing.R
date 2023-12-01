#load libraries
library(Cairo)
library(dplyr)
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
library(invgamma)
library(MotrpacRatTraining6mo)

#### define functions ####
source(file = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/helper_scripts/deg-trait_functions.R")

#general graphical parameters
motrpac_gtex_map = c('t30-blood-rna'='Whole Blood',
                     't52-hippocampus'='Brain - Hippocampus',
                     't53-cortex'='Brain - Cortex',
                     't54-hypothalamus'='Brain - Hypothalamus',
                     't55-gastrocnemius'='Muscle - Skeletal',
                     't56-vastus-lateralis'='Muscle - Skeletal',
                     't58-heart'='Heart - Left Ventricle',
                     't59-kidney'='Kidney - Cortex',
                     't60-adrenal'='Adrenal Gland',
                     't61-colon'='Colon - Transverse',
                     't62-spleen'='Spleen',
                     't63-testes'='Testis',
                     't64-ovaries'='Ovary',
                     't66-lung'='Lung',
                     't67-small-intestine'='Small Intestine - Terminal Ileum',
                     't68-liver'='Liver',
                     't70-white-adipose'='Adipose - Subcutaneous')

cols = list(Tissue=MotrpacRatTraining6moData::TISSUE_COLORS[names(motrpac_gtex_map)], 
            Time=MotrpacRatTraining6moData::GROUP_COLORS[paste0(c(1,2,4,8), "w")],
            Sex=MotrpacRatTraining6moData::SEX_COLORS[c('male','female')])
cols$Tissue[which(is.na(cols$Tissue))] <- '#C0C0C0'

#### figure preprocessing ####
gwas_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/imputed_gwas_hg38_1.1/"
gwas_summary_files <- list.files(gwas_dir)
gwas_summary_files <- gwas_summary_files[-grep(gwas_summary_files, pattern = "README")]
gwas_names <- stringr::str_replace_all(gwas_summary_files, ".txt.gz", "")
coloc_phenotypes <- stringr::str_replace_all(gwas_names, "imputed_", "")
ldsc_output_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/output/"
ldsc_results_paths <- list.files(ldsc_output_dir)
ldsc_log_paths <- paste0(ldsc_output_dir, ldsc_results_paths[grep(ldsc_results_paths, pattern = "log")])
ldsc_results_paths <- paste0(ldsc_output_dir, ldsc_results_paths[grep(ldsc_results_paths, pattern = "results")])
cluster_names <- paste0("Cluster_", 1:15)
ldsc_results_paths <- lapply(gwas_names, function(gwas) ldsc_results_paths[grep(pattern = gwas, ldsc_results_paths)])
ldsc_log_paths <- lapply(gwas_names, function(gwas) ldsc_log_paths[grep(pattern = gwas, ldsc_log_paths)])
names(ldsc_results_paths) <- names(ldsc_log_paths) <- gwas_names


log_files <- as.data.frame(matrix(data = NA, ncol = 5, nrow = length(gwas_names) * length(cluster_names)))
colnames(log_files) <- c("h2", "h2se", "chi2", "gwas", "cluster")
log_files$cluster <- rep(1:length(cluster_names), length(gwas_names))
log_files$gwas <- c(sapply(gwas_names, function(gwas_name) rep(gwas_name, length(cluster_names))))
for(i in 1:nrow(log_files)){
  logfile <- readLines(ldsc_log_paths[[log_files$gwas[i]]][grep(pattern = paste0(cluster_names[log_files$cluster[i]], ".log"), ldsc_log_paths[[log_files$gwas[i]]])][1])
  h2 <- logfile[grep(logfile, pattern = "Total Observed scale h2")]
  log_files$h2[i] <- as.numeric(strsplit(h2, " ")[[1]][5])
  log_files$h2se[i] <- as.numeric(substr(x = strsplit(h2, " ")[[1]][6], start = 2, stop = nchar(strsplit(h2, " ")[[1]][6]) - 1))
  chi2 <- as.numeric(strsplit(logfile[grep(logfile, pattern = "Mean Chi")], split = " ")[[1]][3])
  log_files$chi2[i] <- chi2
}
log_files$gwas <- stringr::str_replace_all(log_files$gwas, "imputed_", "")

#snag estimated heritabilities
log_files$pval <- 1 - pnorm(log_files$h2 / log_files$h2se)
estimated_h2_ldsc <- sapply(setNames(unique(log_files$gwas), unique(log_files$gwas)), function(gwas) mean(log_files$h2[log_files$gwas == gwas]))
estimated_h2_ldsc[estimated_h2_ldsc < 0] <- min(estimated_h2_ldsc[estimated_h2_ldsc > 0])
estimated_h2_ldsc_pval <- sapply(setNames(unique(log_files$gwas), unique(log_files$gwas)), function(gwas) mean(log_files$pval[log_files$gwas == gwas]))

mesc_output_basic <- fread(file = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/mesc_out_basic.txt")
estimated_h2_mesc <- mesc_output_basic[mesc_output_basic$Quantity == "h2",]
estimated_h2_mesc <- sapply(setNames(unique(estimated_h2_mesc$trait), unique(estimated_h2_mesc$trait)), 
                            function(gwas) mean(estimated_h2_mesc$Estimate[estimated_h2_mesc$trait == gwas]))
estimated_h2_mesc[estimated_h2_mesc < 0] <- min(estimated_h2_mesc[estimated_h2_mesc > 0])
estimated_h2_mesc_pval <- mesc_output_basic[mesc_output_basic$Quantity == "h2",]
estimated_h2_mesc_pval <- sapply(setNames(unique(estimated_h2_mesc_pval$trait), unique(estimated_h2_mesc_pval$trait)), 
                                 function(gwas) mean(estimated_h2_mesc_pval$Estimate_pval[estimated_h2_mesc_pval$trait == gwas]))

estimated_h2med_mesc <- mesc_output_basic[mesc_output_basic$Quantity == "h2med",]
estimated_h2med_mesc <- sapply(setNames(unique(estimated_h2med_mesc$trait), unique(estimated_h2med_mesc$trait)), 
                               function(gwas) mean(estimated_h2med_mesc$Estimate[estimated_h2med_mesc$trait == gwas]))
estimated_h2med_mesc[estimated_h2med_mesc < 0] <- min(estimated_h2med_mesc[estimated_h2med_mesc > 0])
estimated_h2med_mesc_pval <- mesc_output_basic[mesc_output_basic$Quantity == "h2med",]
estimated_h2med_mesc_pval <- sapply(setNames(unique(estimated_h2med_mesc_pval$trait), unique(estimated_h2med_mesc_pval$trait)), 
                                    function(gwas) mean(estimated_h2med_mesc_pval$Estimate_pval[estimated_h2med_mesc_pval$trait == gwas]))

estimated_h2med_over_h2_mesc <- mesc_output_basic[mesc_output_basic$Quantity == "h2med",]
estimated_h2med_over_h2_mesc <- sapply(setNames(unique(estimated_h2med_over_h2_mesc$trait), unique(estimated_h2med_over_h2_mesc$trait)), 
                                       function(gwas) mean(estimated_h2med_over_h2_mesc$Estimate_over_h2[estimated_h2med_over_h2_mesc$trait == gwas]))
estimated_h2med_over_h2_mesc[estimated_h2med_over_h2_mesc < 0] <- min(estimated_h2med_over_h2_mesc[estimated_h2med_over_h2_mesc > 0])
estimated_h2med_over_h2_mesc_pval <- mesc_output_basic[mesc_output_basic$Quantity == "h2med",]
estimated_h2med_over_h2_mesc_pval <- sapply(setNames(unique(estimated_h2med_over_h2_mesc_pval$trait), unique(estimated_h2med_over_h2_mesc_pval$trait)), 
                                            function(gwas) mean(estimated_h2med_over_h2_mesc_pval$Estimate_over_h2_pval[estimated_h2med_over_h2_mesc_pval$trait == gwas]))

bonf_p_h2 <- 0.05 / length(estimated_h2med_over_h2_mesc_pval)

# plot(estimated_h2_mesc, estimated_h2_ldsc[names(estimated_h2_mesc)])
# plot(log10(estimated_h2_mesc), log10(estimated_h2_ldsc[names(estimated_h2_mesc)]))
# abline(0,1)

#subset to insteresting trait categories
trait_categories <- read.csv("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/gwas_metadata.csv", header = T)
traitwise_partitions <- trait_categories[,c("Tag", "Category")]
traitname_map <- trait_categories[,c("Tag", "new_Phenotype")]
log_files$trait_category <- traitwise_partitions$Category[match(log_files$gwas, traitwise_partitions$Tag)]
salient.categories <- c("Cardiometabolic", "Aging", "Anthropometric", "Immune", "Psychiatric-neurologic")
salient.categories <- unique(traitwise_partitions$Category)
log_files <- log_files[log_files$trait_category %in% salient.categories,]

total_h2_sigma_thresh <- 7
traits_with_satisfactory_heritaility <- unique(log_files$gwas[log_files$h2 / log_files$h2se > total_h2_sigma_thresh])

#do the gcor
ldsc_directory <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/"

#let's get the gwas summary stats into a format that munge_sumstats.py can process
traitnames <- traits_with_satisfactory_heritaility
gcors <- lapply(traitnames, function(x) integer())
names(gcors) <- traitnames  
for(gi in traitnames){
  logfile <- readLines(paste0(ldsc_directory, "output/signed_gcor/imputed_", gi, "_pairwise-Gcorrs.log"))
  logfile <- logfile[(grep(x = logfile, pattern = "Summary of Genetic Correlation Results")+1):(length(logfile)-3)]
  writeLines(logfile, con = paste0(ldsc_directory, "temp.txt"))
  gcors[[gi]] <- fread(paste0(ldsc_directory, "temp.txt"))
  file.remove(paste0(ldsc_directory, "temp.txt"))
  
  gcors[[gi]]$p1 <- sapply(gcors[[gi]]$p1, function(trait) strsplit(trait, ".txt.gz")[[1]][1])
  gcors[[gi]]$p1 <- sapply(gcors[[gi]]$p1, function(trait) strsplit(trait, "gwas_sumstats/proper_format/signed/imputed_")[[1]][2])
  gcors[[gi]]$p2 <- sapply(gcors[[gi]]$p2, function(trait) strsplit(trait, ".txt.gz")[[1]][1])
  gcors[[gi]]$p2 <- sapply(gcors[[gi]]$p2, function(trait) strsplit(trait, "gwas_sumstats/proper_format/signed/imputed_")[[1]][2])
}

#construct genetic correlation matrix
gcor_mat <- diag(length(traitnames))
colnames(gcor_mat) <- rownames(gcor_mat) <- traitnames
for(ri in rownames(gcor_mat)){
  for(ci in colnames(gcor_mat)){
    gcor_mat[ri, ci] <- gcors[[ri]]$rg[match(ci, gcors[[ri]]$p2)]  
  }
}
diag(gcor_mat) <- rep(1, length(traitnames))

#remove impossible correlations
gcor_mat <- gcor_mat[traits_with_satisfactory_heritaility,traits_with_satisfactory_heritaility]
troublesome_inds_NA_corrs <- as.numeric(names(which(table(which(is.na(gcor_mat), arr.ind = T)) > 1)))
troublesome_inds_NA_corrs_traits <- rownames(gcor_mat)[troublesome_inds_NA_corrs]
if(length(troublesome_inds_NA_corrs) > 0){
  gcor_mat <- gcor_mat[-troublesome_inds_NA_corrs, -troublesome_inds_NA_corrs]
}
troublesome_inds_imposs_corrs <- as.numeric(names(which(table(which(gcor_mat > 1 | gcor_mat < -1, arr.ind = T)) > 1)))
troublesome_inds_imposs_corrs_traits <- rownames(gcor_mat)[troublesome_inds_imposs_corrs]
if(length(troublesome_inds_imposs_corrs) > 0){
  gcor_mat <- gcor_mat[-troublesome_inds_imposs_corrs, -troublesome_inds_imposs_corrs]
}

rownames(gcor_mat) <- colnames(gcor_mat) <- traitname_map$new_Phenotype[match(rownames(gcor_mat), traitname_map$Tag)]

#perform hypothesis testing
gcor_mat_se <- diag(length(traitnames))
colnames(gcor_mat_se) <- rownames(gcor_mat_se) <- traitnames
for(ri in rownames(gcor_mat_se)){
  for(ci in colnames(gcor_mat_se)){
    gcor_mat_se[ri, ci] <- gcors[[ri]]$se[match(ci, gcors[[ri]]$p2)]
  }
}
diag(gcor_mat_se) <- rep(NA, length(traitnames))
rownames(gcor_mat_se) <- colnames(gcor_mat_se) <- traitname_map$new_Phenotype[match(rownames(gcor_mat_se), traitname_map$Tag)]
gcor_mat_se <- gcor_mat_se[rownames(gcor_mat), rownames(gcor_mat)]
gcor_mat_z <- abs(gcor_mat / gcor_mat_se)
zscore_sig_thresh <- abs(qnorm(0.025  / choose(dim(gcor_mat)[1], 2), 0, 1))
gcor_mat_sig <- gcor_mat_z > zscore_sig_thresh
diag(gcor_mat_sig) <- T
sum(gcor_mat_sig, na.rm = T) - nrow(gcor_mat_sig)

#get a few last graphical parameters
unique_trait_categories <- unique(traitwise_partitions$Category)
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

#### fig 1 extra ####

if(figure_id == 1){
  
  #check if grex script has been called
  if(!exists("called_grex_script")){called_grex_script <- F}
  if(!called_grex_script){
    called_grex_script <- T
    source("/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/analyses/analysis_GREx_RelativeEffectSize.R")
  }
  if(!exists("relative_expression_data")){
    load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/relative_expression_motrpac_gtex")
  }
  
  load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/node_metadata_list.RData")
  ensembl_genes <- orig_ensembl_genes <- lapply(split(node_metadata_list$`8w`$human_ensembl_gene[!is.na(node_metadata_list$`8w`$human_ensembl_gene)], 
                                                      node_metadata_list$`8w`$tissue[!is.na(node_metadata_list$`8w`$human_ensembl_gene)]), unique)
  symbol_map <- unique(node_metadata_list$`8w`[,c("human_gene_symbol", "human_ensembl_gene")])
  all_genes <- unlist(orig_ensembl_genes)
  n_tissues_per_gene <- table(all_genes)
  ensembl_genes$THREE <- names(n_tissues_per_gene)[n_tissues_per_gene > 2]
  
  jacmat <- sapply(orig_ensembl_genes, function(x) sapply(orig_ensembl_genes, function(y) jaccard(x,y)))
  jacmat_inds <- order(cmdscale(1-jacmat, k = 1))
  jacmat <- jacmat[jacmat_inds, jacmat_inds]
  n_tissues_per_gene <- table(unlist(orig_ensembl_genes))
  nt2pg <- table(n_tissues_per_gene)
  ensembl_genes_df <- data.frame(gene = unlist(orig_ensembl_genes),
                                 tissue = rep(names(orig_ensembl_genes), sapply(orig_ensembl_genes, length), each = T))
  ensembl_genes_df$n_tissue <- n_tissues_per_gene[match(ensembl_genes_df$gene, names(n_tissues_per_gene))]
  n_per_cat <- lapply(setNames(sort(unique(n_tissues_per_gene)), sort(unique(n_tissues_per_gene))), 
                      function(nt){x <- table(ensembl_genes_df$tissue[ensembl_genes_df$n_tissue == nt]); x})
  prop_per_cat <- lapply(n_per_cat, function(nt) nt / sum(nt))
  
  #now do the Open Targets curves
  if(!exists("n_traits_above_at_least_1")){
    use_indirect <- F
    load(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/open-targets_tissue-x-disease_", 
                ifelse(use_indirect, "indirect", "direct"), 
                "-associations"))
    breakpoints <- 0:100/100
    n_above <- sapply(names(tissue_x_disease), function(tissue) 
      log10((data.frame(cats = cut(tissue_x_disease[[tissue]][tissue_x_disease[[tissue]] > 1E-6], breaks=c(breakpoints, Inf)), ordered_result=TRUE) %>% 
               dplyr::count(cats, .drop=FALSE) %>% arrange(desc(cats)) %>% mutate(cumfreq = cumsum(n)) %>% arrange(cats))$cumfreq))
    prop_above <- log10(10^n_above %*% (diag(1 / sapply(ensembl_genes, length))))
    colnames(prop_above) <- colnames(n_above)
    
    n_above_at_least_1 <- sapply(names(tissue_x_disease), function(tissue) 
      log10((data.frame(cats = cut(apply(tissue_x_disease[[tissue]], 1, max)[apply(tissue_x_disease[[tissue]], 1, max) > 1E-6], breaks=c(breakpoints, Inf)), ordered_result=TRUE) %>% 
               dplyr::count(cats, .drop=FALSE) %>% arrange(desc(cats)) %>% mutate(cumfreq = cumsum(n)) %>% arrange(cats))$cumfreq))
    
    
    n_traits_above_at_least_1 <- sapply(names(tissue_x_disease), function(tissue) 
      log10((data.frame(cats = cut(apply(tissue_x_disease[[tissue]], 2, max)[apply(tissue_x_disease[[tissue]], 2, max) > 1E-6], breaks=c(breakpoints, Inf)), ordered_result=TRUE) %>% 
               dplyr::count(cats, .drop=FALSE) %>% arrange(desc(cats)) %>% mutate(cumfreq = cumsum(n)) %>% arrange(cats))$cumfreq))
  }
  
  tiss_names <- names(tissue_x_disease)
  ntiss <- length(tissue_x_disease)
  rm(map)
  
}

#### fig 2 extra ####

if(figure_id == 2){
  
  source(file = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/analyses/analysis_GREx_RelativeEffectSize.R")
  if(!exists("gcta_output")){
    load(file = paste0(gcta_directory, "gcta_output_GTEx_allTissues_list_IHW.RData"))
    load(paste0(gcta_directory, "gcta_output_GTEx_allTissues.RData"))
  }
  n_h2s <- sum(sapply(gcta_output, function(i) nrow(i)))
  gcta_output_h2s <- lapply(setNames((names(gcta_output)), (names(gcta_output))), function(tissue) gcta_output[[tissue]]$h2)
  if(!exists("h2_freqs")){
    h2_freqs <- lapply(setNames(rev(names(gcta_output)), rev(names(gcta_output))), function(tissue){
      out <- hist(do.call(c, gcta_output_h2s[1:match(tissue, names(gcta_output_h2s))]), breaks = 0:5/5, plot = F)
      out$counts <- out$counts / n_h2s
      out
    })
  }
  tissues <- rev(names(gcta_output))
  
  load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/GREx_sds_expression")
  load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/GREx_invgamma_estimate")
  
  xr <- c(-2,2)
  logfc_dens <- sapply(deg_eqtl_list, function(del)
    density((del$logFC), na.rm = T, from = xr[1], to = xr[2], n = 256)$y)
  rm(map)
  rm(h2)
  
  
}

if(figure_id == "S1"){
  zcor = as.matrix(read.table("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/zcor_transcriptome_pass1b.tsv"))
}
# cols$category <- cols$Tissue[1:length(unique_trait_categories)+1]
# names(cols$category) <- unique_trait_categories
