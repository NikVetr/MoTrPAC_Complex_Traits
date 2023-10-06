# setwd('/oak/stanford/groups/smontgom/nicolerg/MOTRPAC/PASS_ANALYSIS/PASS1B/RNA/eqtl-coloc')

library(ks)
library(data.table)
library(edgeR)
library(limma)
library(DESeq2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(org.Rn.eg.db)
library(biomaRt)
# library(MotrpacBicQC)
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

#gsutil = '~/google-cloud-sdk/bin/gsutil'
#source('/oak/stanford/groups/smontgom/nicolerg/src/MOTRPAC/motrpac-mawg/pass1b-06/integrative/outliers_covariates/pi1_cook_fx.R')
#source('/oak/stanford/groups/smontgom/nicolerg/src/MOTRPAC/motrpac-mawg/pass1b-06/tools/get_fx.R')

#### load data ####
source(file = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/helper_scripts/load_deg-eqtl_merged_file.R")


#### LDSC output ####

gwas_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/imputed_gwas_hg38_1.1/"
gwas_summary_files <- list.files(gwas_dir)
gwas_summary_files <- gwas_summary_files[-grep(gwas_summary_files, pattern = "README")]
gwas_names <- stringr::str_replace_all(gwas_summary_files, ".txt.gz", "")
coloc_phenotypes <- stringr::str_replace_all(gwas_names, "imputed_", "")
ldsc_output_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/output/"
ldsc_results_paths <- list.files(ldsc_output_dir)
ldsc_log_paths <- paste0(ldsc_output_dir, ldsc_results_paths[grep(ldsc_results_paths, pattern = "log")])
ldsc_results_paths <- paste0(ldsc_output_dir, ldsc_results_paths[grep(ldsc_results_paths, pattern = "results")])
# tissue_abbrev <- unique(MotrpacBicQC::bic_animal_tissue_code$abbreviation)
# tissue_abbrev <- tissue_abbrev[!is.na(tissue_abbrev)]
tissue_abbrev <- MotrpacRatTraining6moData::TISSUE_ABBREV
tissue_abbrev <- setdiff(tissue_abbrev, c("PLASMA", "HYPOTH", "TESTES", "OVARY", "VENACV"))
cluster_names <- paste0(tissue_abbrev, "-sex_homogeneous_changing")
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


if(file.exists("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/ldsc_cluster_results.txt")){
  ldsc_results <- as.data.frame(fread("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/ldsc_cluster_results.txt"))
} else{
  ldsc_results <- as.data.frame(matrix(data = NA, ncol = ncol(fread(ldsc_results_paths[[1]][1])), nrow = length(gwas_names) * length(cluster_names)))
  colnames(ldsc_results) <- colnames(fread(ldsc_results_paths[[1]][1]))
  ldsc_results$cluster <- rep(cluster_names, length(gwas_names))
  ldsc_results$gwas <- c(sapply(gwas_names, function(gwas_name) rep(gwas_name, length(cluster_names))))
  for(i in 1:nrow(ldsc_results)){
    output <- fread(ldsc_results_paths[[ldsc_results$gwas[i]]][grep(pattern = paste0(ldsc_results$cluster[i], ".res"), 
                                                                    ldsc_results_paths[[ldsc_results$gwas[i]]])][1])
    ldsc_results[i,1:ncol(output)] <- output[grep(pattern = ldsc_results$cluster[i], output$Category),]
  }
  fwrite(ldsc_results, "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/ldsc_cluster_results.txt")
}

ldsc_results$gwas <- stringr::str_replace_all(ldsc_results$gwas, "imputed_", "")

ldsc_results$logPVal_enrichment <- log10(ldsc_results$Enrichment_p)
# hist(ldsc_results$logPVal_enrichment)
# hist(ldsc_results$logPVal_enrichment + log10(nrow(ldsc_results)))
# abline(v = log10(0.1), col = 2)
# plot((ldsc_results$Prop._h2 - ldsc_results$Prop._SNPs) / ldsc_results$Prop._h2_std_error, ldsc_results$Enrichment_p)
# plot(ldsc_results$Enrichment / ldsc_results$Enrichment_std_error, ldsc_results$Enrichment_p)
# plot(ldsc_results$Prop._h2 / ldsc_results$Prop._SNPs, ldsc_results$Enrichment)

cols = list(Tissue=MotrpacRatTraining6moData::TISSUE_COLORS[tissue_abbrev], 
            Time=MotrpacRatTraining6moData::GROUP_COLORS[paste0(c(1,2,4,8), "w")],
            Sex=MotrpacRatTraining6moData::SEX_COLORS[c('male','female')])
cols$Tissue[which(is.na(cols$Tissue))] <- '#C0C0C0'
names(cols$Tissue)[which(is.na(names(cols$Tissue)))] <- "t1000-gonads"
cols$cluster <- cols$Tissue[1:length(cluster_names)]
names(cols$cluster) <- cluster_names

#read in conditional analyses
if(file.exists("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/ldsc_cluster_results_conditional.txt")){
  ldsc_results_conditional <- as.data.frame(fread("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/ldsc_cluster_results_conditional.txt"))
} else{
  ldsc_results_conditional_paths <- list.files(ldsc_output_dir)
  ldsc_results_conditional_paths <- paste0(ldsc_output_dir, ldsc_results_conditional_paths[grep(ldsc_results_conditional_paths, pattern = "results")])
  # cluster_names <- paste0("Cluster_", 1:15)
  # tissue_abbrev <- unique(MotrpacBicQC::bic_animal_tissue_code$abbreviation)
  # tissue_abbrev <- tissue_abbrev[!is.na(tissue_abbrev)]
  tissue_abbrev <- MotrpacRatTraining6moData::TISSUE_ABBREV
  tissue_abbrev <- setdiff(tissue_abbrev, c("PLASMA", "HYPOTH", "TESTES", "OVARY", "VENACV"))
  cluster_names <- paste0(tissue_abbrev, "-sex_homogeneous_changing")
  ldsc_results_conditional_paths <- lapply(gwas_names, function(gwas) ldsc_results_conditional_paths[grep(pattern = gwas, ldsc_results_conditional_paths)])
  names(ldsc_results_conditional_paths) <-  gwas_names
  ldsc_results_conditional_paths <- lapply(gwas_names, function(gwas) 
    ldsc_results_conditional_paths[[gwas]][grep(pattern = 
                                                  paste0(sort(cluster_names)[c(1, length(cluster_names))], collapse = "-"), ldsc_results_conditional_paths[[gwas]])])
  names(ldsc_results_conditional_paths) <-  gwas_names
  
  ldsc_results_conditional <- as.data.frame(matrix(data = NA, ncol = ncol(fread(ldsc_results_conditional_paths[[1]][1])), nrow = length(gwas_names) * length(cluster_names)))
  colnames(ldsc_results_conditional) <- colnames(fread(ldsc_results_conditional_paths[[1]][1]))
  ldsc_results_conditional$cluster <- rep(cluster_names, length(gwas_names))
  ldsc_results_conditional$gwas <- c(sapply(gwas_names, function(gwas_name) rep(gwas_name, length(cluster_names))))
  for(i in 0:(length(gwas_names)-1) * length(cluster_names) + 1){
    output <- fread(ldsc_results_conditional_paths[[ldsc_results_conditional$gwas[i]]][1])
    output$Category <- stringr::str_replace(pattern = "L2_0", replacement = "", output$Category)
    output <- output[output$Category %in% paste0("Cluster_", cluster_names),]
    reorder_output <- order(match(stringr::str_replace(output$Category, pattern = "Cluster_", replacement = ""), 
                                  ldsc_results_conditional[i:(i+length(cluster_names)-1),"cluster"]))
    ldsc_results_conditional[i:(i+length(cluster_names)-1),1:ncol(output)] <- output[reorder_output,]
  }
  ldsc_results_conditional$gwas <- stringr::str_remove(ldsc_results_conditional$gwas, "imputed_")
  ldsc_results_conditional$logPVal_enrichment <- log10(ldsc_results_conditional$Enrichment_p)
  fwrite(ldsc_results_conditional, "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/ldsc_cluster_results_conditional.txt")
}


#read in ldsc cts results outputted from *old* command
if(file.exists("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/ldsc_cluster_results_cts_alt.txt")){
  ldsc_results_cts_alt <- as.data.frame(fread("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/ldsc_cluster_results_cts_alt.txt"))
} else{
  ldsc_output_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/output/original_command/"
  ldsc_results_cts_alt_paths <- list.files(ldsc_output_dir)
  ldsc_results_cts_alt_paths <- paste0(ldsc_output_dir, ldsc_results_cts_alt_paths[grep(ldsc_results_cts_alt_paths, pattern = "results")])
  # cluster_names <- paste0("Cluster_", 1:15)
  
  # tissue_abbrev <- unique(MotrpacBicQC::bic_animal_tissue_code$abbreviation)
  # tissue_abbrev <- tissue_abbrev[!is.na(tissue_abbrev)]
  tissue_abbrev <- MotrpacRatTraining6moData::TISSUE_ABBREV
  tissue_abbrev <- setdiff(tissue_abbrev, c("PLASMA", "HYPOTH", "TESTES", "OVARY", "VENACV"))
  cluster_names <- paste0(tissue_abbrev, "-sex_homogeneous_changing")
  
  ldsc_results_cts_alt_paths <- lapply(gwas_names, function(gwas) ldsc_results_cts_alt_paths[grep(pattern = gwas, ldsc_results_cts_alt_paths)])
  names(ldsc_results_cts_alt_paths) <-  gwas_names
  
  ldsc_results_cts_alt <- as.data.frame(matrix(data = NA, 
                                               ncol = ncol(fread(ldsc_results_cts_alt_paths[[1]][[1]][1])), 
                                               nrow = length(gwas_names) * length(cluster_names)))
  colnames(ldsc_results_cts_alt) <- colnames(fread(ldsc_results_cts_alt_paths[[1]][[1]][1]))
  ldsc_results_cts_alt$cluster <- rep(cluster_names, length(gwas_names))
  ldsc_results_cts_alt$gwas <- c(sapply(gwas_names, function(gwas_name) rep(gwas_name, length(cluster_names))))
  for(i in 0:(length(gwas_names)-1) * length(cluster_names) + 1){
    paths <- ldsc_results_cts_alt_paths[[ldsc_results_cts_alt$gwas[i]]]
    output <- do.call(rbind, lapply(setNames(cluster_names, cluster_names), function(cli){
      path <- paths[grep(cli, paths)]
      temp_output <- fread(path)
      temp_output <- temp_output[temp_output$Category == "L2_0"]
      temp_output$Category <- cli
      return(temp_output)
    }))
    reorder_output <- order(match(output$Category, ldsc_results_cts_alt[i:(i+length(cluster_names)-1),"cluster"]))
    ldsc_results_cts_alt[i:(i+length(cluster_names)-1),1:ncol(output)] <- output[reorder_output,]
  }
  ldsc_results_cts_alt$gwas <- stringr::str_remove(ldsc_results_cts_alt$gwas, "imputed_")
  ldsc_results_cts_alt$logPVal_enrichment <- log10(ldsc_results_cts_alt$Enrichment_p)
  fwrite(ldsc_results_cts_alt, "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/ldsc_cluster_results_cts_alt.txt")
}
# all(as.numeric(stringr::str_remove_all(ldsc_results_conditional$Category, pattern = "Cluster_")) == ldsc_results_conditional$cluster)


#read in ldsc cts results outputted from *old* command
if(file.exists("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/ldsc_cluster_results_cts_alt_fullyconditional.txt")){
  ldsc_results_cts_alt_fullyconditional <- as.data.frame(fread("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/ldsc_cluster_results_cts_alt_fullyconditional.txt"))
} else{
  ldsc_output_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/output/original_command/"
  ldsc_results_cts_alt_fullyconditional_paths <- list.files(ldsc_output_dir)
  ldsc_results_cts_alt_fullyconditional_paths <- paste0(ldsc_output_dir, ldsc_results_cts_alt_fullyconditional_paths[grep(ldsc_results_cts_alt_fullyconditional_paths, pattern = "results")])
  # cluster_names <- paste0("Cluster_", 1:15)
  
  # tissue_abbrev <- unique(MotrpacBicQC::bic_animal_tissue_code$abbreviation)
  # tissue_abbrev <- tissue_abbrev[!is.na(tissue_abbrev)]
  tissue_abbrev <- MotrpacRatTraining6moData::TISSUE_ABBREV
  tissue_abbrev <- setdiff(tissue_abbrev, c("PLASMA", "HYPOTH", "TESTES", "OVARY", "VENACV"))
  cluster_names <- paste0(tissue_abbrev, "-sex_homogeneous_changing")
  
  ldsc_results_cts_alt_fullyconditional_paths <- lapply(gwas_names, function(gwas) ldsc_results_cts_alt_fullyconditional_paths[grep(pattern = gwas, ldsc_results_cts_alt_fullyconditional_paths)])
  names(ldsc_results_cts_alt_fullyconditional_paths) <-  gwas_names
  
  ldsc_results_cts_alt_fullyconditional <- as.data.frame(matrix(data = NA, 
                                                                ncol = ncol(fread(ldsc_results_cts_alt_fullyconditional_paths[[1]][[1]][1])), 
                                                                nrow = length(gwas_names) * length(cluster_names)))
  colnames(ldsc_results_cts_alt_fullyconditional) <- colnames(fread(ldsc_results_cts_alt_fullyconditional_paths[[1]][[1]][1]))
  ldsc_results_cts_alt_fullyconditional$cluster <- rep(cluster_names, length(gwas_names))
  ldsc_results_cts_alt_fullyconditional$gwas <- c(sapply(gwas_names, function(gwas_name) rep(gwas_name, length(cluster_names))))
  for(i in 0:(length(gwas_names)-1) * length(cluster_names) + 1){
    paths <- ldsc_results_cts_alt_fullyconditional_paths[[ldsc_results_cts_alt_fullyconditional$gwas[i]]]
    output <- do.call(rbind, lapply(setNames(cluster_names, cluster_names), function(cli){
      path <- paths[grep(cli, paths)]
      path <- path[grep("baseline-plus", path)]
      temp_output <- fread(path)
      temp_output <- temp_output[temp_output$Category == "L2_0"]
      temp_output$Category <- cli
      return(temp_output)
    }))
    reorder_output <- order(match(output$Category, ldsc_results_cts_alt_fullyconditional[i:(i+length(cluster_names)-1),"cluster"]))
    ldsc_results_cts_alt_fullyconditional[i:(i+length(cluster_names)-1),1:ncol(output)] <- output[reorder_output,]
  }
  ldsc_results_cts_alt_fullyconditional$gwas <- stringr::str_remove(ldsc_results_cts_alt_fullyconditional$gwas, "imputed_")
  ldsc_results_cts_alt_fullyconditional$logPVal_enrichment <- log10(ldsc_results_cts_alt_fullyconditional$Enrichment_p)
  fwrite(ldsc_results_cts_alt_fullyconditional, "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/ldsc_cluster_results_cts_alt_fullyconditional.txt")
}


#### plot preprocessing + plotting LDSC output ####
#specify graph parameters
max_point_cex <- 3.5
point_cex_power <- 0.35
n_points_for_legend = 7
buffer_min_and_max = 0.05
# minimum_enrichment_logPval <- min(ldsc_results_sub$logPVal_enrichment)
opacity_insig_points <- 0.2
opacity_sig_points <- 0.8

plot_LDSC_comparison = T
reorder_vertical <- T
use_conditional_model = F
use_cts_alt_model = T
use_enrichment = F
partition_by_category <- T
use_heritability = !use_enrichment
fix_axes = F
fix_axes_bounds_enrichment = c(0,5)
fix_axes_bounds_heritability = c(0,1)

#filter by h2 sigma
total_h2_sigma_thresh <- 7
traits_with_satisfactory_heritaility <- unique(log_files$gwas[log_files$h2 / log_files$h2se > total_h2_sigma_thresh])

# ldsc_results_sub <- ldsc_results[ldsc_results$gwas %in% traits_with_satisfactory_heritaility,]
# if(use_conditional_model){
#   ldsc_results_sub <- ldsc_results_conditional[ldsc_results_conditional$gwas %in% traits_with_satisfactory_heritaility,]
# }

trait_categories <- read.csv("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/gwas_metadata.csv", header = T)
traitname_map <- trait_categories[,c("Tag", "new_Phenotype")]
traitwise_partitions <- trait_categories[,c("Tag", "Category")]
categories <- sort(unique(traitwise_partitions$Category))

#filter by category
trait_categories <- read.csv("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/gwas_metadata.csv", header = T)
traitwise_partitions <- trait_categories[,c("Tag", "Category")]
traits_in_focal_categories <- traitwise_partitions$Tag[traitwise_partitions$Category %in% 
                                                         c("Cardiometabolic", "Aging", "Anthropometric", 
                                                           "Immune", "Psychiatric-neurologic")]
traits_in_focal_categories <- traitwise_partitions$Tag

#BH correction significant
use_IHW_ldsc <- T
if(length(ldsc_results$adj_p) == 0){
  if(use_IHW_ldsc){
    ihw_results_ldsc <- IHW::ihw(p ~ trait, data = data.frame(p = ldsc_results$Enrichment_p, trait = as.factor(ldsc_results$gwas)), alpha = 0.05)
    ldsc_results$adj_p <- ihw_results_ldsc@df$adj_pvalue
    
    ihw_results_ldsc <- IHW::ihw(p ~ trait, data = data.frame(p = ldsc_results_conditional$Enrichment_p, trait = as.factor(ldsc_results_conditional$gwas)), alpha = 0.05)
    ldsc_results_conditional$adj_p <- ihw_results_ldsc@df$adj_pvalue
    
    ihw_results_ldsc <- IHW::ihw(p ~ trait, data = data.frame(p = ldsc_results_cts_alt$Enrichment_p, trait = as.factor(ldsc_results_cts_alt$gwas)), alpha = 0.05)
    ldsc_results_cts_alt$adj_p <- ihw_results_ldsc@df$adj_pvalue
    
    ihw_results_ldsc <- IHW::ihw(p ~ trait, data = data.frame(p = ldsc_results_cts_alt_fullyconditional$Enrichment_p, trait = as.factor(ldsc_results_cts_alt_fullyconditional$gwas)), alpha = 0.05)
    ldsc_results_cts_alt_fullyconditional$adj_p <- ihw_results_ldsc@df$adj_pvalue
  } else {
    ldsc_results$adj_p <- p.adjust(ldsc_results$Enrichment_p, "BH")
  }
}

#get final trait list
traits_to_keep <- intersect(traits_with_satisfactory_heritaility, traits_in_focal_categories)
traits_with_significant_hits <- unique(ldsc_results_cts_alt_fullyconditional$gwas[
  ldsc_results_cts_alt_fullyconditional$adj_p < 0.05 & ldsc_results$Enrichment > 0])
traits_to_keep <- intersect(traits_to_keep, traits_with_significant_hits)

coloc_phenotypes_sub <- traits_to_keep

ldsc_results_sub <- ldsc_results[ldsc_results$gwas %in% traits_to_keep,]

if(use_conditional_model){
  ldsc_results_sub <- ldsc_results_conditional[ldsc_results_conditional$gwas %in% traits_to_keep,]
}

if(use_cts_alt_model){
  ldsc_results_sub <- ldsc_results_cts_alt[ldsc_results_cts_alt$gwas %in% traits_to_keep,]
}

use_fullyconditional_cts_alt_model = T
if(use_fullyconditional_cts_alt_model){
  ldsc_results_sub <- ldsc_results_cts_alt_fullyconditional[ldsc_results_cts_alt_fullyconditional$gwas %in% traits_to_keep,]
}

#get minimum p-val
minimum_enrichment_logPval <- log10(min(ldsc_results_sub$adj_p))


#identify traits with negative heritabilities
bad_boys <- sort(unique(ldsc_results_sub$gwas[ldsc_results_sub$Enrichment < 0 | ldsc_results_sub$Prop._h2 < 0]))
bad_boys_cols <- rep(1, length(gwas_names))
# bad_boys_cols[match(bad_boys, coloc_phenotypes_sub)] <- 2

