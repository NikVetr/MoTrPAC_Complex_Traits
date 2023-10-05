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
hist(ldsc_results$logPVal_enrichment)
hist(ldsc_results$logPVal_enrichment + log10(nrow(ldsc_results)))
abline(v = log10(0.1), col = 2)
plot((ldsc_results$Prop._h2 - ldsc_results$Prop._SNPs) / ldsc_results$Prop._h2_std_error, ldsc_results$Enrichment_p)
plot(ldsc_results$Enrichment / ldsc_results$Enrichment_std_error, ldsc_results$Enrichment_p)
plot(ldsc_results$Prop._h2 / ldsc_results$Prop._SNPs, ldsc_results$Enrichment)

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


#### plotting LDSC output ####
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

plot(density(log_files$h2))
segments(x0 = c(log_files$h2[log_files$gwas %in% bad_boys]), x1 = c(log_files$h2[log_files$gwas %in% bad_boys]), y0 = 0, y1 = 5, col = "red", lwd = 0.4)
hist(log_files$h2[log_files$gwas %in% bad_boys])
hist(sapply(log_files$h2[log_files$gwas %in% bad_boys], function(h) mean(h > log_files$h2)), cex.main = 1)
par(mar = c(4,4,4,4))
plot(log_files$h2[log_files$gwas %in% bad_boys], log_files$h2se[log_files$gwas %in% bad_boys], xlim = c(-0.003, 0.006), ylim = c(-0.003, 0.006), 
     xlab = "total heritability", ylab = "total heritability SE")
hist(log_files$h2[log_files$gwas %in% bad_boys] / log_files$h2se[log_files$gwas %in% bad_boys])
hist(log_files$h2 / log_files$h2se, breaks = seq(min(log_files$h2 / log_files$h2se), max(log_files$h2 / log_files$h2se), length.out = 100))
abline(a = 0, b = 1, xpd = F)
hist(log_files$h2, breaks = seq(min(log_files$h2), max(log_files$h2), length.out = 50))
hist(log_files$chi2, breaks = seq(min(log_files$chi2), max(log_files$chi2), length.out = 50))
hist(log_files$chi2[log_files$gwas %in% bad_boys])
log_files$gwas[log_files$h2 > 1]



change_names <- F
change_names_in_plot = T
nicole_mods = T
arnold_mods = F
if(plot_LDSC_comparison){
  
  
  # grDevices::cairo_pdf(filename = paste0("~/Documents/Documents - nikolai/DExEQTL/LDSC_comparison_heritabilities", 
  #                                        ifelse(use_enrichment, "_enrichment", "_heritability"), 
  #                                        ifelse(fix_axes, "_fixedAxes", ""), 
  #                                        ifelse(use_fullyconditional_cts_alt_model, "_fullyConditional", ""), 
  #                                        # ifelse(use_conditional_model, "_conditionalModel", "_unconditionalModel"), 
  #                                        ifelse(use_cts_alt_model, "_cell-type-specific-Model", ""), ".pdf"), 
  #                      width = 1500 / 72, height = 2000 / 72 * 18 / 17 / 114 * length(coloc_phenotypes_sub) * 2, family="Arial Unicode MS")
  
  grDevices::cairo_pdf(filename = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/figures//LDSC_output.pdf"), 
                       width = 1500 / 72, height = 2000 / 72 * 18 / 17 / 114 * length(coloc_phenotypes_sub) * 2, family="Arial Unicode MS")
  
  par(mfrow = c(2, 1), mar = c(3,3,2,3), xpd = NA)
  
  for(type_of_plot in 1:2){
    
    use_enrichment <- c(T, F)[type_of_plot]
    
    max_horiz_axis <- 2.05
    
    
    #first get metainfo
    if(reorder_vertical){
      
      if(use_enrichment){
        order_traits <- order(sapply(coloc_phenotypes_sub, function(coloc_phenotype) mean(ldsc_results_sub[ldsc_results_sub$gwas == coloc_phenotype,]$Enrichment)), decreasing = T)
      } else {
        order_traits <- order(sapply(coloc_phenotypes_sub, function(coloc_phenotype) mean(ldsc_results_sub[ldsc_results_sub$gwas == coloc_phenotype,]$Prop._h2)), decreasing = T)
      }
      
    } else {
      order_traits <- 1:length(coloc_phenotypes_sub)
    }
    
    if(partition_by_category){
      order_traits <- order_traits[order(traitwise_partitions[,2][match(coloc_phenotypes_sub[order_traits], traitwise_partitions[,1])], sort(order_traits))]
      trait_category_counts <- as.data.table(table(traitwise_partitions[,2][match(coloc_phenotypes_sub[order_traits], traitwise_partitions[,1])] ))
      trait_category_counts$cN <- cumsum(trait_category_counts$N)
    }
    
    if(use_enrichment){
      if(fix_axes){
        logprob_range <- fix_axes_bounds_enrichment
      } else {
        logprob_range <- range(ldsc_results_sub$Enrichment)
      }
    } else {
      if(fix_axes){
        logprob_range <- fix_axes_bounds_heritability
      } else {
        logprob_range <- range(ldsc_results_sub$Prop._h2)
      }
    }
    
    logprob_range <- c(logprob_range[1] - diff(logprob_range) * buffer_min_and_max / 2, logprob_range[2] + diff(logprob_range) * buffer_min_and_max / 2)
    logprob_range <- c(0,ifelse(use_enrichment, 
                                max(ldsc_results_sub$Enrichment[ldsc_results_sub$adj_p < 0.05]) * 1.025, 
                                max(ldsc_results_sub$Prop._h2[ldsc_results_sub$adj_p < 0.05])) * 1.025)
    max_prob <- diff(logprob_range)
    
    nameloc <- 0.3 + ifelse(change_names | change_names_in_plot, 0, 0.3)
    
    plot(1,1,xlim = c(0,2.1), ylim = c(0,10), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
    
    ylocs <- seq(10, 0, length.out = length(coloc_phenotypes_sub))
    
    #put boxes around categories, if we're doing that  
    if(partition_by_category){
      segments(x0 = -0.1, x1 = -0.1, y0 = 0+diff(ylocs)[1]/2, y1 = 10-diff(ylocs)[1]/2, col = "grey35")
      for(traitcat in 1:(nrow(trait_category_counts)-1)){
        segments(y0 = ylocs[trait_category_counts$cN[traitcat]] + diff(ylocs)[1]/2, y1 = ylocs[trait_category_counts$cN[traitcat]] + diff(ylocs)[1]/2, x0 = -0.1, x1 = nameloc, col = "grey35")
      }
      
      catylocs <- (ylocs[trait_category_counts$cN] + c(ylocs[1], ylocs[trait_category_counts$cN][-length(trait_category_counts$cN)])) / 2 - diff(ylocs)[1]/1.5
      catnames <- as.vector(trait_category_counts$V1)
      catnames <- sapply(catnames, function(x) strsplit(x, split = c("-"))[[1]][1])
      catnames <- sapply(catnames, function(x) strsplit(x, split = c(" "))[[1]][1])
      
      #manually fudge category locations
      catylocs[catnames == "Hair"] <- catylocs[catnames == "Hair"] + diff(ylocs)[1]/1.5
      catylocs[catnames == "Aging"] <- catylocs[catnames == "Aging"] - diff(ylocs)[1]/3
      catylocs[catnames == "Endocrine"] <- catylocs[catnames == "Endocrine"] - diff(ylocs)[1]/10
      catylocs[catnames == "Allergy"] <- catylocs[catnames == "Allergy"] - diff(ylocs)[1]*1.2
      catylocs[catnames == "Anthropometric"] <- catylocs[catnames == "Anthropometric"] - diff(ylocs)[1]*1.2
      catylocs[catnames == "Cardiometabolic"] <- catylocs[catnames == "Cardiometabolic"] - diff(ylocs)[1]*2
      for(catname in 1:length(catnames)){
        text(x = -0.1, y = catylocs[catname] + ifelse(catnames[catname] == "Endocrine", 0.1, 0), labels = catnames[catname], pos = 2, srt = 90, col = "grey35")
      }
      
    }
    
    #plot names of traits
    if(change_names_in_plot){
      coloc_phenotypes_sub_newnames <- traitname_map[match(coloc_phenotypes_sub, traitname_map[,1]),2]
    } else {
      coloc_phenotypes_sub_newnames <- coloc_phenotypes_sub
    }
    
    if(nicole_mods){
      # text(paste0(1:length(coloc_phenotypes_sub_newnames), ": ", 
      #             coloc_phenotypes_sub_newnames)[order_traits], 
      #      col = bad_boys_cols[order_traits], x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
      text(paste0(coloc_phenotypes_sub_newnames)[order_traits], 
           col = bad_boys_cols[order_traits], x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
    } else {
      text(paste0(1:length(coloc_phenotypes_sub_newnames), ": ", 
                  coloc_phenotypes_sub_newnames)[order_traits], 
           col = pheno_cols, x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
    }
    
    #vert axis label
    text(labels = "Phenotypes", x = nameloc, y = 10.45, pos = 2, cex = 2.5, xpd = NA, family="Courier", col = "grey25")
    #horiz axis label
    text(labels = paste0("Heritability ", ifelse(use_enrichment, "Enrichment", "Proportion"), " Across Tissues"), x = 1.2, y = -0.5, cex = 2.25, pos = 1, xpd = NA)
    
    
    #guiding lines for traits
    if(nicole_mods){
      segments(x0 = nameloc, x1 = max_horiz_axis, y0 = ylocs, y1 = ylocs, col = "grey80", lty = 3, lwd = 0.5)
    }
    
    #axes
    segments(x0 = nameloc, x1 = nameloc, y0 = 10.675, y1 = -0.1, lwd = 2)
    segments(x0 = nameloc, x1 = max_horiz_axis, y0 = -0.1, y1 = -0.1, lwd = 2)
    
    #horiz axis ticks and nums
    # if(use_geometric_mean){
    #   
    #   probs <- ceiling(logprob_range[1]):ceiling(logprob_range[2])
    #   segments(x0 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
    #            x1 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
    #            y0 = -0.1, y1 = -0.15, lwd = 2)
    #   text(x = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
    #        y = -0.135, pos = 1, cex = 1.25,
    #        labels =  c(sapply(probs, function(num) as.expression(bquote(10^.(num))))))
    #   
    # } else {
    #   
    #   segments(x0 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
    #            x1 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
    #            y0 = -0.1, y1 = -0.15, lwd = 2)
    #   text(x = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
    #        y = -0.135, pos = 1, cex = 1.25,
    #        labels =  c(sapply(seq(ceiling(logprob_range[1]), ceiling(logprob_range[2]), by = 1), function(num) as.expression(bquote(10^.(num))))))
    # }
    
    if(use_enrichment){
      probs <- seq(ceiling(logprob_range[1]), floor(logprob_range[2]), by = ifelse(fix_axes, 0.5, 1))  
    } else {
      probs <- seq(ceiling(logprob_range[1]), (logprob_range[2]), by = ifelse(fix_axes, 0.1, 0.1))
    }
    
    segments(x0 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc,
             x1 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc,
             y0 = -0.1, y1 = -0.15, lwd = 2)
    text(x = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc,
         y = -0.135, pos = 1, cex = 1.25,
         labels =  probs)
    
    
    # #let's get minor tick marks in there too
    # if(use_geometric_mean){
    #   minticks <- log10(as.vector(sapply(1:(length(probs)-1), function(p) seq(10^probs[p], 10^probs[p+1], length.out = 10))))
    #   minticks <- (minticks - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc
    #   minticks <- minticks[minticks < max_horiz_axis]
    #   segments(x0 = minticks, x1 =  minticks, y0 = -0.1, y1 = -0.125, lwd = 1)
    #   #hmm being weird
    # } else {
    #   minticks <- log10(as.vector(t(t(1:9)) %*% t(10^seq(0, max_prob, by = 10^round(log(max_prob))))))
    #   minticks <- minticks / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025
    #   minticks <- minticks[minticks < max_horiz_axis]
    #   segments(x0 = minticks, x1 =  minticks, y0 = -0.1, y1 = -0.125, lwd = 1)
    # }
    
    
    #dashed lines for ticks
    segments(x0 = (probs[-1] - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
             x1 = (probs[-1] - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y0 = 10, y1 = -0.075, lwd = 2, lty = 2, col = "grey75")
    segments(x0 = (1 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
             x1 = (1 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
             y0 = 10, y1 = -0.075, lwd = 3, lty = 1, col = "grey75")
    
    
    #helpful hint line about direction
    
    # if(use_geometric_mean){
    #   segments(x0 = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
    #            x1 = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y0 = 10, y1 = 10.4, lwd = 2, lty = 2, col = "grey75")
    #   text(labels = "non-DE Genes", srt = 90, x = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = 10.4, pos = 2, col = "grey75", cex = 0.75)
    #   text(labels = "DE Genes", srt = 270, x = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = 10.4, pos = 4, col = "grey75", cex = 0.75)
    #   text(labels = "Stronger Coloc\nSignal In", srt = 0, x = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = 10.3875, pos = 3, col = "grey75", cex = 0.75)
    # } else {
    #   segments(x0 = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
    #            x1 = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y0 = 10, y1 = 10.4, lwd = 2, lty = 2, col = "grey75")
    #   text(labels = "eQTL > DEG", srt = 90, x = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = 10.4, pos = 2, col = "grey75", cex = 0.75)
    #   text(labels = "DEG > eQTL", srt = 270, x = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = 10.4, pos = 4, col = "grey75", cex = 0.75)
    #   text(labels = "Stronger Coloc\nSignal In", srt = 0, x = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = 10.3875, pos = 3, col = "grey75", cex = 0.75)
    # }
    
    #plot legend for colors etc.
    lower_legend_by <- 1.6
    rect(xleft = 1.97, ybottom = 8.15 - lower_legend_by - 0.05*n_points_for_legend, 
         ytop = 10.1 - lower_legend_by, xright = 2.2, border = NA, col = "white")
    # text(labels = sapply(names(deg_eqtl_list), function(ts) stringr::str_to_title(paste0(strsplit(ts, "-")[[1]][-1], collapse = " "))),
    #      y = seq(9.95,8.75,length.out = length(names(deg_eqtl_list))), x = 2, pos = 4)
    # points(x = rep(1.9925, length(names(deg_eqtl_list))), y = seq(9.955,8.755,length.out = length(names(deg_eqtl_list))), col = cols$Tissue, pch = 15, cex = 1.75)
    # text(labels = stringr::str_replace_all(cluster_names, "_", " "),
    #      y = seq(9.95,8.75,length.out = length(cluster_names)), x = 2, pos = 4)
    text(labels = sapply(cluster_names, function(cli) strsplit(cli, "-")[[1]][1]),
         y = seq(9.95,8.75 - lower_legend_by,length.out = length(cluster_names)) - 0.02, x = 2, pos = 4)
    points(x = rep(1.9925, length(cluster_names)), y = seq(9.955,8.755 - lower_legend_by,length.out = length(cluster_names)), col = cols$cluster, pch = 19, cex = 1.75)
    
    #plot legend for points
    lower_legend_by <- lower_legend_by + 0.3
    pt_loc_expander <- 8
    point_legend_cexes <- seq(from = 0.4, to = max_point_cex, length.out = n_points_for_legend)
    points_legend_pchs <- rep(19, n_points_for_legend)
    point_legend_cexes_log10_pvals <- round((point_legend_cexes / max_point_cex)^(1/point_cex_power) * minimum_enrichment_logPval, 2)
    # points_legend_pchs[point_legend_cexes_log10_pvals < (log10(0.05) - log10(nrow(ldsc_results_sub)))] <- 18
    points_legend_pchs[point_legend_cexes_log10_pvals < (log10(0.05))] <- 18
    
    points(y = 8.55 - lower_legend_by - cumsum(point_legend_cexes + pt_loc_expander) / sum(point_legend_cexes) * 0.05*n_points_for_legend, 
           x = rep(1.9925, n_points_for_legend), cex = point_legend_cexes, col = "grey50", pch = points_legend_pchs)
    text(labels = latex2exp::TeX("$log_{10}$(enrichment p-val)"), y = 8.6 - lower_legend_by, x = 1.975, cex = 1.1, pos = 4,  font = 2)
    # text(labels = paste0("0.05 FDR @ ", round(log10(0.05) - log10(nrow(ldsc_results_sub)),2), "\n           (Bonferroni)"), y = 8.45 , x = 2.1, cex = 0.8, pos = 4,  font = 2)
    text(labels = latex2exp::TeX(paste0("$\\alpha = 0.05$, IHW")), 
         y = 8.35 - lower_legend_by, x = 2.1, cex = 0.8, pos = 4,  font = 2)
    # text(labels = latex2exp::TeX(paste0("(Benjamini-Hochberg)")), 
    #      y = 8.35 - lower_legend_by, x = 2.07, cex = 0.8, pos = 4,  font = 2)
    # text(labels = paste0("< ", round(log10(0.05) - log10(nrow(ldsc_results_sub)),2), "\n> ", round(log10(0.05) - log10(nrow(ldsc_results_sub)),2)), 
    #      y = 8.25 , x = 2.175, cex = 1, pos = 4,  font = 2)
    text(labels = paste0("< ", round(log10(0.05),2), "\n> ", round(log10(0.05),2)), 
         y = 8.25 - lower_legend_by - 0.1 - c(0.15), x = 2.145, cex = 1, pos = 4,  font = 2)
    points(pch = c(19,18), y = c(8.41, 8.2) - lower_legend_by - 0.1 - 0.15, x = rep(2.145,2), cex = c(1.25, 1.75), 
           col = sapply(c(opacity_insig_points, opacity_sig_points), function(opcty) adjustcolor("grey50", alpha.f = opcty)))
    
    
    text(labels = point_legend_cexes_log10_pvals,y = 8.545 - lower_legend_by - cumsum(point_legend_cexes + pt_loc_expander) / sum(point_legend_cexes) * 0.05*n_points_for_legend, 
         x = rep(1.99, n_points_for_legend) + point_legend_cexes / 200, cex = 1, pos = 4, pch = 19)
    
    
    #figure out cex params
    for(cluster in cluster_names){
      
      # horizontal lines for colocalizing traits
      trait_locs <- ylocs
      if(use_enrichment){
        trait_probs <- ldsc_results_sub$Enrichment[ldsc_results_sub$cluster == cluster]
        #for bonferroni
        # trait_logPvals <- (ldsc_results_sub$logPVal_enrichment[ldsc_results_sub$cluster == cluster])[order_traits]
        #for BH adjustment
        trait_logPvals <- log10((ldsc_results_sub$adj_p[ldsc_results_sub$cluster == cluster])[order_traits])
      } else {
        trait_probs <- ldsc_results_sub$Prop._h2[ldsc_results_sub$cluster == cluster]
        trait_logPvals <- (ldsc_results_sub$logPVal_enrichment[ldsc_results_sub$cluster == cluster])[order_traits]
      }
      trait_probs <- trait_probs[order_traits]
      
      point_cex <- (-trait_logPvals / -minimum_enrichment_logPval)^point_cex_power
      point_cex <- point_cex * max_point_cex
      
      good_points <- trait_probs >= logprob_range[1] & trait_probs <= logprob_range[2]
      
      pchs <- rep(19, length(trait_probs))
      opacities <- rep(opacity_insig_points, length(trait_probs))
      
      #mark "significant" points
      # pchs[trait_logPvals < (log10(0.05) - log10(nrow(ldsc_results_sub)))] <- 18
      # opacities[trait_logPvals < (log10(0.05) - log10(nrow(ldsc_results_sub)))] <- opacity_sig_points
      pchs[trait_logPvals < log10(0.05)] <- 23
      opacities[trait_logPvals < log10(0.05)] <- opacity_sig_points
      
      point_cols <- sapply(opacities, function(opcty) grDevices::adjustcolor(cols$Tissue[match(cluster, cluster_names)], alpha.f = opcty))
      point_bgs <- point_cols
      point_cols[trait_logPvals < log10(0.05)] <- 1
      
      
      points(x = ((trait_probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc)[good_points], y = (trait_locs)[good_points], 
             pch = pchs[good_points], bg = point_bgs[good_points], col = point_cols[good_points], cex = point_cex[good_points])
      
      
    }
    
    
    
    if(arnold_mods){
      addImg(obj = arnold_bayes, x = 1.9, y = 1.6775, width = 0.5)
    }
    
  }
  
  dev.off()
}

#quick comparison of conditional vs unconditional models
par(mfrow = c(1,2), mar = c(4,4,4,2))
plot(ldsc_results_conditional$Enrichment,  ldsc_results$Enrichment, main = "Heritability Enrichment",
     xlab = "Conditional Model Enrichments", ylab = "Unconditional Model Enrichments", pch = 19, col = adjustcolor(1, 0.5))
abline(0,1, col = "red", lty = 2)
plot(ldsc_results_conditional$Prop._h2,  ldsc_results$Prop._h2, main = "Proportion Heritability", xpd = NA,
     xlab = latex2exp::TeX("Conditional Model Prop. h^2"), ylab =  latex2exp::TeX("Uncnditional Model Prop. h^2"), pch = 19, col = adjustcolor(1, 0.5))
abline(0,1, col = "red", lty = 2)

plot(ldsc_results$Prop._h2,  ldsc_results_cts_alt$Prop._h2, main = "Proportion Heritability", xpd = NA,
     xlab = latex2exp::TeX("Conditional Model Prop. h^2"), ylab =  latex2exp::TeX("Uncnditional Model Prop. h^2"), pch = 19, col = adjustcolor(1, 0.5))
abline(0,1, col = "red", lty = 2)

plot(ldsc_results$Enrichment,  ldsc_results_cts_alt$Enrichment, main = "Enrichment", xpd = NA,
     xlab = latex2exp::TeX("Conditional Model Prop. h^2"), ylab =  latex2exp::TeX("Uncnditional Model Prop. h^2"), pch = 19, col = adjustcolor(1, 0.5))
abline(0,1, col = "red", lty = 2)

ldsc_results_conditional
ldsc_results_cts_alt
ldsc_results_cts_alt_fullyconditional

#compare fully conditional model to partly conditional model
sub_1 <- ldsc_results_cts_alt_fullyconditional[ldsc_results_cts_alt_fullyconditional$gwas %in% traits_to_keep,]
sub_2 <- ldsc_results_cts_alt[ldsc_results_cts_alt$gwas %in% traits_to_keep,]
tissues <- tissue_order[tissue_order %in% unique(gsub(x = sub_1$cluster, "-sex_homogeneous_changing", ""))]

par(mfrow = c(1,2), mar = c(4,4,4,2))
plot(sub_1$Enrichment,  sub_2$Enrichment, main = latex2exp::TeX("$h^2_{SNP}$ Enrichment"),
     xlab = "Enrichment (Conditional on Tissue Annotation)", ylab = "Enrichment (Unconditional on Tissue Annotation)", pch = 19, 
     col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[gsub(x = sub_1$cluster, "-sex_homogeneous_changing", "")], 0.5))
legend("topleft", pch = 19, col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues], 0.5), legend = tissues, cex = 0.5, bty="n")
legend("bottomright", lty = 2, col = "red", legend = "1-to-1 line", cex = 0.5, bty="n")
abline(0,1, col = "red", lty = 2)

plot(sub_1$Prop._h2,  sub_2$Prop._h2, main = latex2exp::TeX("Proportion $h^2_{SNP}$"), xpd = NA,
     xlab = latex2exp::TeX("Proportion $h^2_{SNP}$ (Conditional on Tissue Annotation)"), 
     ylab =  latex2exp::TeX("Proportion $h^2_{SNP}$ (Unconditional on Tissue Annotation)"), pch = 19, 
     col =  adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[gsub(x = sub_1$cluster, "-sex_homogeneous_changing", "")], 0.5))
legend("topleft", pch = 19, col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues], 0.5), legend = tissues, cex = 0.5, bty="n")
legend("bottomright", lty = 2, col = "red", legend = "1-to-1 line", cex = 0.5, bty="n")
abline(0,1, col = "red", lty = 2)


#### LDSC figure for paper, fig 3 ####


change_names <- F
change_names_in_plot = T
nicole_mods = T
arnold_mods = F

grDevices::cairo_pdf(filename = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/figures/fig_3_LDSC_output.pdf"), 
                     width = 1500 / 72, height = 2000 / 72 * 18 / 17 / 114 * length(coloc_phenotypes_sub), family="Arial Unicode MS")

par(mfrow = c(1, 1), mar = c(3,3,2,3), xpd = NA)

for(type_of_plot in 1:1){
  
  use_enrichment <- c(T, F)[type_of_plot]
  
  max_horiz_axis <- 2.05
  
  
  #first get metainfo
  if(reorder_vertical){
    
    if(use_enrichment){
      order_traits <- order(sapply(coloc_phenotypes_sub, function(coloc_phenotype) mean(ldsc_results_sub[ldsc_results_sub$gwas == coloc_phenotype,]$Enrichment)), decreasing = T)
    } else {
      order_traits <- order(sapply(coloc_phenotypes_sub, function(coloc_phenotype) mean(ldsc_results_sub[ldsc_results_sub$gwas == coloc_phenotype,]$Prop._h2)), decreasing = T)
    }
    
  } else {
    order_traits <- 1:length(coloc_phenotypes_sub)
  }
  
  if(partition_by_category){
    order_traits <- order_traits[order(traitwise_partitions[,2][match(coloc_phenotypes_sub[order_traits], traitwise_partitions[,1])], sort(order_traits))]
    trait_category_counts <- as.data.table(table(traitwise_partitions[,2][match(coloc_phenotypes_sub[order_traits], traitwise_partitions[,1])] ))
    trait_category_counts$cN <- cumsum(trait_category_counts$N)
  }
  
  if(use_enrichment){
    if(fix_axes){
      logprob_range <- fix_axes_bounds_enrichment
    } else {
      logprob_range <- range(ldsc_results_sub$Enrichment)
    }
  } else {
    if(fix_axes){
      logprob_range <- fix_axes_bounds_heritability
    } else {
      logprob_range <- range(ldsc_results_sub$Prop._h2)
    }
  }
  
  logprob_range <- c(logprob_range[1] - diff(logprob_range) * buffer_min_and_max / 2, logprob_range[2] + diff(logprob_range) * buffer_min_and_max / 2)
  logprob_range <- c(0,ifelse(use_enrichment, 
                              max(ldsc_results_sub$Enrichment[ldsc_results_sub$adj_p < 0.05]) * 1.025, 
                              max(ldsc_results_sub$Prop._h2[ldsc_results_sub$adj_p < 0.05])) * 1.025)
  max_prob <- diff(logprob_range)
  
  nameloc <- 0.3 + ifelse(change_names | change_names_in_plot, 0, 0.3)
  
  plot(1,1,xlim = c(0,2.1), ylim = c(0,10), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
  
  ylocs <- seq(10, 0, length.out = length(coloc_phenotypes_sub))
  
  #put boxes around categories, if we're doing that  
  if(partition_by_category){
    segments(x0 = -0.1, x1 = -0.1, y0 = 0+diff(ylocs)[1]/2, y1 = 10-diff(ylocs)[1]/2, col = "grey35")
    for(traitcat in 1:(nrow(trait_category_counts)-1)){
      segments(y0 = ylocs[trait_category_counts$cN[traitcat]] + diff(ylocs)[1]/2, y1 = ylocs[trait_category_counts$cN[traitcat]] + diff(ylocs)[1]/2, x0 = -0.1, x1 = nameloc, col = "grey35")
    }
    
    catylocs <- (ylocs[trait_category_counts$cN] + c(ylocs[1], ylocs[trait_category_counts$cN][-length(trait_category_counts$cN)])) / 2 - diff(ylocs)[1]/1.5
    catnames <- as.vector(trait_category_counts$V1)
    catnames <- sapply(catnames, function(x) strsplit(x, split = c("-"))[[1]][1])
    catnames <- sapply(catnames, function(x) strsplit(x, split = c(" "))[[1]][1])
    
    #manually fudge category locations
    catylocs[catnames == "Hair"] <- catylocs[catnames == "Hair"] + diff(ylocs)[1]/1.5
    catylocs[catnames == "Aging"] <- catylocs[catnames == "Aging"] - diff(ylocs)[1]/3
    catylocs[catnames == "Endocrine"] <- catylocs[catnames == "Endocrine"] - diff(ylocs)[1]/10
    catylocs[catnames == "Allergy"] <- catylocs[catnames == "Allergy"] - diff(ylocs)[1]*1.2
    catylocs[catnames == "Anthropometric"] <- catylocs[catnames == "Anthropometric"] - diff(ylocs)[1]*1.2
    catylocs[catnames == "Cardiometabolic"] <- catylocs[catnames == "Cardiometabolic"] - diff(ylocs)[1]*2
    for(catname in 1:length(catnames)){
      text(x = -0.1, y = catylocs[catname] + ifelse(catnames[catname] == "Endocrine", 0.1, 0), labels = catnames[catname], pos = 2, srt = 90, col = "grey35")
    }
    
  }
  
  #plot names of traits
  if(change_names_in_plot){
    coloc_phenotypes_sub_newnames <- traitname_map[match(coloc_phenotypes_sub, traitname_map[,1]),2]
  } else {
    coloc_phenotypes_sub_newnames <- coloc_phenotypes_sub
  }
  
  if(nicole_mods){
    # text(paste0(1:length(coloc_phenotypes_sub_newnames), ": ", 
    #             coloc_phenotypes_sub_newnames)[order_traits], 
    #      col = bad_boys_cols[order_traits], x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
    text(paste0(coloc_phenotypes_sub_newnames)[order_traits], 
         col = bad_boys_cols[order_traits], x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
  } else {
    text(paste0(1:length(coloc_phenotypes_sub_newnames), ": ", 
                coloc_phenotypes_sub_newnames)[order_traits], 
         col = pheno_cols, x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
  }
  
  #vert axis label
  text(labels = "Phenotypes", x = nameloc, y = 10.45, pos = 2, cex = 2.5, xpd = NA, family="Courier", col = "grey25")
  #horiz axis label
  text(labels = paste0("Heritability ", ifelse(use_enrichment, "Enrichment", "Proportion"), " Across Tissues"), x = 1.2, y = -0.5, cex = 2.25, pos = 1, xpd = NA)
  
  
  #guiding lines for traits
  if(nicole_mods){
    segments(x0 = nameloc, x1 = max_horiz_axis, y0 = ylocs, y1 = ylocs, col = "grey80", lty = 3, lwd = 0.5)
  }
  
  #axes
  segments(x0 = nameloc, x1 = nameloc, y0 = 10.675, y1 = -0.1, lwd = 2)
  segments(x0 = nameloc, x1 = max_horiz_axis, y0 = -0.1, y1 = -0.1, lwd = 2)
  
  #horiz axis ticks and nums
  # if(use_geometric_mean){
  #   
  #   probs <- ceiling(logprob_range[1]):ceiling(logprob_range[2])
  #   segments(x0 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
  #            x1 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
  #            y0 = -0.1, y1 = -0.15, lwd = 2)
  #   text(x = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
  #        y = -0.135, pos = 1, cex = 1.25,
  #        labels =  c(sapply(probs, function(num) as.expression(bquote(10^.(num))))))
  #   
  # } else {
  #   
  #   segments(x0 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
  #            x1 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
  #            y0 = -0.1, y1 = -0.15, lwd = 2)
  #   text(x = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
  #        y = -0.135, pos = 1, cex = 1.25,
  #        labels =  c(sapply(seq(ceiling(logprob_range[1]), ceiling(logprob_range[2]), by = 1), function(num) as.expression(bquote(10^.(num))))))
  # }
  
  if(use_enrichment){
    probs <- seq(ceiling(logprob_range[1]), floor(logprob_range[2]), by = ifelse(fix_axes, 0.5, 1))  
  } else {
    probs <- seq(ceiling(logprob_range[1]), (logprob_range[2]), by = ifelse(fix_axes, 0.1, 0.1))
  }
  
  segments(x0 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc,
           x1 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc,
           y0 = -0.1, y1 = -0.15, lwd = 2)
  text(x = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc,
       y = -0.135, pos = 1, cex = 1.25,
       labels =  probs)
  
  
  # #let's get minor tick marks in there too
  # if(use_geometric_mean){
  #   minticks <- log10(as.vector(sapply(1:(length(probs)-1), function(p) seq(10^probs[p], 10^probs[p+1], length.out = 10))))
  #   minticks <- (minticks - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc
  #   minticks <- minticks[minticks < max_horiz_axis]
  #   segments(x0 = minticks, x1 =  minticks, y0 = -0.1, y1 = -0.125, lwd = 1)
  #   #hmm being weird
  # } else {
  #   minticks <- log10(as.vector(t(t(1:9)) %*% t(10^seq(0, max_prob, by = 10^round(log(max_prob))))))
  #   minticks <- minticks / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025
  #   minticks <- minticks[minticks < max_horiz_axis]
  #   segments(x0 = minticks, x1 =  minticks, y0 = -0.1, y1 = -0.125, lwd = 1)
  # }
  
  
  #dashed lines for ticks
  segments(x0 = (probs[-1] - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
           x1 = (probs[-1] - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y0 = 10, y1 = -0.075, lwd = 2, lty = 2, col = "grey75")
  segments(x0 = (1 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
           x1 = (1 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
           y0 = 10, y1 = -0.075, lwd = 3, lty = 1, col = "grey75")
  
  
  #helpful hint line about direction
  
  # if(use_geometric_mean){
  #   segments(x0 = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
  #            x1 = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y0 = 10, y1 = 10.4, lwd = 2, lty = 2, col = "grey75")
  #   text(labels = "non-DE Genes", srt = 90, x = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = 10.4, pos = 2, col = "grey75", cex = 0.75)
  #   text(labels = "DE Genes", srt = 270, x = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = 10.4, pos = 4, col = "grey75", cex = 0.75)
  #   text(labels = "Stronger Coloc\nSignal In", srt = 0, x = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = 10.3875, pos = 3, col = "grey75", cex = 0.75)
  # } else {
  #   segments(x0 = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
  #            x1 = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y0 = 10, y1 = 10.4, lwd = 2, lty = 2, col = "grey75")
  #   text(labels = "eQTL > DEG", srt = 90, x = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = 10.4, pos = 2, col = "grey75", cex = 0.75)
  #   text(labels = "DEG > eQTL", srt = 270, x = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = 10.4, pos = 4, col = "grey75", cex = 0.75)
  #   text(labels = "Stronger Coloc\nSignal In", srt = 0, x = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = 10.3875, pos = 3, col = "grey75", cex = 0.75)
  # }
  
  #plot legend for colors etc.
  lower_legend_by <- 1.6
  rect(xleft = 1.97, ybottom = 8.15 - lower_legend_by - 0.05*n_points_for_legend, 
       ytop = 10.1 - lower_legend_by, xright = 2.2, border = NA, col = "white")
  # text(labels = sapply(names(deg_eqtl_list), function(ts) stringr::str_to_title(paste0(strsplit(ts, "-")[[1]][-1], collapse = " "))),
  #      y = seq(9.95,8.75,length.out = length(names(deg_eqtl_list))), x = 2, pos = 4)
  # points(x = rep(1.9925, length(names(deg_eqtl_list))), y = seq(9.955,8.755,length.out = length(names(deg_eqtl_list))), col = cols$Tissue, pch = 15, cex = 1.75)
  # text(labels = stringr::str_replace_all(cluster_names, "_", " "),
  #      y = seq(9.95,8.75,length.out = length(cluster_names)), x = 2, pos = 4)
  text(labels = sapply(cluster_names, function(cli) strsplit(cli, "-")[[1]][1]),
       y = seq(9.95,8.75 - lower_legend_by,length.out = length(cluster_names)) - 0.02, x = 2, pos = 4)
  points(x = rep(1.9925, length(cluster_names)), y = seq(9.955,8.755 - lower_legend_by,length.out = length(cluster_names)), col = cols$cluster, pch = 19, cex = 1.75)
  
  #plot legend for points
  lower_legend_by <- lower_legend_by + 0.3
  pt_loc_expander <- 8
  point_legend_cexes <- seq(from = 0.4, to = max_point_cex, length.out = n_points_for_legend)
  points_legend_pchs <- rep(19, n_points_for_legend)
  point_legend_cexes_log10_pvals <- round((point_legend_cexes / max_point_cex)^(1/point_cex_power) * minimum_enrichment_logPval, 2)
  # points_legend_pchs[point_legend_cexes_log10_pvals < (log10(0.05) - log10(nrow(ldsc_results_sub)))] <- 18
  points_legend_pchs[point_legend_cexes_log10_pvals < (log10(0.05))] <- 18
  
  points(y = 8.55 - lower_legend_by - cumsum(point_legend_cexes + pt_loc_expander) / sum(point_legend_cexes) * 0.05*n_points_for_legend, 
         x = rep(1.9925, n_points_for_legend), cex = point_legend_cexes, col = "grey50", pch = points_legend_pchs)
  text(labels = latex2exp::TeX("$log_{10}$(enrichment p-val)"), y = 8.6 - lower_legend_by, x = 1.975, cex = 1.1, pos = 4,  font = 2)
  # text(labels = paste0("0.05 FDR @ ", round(log10(0.05) - log10(nrow(ldsc_results_sub)),2), "\n           (Bonferroni)"), y = 8.45 , x = 2.1, cex = 0.8, pos = 4,  font = 2)
  text(labels = latex2exp::TeX(paste0("$\\alpha = 0.05$, IHW")), 
       y = 8.35 - lower_legend_by, x = 2.1, cex = 0.8, pos = 4,  font = 2)
  # text(labels = latex2exp::TeX(paste0("(Benjamini-Hochberg)")), 
  #      y = 8.35 - lower_legend_by, x = 2.07, cex = 0.8, pos = 4,  font = 2)
  # text(labels = paste0("< ", round(log10(0.05) - log10(nrow(ldsc_results_sub)),2), "\n> ", round(log10(0.05) - log10(nrow(ldsc_results_sub)),2)), 
  #      y = 8.25 , x = 2.175, cex = 1, pos = 4,  font = 2)
  text(labels = paste0("< ", round(log10(0.05),2), "\n> ", round(log10(0.05),2)), 
       y = 8.25 - lower_legend_by - 0.1 - c(0.15), x = 2.145, cex = 1, pos = 4,  font = 2)
  points(pch = c(19,18), y = c(8.41, 8.2) - lower_legend_by - 0.1 - 0.15, x = rep(2.145,2), cex = c(1.25, 1.75), 
         col = sapply(c(opacity_insig_points, opacity_sig_points), function(opcty) adjustcolor("grey50", alpha.f = opcty)))
  
  
  text(labels = point_legend_cexes_log10_pvals,y = 8.545 - lower_legend_by - cumsum(point_legend_cexes + pt_loc_expander) / sum(point_legend_cexes) * 0.05*n_points_for_legend, 
       x = rep(1.99, n_points_for_legend) + point_legend_cexes / 200, cex = 1, pos = 4, pch = 19)
  
  
  #figure out cex params
  for(cluster in cluster_names){
    
    # horizontal lines for colocalizing traits
    trait_locs <- ylocs
    if(use_enrichment){
      trait_probs <- ldsc_results_sub$Enrichment[ldsc_results_sub$cluster == cluster]
      #for bonferroni
      # trait_logPvals <- (ldsc_results_sub$logPVal_enrichment[ldsc_results_sub$cluster == cluster])[order_traits]
      #for BH adjustment
      trait_logPvals <- log10((ldsc_results_sub$adj_p[ldsc_results_sub$cluster == cluster])[order_traits])
    } else {
      trait_probs <- ldsc_results_sub$Prop._h2[ldsc_results_sub$cluster == cluster]
      trait_logPvals <- (ldsc_results_sub$logPVal_enrichment[ldsc_results_sub$cluster == cluster])[order_traits]
    }
    trait_probs <- trait_probs[order_traits]
    
    point_cex <- (-trait_logPvals / -minimum_enrichment_logPval)^point_cex_power
    point_cex <- point_cex * max_point_cex
    
    good_points <- trait_probs >= logprob_range[1] & trait_probs <= logprob_range[2]
    
    pchs <- rep(19, length(trait_probs))
    opacities <- rep(opacity_insig_points, length(trait_probs))
    
    #mark "significant" points
    # pchs[trait_logPvals < (log10(0.05) - log10(nrow(ldsc_results_sub)))] <- 18
    # opacities[trait_logPvals < (log10(0.05) - log10(nrow(ldsc_results_sub)))] <- opacity_sig_points
    pchs[trait_logPvals < log10(0.05)] <- 23
    opacities[trait_logPvals < log10(0.05)] <- opacity_sig_points
    
    point_cols <- sapply(opacities, function(opcty) grDevices::adjustcolor(cols$Tissue[match(cluster, cluster_names)], alpha.f = opcty))
    point_bgs <- point_cols
    point_cols[trait_logPvals < log10(0.05)] <- 1
    
    
    points(x = ((trait_probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc)[good_points], y = (trait_locs)[good_points], 
           pch = pchs[good_points], bg = point_bgs[good_points], col = point_cols[good_points], cex = point_cex[good_points])
    
    
  }
  
  
  
  if(arnold_mods){
    addImg(obj = arnold_bayes, x = 1.9, y = 1.6775, width = 0.5)
  }
  
}

dev.off()




grDevices::cairo_pdf(filename = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/figures/fig_S3_LDSC_output.pdf"), 
                     width = 1500 / 72, height = 2000 / 72 * 18 / 17 / 114 * length(coloc_phenotypes_sub) * 2, family="Arial Unicode MS", pointsize = 16)

par(mar = c(3,3,4,3), xpd = NA)
layout(matrix(c(1,2,1,3), 2, 2), heights = c(1,0.75))

for(type_of_plot in 2:2){
  
  use_enrichment <- c(T, F)[type_of_plot]
  
  max_horiz_axis <- 2.05
  
  
  #first get metainfo
  if(reorder_vertical){
    
    if(use_enrichment){
      order_traits <- order(sapply(coloc_phenotypes_sub, function(coloc_phenotype) mean(ldsc_results_sub[ldsc_results_sub$gwas == coloc_phenotype,]$Enrichment)), decreasing = T)
    } else {
      order_traits <- order(sapply(coloc_phenotypes_sub, function(coloc_phenotype) mean(ldsc_results_sub[ldsc_results_sub$gwas == coloc_phenotype,]$Prop._h2)), decreasing = T)
    }
    
  } else {
    order_traits <- 1:length(coloc_phenotypes_sub)
  }
  
  if(partition_by_category){
    order_traits <- order_traits[order(traitwise_partitions[,2][match(coloc_phenotypes_sub[order_traits], traitwise_partitions[,1])], sort(order_traits))]
    trait_category_counts <- as.data.table(table(traitwise_partitions[,2][match(coloc_phenotypes_sub[order_traits], traitwise_partitions[,1])] ))
    trait_category_counts$cN <- cumsum(trait_category_counts$N)
  }
  
  if(use_enrichment){
    if(fix_axes){
      logprob_range <- fix_axes_bounds_enrichment
    } else {
      logprob_range <- range(ldsc_results_sub$Enrichment)
    }
  } else {
    if(fix_axes){
      logprob_range <- fix_axes_bounds_heritability
    } else {
      logprob_range <- range(ldsc_results_sub$Prop._h2)
    }
  }
  
  logprob_range <- c(logprob_range[1] - diff(logprob_range) * buffer_min_and_max / 2, logprob_range[2] + diff(logprob_range) * buffer_min_and_max / 2)
  logprob_range <- c(0,ifelse(use_enrichment, 
                              max(ldsc_results_sub$Enrichment[ldsc_results_sub$adj_p < 0.05]) * 1.025, 
                              max(ldsc_results_sub$Prop._h2[ldsc_results_sub$adj_p < 0.05])) * 1.025)
  max_prob <- diff(logprob_range)
  
  nameloc <- 0.3 + ifelse(change_names | change_names_in_plot, 0, 0.3)
  
  plot(1,1,xlim = c(0,2.1), ylim = c(0,10), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
  
  ylocs <- seq(10, 0, length.out = length(coloc_phenotypes_sub))
  
  #put boxes around categories, if we're doing that  
  if(partition_by_category){
    segments(x0 = -0.1, x1 = -0.1, y0 = 0+diff(ylocs)[1]/2, y1 = 10-diff(ylocs)[1]/2, col = "grey35")
    for(traitcat in 1:(nrow(trait_category_counts)-1)){
      segments(y0 = ylocs[trait_category_counts$cN[traitcat]] + diff(ylocs)[1]/2, y1 = ylocs[trait_category_counts$cN[traitcat]] + diff(ylocs)[1]/2, x0 = -0.1, x1 = nameloc, col = "grey35")
    }
    
    catylocs <- (ylocs[trait_category_counts$cN] + c(ylocs[1], ylocs[trait_category_counts$cN][-length(trait_category_counts$cN)])) / 2 - diff(ylocs)[1]/1.5
    catnames <- as.vector(trait_category_counts$V1)
    catnames <- sapply(catnames, function(x) strsplit(x, split = c("-"))[[1]][1])
    catnames <- sapply(catnames, function(x) strsplit(x, split = c(" "))[[1]][1])
    
    #manually fudge category locations
    catylocs[catnames == "Hair"] <- catylocs[catnames == "Hair"] + diff(ylocs)[1]/1.5
    catylocs[catnames == "Aging"] <- catylocs[catnames == "Aging"] - diff(ylocs)[1]/3
    catylocs[catnames == "Endocrine"] <- catylocs[catnames == "Endocrine"] - diff(ylocs)[1]/10
    catylocs[catnames == "Allergy"] <- catylocs[catnames == "Allergy"] - diff(ylocs)[1]*1.2
    catylocs[catnames == "Anthropometric"] <- catylocs[catnames == "Anthropometric"] - diff(ylocs)[1]*1.2
    catylocs[catnames == "Cardiometabolic"] <- catylocs[catnames == "Cardiometabolic"] - diff(ylocs)[1]*2
    for(catname in 1:length(catnames)){
      text(x = -0.1, y = catylocs[catname] + ifelse(catnames[catname] == "Endocrine", 0.1, 0), labels = catnames[catname], pos = 2, srt = 90, col = "grey35")
    }
    
  }
  
  #plot names of traits
  if(change_names_in_plot){
    coloc_phenotypes_sub_newnames <- traitname_map[match(coloc_phenotypes_sub, traitname_map[,1]),2]
  } else {
    coloc_phenotypes_sub_newnames <- coloc_phenotypes_sub
  }
  
  if(nicole_mods){
    # text(paste0(1:length(coloc_phenotypes_sub_newnames), ": ", 
    #             coloc_phenotypes_sub_newnames)[order_traits], 
    #      col = bad_boys_cols[order_traits], x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
    text(paste0(coloc_phenotypes_sub_newnames)[order_traits], 
         col = bad_boys_cols[order_traits], x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
  } else {
    text(paste0(1:length(coloc_phenotypes_sub_newnames), ": ", 
                coloc_phenotypes_sub_newnames)[order_traits], 
         col = pheno_cols, x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
  }
  
  #vert axis label
  text(labels = "Phenotypes", x = nameloc, y = 10.45, pos = 2, cex = 2.5, xpd = NA, family="Courier", col = "grey25")
  #horiz axis label
  text(labels = paste0("Heritability ", ifelse(use_enrichment, "Enrichment", "Proportion"), " Across Tissues"), 
       x = 1.2, y = -0.675, cex = 2.25, pos = 1, xpd = NA)
  
  
  #guiding lines for traits
  if(nicole_mods){
    segments(x0 = nameloc, x1 = max_horiz_axis, y0 = ylocs, y1 = ylocs, col = "grey80", lty = 3, lwd = 0.5)
  }
  
  #axes
  segments(x0 = nameloc, x1 = nameloc, y0 = 10.675, y1 = -0.1, lwd = 2)
  segments(x0 = nameloc, x1 = max_horiz_axis, y0 = -0.1, y1 = -0.1, lwd = 2)
  
  #horiz axis ticks and nums
  # if(use_geometric_mean){
  #   
  #   probs <- ceiling(logprob_range[1]):ceiling(logprob_range[2])
  #   segments(x0 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
  #            x1 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
  #            y0 = -0.1, y1 = -0.15, lwd = 2)
  #   text(x = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
  #        y = -0.135, pos = 1, cex = 1.25,
  #        labels =  c(sapply(probs, function(num) as.expression(bquote(10^.(num))))))
  #   
  # } else {
  #   
  #   segments(x0 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
  #            x1 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
  #            y0 = -0.1, y1 = -0.15, lwd = 2)
  #   text(x = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
  #        y = -0.135, pos = 1, cex = 1.25,
  #        labels =  c(sapply(seq(ceiling(logprob_range[1]), ceiling(logprob_range[2]), by = 1), function(num) as.expression(bquote(10^.(num))))))
  # }
  
  if(use_enrichment){
    probs <- seq(ceiling(logprob_range[1]), floor(logprob_range[2]), by = ifelse(fix_axes, 0.5, 1))  
  } else {
    probs <- seq(ceiling(logprob_range[1]), (logprob_range[2]), by = ifelse(fix_axes, 0.1, 0.1))
  }
  
  segments(x0 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc,
           x1 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc,
           y0 = -0.1, y1 = -0.25, lwd = 2)
  text(x = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc,
       y = -0.25, pos = 1, cex = 1.5,
       labels =  probs)
  
  
  # #let's get minor tick marks in there too
  # if(use_geometric_mean){
  #   minticks <- log10(as.vector(sapply(1:(length(probs)-1), function(p) seq(10^probs[p], 10^probs[p+1], length.out = 10))))
  #   minticks <- (minticks - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc
  #   minticks <- minticks[minticks < max_horiz_axis]
  #   segments(x0 = minticks, x1 =  minticks, y0 = -0.1, y1 = -0.125, lwd = 1)
  #   #hmm being weird
  # } else {
  #   minticks <- log10(as.vector(t(t(1:9)) %*% t(10^seq(0, max_prob, by = 10^round(log(max_prob))))))
  #   minticks <- minticks / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025
  #   minticks <- minticks[minticks < max_horiz_axis]
  #   segments(x0 = minticks, x1 =  minticks, y0 = -0.1, y1 = -0.125, lwd = 1)
  # }
  
  
  #dashed lines for ticks
  segments(x0 = (probs[-1] - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
           x1 = (probs[-1] - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y0 = 10, y1 = -0.075, lwd = 2, lty = 2, col = "grey75")
  segments(x0 = (1 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
           x1 = (1 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
           y0 = 10, y1 = -0.075, lwd = 3, lty = 1, col = "grey75")
  
  
  #helpful hint line about direction
  
  # if(use_geometric_mean){
  #   segments(x0 = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
  #            x1 = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y0 = 10, y1 = 10.4, lwd = 2, lty = 2, col = "grey75")
  #   text(labels = "non-DE Genes", srt = 90, x = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = 10.4, pos = 2, col = "grey75", cex = 0.75)
  #   text(labels = "DE Genes", srt = 270, x = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = 10.4, pos = 4, col = "grey75", cex = 0.75)
  #   text(labels = "Stronger Coloc\nSignal In", srt = 0, x = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = 10.3875, pos = 3, col = "grey75", cex = 0.75)
  # } else {
  #   segments(x0 = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
  #            x1 = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y0 = 10, y1 = 10.4, lwd = 2, lty = 2, col = "grey75")
  #   text(labels = "eQTL > DEG", srt = 90, x = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = 10.4, pos = 2, col = "grey75", cex = 0.75)
  #   text(labels = "DEG > eQTL", srt = 270, x = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = 10.4, pos = 4, col = "grey75", cex = 0.75)
  #   text(labels = "Stronger Coloc\nSignal In", srt = 0, x = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = 10.3875, pos = 3, col = "grey75", cex = 0.75)
  # }
  
  #plot legend for colors etc.
  lower_legend_by <- 1.6
  rect(xleft = 1.97, ybottom = 8.15 - lower_legend_by - 0.05*n_points_for_legend, 
       ytop = 10.1 - lower_legend_by, xright = 2.2, border = NA, col = "white")
  # text(labels = sapply(names(deg_eqtl_list), function(ts) stringr::str_to_title(paste0(strsplit(ts, "-")[[1]][-1], collapse = " "))),
  #      y = seq(9.95,8.75,length.out = length(names(deg_eqtl_list))), x = 2, pos = 4)
  # points(x = rep(1.9925, length(names(deg_eqtl_list))), y = seq(9.955,8.755,length.out = length(names(deg_eqtl_list))), col = cols$Tissue, pch = 15, cex = 1.75)
  # text(labels = stringr::str_replace_all(cluster_names, "_", " "),
  #      y = seq(9.95,8.75,length.out = length(cluster_names)), x = 2, pos = 4)
  text(labels = sapply(cluster_names, function(cli) strsplit(cli, "-")[[1]][1]),
       y = seq(9.95,8.75 - lower_legend_by,length.out = length(cluster_names)) - 0.02, x = 2, pos = 4)
  points(x = rep(1.9925, length(cluster_names)), y = seq(9.955,8.755 - lower_legend_by,length.out = length(cluster_names)), col = cols$cluster, pch = 19, cex = 1.75)
  
  #plot legend for points
  lower_legend_by <- lower_legend_by + 0.3
  pt_loc_expander <- 8
  point_legend_cexes <- seq(from = 0.4, to = max_point_cex, length.out = n_points_for_legend)
  points_legend_pchs <- rep(19, n_points_for_legend)
  point_legend_cexes_log10_pvals <- round((point_legend_cexes / max_point_cex)^(1/point_cex_power) * minimum_enrichment_logPval, 2)
  # points_legend_pchs[point_legend_cexes_log10_pvals < (log10(0.05) - log10(nrow(ldsc_results_sub)))] <- 18
  points_legend_pchs[point_legend_cexes_log10_pvals < (log10(0.05))] <- 18
  
  points(y = 8.55 - lower_legend_by - cumsum(point_legend_cexes + pt_loc_expander) / sum(point_legend_cexes) * 0.05*n_points_for_legend, 
         x = rep(1.9925, n_points_for_legend), cex = point_legend_cexes, col = "grey50", pch = points_legend_pchs)
  text(labels = latex2exp::TeX("$log_{10}$(enrichment p-val)"), y = 8.6 - lower_legend_by, x = 1.975, cex = 1.1, pos = 4,  font = 2)
  # text(labels = paste0("0.05 FDR @ ", round(log10(0.05) - log10(nrow(ldsc_results_sub)),2), "\n           (Bonferroni)"), y = 8.45 , x = 2.1, cex = 0.8, pos = 4,  font = 2)
  text(labels = latex2exp::TeX(paste0("$\\alpha = 0.05$, IHW")), 
       y = 8.35 - lower_legend_by, x = 2.1, cex = 0.8, pos = 4,  font = 2)
  # text(labels = latex2exp::TeX(paste0("(Benjamini-Hochberg)")), 
  #      y = 8.35 - lower_legend_by, x = 2.07, cex = 0.8, pos = 4,  font = 2)
  # text(labels = paste0("< ", round(log10(0.05) - log10(nrow(ldsc_results_sub)),2), "\n> ", round(log10(0.05) - log10(nrow(ldsc_results_sub)),2)), 
  #      y = 8.25 , x = 2.175, cex = 1, pos = 4,  font = 2)
  text(labels = paste0("< ", round(log10(0.05),2), "\n> ", round(log10(0.05),2)), 
       y = 8.25 - lower_legend_by - 0.1 - c(0.15), x = 2.145, cex = 1, pos = 4,  font = 2)
  points(pch = c(19,18), y = c(8.41, 8.2) - lower_legend_by - 0.1 - 0.15, x = rep(2.145,2), cex = c(1.25, 1.75), 
         col = sapply(c(opacity_insig_points, opacity_sig_points), function(opcty) adjustcolor("grey50", alpha.f = opcty)))
  
  
  text(labels = point_legend_cexes_log10_pvals,y = 8.545 - lower_legend_by - cumsum(point_legend_cexes + pt_loc_expander) / sum(point_legend_cexes) * 0.05*n_points_for_legend, 
       x = rep(1.99, n_points_for_legend) + point_legend_cexes / 200, cex = 1, pos = 4, pch = 19)
  
  
  #figure out cex params
  for(cluster in cluster_names){
    
    # horizontal lines for colocalizing traits
    trait_locs <- ylocs
    if(use_enrichment){
      trait_probs <- ldsc_results_sub$Enrichment[ldsc_results_sub$cluster == cluster]
      #for bonferroni
      # trait_logPvals <- (ldsc_results_sub$logPVal_enrichment[ldsc_results_sub$cluster == cluster])[order_traits]
      #for BH adjustment
      trait_logPvals <- log10((ldsc_results_sub$adj_p[ldsc_results_sub$cluster == cluster])[order_traits])
    } else {
      trait_probs <- ldsc_results_sub$Prop._h2[ldsc_results_sub$cluster == cluster]
      trait_logPvals <- (ldsc_results_sub$logPVal_enrichment[ldsc_results_sub$cluster == cluster])[order_traits]
    }
    trait_probs <- trait_probs[order_traits]
    
    point_cex <- (-trait_logPvals / -minimum_enrichment_logPval)^point_cex_power
    point_cex <- point_cex * max_point_cex
    
    good_points <- trait_probs >= logprob_range[1] & trait_probs <= logprob_range[2]
    
    pchs <- rep(19, length(trait_probs))
    opacities <- rep(opacity_insig_points, length(trait_probs))
    
    #mark "significant" points
    # pchs[trait_logPvals < (log10(0.05) - log10(nrow(ldsc_results_sub)))] <- 18
    # opacities[trait_logPvals < (log10(0.05) - log10(nrow(ldsc_results_sub)))] <- opacity_sig_points
    pchs[trait_logPvals < log10(0.05)] <- 23
    opacities[trait_logPvals < log10(0.05)] <- opacity_sig_points
    
    point_cols <- sapply(opacities, function(opcty) grDevices::adjustcolor(cols$Tissue[match(cluster, cluster_names)], alpha.f = opcty))
    point_bgs <- point_cols
    point_cols[trait_logPvals < log10(0.05)] <- 1
    
    
    points(x = ((trait_probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc)[good_points], y = (trait_locs)[good_points], 
           pch = pchs[good_points], bg = point_bgs[good_points], col = point_cols[good_points], cex = point_cex[good_points])
    
    
  }
  
  
  
  if(arnold_mods){
    addImg(obj = arnold_bayes, x = 1.9, y = 1.6775, width = 0.5)
  }
  
}


sub_1 <- ldsc_results_cts_alt_fullyconditional[ldsc_results_cts_alt_fullyconditional$gwas %in% traits_to_keep,]
sub_2 <- ldsc_results_cts_alt[ldsc_results_cts_alt$gwas %in% traits_to_keep,]
tissues <- tissue_order[tissue_order %in% unique(gsub(x = sub_1$cluster, "-sex_homogeneous_changing", ""))]


par(mar = c(4.25,4.25,4,2) * 1.5, xpd = F)
xlims <- range(sub_1$Enrichment)
ylims <- range(sub_2$Enrichment)
plot(sub_1$Enrichment,  sub_2$Enrichment, main = latex2exp::TeX("$h^2_{SNP}$ Enrichment"),
     xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "", pch = 19, cex.main = 2.5,
     col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[gsub(x = sub_1$cluster, "-sex_homogeneous_changing", "")], 0.5),
     xlim = xlims, ylim = ylims, cex = 2)
rect(xleft = par("usr")[1], ybottom = par("usr")[3], xright = par("usr")[2], ytop = par("usr")[4], lwd = 2, xpd = NA)
text(x = mean(xlims), y = ylims[1] - diff(ylims) / 7, labels = "Enrichment (Conditional on Tissue Annotation)", cex = 1.75, pos = 1, xpd = NA)
text(x = xlims[1] - diff(xlims) / 6, y = mean(ylims), labels = "Enrichment (Unconditional on Tissue Annotation)", cex = 1.75, srt = 90, xpd = NA)

#axes
xvals <- round(seq(ceiling(xlims[1]), floor(xlims[2]), length.out = 5))
xvals <- seq(ceiling(xlims[1]), floor(xlims[2]), by = diff(xvals)[1])
segments(x0 = xvals, x1 = xvals, y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4]) / 50, lwd = 2, xpd = NA)
text(xvals, y = par("usr")[3] - diff(par("usr")[3:4]) / 40, labels = xvals, cex = 1.75, pos = 1, xpd = NA)

yvals <- round(seq(ceiling(ylims[1]), floor(ylims[2]), length.out = 5))
yvals <- seq(ceiling(ylims[1]), floor(ylims[2]), by = diff(yvals)[1])
segments(y0 = yvals, y1 = yvals, x0 = par("usr")[1], x1 = par("usr")[1] - diff(par("usr")[1:2]) / 50, lwd = 2, xpd = NA)
text(yvals, x = par("usr")[1] - diff(par("usr")[1:2]) / 50, labels = yvals, cex = 1.75, pos = 2, xpd = NA)


legend("topleft", pch = 19, col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues], 0.5), legend = tissues, cex = 1.25, pt.cex = 2.5, bty="n")
legend("bottomright", lty = 2, col = "red", legend = "1-to-1 line", bty="n", lwd = 3, cex = 1.5)
abline(0,1, col = "red", lty = 2, lwd = 4)


#second figure, prop h2

xlims <- range(sub_1$Prop._h2)
ylims <- range(sub_2$Prop._h2)
plot(sub_1$Prop._h2,  sub_2$Prop._h2, main = latex2exp::TeX("Proportion $h^2_{SNP}$"),
     xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "", pch = 19, cex.main = 2.5,
     col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[gsub(x = sub_1$cluster, "-sex_homogeneous_changing", "")], 0.5),
     xlim = xlims, ylim = ylims, cex = 2)
rect(xleft = par("usr")[1], ybottom = par("usr")[3], xright = par("usr")[2], ytop = par("usr")[4], lwd = 2, xpd = NA)
text(x = mean(xlims), y = ylims[1] - diff(ylims) / 7, labels = latex2exp::TeX("Proportion $h^2_{SNP}$ (Conditional on Tissue Annotation)"), cex = 1.75, pos = 1, xpd = NA)
text(x = xlims[1] - diff(xlims) / 4.5, y = mean(ylims), labels = latex2exp::TeX("Proportion $h^2_{SNP}$ (Unconditional on Tissue Annotation)"), cex = 1.75, srt = 90, xpd = NA)

#axes
xvals <- round(seq((xlims[1]), (xlims[2]), length.out = 5), 2)
segments(x0 = xvals, x1 = xvals, y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4]) / 50, lwd = 2, xpd = NA)
text(xvals, y = par("usr")[3] - diff(par("usr")[3:4]) / 40, labels = xvals, cex = 1.75, pos = 1, xpd = NA)

yvals <- round(seq((ylims[1]), (ylims[2]), length.out = 5), 2)
segments(y0 = yvals, y1 = yvals, x0 = par("usr")[1], x1 = par("usr")[1] - diff(par("usr")[1:2]) / 50, lwd = 2, xpd = NA)
text(yvals, x = par("usr")[1] - diff(par("usr")[1:2]) / 50, labels = yvals, cex = 1.75, pos = 2, xpd = NA)

legend("topleft", pch = 19, col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues], 0.5), legend = tissues, cex = 1.25, pt.cex = 2.5, bty="n")
legend("bottomright", lty = 2, col = "red", legend = "1-to-1 line", bty="n", lwd = 3, cex = 1.5)
abline(0,1, col = "red", lty = 2, lwd = 4)

fig_label("c)", xpd = NA, cex = 3, shrinkX = 0.75)
fig_label("b)", xpd = NA, cex = 3, shrinkX = 7.5)
fig_label("a)", xpd = NA, cex = 3, shrinkX = 7.75, shrinkY = 2.65)

dev.off()
