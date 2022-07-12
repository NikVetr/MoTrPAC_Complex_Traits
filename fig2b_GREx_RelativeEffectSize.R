#### GREX, or just expression data works too ####
# load libraries
library(ks)
library(arrow)
library(data.table)
library(edgeR)
library(limma)
library(DESeq2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(org.Rn.eg.db)
library(biomaRt)
library(MotrpacBicQC)
library(plotrix)
library(ggplot2)
library(testit)
library(circlize)
library(jpeg)
library(doParallel)
library(pracma)
library(data.table)
library(foreach)
library(doParallel)
library(parallel)

#### define functions ####
source(file = "~/scripts/montgomery_lab/deg-trait_functions.R")

#### load deg-eqtl data ####
source(file = "~/scripts/montgomery_lab/load_deg-eqtl_merged_file.R")

#initialize parallelization
foreach_parallel <- F
if(foreach_parallel){
  if(!exists("cl")){
    cl <- makeCluster(3, outfile="")
    registerDoParallel(cl)
  }
  getDoParWorkers()
}

#packages
library(cmdstanr)
library(posterior)
library(caret)
library(MASS)
library(invgamma)
library(fitdistrplus)

#stan model for sd
# stan_program <- '
# data {
#     int<lower=0> n_gene;
#     int<lower=0> n_male;
#     int<lower=0> n_female;
#     vector[n_gene * n_male] male;
#     vector[n_gene * n_female] female;
#     int<lower=1, upper=n_gene> male_gene_index[n_gene * n_male];
#     int<lower=1, upper=n_gene> female_gene_index[n_gene * n_female];
#     
# }
# parameters {
#     real mean_log_sd;
#     real<lower=0> sd_logsd;
#     vector[n_gene] raw_gene_logsd;
#     real<lower=0> sd_gene_logsd;
#     vector[n_gene] raw_male_gene_logsd;
#     vector[n_gene] raw_female_gene_logsd;
# }
# transformed parameters {
#     vector[n_gene] gene_logsd = raw_gene_logsd * sd_logsd + mean_log_sd;
#     vector[n_gene] male_gene_logsd = raw_male_gene_logsd * sd_gene_logsd + gene_logsd;
#     vector[n_gene] female_gene_logsd = raw_female_gene_logsd * sd_gene_logsd + gene_logsd;
# }
# model {
#     //priors
#     mean_log_sd ~ normal(0,3);
#     sd_logsd ~ std_normal();
#     raw_gene_logsd ~ std_normal();
#     sd_gene_logsd ~ std_normal();
#     raw_male_gene_logsd ~ std_normal();
#     raw_female_gene_logsd ~ std_normal();
#     
#     //likelihood
#     male ~ normal(0, exp(male_gene_logsd)[male_gene_index]);
#     female ~ normal(0, exp(female_gene_logsd)[female_gene_index]);
# }
# '
# 
# 
# if(!exists("curr_stan_program") || stan_program != curr_stan_program){
#   curr_stan_program <- stan_program
#   f <- write_stan_file(stan_program)
# }
# mod <- cmdstan_model(f)


#### estimate variances ####
tissues <- names(motrpac_gtex_map)

run_sds_estimation <- F
if(run_sds_estimation){
  sds_expression <- list()
  invgamma_estimates <- list()
  n_peer_factors <- data.frame(n_factors = c(15,30,45,60), N = c(0,150,250,350))
  
  for(tissue in tissues){
    print(tissue)
    tissue_filename_keyword <- gsub(" ", "_", gsub(" - ", "_", motrpac_gtex_map[tissue]))
    
    #read in new and old expression files
    d <- fread(file = paste0("/Volumes/SSD500GB/gtex-pipeline/log2-normalized-expression/log2-normalized-expression_",tissue,".expression.bed.gz"))
    d_orig_genes <- gsub('\\..*','',fread(file = paste0("/Volumes/SSD500GB/gtex-pipeline/GTEx_Analysis_v8_eQTL_expression_matrices/",
                                                        tissue_filename_keyword,".v8.normalized_expression.bed.gz"))$gene_id)
    d$gene_id <- gsub('\\..*','',d$gene_id)
    d <- d[d$gene_id %in% d_orig_genes,]
    
    #get useful informations
    expr_cols <- grep("GTEX", colnames(d))
    n_indiv <- length(expr_cols)
    n_peer_factors_to_use <- n_peer_factors$n_factors[sum(n_peer_factors$N <= n_indiv)]
    
    #read in covariates and subset to the correct # of peer factors
    covariates <- as.data.frame(fread(file = paste0("/Volumes/SSD500GB/gtex-pipeline/GTEx_Analysis_v8_eQTL_covariates/", tissue_filename_keyword,".v8.covariates.txt")))
    covariates_to_include <- match(setdiff(covariates$ID, covariates$ID[grep("Inferr", covariates$ID)][-(1:n_peer_factors_to_use)]), covariates$ID)
    covariates <- covariates[covariates_to_include,]
    
    #compute residuals of genes from covariates file
    d_resid <- do.call(rbind, mclapply(1:nrow(d), function(ri){
      # if(ri%%1000==0){system(sprintf('echo "%s"', paste0(ri, " ", collapse="")))}
      lm(as.numeric(d[ri,..expr_cols]) ~ t(as.matrix(covariates[,-1])))$resid
    }, mc.cores = 12))
    colnames(d_resid) <- colnames(d)[expr_cols]
    
    # ec <- fread(file = paste0("~/repos/gtex-pipeline/expression_data/GTEx_Analysis_v8_",tissue,"_expected_count.gct.gz"))
    # ec$gene_id <- gsub('\\..*','',ec$gene_id)
    # 
    # tpm <- fread(file = paste0("~/repos/gtex-pipeline/expression_data/GTEx_Analysis_v8_",tissue,"_tpm.gct.gz"))
    # tpm$gene_id <- gsub('\\..*','',tpm$gene_id)
    # 
    # diff(as.integer(sapply(colnames(d), function(name) grep(name, colnames(ec)))))
    # diff(as.integer(sapply(colnames(d), function(name) grep(name, colnames(tpm)))))
    # 
    # plot(as.numeric(d[1,-c(1:4)]), log2(as.numeric(ec[match(d$gene_id[1], ec$gene_id),-1])+1))
    # plot(as.numeric(d[3,-c(1:4)]), log2(as.numeric(tpm[match(d$gene_id[3], tpm$gene_id),-1])+1))
    # plot(as.numeric(tpm[1,-1]), as.numeric(ec[1,-1]))
    
    #get sex information
    #in GTEx, 1=Male	2=Female, https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_Analysis_v7_Annotations_SubjectPhenotypesDD.xlsx
    gtex_sex <- covariates[covariates$ID == "sex",] 
    males <- colnames(gtex_sex)[gtex_sex == 1]
    females <- colnames(gtex_sex)[gtex_sex == 2]
    
    # males <- as.matrix(d[,..males]) - apply(as.matrix(d[,..males]), 1, mean)
    # females <- as.matrix(d[,..females]) - apply(as.matrix(d[,..females]), 1, mean)
    # mean_centered <- cbind(males, females)
    # plot(setNames(apply(females, 1, sd), d$gene_id), setNames(apply(mean_centered, 1, sd), d$gene_id))
    
    #compute sample variance and estimate invgamma dist hyperparameters
    sample_pooled_vars <- apply(d_resid, 1, var)
    sample_pooled_vars <- sample_pooled_vars[sample_pooled_vars != 0]
    invgamma_estimate <- fitdist(sample_pooled_vars, "invgamma", method = "mge")$estimate
    
    # print(invgamma_estimate)
    # hist(sample_pooled_vars, probability = T, breaks = 100)
    # lines(0:1000/100, dinvgamma(0:1000/100, shape = invgamma_estimate[1], scale = invgamma_estimate[2]))
    
    #estimate posterior mean SDs
    posterior_mean_sds <- as.data.frame(do.call(rbind, lapply(1:nrow(d_resid), function(i){
      
      posterior_params <- c(alpha = invgamma_estimate[1] + ncol(d_resid),
                            beta = invgamma_estimate[2] + sum((d_resid[i,])^2) / 2)
      est_var <- as.numeric(posterior_params["beta.scale"] / (posterior_params["alpha.shape"]-1))
      
      posterior_params_m <- c(alpha = invgamma_estimate[1] + length(males),
                              beta = invgamma_estimate[2] + sum((d_resid[i,males])^2) / 2)
      est_var_m <- as.numeric(posterior_params_m["beta.scale"] / (posterior_params_m["alpha.shape"]-1))
      
      posterior_params_f <- c(alpha = invgamma_estimate[1] + length(females),
                              beta = invgamma_estimate[2] + sum((d_resid[i,females])^2) / 2)
      est_var_f <- as.numeric(posterior_params_f["beta.scale"] / (posterior_params_f["alpha.shape"]-1))
      
      return(c(est_sd = sqrt(est_var), est_sd_m = sqrt(est_var_m), est_sd_f = sqrt(est_var_f)))
    })))
    rownames(posterior_mean_sds) <- d$gene_id
    posterior_mean_sds$tissue <- tissue
    
    # plot(apply(d_resid, 1, sd), posterior_mean_sds$est_sd); abline(0,1,col=2,lwd=2)
    # plot(posterior_mean_sds$est_sd_m, posterior_mean_sds$est_sd_f); abline(0,1,col=2,lwd=2)
    
    
    # data <- list(
    #   n_gene = length(d$gene_id),
    #   n_male = length(males),
    #   n_female = length(females),
    #   male = c(as.matrix(d[,..males]) - apply(as.matrix(d[,..males]), 1, mean)),
    #   female = c(as.matrix(d[,..females] - apply(as.matrix(d[,..females]), 1, mean))),
    #   male_gene_index = rep(1:length(d$gene_id), length(males)),
    #   female_gene_index = rep(1:length(d$gene_id), length(females))
    # )
    #   
    # out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3,
    #                   data = data,
    #                   parallel_chains = 4,
    #                   adapt_delta = 0.85, refresh = 10, init = 0.1, max_treedepth = 15, thin = 5)
    # samps <- data.frame(as_draws_df(out$draws()))
    # hist(samps$gene_logsd.3., breaks =100)
    
    # plot(setNames(apply(d[,..males], 1, sd), d$gene_id), setNames(apply(d[,..females], 1, sd), d$gene_id))
    
    #write to object
    invgamma_estimates[[tissue]] <- invgamma_estimate
    sds_expression[[tissue]] <- posterior_mean_sds
  }
  
  save(x = sds_expression, file = "~/data/smontgom/GREx_sds_expression")
  save(x = invgamma_estimates, file = "~/data/smontgom/GREx_invgamma_estimate")
} else {
  load("~/data/smontgom/GREx_sds_expression")
  load("~/data/smontgom/GREx_invgamma_estimate")
}



# sapply(sds_expression, function(x) length(x))
# if(all(sapply(tissues, function(ti) all(names(sds_expression[[1]]) == names(sds_expression[[ti]]))))){
#   sds_expression <- do.call(rbind, sds_expression)
# }
# rownames(sds_expression) <- tissues


GTEx_SampleSize <- read.csv("~/data/smontgom/GTEx_SampleSizes.csv", header = T)
colnames(GTEx_SampleSize) <- gsub("X..", "", colnames(GTEx_SampleSize))
GTEx_SampleSize <- data.frame(tissue = names(motrpac_gtex_map), sample_size = GTEx_SampleSize$RNASeq.and.Genotyped.samples[match(motrpac_gtex_map, GTEx_SampleSize$Tissue)])
# dev.off()
# hist(sds_expression)


#### process GWAS ####
ldsc_directory <- "~/repos/ldsc/"

process_GWASy_sumstats <- F
if(process_GWASy_sumstats){
  foreach(tissue=tissues, .packages = c("data.table", "arrow")) %dopar% {
    # for(tissue in tissues){
    print(tissue)
    if(!dir.exists(paste0(ldsc_directory, "/GTEx_v8_log2norm_sumstats/", tissue))){dir.create(paste0(ldsc_directory, "/GTEx_v8_log2norm_sumstats/", tissue))}
    
    for(cri in c(1:22,"X")){
      
      cat(paste0(" ", cri))
      if(!dir.exists(paste0(ldsc_directory, "/GTEx_v8_log2norm_sumstats/", tissue, "/", cri))){dir.create(paste0(ldsc_directory, "/GTEx_v8_log2norm_sumstats/", tissue, "/", cri))}
      
      RSID_POS_MAP <- fread(paste0("~/data/smontgom/RSID_POS_MAP_",cri,".txt"))
      
      eQTL_sumstats <- read_parquet(paste0("~/repos/gtex-pipeline/tensorQTL_output/log2-normalized-expression_", tissue,".cis_qtl_pairs.chr", cri, ".parquet"))
      eQTL_sumstats$ENSG <- gsub('\\..*','',eQTL_sumstats$phenotype_id)
      
      vids <- eQTL_sumstats$variant_id
      vids <- substr(vids, start = 4, nchar(vids))
      vids <- strsplit(vids, "_", T)
      vids <- as.data.table(data.frame(data.table::transpose(vids)))
      colnames(vids) <- c("CHROM", "POS", "REF", "ALT", "BUILD")
      vids$RSID <- RSID_POS_MAP$ID[match(vids$POS, RSID_POS_MAP$POS)]
      eQTL_sumstats <- cbind(eQTL_sumstats, vids)
      eQTL_sumstats <- eQTL_sumstats[!is.na(eQTL_sumstats$RSID),]
      eQTL_sumstats$Z <- round(eQTL_sumstats$slope / eQTL_sumstats$slope_se, 5)
      eQTL_sumstats$N <-  paste0(GTEx_SampleSize$sample_size[GTEx_SampleSize$tissue == tissue], ".00")
      
      genes <- unique(eQTL_sumstats$ENSG)
      eQTL_sumstats <- split(eQTL_sumstats , f = eQTL_sumstats$ENSG)
      eQTL_sumstats <- lapply(eQTL_sumstats, function(gid){
        gene_info <- gid[,c("RSID", "REF", "ALT", "Z", "N")]
        colnames(gene_info) <- c("SNP", "A1", "A2", "Z", "N")
        return(gene_info)
      })
      
      for(gene_name in names(eQTL_sumstats)){
        
        fwrite(eQTL_sumstats[[gene_name]], file = paste0(ldsc_directory, "/GTEx_v8_log2norm_sumstats/", tissue, "/", cri, "/", gene_name, ".sumstats.gz"), sep = "\t", append = F, col.names = T)
      }
      
    }
  }
}

#compare to heritability estimates from https://github.com/WheelerLab/GenArchDB
library(data.table)
h2dir <- "~/repos/GenArchDB/GenArch_reml-no-constrain_h2/"
h2_files <- list.files(h2dir)
h2_files <- h2_files[intersect(grep("GTEx", h2_files), grep(".TS.", h2_files))]
h2_tissue_to_file <- sapply(motrpac_gtex_map, function(tissue) h2_files[intersect_rec(sapply(strsplit(tissue, "_")[[1]], function(x) grep(x, h2_files)))][1])
tissue = names(motrpac_gtex_map)[7]
h2 <- fread(paste0(h2dir, h2_tissue_to_file[tissue]))
h2$gene_id <- gsub('\\..*','',h2$ensid)

#compare to gcta heritability estimates
gcta_directory <- "/Volumes/SSD500GB/gcta_1.93.2beta_mac/"
load(paste0(gcta_directory, "gcta_output_GTEx_allTissues.RData")) #gcta_output

do_EDA <- F
if(do_EDA){
  par(mfrow = c(4,4))
  for(tissue in setdiff(names(motrpac_gtex_map), "t59-kidney")){
    h2 <- fread(paste0(h2dir, h2_tissue_to_file[tissue]))
    h2$gene_id <- gsub('\\..*','',h2$ensid)
    h2s <- h2$global.h2[match(deg_eqtl_list[[tissue]]$human_ensembl_gene.x, h2$gene_id)]
    prop_pos <- round(sum(h2s > 0, na.rm = T) / sum(!is.na(h2s)) * 100, 1)
    hist(h2s, main= paste0(tissue, " (", prop_pos, "% positive)"), xlab = "heritability")
  }
  for(tissue in setdiff(names(motrpac_gtex_map), "t59-kidney")){
    h2 <- as.data.frame(fread(paste0(h2dir, h2_tissue_to_file[tissue])))
    h2$gene_id <- gsub('\\..*','',h2$ensid)
    h2$myh2 <- gcta_output[[tissue]]$h2[match(h2$gene_id, gcta_output[[tissue]]$ENSG)]
    h2s <- h2[match(deg_eqtl_list[[tissue]]$human_ensembl_gene.x[deg_eqtl_list[[tissue]]$comparison_group == "1w" & deg_eqtl_list[[tissue]]$sex == deg_eqtl_list[[tissue]]$sex[1]], 
                    h2$gene_id), c("local.h2", "myh2")]
    plot(h2s[complete.cases(h2s),], main= paste0(tissue, " (r = ", round(cor(h2s[complete.cases(h2s),])[1,2], 3), ")"), 
         xlab = "wheeler h2", ylab = "my gcta h2", pch = 19, col = adjustcolor(1, 0.1))
  }
}

#find where GCTA reliably estimated h2
do_ihw <- F
if(do_ihw){
  load(paste0(gcta_directory, "gcta_output_GTEx_allTissues.RData")) #gcta_output
  gcta_output_df <- do.call(rbind, gcta_output)
  gcta_output_df$tissue <- as.factor(gcta_output_df$tissue)
  gcta_output_df$p <- gcta_output_df$p*2
  ihw_results_gcta <- IHW::ihw(p ~ tissue, data = gcta_output_df, alpha = 0.1)
  gcta_alpha <- 0.1
  sum(ihw_results_gcta@df$adj_pvalue < gcta_alpha)
  # hist(ihw_results_gcta@df$adj_pvalue, breaks = 20)
  gcta_output_df$adj_p <- ihw_results_gcta@df$adj_pvalue
  gcta_output <- lapply(setNames(tissues, tissues), function(tiss) gcta_output_df[gcta_output_df$tissue == tiss & gcta_output_df$adj_p < gcta_alpha,])
  save(gcta_output, file = paste0(gcta_directory, "gcta_output_GTEx_allTissues_list_IHW.RData"))
} else {
  load(paste0(gcta_directory, "gcta_output_GTEx_allTissues_list_IHW.RData"))
}



#### compute relative exercise effect z-scores ####
load(paste0(gcta_directory, "gcta_output_GTEx_allTissues_list_IHW.RData"))
for(tissue in tissues){
  print(tissue)
  male_inds <- which(deg_eqtl_list[[tissue]]$sex == "male")
  female_inds <- which(deg_eqtl_list[[tissue]]$sex == "female")
  
  deg_eqtl_list[[tissue]]$phenotypic_expression_Z <- 0
  deg_eqtl_list[[tissue]]$phenotypic_expression_Z[male_inds] <- deg_eqtl_list[[tissue]]$logFC[male_inds] / 
    sds_expression[[tissue]]$est_sd_m[match(deg_eqtl_list[[tissue]]$human_ensembl_gene.x[male_inds], rownames(sds_expression[[tissue]]))]
  deg_eqtl_list[[tissue]]$phenotypic_expression_Z[female_inds] <- deg_eqtl_list[[tissue]]$logFC[female_inds] / 
    sds_expression[[tissue]]$est_sd_f[match(deg_eqtl_list[[tissue]]$human_ensembl_gene.x[female_inds], rownames(sds_expression[[tissue]]))]
  
  deg_eqtl_list[[tissue]]$genetic_expression_Z <- 0
  deg_eqtl_list[[tissue]]$genetic_expression_Z[male_inds] <- deg_eqtl_list[[tissue]]$logFC[male_inds] / 
    sds_expression[[tissue]]$est_sd_m[match(deg_eqtl_list[[tissue]]$human_ensembl_gene.x[male_inds], rownames(sds_expression[[tissue]]))] / 
    sqrt(gcta_output[[tissue]]$h2[match(deg_eqtl_list[[tissue]]$human_ensembl_gene.x[male_inds], gcta_output[[tissue]]$ENSG)])
  deg_eqtl_list[[tissue]]$genetic_expression_Z[female_inds] <- deg_eqtl_list[[tissue]]$logFC[female_inds] / 
    sds_expression[[tissue]]$est_sd_f[match(deg_eqtl_list[[tissue]]$human_ensembl_gene.x[female_inds], rownames(sds_expression[[tissue]]))] / 
    sqrt(gcta_output[[tissue]]$h2[match(deg_eqtl_list[[tissue]]$human_ensembl_gene.x[female_inds], gcta_output[[tissue]]$ENSG)])
  
  #remove unexpressed genes  
  bad_inds <- which(deg_eqtl_list[[tissue]]$reference_average_intensity == 0 | deg_eqtl_list[[tissue]]$comparison_average_intensity == 0)
  deg_eqtl_list[[tissue]]$genetic_expression_Z[bad_inds] <- NA
  deg_eqtl_list[[tissue]]$phenotypic_expression_Z[bad_inds] <- NA
  
  # deg_eqtl_list[[tissue]]$genetic_expression_plus2SE_Z <- deg_eqtl_list[[tissue]]$logFC / 
  #   (sds_expression[tissue, match(deg_eqtl_list[[tissue]]$human_ensembl_gene.x, colnames(sds_expression))] * 
  #      sqrt(gcta_output[[tissue]]$h2[match(deg_eqtl_list[[tissue]]$human_ensembl_gene.x, gcta_output[[tissue]]$ENSG)] + 
  #             2 * gcta_output[[tissue]]$SE[match(deg_eqtl_list[[tissue]]$human_ensembl_gene.x, gcta_output[[tissue]]$ENSG)]))
}
save(deg_eqtl_list, file = "~/data/smontgom/relative_effect_sizes_deg_eqtl_list.RData")
load("~/data/smontgom/relative_effect_sizes_deg_eqtl_list.RData")

#remove unexpressed genes
# where either "reference_average_intensity" or "comparison_timewise_intensity" in the timewise-dea tables is 0
# load('~/data/smontgom/dea/transcript_rna_seq_20210804.RData')
# sum(transcript_rna_seq$timewise_dea$reference_average_intensity == 0 | transcript_rna_seq$timewise_dea$comparison_average_intensity == 0)
# bad_genes <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$reference_average_intensity == 0 | 
#                                              transcript_rna_seq$timewise_dea$comparison_average_intensity == 0,
#                                              c("feature_ID", "tissue", "sex", "comparison_group")]
# bad_genes <- lapply(setNames(tissues, tissues), function(tiss) bad_genes[bad_genes$tissue == tiss,])

#### compute quantile objects ####
qs2use <- 1:9999/10000
EZ_PZ <- lapply(setNames(paste0(2^(0:3), "w"),paste0(2^(0:3), "w")), function(ti) 
  sapply(tissues, function(tissue) quantile(
    x = deg_eqtl_list[[tissue]]$phenotypic_expression_Z[deg_eqtl_list[[tissue]]$comparison_group == ti], 
    probs = qs2use, na.rm = T)))
EZ_PZ <- lapply(setNames(paste0(2^(0:3), "w"),paste0(2^(0:3), "w")), function(ti) 
  sapply(tissues, function(tissue) quantile(
    x = deg_eqtl_list[[tissue]]$genetic_expression_Z[deg_eqtl_list[[tissue]]$comparison_group == ti], 
    probs = qs2use, na.rm = T)))


load(file = "~/data/smontgom/node_metadata_list.RData")
relative_expression_data <- 
  lapply(setNames(c("male", "female"),c("male", "female")), function(sex_i) list(
    
    phenotypic_expression = lapply(setNames(paste0(2^(0:3), "w"),paste0(2^(0:3), "w")), function(ti) 
      sapply(tissues, function(tissue) quantile(
        x = deg_eqtl_list[[tissue]]$phenotypic_expression_Z[deg_eqtl_list[[tissue]]$comparison_group == ti & deg_eqtl_list[[tissue]]$sex == sex_i], 
        probs = qs2use, na.rm = T))),
    
    genetic_expression = lapply(setNames(paste0(2^(0:3), "w"),paste0(2^(0:3), "w")), function(ti) 
      sapply(tissues, function(tissue) quantile(
        x = deg_eqtl_list[[tissue]]$genetic_expression_Z[deg_eqtl_list[[tissue]]$comparison_group == ti & deg_eqtl_list[[tissue]]$sex == sex_i], 
        probs = qs2use, na.rm = T))),
    
    phenotypic_expression_sexhomo = lapply(setNames(paste0(2^(0:3), "w"),paste0(2^(0:3), "w")), function(ti) 
      sapply(setdiff(tissues, c("t64-ovaries", "t63-testes", "t59-kidney")), function(tissue) quantile(
        x = deg_eqtl_list[[tissue]]$phenotypic_expression_Z[deg_eqtl_list[[tissue]]$comparison_group == ti & deg_eqtl_list[[tissue]]$sex == sex_i & 
                                                              deg_eqtl_list[[tissue]]$human_ensembl_gene.x %in% node_metadata_list[[ti]]$human_ensembl_gene], 
        probs = qs2use, na.rm = T))),
    
    genetic_expression_sexhomo = lapply(setNames(paste0(2^(0:3), "w"),paste0(2^(0:3), "w")), function(ti) 
      sapply(setdiff(tissues, c("t64-ovaries", "t63-testes", "t59-kidney")), function(tissue) quantile(
        x = deg_eqtl_list[[tissue]]$genetic_expression_Z[deg_eqtl_list[[tissue]]$comparison_group == ti & deg_eqtl_list[[tissue]]$sex == sex_i & 
                                                           deg_eqtl_list[[tissue]]$human_ensembl_gene.x %in% node_metadata_list[[ti]]$human_ensembl_gene], 
        probs = qs2use, na.rm = T)))
    
  ))
save(file = "~/data/smontgom/relative_expression_motrpac_gtex", relative_expression_data)

# tissue <- tissues[3]
# dev.off()
# plot(x = deg_eqtl_list[[tissue]]$logFC, y = deg_eqtl_list[[tissue]]$genetic_expression_Z)

#### just the plotting ####

load("~/data/smontgom/relative_expression_motrpac_gtex")

EZ_PZ <- relative_expression_data[["male"]][[1]]

ti = "8w"
f_p <- 0.4
f_x <- 2
ylims = c(min(sort(EZ_PZ[[ti]])[sort(EZ_PZ[[ti]]) != -Inf]),max(sort(EZ_PZ[[ti]],T)[sort(EZ_PZ[[ti]],T) != Inf]))
ylims <- squish_middle_x(ylims, f_x)

plot(100,100,xlim = c(0,1.25), ylim = ylims, xpd=NA, ylab = "Relative Effect Size",
     main = latex2exp::TeX("Ratio of Exercise DE to \\sqrt{Variance in log_2(Gene Expression)}"), xlab = "Quantile", xaxt = "n", yaxt = "n", bty="n", cex.lab = 1.25, cex.axis = 1.25)

xylocs_tissue_names <- cbind(rep(1.05, ncol(EZ_PZ[[ti]])), redistribute(as.numeric(squish_middle_x(tail(EZ_PZ[[ti]], 1), f_x)), diff(ylims) / 30))
rownames(xylocs_tissue_names) <- colnames(EZ_PZ[[ti]]); colnames(xylocs_tissue_names) <- c("x", "y")

#horiz axis
segments(x0 = 0, y0 = ylims[1] - diff(ylims) / 100, x1 = 1, y1 = ylims[1] - diff(ylims) / 100, lwd = 2)
segments(x0 = 0:10/10, y0 = ylims[1] - diff(ylims) / 100, x1 = 0:10/10, y1 = ylims[1] - diff(ylims) / 50, lwd = 2, xpd = NA)
horiz_axis_labels <- round(unsquish_middle_p(0:10/10, f_p), 3)
horiz_axis_labels[1] <- 0; horiz_axis_labels[length(horiz_axis_labels)] <- 1;
text(labels = horiz_axis_labels, x = 0:10/10, y = rep(ylims[1] - diff(ylims) / 50, 10), pos = 1, xpd = NA)
segments(x0 = 0.5, y0 = ylims[1], x1 = 0.5, y1 = ylims[2], lwd = 2, lty = 2, col = "grey50")
segments(x0 = 0:10/10, y0 = ylims[1], x1 = 0:10/10, y1 = ylims[2], lwd = 1, lty = 3, col = "grey75")

#vert axis
segments(x0 = -1/100, y0 = ylims[1] + diff(ylims)/100, x1 = -1/100, y1 = ylims[2], lwd = 2)
segments(x0 = -1/100, y0 = seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), x1 = -1/50,
         y1 = seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), lwd = 2)
text(labels = round(unsquish_middle_x(seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), f_x), 2),
     x = -1/50, y = seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), pos = 2, xpd = NA)
segments(x0 = 0, y0 = 0, x1 = 1, y1 = 0, lwd = 2, lty = 2, col = "grey50")
segments(x0 = 0, y0 = seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), x1 = 1, 
         y1 = seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), lwd = 1, lty = 3, col = "grey75")


for(tissue in colnames(EZ_PZ[[ti]])){
  lines(squish_middle_p(qs2use, f_p), squish_middle_x(EZ_PZ[[ti]][,tissue], f_x), col = cols$Tissue[tissue])
  text(tissue, x = xylocs_tissue_names[tissue,"x"], y = xylocs_tissue_names[tissue,"y"], pos = 4, xpd = NA,
       col = cols$Tissue[tissue])
  segments(x0 = squish_middle_p(tail(qs2use, 1), f_p), x1 = xylocs_tissue_names[tissue,"x"]+0.01,
           y0 = squish_middle_x(tail(EZ_PZ[[ti]][,tissue], 1), f_x), y1 = xylocs_tissue_names[tissue,"y"],
           lty = 3, col = cols$Tissue[tissue])
}
shadowtext(x = xylocs_tissue_names[1,"x"], y = min(xylocs_tissue_names[,"y"]) - diff(ylims) / 7.5, labels = ti, 
           cex = 5, col = cols$Time[ti], pos = 4, r = 0.2) 


#### composite plot that tells the whole story ####

grDevices::cairo_pdf(filename = paste0("~/Documents/Documents - nikolai/motrpac_companion/figures/figure2_relative-expression_composite.pdf"), 
                     width = 1400 / 72, height = 600 / 72, family="Arial Unicode MS", pointsize = 20)

#layout matrix
layout(rbind(
  c(11,8,1,1,1,9,10,10,4,4,4,11,6,6,6),
  c(11,8,1,1,1,9,10,10,4,4,4,11,6,6,6),
  c(11,2,2,11,3,3,10,10,5,5,5,11,7,7,7),
  c(11,2,2,11,3,3,10,10,5,5,5,11,7,7,7)
))



par(mar = c(4,2,2,2)+0.5)
#logFCs
xr <- c(-2,2)
x_logfc <- seq(from = xr[1], to = xr[2], length.out = 256)
logfc_dens <- sapply(deg_eqtl_list, function(del)
  density((del$logFC), na.rm = T, from = xr[1], to = xr[2], n = 256)$y)
plot(NA,NA, xlim = xr, ylim = c(0, quantile(logfc_dens, 1)), xlab = "", ylab = "",
     xaxt = "n", xpd = NA)
xvals <- seq(xr[1], xr[2], length.out = 5)
segments(y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4])/50, x0 = xvals, x1 = xvals, xpd = NA)
text(y = par("usr")[3] - diff(par("usr")[3:4])/50, pos = 1, x = xvals, labels = round((xvals), 2), xpd = NA)
text(x = par("usr")[1] - diff(par("usr")[1:2])/4.5, y = mean(par("usr")[3:4]), labels = "", srt = 90, xpd = NA)
text(x = par("usr")[1] - diff(par("usr")[1:2])/4.5, y = mean(par("usr")[3:4]), labels = "Density", srt = 90, xpd = NA)
text(x = mean(par("usr")[1:2]), y = par("usr")[3] - diff(par("usr")[3:4])/4.5, labels = latex2exp::TeX("Exercise-induced $log_2FC$"), srt = 0, xpd = NA)

for(tissue in colnames(logfc_dens)){
  polygon(c(x_logfc, rev(x_logfc)), c(logfc_dens[,tissue], rep(0,length(x_logfc))), col = adjustcolor(cols$Tissue[tissue], 0.5))  
}

#surrounding text
# text(x = -6, y = -6.5, labels = "(", cex = 16, xpd = NA, family = "Helvetica Narrow", srt = 0, col = 1)
par(xpd = NA)
arc(t1 = pi/2, t2 = 3*pi/2, r1 = 3, r2 = 3, center = c(-5.8, -7), adjx = 0.25, lwd = 3)
arc(t1 = pi/2-1E-5, t2 = 3*pi/2, r1 = 3, r2 = 3, center = c(5.3, -7), adjx = 0.25, lwd = 3)
text(x = 6, -3.75, labels = "1/2", cex = 2)
segments(x0 = -5.5, x1 = 5, y0 = -2.75, y1 = -2.75, lwd = 5, xpd = NA)
segments(x0 = 6, x1 = 7, y0 = -2.5, y1 = -2.5, lwd = 10, xpd = NA, col = 2)
points(7.25, -2.5, pch = -9658, cex = 5, xpd = NA, col = 2)
par(xpd = T)

#snp heritability
par(mar = c(5,0,2,0))
if(!exists("gcta_output")){
  load(file = "gcta_output_GTEx_allTissues_list_IHW.RData")
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
plot(h2_freqs[[tissues[1]]], col = adjustcolor(cols$Tissue[tissues[1]], 0.5), xlab = "", ylab = "Density",
     main = "")
for(tissue in tissues[-1]){
  plot(h2_freqs[[tissue]], col = adjustcolor(cols$Tissue[tissue], 0.5), add = T)
}
text(x = par("usr")[1] - diff(par("usr")[1:2])/4, y = mean(par("usr")[3:4]), labels = "Density", srt = 90, xpd = NA)
text(x = mean(par("usr")[1:2]), y = par("usr")[3] - diff(par("usr")[3:4])/4.5, labels = latex2exp::TeX("Estimated $h^2_{SNP}$"), srt = 0, xpd = NA)
box("plot")

#invgamma hyperpriors
x_var <- 1:500/1000
var_dens <- sapply(invgamma_estimates, function(invgamma_estimate)
  dinvgamma(x_var, shape = invgamma_estimate[1], scale = invgamma_estimate[2]))
colnames(var_dens) <- names(invgamma_estimates)
plot(NA,NA, xlim = range(x_var), ylim = c(0, quantile(var_dens, 1)), xlab = "", ylab = "Density")
for(tissue in colnames(var_dens)){
  polygon(c(x_var, rev(x_var)), c(var_dens[,tissue], rep(0,length(x_var))), col = adjustcolor(cols$Tissue[tissue], 0.5))  
}
text(x = par("usr")[1] - diff(par("usr")[1:2])/4, y = mean(par("usr")[3:4]), labels = "Density", srt = 90, xpd = NA)
text(x = par("usr")[1] - diff(par("usr")[1:2])/2.65, y = mean(par("usr")[3:4]), labels = "•", srt = 90, xpd = NA, cex = 4)
text(x = mean(par("usr")[1:2]), y = par("usr")[3] - diff(par("usr")[3:4])/4.5, labels = latex2exp::TeX("Expression Variance"), srt = 0, xpd = NA)

#quantile plots
par(mar = c(4,1,2,1), xpd = NA)
nnmap <- as.data.frame(bic_animal_tissue_code[,4:5])
qs2use <- round(invlogit(seq(-9,9, length.out = 75)), 4)


for(sex_i in c("male", "female")){
for(type in 3:4){
    
  EZ_PZ <- relative_expression_data[[sex_i]][[type]]
  EZ_PZ <- lapply(EZ_PZ, function(x) x[match(paste0(sprintf("%.2f", qs2use*100), "%"), paste0(sprintf("%.2f", as.numeric(gsub("%", "", rownames(x)))), "%")),])
  
  ti = "8w"
  f_p <- 0.4
  f_x <- 2
  ylims = c(min(sort(EZ_PZ[[ti]])[sort(EZ_PZ[[ti]]) != -Inf]),max(sort(EZ_PZ[[ti]],T)[sort(EZ_PZ[[ti]],T) != Inf]))
  ylims <- squish_middle_x(ylims, f_x)
  
  plot(100,100,xlim = c(0,1.25), ylim = ylims, xpd=NA, ylab = "", xlab = "", xaxt = "n", yaxt = "n", bty="n", cex.lab = 1.25, cex.axis = 1.25)
  text(x = par("usr")[1] - diff(par("usr")[1:2])/4.5, y = mean(par("usr")[3:4]), labels = "Standardized Effect Size (SD)", srt = 90, xpd = NA)
  
  text("Quantile", x = 0.5, y = ylims[1] - diff(ylims)/5, pos = 1, cex = 1)
  if(ti == "2w"){
    text(latex2exp::TeX(paste0("Ratio of Exercise DE to \\sqrt{", ifelse(type == 1, "Phenotypic", "Genetic"), " Variance in $log_2$(Gene Expression)}")), x = 1.35, y = ylims[2] + diff(ylims) / 15, pos = 3, cex = 3)}
  if(ti == "1w"){
    fig_label(ifelse(type == 1, "a)", "b)"), region = "figure", pos = "topleft", cex = 3.5)
  }
  
  xylocs_tissue_names <- cbind(rep(1.05, ncol(EZ_PZ[[ti]])), redistribute(as.numeric(squish_middle_x(tail(EZ_PZ[[ti]], 1), f_x)), diff(ylims) / 20))
  rownames(xylocs_tissue_names) <- colnames(EZ_PZ[[ti]]); colnames(xylocs_tissue_names) <- c("x", "y")
  
  #horiz axis
  segments(x0 = 0, y0 = ylims[1] - diff(ylims) / 100, x1 = 1, y1 = ylims[1] - diff(ylims) / 100, lwd = 2)
  segments(x0 = 0:10/10, y0 = ylims[1] - diff(ylims) / 100, x1 = 0:10/10, y1 = ylims[1] - diff(ylims) / 50, lwd = 2, xpd = NA)
  horiz_axis_labels <- round(unsquish_middle_p(0:10/10, f_p), 3)
  horiz_axis_labels[1] <- 0; horiz_axis_labels[length(horiz_axis_labels)] <- 1;
  text(labels = horiz_axis_labels, x = 0:10/10 + 0.035, y = rep(ylims[1] - diff(ylims) / 25, 10), pos = 2, xpd = NA, srt = 90, cex = 0.75)
  segments(x0 = 0.5, y0 = ylims[1], x1 = 0.5, y1 = ylims[2], lwd = 2, lty = 2, col = "grey50")
  segments(x0 = 0:10/10, y0 = ylims[1], x1 = 0:10/10, y1 = ylims[2], lwd = 1, lty = 3, col = "grey75")
  
  #vert axis
  segments(x0 = -1/100, y0 = ylims[1] + diff(ylims)/100, x1 = -1/100, y1 = ylims[2], lwd = 2)
  segments(x0 = -1/100, y0 = seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), x1 = -1/50,
           y1 = seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), lwd = 2)
  text(labels = round(unsquish_middle_x(seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), f_x), 2),
       x = -1/50, y = seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), pos = 2, xpd = NA, cex = 0.75)
  segments(x0 = 0, y0 = 0, x1 = 1, y1 = 0, lwd = 2, lty = 2, col = "grey50")
  segments(x0 = 0, y0 = seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), x1 = 1, 
           y1 = seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), lwd = 1, lty = 3, col = "grey75")
  
  
  for(tissue in colnames(EZ_PZ[[ti]])){
    lines(squish_middle_p(qs2use, f_p), squish_middle_x(EZ_PZ[[ti]][,tissue], f_x), col = cols$Tissue[tissue])
    text(nnmap$abbreviation[match(tissue, nnmap$tissue_name_release)], x = xylocs_tissue_names[tissue,"x"], y = xylocs_tissue_names[tissue,"y"], pos = 4, xpd = NA,
         col = cols$Tissue[tissue], cex = 0.75)
    segments(x0 = squish_middle_p(tail(qs2use, 1), f_p), x1 = xylocs_tissue_names[tissue,"x"]+0.01,
             y0 = squish_middle_x(tail(EZ_PZ[[ti]][,tissue], 1), f_x), y1 = xylocs_tissue_names[tissue,"y"],
             lty = 3, col = cols$Tissue[tissue])
  }
  text(x = -0.01, y = ylims[2] - diff(ylims) / 8, labels = ti, 
             cex = 3, col = cols$Time[ti], pos = 4)
  text(x = 0.275, y = ylims[2] - diff(ylims) / 8, labels = ",", 
       cex = 3, col = cols$Time[ti], pos = 4)
  text(labels = c(female = "\u2640", male = "\u2642")[sex_i], x = 0.35, y = ylims[2] - diff(ylims) / 8, pos = 4, cex = 3, 
       col = sex_cols[sex_i], family = "Arial Unicode MS")
  
  if(sex_i == "male"){
    text(latex2exp::TeX(paste0("Ratio of Exercise DE to \\sqrt{", 
                               ifelse(type == 3, "Phenotypic", "Genetic"), " Variance in $log_2$(Gene Expression)}")), 
         x = 1.5, y = ylims[2] + diff(ylims) / 50, pos = 3, cex = 1)
  }
}
  
  
}



dev.off()


#### histograms ####
# tissues <- colnames(relative_expression_data[["male"]][[3]][["8w"]])
# xr <- c(-4, 4)
# x_gen <- seq(from = xr[1], to = xr[2], length.out = 256)
# gen_dens <- lapply(setNames(c("male", "female"),c("male", "female")), function(sex_i)
#   sapply(setNames(tissues, tissues), function(tissue){
#     dens_input <- deg_eqtl_list[[tissue]]$genetic_expression_Z[deg_eqtl_list[[tissue]]$comparison_group == "8w" &
#                                                                  deg_eqtl_list[[tissue]]$sex == sex_i]
#     density(dens_input, na.rm = T, from = xr[1], to = xr[2], n = 256)$y
#   }
#   ))
# 
# lapply(setNames(c("male", "female"),c("male", "female")), function(sex_i)
#   sapply(setNames(tissues, tissues), function(tissue){
#     dens_input <- deg_eqtl_list[[tissue]]$genetic_expression_Z[deg_eqtl_list[[tissue]]$comparison_group == "8w" &
#                                                                  deg_eqtl_list[[tissue]]$sex == sex_i]
#     quantile(dens_input, na.rm = T)
#   }
#   ))
# 
# for(sex in c("male", "female")){
#   plot(NA,NA, xlim = xr, ylim = c(0, quantile(gen_dens[[sex]], 1)), xlab = "", ylab = "",
#        xaxt = "n", xpd = NA)
#   xvals <- seq(xr[1], xr[2], length.out = 5)
#   segments(y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4])/50, x0 = xvals, x1 = xvals, xpd = NA)
#   text(y = par("usr")[3] - diff(par("usr")[3:4])/50, pos = 1, x = xvals, labels = round((xvals), 2), xpd = NA)
#   text(x = par("usr")[1] - diff(par("usr")[1:2])/4.5, y = mean(par("usr")[3:4]), labels = "", srt = 90, xpd = NA)
#   text(x = par("usr")[1] - diff(par("usr")[1:2])/4.5, y = mean(par("usr")[3:4]), labels = "Density", srt = 90, xpd = NA)
#   text(x = mean(par("usr")[1:2]), y = par("usr")[3] - diff(par("usr")[3:4])/4.5, labels = latex2exp::TeX("Standardized Effect Size (SD)"), srt = 0, xpd = NA)
#   
#   for(tissue in colnames(gen_dens[[sex]])){
#     polygon(c(x_gen, rev(x_gen)), c(gen_dens[[sex]][,tissue], rep(0,length(x_gen))), col = adjustcolor(cols$Tissue[tissue], 0.5))
#   }
# }

# deg_eqtl_list$`t56-vastus-lateralis`$phenotypic_expression_Z
# hist(asinh(asinh(deg_eqtl_list$`t56-vastus-lateralis`$genetic_expression_Z)))
# 
# #figure out height of histogram
# histogram_range_rescaler <- 1.4 * max(1, max(hist_corrs_prots.trans$counts[1:5]) /  max(hist_corrs_prots.trans$counts) * 2)
# range_hist_corrs_prots.trans <- c(0,max(hist_corrs_prots.trans$counts)) * histogram_range_rescaler
# range_hist_corrs_trans_logFC.zscores <- c(0,max(hist_corrs_trans_logFC.zscores$counts)) * histogram_range_rescaler
# 
# 
# #plot histogram blocks
# for(bin in 1:(length(hist_corrs_prots.trans$breaks) - 1)){
#   rect(xleft = (hist_corrs_prots.trans$breaks[bin] + 1)*(1-buffer_hist)+buffer_hist,
#        xright = (hist_corrs_prots.trans$breaks[bin+1] + 1)*(1-buffer_hist)+buffer_hist,
#        # ybottom = (-lower_hist_bin_by/range_hist_corrs_prots.trans[2])*(1-buffer_hist)+buffer_hist/2,
#        ybottom = 0,
#        ytop = (hist_corrs_prots.trans$counts[bin] / range_hist_corrs_prots.trans[2])*(1-buffer_hist)+buffer_hist/2,
#        col = adjustcolor(ome_cols[5], 0.75),
#        border = NA
#   )
# }
# 
# for(gbi in 1:length(genes_in_bins)){
#   
#   #write number of genes of each sex above rectangle
#   n_male <- length(grep(pattern = "\u2642", genes_in_bins[[gbi]]))
#   n_female <- length(grep(pattern = "\u2640", genes_in_bins[[gbi]]))
#   sex_summary_text <- paste0(n_male, "\u2642,", n_female, "\u2640")
#   sex_summary_text_cols <- c(rep(1,nchar(n_male)), cols$Sex["male"], rep(1,nchar(n_female)+1), cols$Sex["female"])
#   text_cols(string = sex_summary_text, pos = 3, cex = 0.75, cols = sex_summary_text_cols,
#             x = (hist_corrs_prots.trans$mids[gbi] + 1)*(1-buffer_hist)+buffer_hist, 
#             y = (hist_corrs_prots.trans$counts[gbi] / range_hist_corrs_prots.trans[2])*(1-buffer_hist)+buffer_hist/2 - 0.05)
#   
#   if(length(genes_in_bins[[gbi]]) == 0){
#     next()
#   }
# 
#   txt <- genes_in_bins[[gbi]]
# 
#   #truncate long gene names
#   max_gene_name_length <- 10
#   short_txt <- sapply(txt, function(stxt) strsplit(stxt, "-")[[1]][-length(strsplit(stxt, "-")[[1]])])
#   short_txt <- sapply(short_txt, function(stxt) paste0(substr(stxt, 1, max_gene_name_length), ifelse(nchar(stxt) > max_gene_name_length, ".", "")))
#   txt <- sapply(1:length(txt), function(txti) paste0(short_txt[txti], "-", strsplit(txt[txti], "-")[[1]][length(strsplit(txt[txti], "-")[[1]])]))
#   
#   rect_coords <- list(x0 = (hist_corrs_prots.trans$breaks[gbi] + 1)*(1-buffer_hist)+buffer_hist,
#                       x1 = (hist_corrs_prots.trans$breaks[gbi+1] + 1)*(1-buffer_hist)+buffer_hist,
#                       y0 = (hist_corrs_prots.trans$counts[gbi] / range_hist_corrs_prots.trans[2])*(1-buffer_hist)+buffer_hist/2 - 0.01,
#                       y1 = 0*(1-buffer_hist)+buffer_hist/2)
#   optimal_word_placement_inf <- find_optimal_cex_and_lines(txt = txt, rect_coords = rect_coords, rect_rescaling_ratio = 0.925)
#   #1 is male, 2 female, 3 both
#   txt_sexes <- as.numeric(nchar(txt) - nchar(sapply(txt, function(stxt) strsplit(stxt, "-")[[1]][1])) == 3) + 2
#   txt_sexes[txt_sexes == 2] <- as.numeric(sapply(txt, function(stxt) strsplit(stxt, "-")[[1]][2])[txt_sexes == 2] == "\u2640" ) + 1
#   sex_cols_map <- list(male = cols$Sex[c("male")], female = cols$Sex[c("female")], both = cols$Sex[c("male", "female")])
#   cols_list <- lapply(1:length(txt), function(word_i) as.vector(c(rep(0, nchar(txt[word_i]) - 1 - floor(txt_sexes[word_i]/3)), unlist(sex_cols_map[txt_sexes[word_i]]))))
#   text_wrapped_words(txt, rect_coords, optimal_word_placement_inf, justified = T, col = "white", str_width_lefter_start_ratio = 0.15, rect_rescaling_ratio = 0.925,
#                      cols_list = cols_list, multicolor_words = T)
# }