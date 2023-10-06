#### compile data ####
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
# plot(data$count_compl / data$total_compl, data$count_inters / data$total_inters,
#      main = paste0("prop above line = ", round(mean((data$count_inters / data$total_inters > data$count_compl / data$total_compl), na.rm = T), 2)),
#      cex = data$total_inters / max(data$total_inters) * 5, pch = 19, col = adjustcolor(1,0.5))
# abline(0,1,lwd=2,lty=2,col=adjustcolor(2,0.5))
# abline(h=0.5,lwd=2,lty=2,col=adjustcolor(1,0.5))
# plot(data$count_all / data$total_all, data$count_inters / data$total_inters,
#      main = paste0("prop above line = ", round(mean((data$count_inters / data$total_inters > data$count_all / data$total_all), na.rm = T), 2)),
#      cex = data$total_inters / max(data$total_inters) * 5, pch = 19, col = adjustcolor(1,0.5))
# abline(0,1,lwd=2,lty=2,col=adjustcolor(2,0.5))
# abline(h=0.5,lwd=2,lty=2,col=adjustcolor(1,0.5))

#make slightly fancier plot
category_shapes <- setNames(15:19, c("Cardiometabolic", "Aging", "Anthropometric", 
                                     "Immune", "Psychiatric-neurologic")
)
focal_traitcats <- c("Cardiometabolic", "Immune", "Endocrine system", "Anthropometric", "Allergy", "Aging")
category_shapes <- setNames(c(24,19,15,43,23,25)[1:length(focal_traitcats)], focal_traitcats)

pchs <- category_shapes[traitwise_partitions$Category[match(data$trait, traitwise_partitions$Tag)]]
pchs[is.na(pchs)] <- 1
# plot(data$count_all / data$total_all, data$count_inters / data$total_inters,
#      main = "", 
#      pch = pchs,
#      cex = data$total_inters / max(data$total_inters) * 5, col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[data$tissue],0.5))
# abline(h=0.5,lwd=2,lty=2,col=adjustcolor(1,0.5))
# abline(v=0.5,lwd=2,lty=2,col=adjustcolor(1,0.5))


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

# dev.off()
# hist(apply(subset_samps("cell_total_prob_bias", c("raw", "sd"), samps = samps), 2, mean))
# hist(apply(subset_samps("col_bias", c("raw", "sd"), samps = samps), 2, mean), breaks = 10)
# hist(samps$col_sd_comp)
# hist(invlogit(apply(subset_samps("col_means_compl", c("raw", "sd"), samps = samps), 2, mean)))
# hist(invlogit(apply(subset_samps("logodds_compl", c("raw", "sd"), samps = samps), 2, mean)))
# hist(invlogit(apply(subset_samps("logodds_focal", c("raw", "sd"), samps = samps), 2, mean)))
# hist(invlogit(apply(subset_samps("logodds", c("raw", "sd"), samps = samps), 2, mean)))

# hist(apply(subset_samps("logodds", c("raw", "compl", "sd", "mean", "bias"), samps = samps), 2, mean))
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
# plot(d$count / d$total, apply(subset_samps("cell_total_prob_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0))
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
# hist(samps$prop_gcor, probability = T)
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


#### evaluate directional proportion ####
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


#### load graphical parameters ####

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

