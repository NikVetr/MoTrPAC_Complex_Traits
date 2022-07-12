library(MotrpacBicQC)
library(data.table)
library(EnsDb.Hsapiens.v79)

#read in rat-human mapping
rgd_orthologs <- fread("~/data/smontgom/RGD_ORTHOLOGS_20201001.txt", header = T, sep = "\t")
gencode_gene_map <- rdg_mapping <- fread("~/data/smontgom/gencode.v39.RGD.20201001.human.rat.gene.ids.txt")
gencode_gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID <- gsub(gencode_gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID, pattern = "\\..*", replacement = "")
gene_map <- gencode_gene_map
# feature_to_gene_map <- fread("~/data/smontgom/motrpac-mappings-master_feature_to_gene.txt", header = T, sep = "\t", fill = T)

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

#read in human sex-biased genes from gtex
sexbg <- read.table(file = "~/data/smontgom/GTEx_Analysis_v8_sbgenes/signif.sbgenes.txt", header = T)
sexbg$ENSG <- gsub('\\..*','',sexbg$gene)

sex_DE_humans <- lapply(setNames(setdiff(motrpac_gtex_map, c("Testis", "Ovary")), 
                               setdiff(motrpac_gtex_map, c("Testis", "Ovary"))), function(x) NULL)
for(tissue in names(sex_DE_humans)){
  print(tissue)
  sex_DE_humans[[tissue]] <- sexbg[sexbg$tissue == tissue,c("ENSG","effsize","lfsr")]
  # map <- ensembldb::select(EnsDb.Hsapiens.v79, 
  #                          keys = sex_DE_humans[[tissue]]$ENSG, 
  #                          keytype = "GENEID", columns = c("GENEID", "SYMBOL"))
  # sex_DE_humans[[tissue]]$HUMAN_ORTHOLOG_SYMBOL <- map$SYMBOL[match(sex_DE_humans[[tissue]]$ENSG, map$GENEID)]
  sex_DE_humans[[tissue]]$HUMAN_ORTHOLOG_SYMBOL <- gene_map$HUMAN_ORTHOLOG_SYMBOL[match(sex_DE_humans[[tissue]]$ENSG, gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID)]
}
names(sex_DE_humans) <- names(motrpac_gtex_map)[match(names(sex_DE_humans), motrpac_gtex_map)]
sex_DE_humans[[length(sex_DE_humans) + 1]] <- sex_DE_humans$`t55-gastrocnemius`
names(sex_DE_humans)[length(sex_DE_humans)] <- "t56-vastus-lateralis"

#read in motrpac data
bic_animal_tissue_code = as.data.table(MotrpacBicQC::bic_animal_tissue_code)
bic_animal_tissue_code = bic_animal_tissue_code[tissue_name_release!='']
bic_animal_tissue_code[,my_tissue := tolower(gsub(" Powder", "", bic_tissue_name))]
bic_animal_tissue_code[,my_tissue := gsub(' ','_',my_tissue)]
bic_animal_tissue_code[my_tissue == 'blood_rna', my_tissue := 'paxgene_rna']
tissue_codes = bic_animal_tissue_code[, tissue_name_release]
tissue_codes = tissue_codes[tissue_codes != '']

sex_DE_rats <- lapply(setNames(setdiff(names(motrpac_gtex_map), c("t63-testes", "t64-ovaries")), 
                               setdiff(names(motrpac_gtex_map), c("t63-testes", "t64-ovaries"))), function(x) NULL)
for(tissue in setdiff(names(motrpac_gtex_map), c("t63-testes", "t64-ovaries"))){
  print(tissue)
  load(paste0("~/data/smontgom/old_dea_deseq_20201121/", tissue, "_training-dea_20201121.RData"))
  p_values <- pnorm(dds@rowRanges@elementMetadata@listData$sex_male_vs_female / 
                    dds@rowRanges@elementMetadata@listData$SE_sex_male_vs_female)
  p_values[p_values > 0.5] <- 1-p_values[p_values > 0.5]
  p_values <- p_values*2
  sex_DE_rats[[tissue]] <-data.table(ENSRNOG = rownames(dds),
                                     p_value = p_values,
                                     effect_size = dds@rowRanges@elementMetadata@listData$sex_male_vs_female)
}

n_entries <- sapply(names(sex_DE_rats), function(tissue) nrow(sex_DE_rats[[tissue]]))
n_comparisons <- sum(n_entries)
fdr_pvals <- p.adjust(unlist(sapply(names(sex_DE_rats), function(tissue) sex_DE_rats[[tissue]]$p_value)), "fdr")

ihw_results <- IHW::ihw(pval ~ tissue,
                        data = data.frame(pval = unlist(sapply(names(sex_DE_rats), function(tissue) sex_DE_rats[[tissue]]$p_value)),
                                          tissue = factor(unlist(sapply(names(sex_DE_rats), function(tissue) rep(tissue, length(sex_DE_rats[[tissue]]$p_value)))))),
                        alpha = 0.05)
fdr_pvals <- ihw_results@df$adj_pvalue

alpha = 0.05
for(tissue in names(sex_DE_rats)){
  print(tissue)
  
  #snag p-vals
  sex_DE_rats[[tissue]]$bonferroni_significant <- log(sex_DE_rats[[tissue]]$p_value) < (log(alpha) - log(n_comparisons))
  sex_DE_rats[[tissue]]$fdr_significant <- fdr_pvals[(cumsum(c(0,n_entries))[match(tissue, names(sex_DE_rats))] + 1):(cumsum(c(0,n_entries))[match(tissue, names(sex_DE_rats)) + 1])] < alpha
  
  #map names and subset
  sex_DE_rats[[tissue]]$RAT_GENE_SYMBOL <- feature_to_gene_map$gene_symbol[match(sex_DE_rats[[tissue]]$ENSRNOG, feature_to_gene_map$ensembl_gene)]
  sex_DE_rats[[tissue]]$HUMAN_ORTHOLOG_SYMBOL <- rgd_orthologs$HUMAN_ORTHOLOG_SYMBOL[match(sex_DE_rats[[tissue]]$RAT_GENE_SYMBOL, rgd_orthologs$RAT_GENE_SYMBOL)]
  # sex_DE_rats[[tissue]] <- sex_DE_rats[[tissue]][sex_DE_rats[[tissue]]$bonferroni_significant,]
  sex_DE_rats[[tissue]] <- sex_DE_rats[[tissue]][sex_DE_rats[[tissue]]$fdr_significant,]
}
sapply(names(sex_DE_rats), function(tissue) sum(sex_DE_rats[[tissue]]$bonferroni_significant))
sapply(names(sex_DE_rats), function(tissue) sum(sex_DE_rats[[tissue]]$fdr_significant))

#### do plotting ####

grDevices::cairo_pdf(filename = paste0("~/Documents/Documents - nikolai/motrpac_companion/figures/rat-man_sex_comparison.pdf"), 
                     width = 900 / 72, height = 1050 / 72, family="Arial Unicode MS", pointsize = 19)

layout(rbind(matrix(1:15, 3, 5, T), matrix(16, 3, 5, T)))

for(tissue in setdiff(names(motrpac_gtex_map), c("t63-testes", "t64-ovaries"))){
  
  if(match(tissue, setdiff(names(motrpac_gtex_map), c("t63-testes", "t64-ovaries"))) %% 5 == 1){
    par(mar = c(2,4,2,0))
  } else if(match(tissue, setdiff(names(motrpac_gtex_map), c("t63-testes", "t64-ovaries"))) %% 5 == 2){
    par(mar = c(2,3.75,2,0.25))  
  } else if(match(tissue, setdiff(names(motrpac_gtex_map), c("t63-testes", "t64-ovaries"))) %% 5 == 3){
    par(mar = c(2,3.5,2,0.5))  
  } else if(match(tissue, setdiff(names(motrpac_gtex_map), c("t63-testes", "t64-ovaries"))) %% 5 == 4){
    par(mar = c(2,3.25,2,0.75))  
  } else if(match(tissue, setdiff(names(motrpac_gtex_map), c("t63-testes", "t64-ovaries"))) %% 5 == 0){
    par(mar = c(2,3,2,1))  
  }
  
  # tissue <- names(motrpac_gtex_map)[2]
  gene_overlap <- intersect(sex_DE_rats[[tissue]]$HUMAN_ORTHOLOG_SYMBOL, sex_DE_humans[[tissue]]$HUMAN_ORTHOLOG_SYMBOL)
  if(length(gene_overlap) == 0){
    plot(0,0,col="white", main = tissue, xlab = "rat DE", ylab = "human DE")
  } else {
    xlims = range(sex_DE_rats[[tissue]]$effect_size[match(gene_overlap, sex_DE_rats[[tissue]]$HUMAN_ORTHOLOG_SYMBOL)])
    if(xlims[1] > 0){xlims[1] <- -1}
    if(xlims[2] < 0){xlims[1] <- 1}
    ylims = range(-sex_DE_humans[[tissue]]$effsize[match(gene_overlap, sex_DE_humans[[tissue]]$HUMAN_ORTHOLOG_SYMBOL)])
    if(ylims[1] > 0){ylims[1] <- -1}
    if(ylims[2] < 0){ylims[1] <- 1}
    #in GTEx, 1=Male	2=Female, https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_Analysis_v7_Annotations_SubjectPhenotypesDD.xlsx
    #so it's the opposite of the MoTrPAC encoding
    hax <- sex_DE_rats[[tissue]]$effect_size[match(gene_overlap, sex_DE_rats[[tissue]]$HUMAN_ORTHOLOG_SYMBOL)]
    vax <- -sex_DE_humans[[tissue]]$effsize[match(gene_overlap, sex_DE_humans[[tissue]]$HUMAN_ORTHOLOG_SYMBOL)]
    plot(hax, vax, xpd = F,
       main = "", xlab = "", ylab = "", pch = 19, col = adjustcolor(1, 0.4),
       xlim = xlims, ylim = ylims)
    if(tissue == "t56-vastus-lateralis"){ 
      text(latex2exp::TeX("Human Differential Expression (inverse normal transform)"), 
           x = par("usr")[1] - diff(range(par("usr")[1:2])) / 2.25, 
           y = mean(par("usr")[3:4]),
           cex = 1.5, pos = 1, xpd = NA, srt = 90)
      }
    if(tissue == "t67-small-intestine"){
      text(latex2exp::TeX("Rat Differential Expression ($log_2FC$)"), x = mean(par("usr")[1:2]), 
           y = par("usr")[3] - diff(range(par("usr")[3:4])) / 3,
           cex = 1.5, pos = 1, xpd = NA)
    }
    title(tissue, line = 0.5, xpd = NA)
    # x1, x2, y1, y2
    xdisp <- diff(range(hax)) / 10
    ydisp <- diff(range(vax)) / 10
    text(par("usr")[1] + xdisp, par("usr")[3] + ydisp, round(sum(hax < 0 & vax < 0) / length(hax) * 100), 
         col = adjustcolor("darkgreen", 0.8), cex = 1.5, font = 2) #bl
    text(par("usr")[1] + xdisp, par("usr")[4] - ydisp, round(sum(hax < 0 & vax > 0) / length(hax) * 100), 
         col = adjustcolor("darkred", 0.8), cex = 1.5, font = 2) #tl
    text(par("usr")[2] - xdisp, par("usr")[3] + ydisp, round(sum(hax > 0 & vax < 0) / length(hax) * 100), 
         col = adjustcolor("darkred", 0.8), cex = 1.5, font = 2) #br
    text(par("usr")[2] - xdisp, par("usr")[4] - ydisp, round(sum(hax > 0 & vax > 0) / length(hax) * 100), 
         col = adjustcolor("darkgreen", 0.8), cex = 1.5, font = 2) #tr
    
    abline(h = 0, col = adjustcolor("gray50", 0.5), lty = 2, xpd = F)
    abline(v = 0, col = adjustcolor("gray50", 0.5), lty = 2, xpd = F)
  }
}

#do boxplot

prop_interval <- .95
prop_same_dir_flatbeta <- sapply(setdiff(names(motrpac_gtex_map), c("t63-testes", "t64-ovaries")), function(tissue) {
  gene_overlap <- intersect(sex_DE_rats[[tissue]]$HUMAN_ORTHOLOG_SYMBOL, sex_DE_humans[[tissue]]$HUMAN_ORTHOLOG_SYMBOL)
  hax <- sex_DE_rats[[tissue]]$effect_size[match(gene_overlap, sex_DE_rats[[tissue]]$HUMAN_ORTHOLOG_SYMBOL)]
  vax <- -sex_DE_humans[[tissue]]$effsize[match(gene_overlap, sex_DE_humans[[tissue]]$HUMAN_ORTHOLOG_SYMBOL)]
  c((sum(hax < 0 & vax < 0) + sum(hax > 0 & vax > 0)) / length(hax),
    qbeta(p = c((1-prop_interval)/2, 1-(1-prop_interval)/2), 
          shape1 = 1 + sum(hax < 0 & vax < 0) + sum(hax > 0 & vax > 0), 
          shape2 = 1 + (length(hax) - (sum(hax < 0 & vax < 0) + sum(hax > 0 & vax > 0)))))
})

same_direction <- data.frame(t(sapply(setdiff(names(motrpac_gtex_map), c("t63-testes", "t64-ovaries")), function(tissue) {
  gene_overlap <- intersect(sex_DE_rats[[tissue]]$HUMAN_ORTHOLOG_SYMBOL, sex_DE_humans[[tissue]]$HUMAN_ORTHOLOG_SYMBOL)
  hax <- sex_DE_rats[[tissue]]$effect_size[match(gene_overlap, sex_DE_rats[[tissue]]$HUMAN_ORTHOLOG_SYMBOL)]
  vax <- -sex_DE_humans[[tissue]]$effsize[match(gene_overlap, sex_DE_humans[[tissue]]$HUMAN_ORTHOLOG_SYMBOL)]
  c(n_same_direction = (sum(hax < 0 & vax < 0) + sum(hax > 0 & vax > 0)), n_genes = length(hax))
})))

#what about a hiearchical model?
library(cmdstanr)
library(posterior)
library(caret)
library(MASS)

stan_program <- "
data { 
  int<lower=0> n_tissues;
  int<lower=0> n_genes[n_tissues];
  int<lower=0> n_same_direction[n_tissues];
} 
parameters { 
  real<lower=0, upper=1> mean_prob;
  real concentration_offset_log;
  vector<lower=0, upper=1>[n_tissues] prob_same_direction;
} 
transformed parameters {
  real<lower=1> concentration = pow(10, concentration_offset_log) + 1;
}
model { 
  mean_prob ~ beta(1,1);
  concentration_offset_log ~ normal(0, 1);
  prob_same_direction ~ beta(mean_prob * concentration, (1 - mean_prob) * concentration);
  n_same_direction ~ binomial(n_genes, prob_same_direction);
} 
"

d <- list(n_tissues = nrow(same_direction),
          n_genes = same_direction$n_genes,
          n_same_direction = same_direction$n_same_direction)
if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)
}
mod <- cmdstan_model(f)

#fit model
# out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = d, parallel_chains = 4, adapt_delta = 0.95, refresh = 10, init = lapply(1:4, function(x) list(sigma = rep(3,n_dim))))
out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = d, parallel_chains = 4, adapt_delta = 0.9, refresh = 1000, init = 0.1, max_treedepth = 20)
# out <- mod$variational(data = d)
summ <- out$summary()
summ[order(summ$ess_bulk),]
samps <- data.frame(as_draws_df(out$draws()))
relev_samps <- samps[,grep("prob_same_direction", colnames(samps))]
prop_same_dir_hier <- rbind(mean = apply(relev_samps, 2, mean),
                                   apply(relev_samps, 2, quantile, 
                                    p = c((1-prop_interval)/2, 1-(1-prop_interval)/2)))
colnames(prop_same_dir_flatbeta) <- bic_animal_tissue_code$abbreviation[match(colnames(prop_same_dir_flatbeta), 
                                                                              bic_animal_tissue_code$tissue_name_release)]
colnames(prop_same_dir_hier) <- colnames(relev_samps) <- colnames(prop_same_dir_flatbeta)
write.table(relev_samps, "~/data/smontgom/ratman_sex_comparison_samps.txt")

# prop_same_dir <- prop_same_dir_flatbeta
prop_same_dir <- prop_same_dir_hier
par(mar = c(8,5,5,1))
bp <- barplot(prop_same_dir[1,], las = 2, ylim = c(0,1), add = F, cex.names = 1.5, cex.axis = 1.2,
        ylab = "Proportion Same Direction", cex.lab = 1.25, col = "white", border = "white", cex.lab = 1.5)
bp <- c(bp)
text(x = length(prop_same_dir) / 2, y = -0.35, labels = "Tissue", xpd = NA, pos = 4, cex = 1.25)
abline(h = 0.5, lty = 1, col = adjustcolor("darkgreen", 1), lwd = 2, xpd = F)
# barplot(prop_same_dir[1,], las = 2, ylim = c(0,1), add = T, ylab = "Proportion Positive Direction")
subset_tissues <- setdiff(names(motrpac_gtex_map), c("t63-testes", "t64-ovaries"))
for(i in 1:length(subset_tissues)){
  rect(xleft = bp[i] - 0.5, xright = bp[i] + 0.5, ybottom = prop_same_dir[2,i], ytop = prop_same_dir[3,i], col = "lightgrey", lwd = 2.5)
  segments(x0 = bp[i] - 0.5, x1 = bp[i] + 0.5, y0 = prop_same_dir[1,i], y1 = prop_same_dir[1,i], lwd = 4)
  segments(x0 = bp[i], x1 = bp[i], y0 = -0.03, y1 = prop_same_dir[2,i], lwd = 1.5, lty = 2, col = "grey70", xpd = NA)
}
abline(h = 0.5, lty = 3, col = adjustcolor("darkgreen", 1), lwd = 1.5, xpd = F)
title("Boxes Represent Posterior Means + 95% Credible Intervals", cex.main = 2, line = 0)

dev.off()
