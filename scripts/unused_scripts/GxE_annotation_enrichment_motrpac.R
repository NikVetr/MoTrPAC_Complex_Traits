library(data.table)
library(ivs)
library(foreach)
library(doParallel)
library(parallel)
library(MotrpacBicQC)

#initialize parallelization
if(!exists("cl")){
  cl <- makeCluster(8, outfile="")
  registerDoParallel(cl)
}
getDoParWorkers()

#load all genes tested
source(file = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/load_deg-eqtl_merged_file.R")
tissues <- setNames(names(deg_eqtl_list), names(deg_eqtl_list))
tissue_abbr <- MotrpacBicQC::tissue_abbr
tissue_abbr_rev <- setNames(names(tissue_abbr[tissues]), tissue_abbr[tissues])
motrpac_expressed_genes <- lapply(tissues, function(ti) unique(deg_eqtl_list[[ti]]$human_ensembl_gene))

#load gene map
gene_map <- fread("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/gencode.v39.RGD.20201001.human.rat.gene.ids.txt")
gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID <- gsub(gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID, pattern = "\\..*", replacement = "")

#load gene sets
if(!exists("cluster_membership") | !exists("node_metadata")){
  # load("~/data/dea_clustering_0.1-FDR-ftest_kmeans-15.RData")
  # cluster_membership <- cluster_membership[cluster_membership$cluster %in% c(1,3,7,15) & cluster_membership$assay == "TRNSCRPT",]
  # map = fread('/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/pass1b-06_transcript-rna-seq_feature-mapping_20210721.txt', sep='\t', header=T)
  # cluster_genes <- map$human_gene_symbol[match(cluster_membership$feature_ID, map$feature_ID)]
  # cluster_genes <- cluster_genes[!is.na(cluster_genes)]
  
  load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/graphical_analysis_results_20211220.RData")
  nodes_to_look_at_list <- list(c("1w_F1_M1", "1w_F-1_M-1"),
                                c("2w_F1_M1", "2w_F-1_M-1"),
                                c("4w_F1_M1", "4w_F-1_M-1"),
                                c("8w_F1_M1", "8w_F-1_M-1"))
  
  node_metadata_list <- lapply(setNames(nodes_to_look_at_list, paste0(2^(0:3), "w")), function(nodes_to_look_at){
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
}

#load other info
ENSG_coord <- read.table(file = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/ENSG_coord.txt", header = T)
rownames(ENSG_coord) <- ENSG_coord$GENE
RSID_POS_MAP <- fread("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/RSID_POS_MAP.txt")
RSID_POS_MAP <- split(RSID_POS_MAP, RSID_POS_MAP$`#CHROM`)

#load GxE results
GxE <- fread("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/gwas.chrALL.gxe.bmi.all.glm.linear.ADDxPA.hg38.txt")
rownames(GxE) <- GxE$V3
GxE$V1 <- gsub(GxE$V1, pattern = "_.*", replacement = "")
GxE <- split(GxE, GxE$V1)

SNP_frx <- fread("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ukbb.ALL.impute.frqx")
SNP_frx <- split(SNP_frx, SNP_frx$CHR)

#subset to focal genes
tissue_nodes <- split(node_metadata_list$`8w`, node_metadata_list$`8w`$tissue)
tissue_effects_list <- mclapply(tissue_abbr[tissues], function(tissue){
  
  focal_genes <- tissue_nodes[[tissue]]$human_ensembl_gene[tissue_nodes[[tissue]]$human_ensembl_gene %in% ENSG_coord$GENE]
  focal_genes <- focal_genes[!is.na(focal_genes)]
  alt_genes <- intersect(ENSG_coord$GENE, setdiff(motrpac_expressed_genes[[tissue_abbr_rev[tissue]]], focal_genes))
  
  #sanity check null model
  # focal_genes <- sample(c(focal_genes, alt_genes), size = length(focal_genes), replace = F)
  # alt_genes <- intersect(ENSG_coord$GENE, setdiff(motrpac_expressed_genes[[tissue_abbr_rev[tissue]]], focal_genes))
  
  focal_gene_coords <- ENSG_coord[focal_genes,]
  focal_gene_coords$CHR <- as.numeric(gsub("chr", "", focal_gene_coords$CHR))
  focal_gene_coords <- focal_gene_coords[order(focal_gene_coords$CHR),]
  focal_gene_coords <- split(focal_gene_coords, focal_gene_coords$CHR)
  
  alt_gene_coords <- ENSG_coord[alt_genes,]
  alt_gene_coords$CHR <- as.numeric(gsub("chr", "", alt_gene_coords$CHR))
  alt_gene_coords <- alt_gene_coords[order(alt_gene_coords$CHR),]
  alt_gene_coords <- split(alt_gene_coords, alt_gene_coords$CHR)
  
  window_size <- 1E5
  
  cris <- setNames(names(focal_gene_coords),names(focal_gene_coords))
  effects_list <- lapply(cris, function(cri){
    
    system(sprintf('echo "%s"', paste0(paste0("(", tissue, ", ", cri, ") "), collapse="")))
    
    focal_locs <- focal_gene_coords[[cri]]
    focal_locs$START <- focal_locs$START - window_size
    focal_locs$END <- focal_locs$END - window_size
    
    alt_locs <- alt_gene_coords[[cri]]
    alt_locs$START <- alt_locs$START - window_size
    alt_locs$END <- alt_locs$END - window_size
    
    #which is faster?
    # sum(RSID_POS_MAP[[cri]]$POS %inrange% data.table(focal_locs[,c("START", "END")]))
    # sum(iv_between(RSID_POS_MAP[[cri]]$POS, iv(focal_locs$START, focal_locs$END)))
    # 
    # sum(RSID_POS_MAP[[cri]]$POS %inrange% data.table(alt_locs[,c("START", "END")]))
    # sum(iv_between(RSID_POS_MAP[[cri]]$POS, iv(alt_locs$START, alt_locs$END)))
    
    focal_rsids <- RSID_POS_MAP[[cri]]$ID[RSID_POS_MAP[[cri]]$POS %inrange% data.table(focal_locs[,c("START", "END")])]
    alt_rsids <- RSID_POS_MAP[[cri]]$ID[RSID_POS_MAP[[cri]]$POS %inrange% data.table(alt_locs[,c("START", "END")])]
    alt_rsids <- setdiff(alt_rsids, focal_rsids)
    
    focal_effect_inds <- match(focal_rsids, GxE[[cri]]$V3)
    focal_effect_inds <- focal_effect_inds[!is.na(focal_effect_inds)]
    focal_effects <- GxE[[cri]]$V9[focal_effect_inds]
    focal_effect_freqs <- SNP_frx[[cri]]$MAF[match(GxE[[cri]]$V3[focal_effect_inds], SNP_frx[[cri]]$SNP)]
    focal_effects_sd <- sqrt(2 * focal_effect_freqs * (1-focal_effect_freqs))
    focal_effects <- focal_effects / focal_effects_sd
    
    alt_effect_inds <- match(alt_rsids, GxE[[cri]]$V3)
    alt_effect_inds <- alt_effect_inds[!is.na(alt_effect_inds)]
    alt_effects <- GxE[[cri]]$V9[alt_effect_inds]
    alt_effect_freqs <- SNP_frx[[cri]]$MAF[match(GxE[[cri]]$V3[alt_effect_inds], SNP_frx[[cri]]$SNP)]
    alt_effects_sd <- sqrt(2 * alt_effect_freqs * (1-alt_effect_freqs))
    alt_effects <- alt_effects / alt_effects_sd
    
    data.frame(effect = c(focal_effects, alt_effects), group = as.factor(c(rep(1, length(focal_effects)), rep(2, length(alt_effects)))), cri = cri)
    
  })
  
  effects <- do.call(rbind, effects_list)
  # effects <- split(effects, effects$group)
  return(effects)
  
}, mc.cores = 4)

names(tissue_effects_list) <- tissue_abbr[names(tissue_effects_list)]

ftest_results <- lapply(setNames(tissue_abbr[tissues], tissue_abbr[tissues]), function(tissue){
  if(is.null(tissue_effects_list[[tissue]])){return(NULL)}
  print(tissue)
  
  focal_effects <- tissue_effects_list[[tissue]]$effect[tissue_effects_list[[tissue]]$group == 1]
  alt_effects <- tissue_effects_list[[tissue]]$effect[tissue_effects_list[[tissue]]$group == 2]
  
  #sanity test
  # focal_effects_inds <- sample(1:length(c(focal_effects, alt_effects)), size = length(focal_effects), replace = F)
  # focal_effects <- c(focal_effects, alt_effects)[focal_effects_inds]
  # alt_effects <- c(focal_effects, alt_effects)[-focal_effects_inds]
  
  out_ftest <- var.test(focal_effects, alt_effects, alternative = "two.sided")
  # out_levene <- car::leveneTest(effect ~ group, data = data.frame(effect = c(focal_effects, alt_effects), group = as.factor(c(rep(1, length(focal_effects)), rep(2, length(alt_effects))))))
  
  log10pval <- (pf(q = out_ftest$statistic, df1 = out_ftest$parameter[1], df2 = out_ftest$parameter[2], log.p = T, 
                   lower.tail = ifelse(out_ftest$statistic > 1, F, T)) + log(2)) / log(10)
  
  #what if we set df to the # genes?
  focal_genes <- tissue_nodes[[tissue]]$human_ensembl_gene[tissue_nodes[[tissue]]$human_ensembl_gene %in% ENSG_coord$GENE]
  focal_genes <- focal_genes[!is.na(focal_genes)]
  alt_genes <- intersect(ENSG_coord$GENE, setdiff(motrpac_expressed_genes[[tissue_abbr_rev[tissue]]], focal_genes))
  log10pval_adjDF <- (pf(q = out_ftest$statistic, 
                         df1 = length(focal_genes), 
                         df2 = length(alt_genes), 
                         log.p = T, 
                   lower.tail = ifelse(out_ftest$statistic > 1, F, T)) + log(2)) / log(10)
  
  x <- rnorm(length(focal_genes)); x <- x / sd(x) * sd(focal_effects)
  y <- rnorm(length(alt_genes)); y <- y / sd(y) * sd(alt_effects)
  fakeout_ftest <- var.test(x, y, alternative = "two.sided")
  
  #other tests
  out_ttest <- t.test(x = abs(focal_effects), y = abs(alt_effects))
  out_kstest <- ks.test(x = focal_effects, y = alt_effects)
  
  data.frame(pval = out_ftest$p.value, log10pval = log10pval, log10pval_adjDF = log10pval_adjDF, ratio = out_ftest$estimate, lower95 = out_ftest$conf.int[1],  upper95 = out_ftest$conf.int[2], tissue = tissue,
             ttest.pval = out_ttest$p.value, ttest.relative_effect = -diff(out_ttest$estimate) / out_ttest$estimate[2], kstest.pval = out_kstest$p.value,
             lower95_adj = fakeout_ftest$conf.int[1], upper95_adj = fakeout_ftest$conf.int[2])
  
})

ftest_results <- do.call(rbind, ftest_results)
ftest_results <- ftest_results[order(ftest_results$ratio),]

plot(ftest_results$log10pval, log10(ftest_results$ttest.pval))
abline(0,1)

#### plot the results ####
# par(mfrow = c(1,2))
max_ptcex <- 3
pow_ptcex <- 0.65
assume_SNP_independence <- F

if(assume_SNP_independence){
  ptcex <- -ftest_results$log10pval
  ptcex <- ptcex - min(-ftest_results$log10pval)
  ptcex <- ptcex / max(-ftest_results$log10pval)
} else {
  ptcex <- -ftest_results$log10pval_adjDF
  ptcex <- ptcex - min(-ftest_results$log10pval_adjDF)
  ptcex <- ptcex / max(-ftest_results$log10pval_adjDF)
}
ptcex <- ptcex ^ pow_ptcex
ptcex <- ptcex * (max_ptcex-1) + 1

if(assume_SNP_independence){
  legend_vals <- 2^(0:ceiling(log2(floor(max(-ftest_results$log10pval)))))
  legend_cex <- legend_vals
  legend_cex <- legend_cex - min(-ftest_results$log10pval)
  legend_cex <- legend_cex / max(-ftest_results$log10pval)
} else {
  legend_vals <- 2^(0:ceiling(log2(floor(max(-ftest_results$log10pval_adjDF)))))
  legend_cex <- legend_vals
  legend_cex <- legend_cex - min(-ftest_results$log10pval_adjDF)
  legend_cex <- legend_cex / max(-ftest_results$log10pval_adjDF)
}
legend_cex <- legend_cex ^ pow_ptcex
legend_cex <- legend_cex * (max_ptcex-1) + 1


if(assume_SNP_independence){
  ylims <- range(c(ftest_results$lower95, ftest_results$upper95))
} else {
  ylims <- range(c(ftest_results$lower95_adj, ftest_results$upper95_adj))
}

par(mar = c(6,5,2,3))

plot(NULL, xlim = c(1, nrow(ftest_results)), ylim = ylims, xlab = "", ylab = "Ratio of Focal Variance / Alt Variance", xaxt = "n",
     main = ifelse(assume_SNP_independence, "each SNP is assumed independent ('df's = # SNPs)", "F-test 'df's are set to the # of genes"))
abline(h = 1, col = "grey50", lty = 2, lwd = 2)
abline(v = 1:nrow(ftest_results), col = "grey50", lty = 3, lwd = 1)

points(ftest_results$ratio, col = tissue_cols[ftest_results$tissue], pch = 19, cex = ptcex)
if(assume_SNP_independence){
  segments(x0 = 1:nrow(ftest_results), x1 = 1:nrow(ftest_results), y0 = ftest_results$lower95, y1 = ftest_results$upper95, col = tissue_cols[ftest_results$tissue], lwd = 2)
} else {
  segments(x0 = 1:nrow(ftest_results), x1 = 1:nrow(ftest_results), y0 = ftest_results$lower95_adj, y1 = ftest_results$upper95_adj, col = tissue_cols[ftest_results$tissue], lwd = 2)
}

legend_labs <- sapply(-legend_vals, function(lv) latex2exp::TeX(paste0("$10^{", lv, "}$")))
pval_legend <- legend(x = "topleft", legend = legend_labs, pch = 19, pt.cex = legend_cex, y.intersp = 1.2, 
                      col = "grey50", box.lty = 2, title = "p-value", plot = FALSE)
legend(x = "topleft", legend = legend_labs, pch = 19, pt.cex = legend_cex, y.intersp = 1.2, 
       col = "grey50", box.lty = 2, title = "p-value")

legend(x = pval_legend$rect$left + pval_legend$rect$w, y = par("usr")[4], 
       legend = c("point estimate", "95% CI"), lwd = c(NA, 2), pch = c(19, NA), col = "grey50", box.lty = 2)
# text(x = par("usr")[2], y = mean(ylims) + cumsum(legend_cex) / max(cumsum(legend_cex)) * diff(ylims) / 2, latex2exp::TeX(paste0("$10^{", -legend_vals, "}$")), pos = 4, xpd = NA)

# hax
segments(x0 = 1:nrow(ftest_results), x1 = 1:nrow(ftest_results), y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4]) / 50, col = 1, xpd = NA)
text(x = 1:nrow(ftest_results) + 0.2, y = par("usr")[3] - diff(par("usr")[3:4]) / 30, labels = ftest_results$tissue, xpd = NA, pos = 2, srt = 45)

#PROBLEM: tests are very miscalibrated in the presence of LD

#### let's try running an ldsc h^2_SNP enrichment analysis ####

#get sumstants in appropriate format for munge_sumstats.py
ldsc_directory <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/"
geneset_directory <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets/"
command_changedir <- paste0("source ~/.bash_profile; cd ", ldsc_directory, "; ")
command_initiateLDSC <- "source activate ldsc; "

GxE <- fread("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/gwas.chrALL.gxe.bmi.all.glm.linear.ADDxPA.hg38.txt")
SNP_frx <- fread("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ukbb.ALL.impute.frqx")
GxE$MAF <- SNP_frx$MAF[match(GxE$V3, SNP_frx$SNP)]

sumstats_file_subset <- GxE[,c("V3", "V4", "V5", "MAF", "V12", "V8")]
hist(log10(sumstats_file_subset$p), xlim = c(-10,0)); abline(v = log10(0.05 / length(sumstats_file_subset$p)), lwd = 2, col = 2)
colnames(sumstats_file_subset) <- c("MarkerName", "Allele1", "Allele2", "Freq.Allele1.HapMapCEU", "p", "N")
fwrite(sumstats_file_subset, file = paste0(ldsc_directory, "gwas_sumstats/", "BMI_GxE.txt.gz"), sep = " ")

command_ldsc <- paste0("python2 munge_sumstats.py --sumstats /Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/gwas_sumstats/", "BMI_GxE.txt.gz", 
                       " --merge-alleles w_hm3.snplist --out gwas_sumstats/proper_format/", "BMI_GxE.txt.gz"," --a1-inc")
command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
system(command)

#### prepare alternate geneset annotations ####
nDE_genes <- lapply(tissue_abbr[tissues], function(tissue){
  focal_genes <- tissue_nodes[[tissue]]$human_ensembl_gene[tissue_nodes[[tissue]]$human_ensembl_gene]
  focal_genes <- focal_genes[!is.na(focal_genes)]
  alt_genes <- setdiff(motrpac_expressed_genes[[tissue_abbr_rev[tissue]]], focal_genes)
  return(alt_genes)
})
names(nDE_genes) <- paste0("not-", names(nDE_genes))

cluster_ids <- names(nDE_genes)
for(cli in cluster_ids){
  sink(paste0(geneset_directory, "/cluster_", cli, ".GeneSet"))
  cat(paste0(unique(nDE_genes[[cli]]), "\n"), sep = "")
  sink()
}

#construct annotation files
windowsize <- "100000"
n_chromosomes <- 22

#make .annot files
foreach(cri=1:n_chromosomes) %dopar% {
  # for(cri in 1:n_chromosomes){
  for(cli in cluster_ids){
    cat(paste0("(chr: ", cri, ", clus: ", cli, ") "))
    command_ldsc <- paste0("python2 make_annot.py --gene-set-file custom_genesets/cluster_", cli, ".GeneSet --gene-coord-file ENSG_coord.txt --windowsize ", windowsize, " --bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri,".bim --annot-file ", geneset_directory, "cluster_", cli, ".chr_", cri,".annot")
    command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
    system(command, ignore.stdout = T)
  }
}

#thicken up those .annots
subset_baseline_to_BIM <- T
# foreach(cri=1:n_chromosomes) %dopar% {
foreach(cri=1:n_chromosomes, .packages = c("data.table")) %dopar% {
  # for(cri in 1:n_chromosomes){
  cat(paste0("\n(chr: ", cri, "): "))
  baseline_annot <- fread(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/baseline/baseline.", cri, ".annot"))
  bim_file <- fread(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri, ".bim"), header = F)
  colnames(bim_file) <- c("CHR", "SNP", "CM", "BP", "minor_allele", "major_allele")
  bim_file <- bim_file[,-c("minor_allele", "major_allele")]
  
  if(subset_baseline_to_BIM){
    SNP_inds <- match(baseline_annot$SNP, bim_file$SNP)
    baseline_present_SNPs <- SNP_inds[!is.na(SNP_inds)]
    baseline_absent_SNPs <- which(is.na(match(bim_file$SNP, baseline_annot$SNP)))
  } else {
    SNP_inds <- match(bim_file$SNP, baseline_annot$SNP)
  }
  
  missing_SNP_inds <- which(is.na(SNP_inds))
  
  # intersect(baseline_annot$SNP, bim_file$SNP)
  for(cli in cluster_ids){
    cat(paste0(cli, " "))
    
    cluster_code <- fread(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets/cluster_", cli, ".chr_", cri, ".annot"))
    cluster_annot <- cbind(bim_file, cluster_code)
    # cluster_annot <- merge(baseline_annot, cluster_annot, by = "SNP")
    # sum(cluster_annot$BP.x != cluster_annot$BP.y)
    # cor(cluster_annot$CM.x, cluster_annot$CM.y)
    
    if(subset_baseline_to_BIM){
      
      merged_annot <- cluster_annot
      merged_annot <- cbind(merged_annot, base = 1)
      merged_annot <- cbind(merged_annot, matrix(0, nrow = nrow(merged_annot), ncol = ncol(baseline_annot) - 5))
      merged_annot[baseline_present_SNPs,7:ncol(merged_annot)] <- baseline_annot[!is.na(SNP_inds), 6:ncol(baseline_annot)]
      
      colnames(merged_annot)[colnames(merged_annot) == "ANNOT"] <- paste0("Cluster_", cli)
      colnames(merged_annot)[startsWith(colnames(merged_annot), "V")] <- colnames(baseline_annot)[6:ncol(baseline_annot)]
      # tail(merged_annot, n = 100)
      
    } else {
      
      merged_annot <- baseline_annot
      merged_annot$cluster_TEMP <- 0
      merged_annot$cluster_TEMP[SNP_inds[!is.na(SNP_inds)]] <- cluster_annot$ANNOT[!is.na(SNP_inds)]
      
      missing_SNPS <- cluster_annot[missing_SNP_inds,]
      missing_SNPS <- cbind(missing_SNPS[,-"ANNOT"], base = rep(1, nrow(missing_SNPS)), 
                            matrix(0, nrow = nrow(missing_SNPS), ncol = ncol(baseline_annot) - 5), missing_SNPS[,"ANNOT"])
      merged_annot <- rbind(merged_annot, missing_SNPS, use.names = F)
      colnames(merged_annot)[colnames(merged_annot) == "cluster_TEMP"] <- paste0("Cluster_", cli)
    }
    
    merged_annot$CM <- as.character(merged_annot$CM)
    merged_annot$CM[merged_annot$CM == "0"] <- "0.0"
    
    #write to disk
    fwrite(merged_annot, paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets/thickCluster_", cli, ".chr_", cri, ".annot"), sep = "\t")
    
  }
}

#process ld files
if(exists("cl")){
  stopCluster(cl)
  remove(cl)
  cl <- makeCluster(2, outfile="")
  registerDoParallel(cl)
}

# detach(“package:data.table”, unload=TRUE) 
foreach(cri=1:n_chromosomes) %dopar% {
  # for(cri in 1:n_chromosomes){
  for(cli in cluster_ids){ #redundant across clusters
    print(paste0("now processing cluster ", cli))
    # command_ldsc <- paste0("python2 ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri," --ld-wind-cm 1 --annot ", geneset_directory, "thickCluster_", cli, ".chr_", cri,".annot --thin-annot --out ", geneset_directory, "cluster_", cli, " --print-snps hapmap3_snps/hm.", cri,".snp")
    command_ldsc <- paste0("python2 ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri," --ld-wind-cm 1 --annot ", 
                           geneset_directory, "thickCluster_", cli, ".chr_", cri,".annot --out ", geneset_directory, 
                           "cluster_", cli, ".chr_", cri, " --print-snps hapmap3_snps/hm.", cri,".snp")
    command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
    system(command)
  }
}  

# for(cri in 1:n_chromosomes){
#   for(cli in cluster_ids){ 
#     file.rename(from = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets/cluster_", cli, ".chr_", cri, ".annot"), 
#                 to = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets/thincluster_", cli, ".chr_", cri, ".annot"))
#     file.rename(from = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets/cluster_", cli, ".chr_", cri, ".annot.gz"), 
#                 to = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets/thincluster_", cli, ".chr_", cri, ".annot.gz"))
#     file.rename(from = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets/thickCluster_", cli, ".chr_", cri, ".annot"), 
#                 to = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets/cluster_", cli, ".chr_", cri, ".annot"))
#   }
# }

if(exists("cl")){
  stopCluster(cl)
  remove(cl)
  cl <- makeCluster(8, outfile="")
  registerDoParallel(cl)
}

foreach(cli=cluster_ids) %dopar% {
  command_ldsc <- paste0("python2 ldsc.py ",
                         "--n-blocks 1000 ",
                         "--h2 ", "gwas_sumstats/proper_format/", "BMI_GxE.txt.gz", ".sumstats.gz ",
                         "--ref-ld-chr custom_genesets/cluster_", cli, ".chr_ ",
                         "--w-ld-chr weights_hm3_no_hla/weights. ",
                         "--overlap-annot ",
                         "--frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. ",
                         "--out output/original_command/", 
                         "BMI_GxE", "_cluster_", cli, " --print-coefficients")
  
  command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
  system(command)
}


foreach(cli=paste0(tissue_abbr[tissues], "-sex_homogeneous_changing")) %dopar% {
  command_ldsc <- paste0("python2 ldsc.py ",
                         "--n-blocks 1000 ",
                         "--h2 ", "gwas_sumstats/proper_format/", "BMI_GxE.txt.gz", ".sumstats.gz ",
                         "--ref-ld-chr custom_genesets/cluster_", cli, ".chr_ ",
                         "--w-ld-chr weights_hm3_no_hla/weights. ",
                         "--overlap-annot ",
                         "--frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. ",
                         "--out output/original_command/", 
                         "BMI_GxE", "_cluster_", cli, " --print-coefficients")
  
  command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
  system(command)
}

#read the results back in
gwas_names <- "BMI_GxE"
ldsc_output_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/output/original_command/"
ldsc_results_paths <- list.files(ldsc_output_dir)
ldsc_log_paths <- paste0(ldsc_output_dir, ldsc_results_paths[grep(ldsc_results_paths, pattern = "log")])
ldsc_results_paths <- paste0(ldsc_output_dir, ldsc_results_paths[grep(ldsc_results_paths, pattern = "results")])
tissue_abbrev <- unique(MotrpacBicQC::bic_animal_tissue_code$abbreviation)
tissue_abbrev <- tissue_abbrev[!is.na(tissue_abbrev)]
tissue_abbrev <- setdiff(tissue_abbrev, c("PLASMA", "HYPOTH", "TESTES", "OVARY", "VENACV"))
cluster_ids <- unique(c(cluster_ids, paste0(tissue_abbrev, "-sex_homogeneous_changing")))
ldsc_results_paths <- lapply(gwas_names, function(gwas) ldsc_results_paths[grep(pattern = gwas, ldsc_results_paths)])
ldsc_log_paths <- lapply(gwas_names, function(gwas) ldsc_log_paths[grep(pattern = gwas, ldsc_log_paths)])
names(ldsc_results_paths) <- names(ldsc_log_paths) <- gwas_names

log_files <- as.data.frame(matrix(data = NA, ncol = 5, nrow = length(gwas_names) * length(cluster_ids)))
colnames(log_files) <- c("h2", "h2se", "chi2", "gwas", "cluster")
log_files$cluster <- rep(cluster_ids, length(gwas_names))
log_files$gwas <- c(sapply(gwas_names, function(gwas_name) rep(gwas_name, length(cluster_ids))))
for(i in 1:nrow(log_files)){
  if(!file.exists(paste0(ldsc_output_dir, log_files$gwas[i], "_cluster_", log_files$cluster[i], ".log"))){next()}
  logfile <- readLines(paste0(ldsc_output_dir, log_files$gwas[i], "_cluster_", log_files$cluster[i], ".log"))
  h2 <- logfile[grep(logfile, pattern = "Total Observed scale h2")]
  log_files$h2[i] <- as.numeric(strsplit(h2, " ")[[1]][5])
  log_files$h2se[i] <- as.numeric(substr(x = strsplit(h2, " ")[[1]][6], start = 2, stop = nchar(strsplit(h2, " ")[[1]][6]) - 1))
  chi2 <- as.numeric(strsplit(logfile[grep(logfile, pattern = "Mean Chi")], split = " ")[[1]][3])
  log_files$chi2[i] <- chi2
}
log_files$gwas <- stringr::str_replace_all(log_files$gwas, "imputed_", "")


ldsc_results <- as.data.frame(matrix(data = NA, ncol = ncol(fread(ldsc_results_paths[[1]][1])), nrow = length(gwas_names) * length(cluster_ids)))
colnames(ldsc_results) <- colnames(fread(ldsc_results_paths[[1]][1]))
ldsc_results$cluster <- rep(cluster_ids, length(gwas_names))
ldsc_results$gwas <- c(sapply(gwas_names, function(gwas_name) rep(gwas_name, length(cluster_ids))))
for(i in 1:nrow(ldsc_results)){
  if(!file.exists(paste0(ldsc_output_dir, log_files$gwas[i], "_cluster_", log_files$cluster[i], ".results"))){next()}
  output <- fread(paste0(ldsc_output_dir, log_files$gwas[i], "_cluster_", log_files$cluster[i], ".results"))
  ldsc_results[i,1:ncol(output)] <- output[grep(pattern = ldsc_results$cluster[i], output$Category),]
}


ldsc_results$logPVal_enrichment <- log10(ldsc_results$Enrichment_p)
hist(ldsc_results$Enrichment_p[!grepl("not", ldsc_results$cluster)], breaks = 100, xlim = c(0,1),
     main = "enrichment p-values for tissue-specific DE genes", xlab = "p-value")
abline(v = 0.05 / sum(!grepl("not", ldsc_results$cluster)), col = 2, lwd = 2)
text(x = 0.05 / sum(!grepl("not", ldsc_results$cluster)), y = par('usr')[4], pos = 3, labels = "bonferroni alpha = 0.05", xpd = NA, col = 2)


##### try now using the varexp package
library(VarExp)

RSID_POS_MAP <- fread("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/RSID_POS_MAP.txt")
RSID_POS_MAP <- split(RSID_POS_MAP, RSID_POS_MAP$`#CHROM`)
GxE <- fread("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/gwas.chrALL.gxe.bmi.all.glm.linear.ADDxPA.hg38.txt")
GxE$V1 <- gsub(GxE$V1, pattern = "_.*", replacement = "")
SNP_frx <- fread("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ukbb.ALL.impute.frqx")
GxE$MAF <- SNP_frx$MAF[match(GxE$V3, SNP_frx$SNP)]
GxE$na <- NA

GWAS <- GxE[,c("V3", "V1", "na", "V4", "MAF", "na", "V9")]
colnames(GWAS) <- c("RSID", "CHR", "POS", "A0", "FREQ_A0", "MAIN_EFFECT", "INT_EFFECT")
GWAS <- do.call(rbind, lapply(split(GWAS, GWAS$CHR), function(gwas_i){
  chr <- gwas_i$CHR[1]
  cat(paste0(chr, " "))
  gwas_i$POS <- RSID_POS_MAP[[chr]]$POS[match(gwas_i$RSID, RSID_POS_MAP[[chr]]$ID)]
  gwas_i <- gwas_i[!is.na(gwas_i$POS),]
  return(gwas_i)
}))

C <- getGenoCorMatrix(GWAS$RSID, GWAS$CHR, GWAS$POS, GWAS$A0, "EUR", pruning = FALSE)


  
