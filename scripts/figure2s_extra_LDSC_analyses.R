# libraries
library(foreach)
library(doParallel)
library(parallel)
library(EnsDb.Hsapiens.v79)
library(MotrpacRatTraining6mo) # v1.6.0

# initialize parallelization
if(!exists("cl")){
  cl <- makeCluster(8, outfile="")
  registerDoParallel(cl)
}
getDoParWorkers()

library(data.table)

#### estimate h2_SNP in all genes and in all rat orthologs ####

#load tested genes
genes_tested_in_transcriptome_DEA <- unique(unlist(MotrpacRatTraining6moData::GENE_UNIVERSES$ensembl_gene$TRNSCRPT))

# load mapping of rat -> human genes
gene_map <- MotrpacRatTraining6moData::RAT_TO_HUMAN_GENE
human_ensembl_genes_in_DEA <- gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID[match(genes_tested_in_transcriptome_DEA, gene_map$RAT_ENSEMBL_ID)]
human_ensembl_genes_in_DEA <- gsub(human_ensembl_genes_in_DEA, pattern = "\\..*", replacement = "")
human_ensembl_genes_in_DEA <- human_ensembl_genes_in_DEA[-which(is.na(human_ensembl_genes_in_DEA))]

#confirm ENSG coordinate mapping
ENSG_coord <- read.table(file = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/ENSG_coord.txt", header = T)
100-length(setdiff(human_ensembl_genes_in_DEA, ENSG_coord$GENE)) / length(human_ensembl_genes_in_DEA) * 100

#specify gene sets
genesets <- list(all_human_genes = ENSG_coord$GENE, all_rat_human_orthologs = human_ensembl_genes_in_DEA)
geneset_ids <- names(genesets)

#write gene sets to file
ldsc_directory <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/"
geneset_directory <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets/"

if(!dir.exists(geneset_directory)){dir.create(geneset_directory)}
for(cli in geneset_ids){
  sink(paste0(geneset_directory, "/geneset_", cli, ".GeneSet"))
  cat(paste0(unique(genesets[[cli]]), "\n"), sep = "")
  sink()
}

#write gene sets to file
ldsc_directory <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/"
geneset_directory <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets/"
sink(paste0(geneset_directory, "/Cluster_Control.GeneSet"))
cat(paste0(human_ensembl_genes_in_DEA, "\n"), sep = "")
sink()

#generic commands
command_changedir <- paste0("source ~/.bash_profile; cd ", ldsc_directory, "; ")
command_initiateLDSC <- "source activate ldsc; "

#construct annotation files
windowsize <- "100000"
n_chromosomes <- 22

#make .annot files
foreach(cri=1:n_chromosomes) %dopar% {
  for(cli in geneset_ids){
    cat(paste0("(chr: ", cri, ", clus: ", cli, ") "))
    command_ldsc <- paste0("python2 make_annot.py --gene-set-file custom_genesets/geneset_", cli, ".GeneSet --gene-coord-file ENSG_coord.txt --windowsize ", windowsize, " --bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri,".bim --annot-file ", geneset_directory, "thingeneset_", cli, ".chr_", cri,".annot")
    command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
    system(command, ignore.stdout = T)
  }
}

#thicken up those .annots
#add both all human and all rat genes in
add_all_human_genes_to_rat <- F
foreach(cri=1:n_chromosomes, .packages = c("data.table")) %dopar% {
  
  if(readLines("~/foreach_continue.txt") != "T"){
    next()
  }
  
  cat(paste0("\n(chr: ", cri, "): "))
  baseline_annot <- fread(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/baseline/baseline.", cri, ".annot"))
  bim_file <- fread(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri, ".bim"), header = F)
  colnames(bim_file) <- c("CHR", "SNP", "CM", "BP", "minor_allele", "major_allele")
  bim_file <- bim_file[,-c("minor_allele", "major_allele")]
  
  #subset baseline to bim
  SNP_inds <- match(baseline_annot$SNP, bim_file$SNP)
  baseline_present_SNPs <- SNP_inds[!is.na(SNP_inds)]
  baseline_absent_SNPs <- which(is.na(match(bim_file$SNP, baseline_annot$SNP)))
  missing_SNP_inds <- which(is.na(SNP_inds))
  
  # merge geneset .annots in and write to disk
  for(cli in geneset_ids){
    
    cat(paste0(cli, " "))
    
    geneset_code <- fread(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets/thingeneset_", cli, ".chr_", cri, ".annot"))
    geneset_annot <- cbind(bim_file, geneset_code)
    
    merged_annot <- geneset_annot
    merged_annot <- cbind(merged_annot, base = 1)
    merged_annot <- cbind(merged_annot, matrix(0, nrow = nrow(merged_annot), ncol = ncol(baseline_annot) - 5))
    merged_annot[baseline_present_SNPs,7:ncol(merged_annot)] <- baseline_annot[!is.na(SNP_inds), 6:ncol(baseline_annot)]
    
    colnames(merged_annot)[colnames(merged_annot) == "ANNOT"] <- paste0("geneset_", cli)
    colnames(merged_annot)[startsWith(colnames(merged_annot), "V")] <- colnames(baseline_annot)[6:ncol(baseline_annot)]
    
    #make some minor lexical / syntactical adjustments
    merged_annot$CM <- as.character(merged_annot$CM)
    merged_annot$CM[merged_annot$CM == "0"] <- "0.0"
    
    if(cli != "all_human_genes"){
      if(add_all_human_genes_to_rat){
        extra_baseline <- fread(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets/thingeneset_", "all_human_genes", ".chr_", cri, ".annot"))
        merged_annot <- cbind(merged_annot, all_human_genes = extra_baseline)
      } 
    } else {
      if(add_all_human_genes_to_rat){
        
      #let's test if subsetting human genes helps to make things less certain (remove for real analysis)
      all_human_subset <- fread(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets/thingeneset_", cli, ".chr_", cri, ".annot"))
      hits <- which(all_human_subset == 1, arr.ind = T)[,1]
      subhits <- sample(hits, round(length(hits) / 50))
      all_human_subset[subhits, 1] <- 0
      merged_annot <- cbind(merged_annot, all_human_subset = all_human_subset)
      }
    }
    
    #write to disk
    fwrite(merged_annot, paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets/geneset_", cli, ".chr_", cri, ".annot"), sep = "\t")
    
    #also write thin annot with thicker columns to disk
    merged_annot_semithin <- merged_annot[,1:5]
    fwrite(merged_annot_semithin, 
           paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets/cts_thin_annot/thingeneset_", cli, 
                  ".chr_", cri, ".annot"), sep = "\t")
    
  }
}

#process ld files
if(exists("cl")){
  stopCluster(cl)
  remove(cl)
  cl <- makeCluster(2, outfile="")
  registerDoParallel(cl)
}
getDoParWorkers()

foreach(cri=1:n_chromosomes) %dopar% {
  for(cli in geneset_ids){
    print(paste0("now processing geneset ", cli))
    command_ldsc <- paste0("python2 ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri,
                           " --ld-wind-cm 1 --annot ", geneset_directory, "geneset_", cli, 
                           ".chr_", cri,".annot --out ", geneset_directory, "geneset_", cli, ".chr_", cri, " --print-snps hapmap3_snps/hm.", cri,".snp")
    command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
    system(command)
  }
}  

#alternatively do it for the thinned annots
geneset_directory_thin <- paste0(geneset_directory, "cts_thin_annot/")
foreach(cri=1:n_chromosomes) %dopar% {
  
  if(readLines("~/foreach_continue.txt") != "T"){
    next()
  }
  
  for(cli in geneset_ids){
    print(paste0("now processing geneset ", cli))
    command_ldsc <- paste0("python2 ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri,
                           " --ld-wind-cm 1 --annot ", geneset_directory_thin, "thingeneset_", cli, 
                           ".chr_", cri,".annot --out ", geneset_directory_thin, "thingeneset_", cli, ".chr_", cri, 
                           " --print-snps hapmap3_snps/hm.", cri,".snp")
    command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
    system(command)
  }
}  

# do some renaming for later convenience (note to self -- don't do this in the future, it's very hacky)
#it's because input above is thickgeneset_ and output is geneset_, so we gotta respect the names
# for(cri in 1:n_chromosomes){
#   for(cli in geneset_ids){ 
#     file.rename(from = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets/geneset_", cli, ".chr_", cri, ".annot"), 
#                 to = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets/thingeneset_", cli, ".chr_", cri, ".annot"))
#     file.rename(from = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets/geneset_", cli, ".chr_", cri, ".annot.gz"), 
#                 to = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets/thingeneset_", cli, ".chr_", cri, ".annot.gz"))
#     file.rename(from = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets/thickgeneset_", cli, ".chr_", cri, ".annot"), 
#                 to = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets/geneset_", cli, ".chr_", cri, ".annot"))
#   }
# }

# make correction to geneset annotation file
foreach(cri=1:n_chromosomes, .packages = c("data.table")) %dopar% {
  # for(cri in 1:n_chromosomes){
  for(cli in geneset_ids){
    cat(paste0("(", cri, ", ", cli, ")"))
    #make geneset annotation thin
    thick_annot <- fread(paste0(ldsc_directory, "custom_genesets/geneset_", cli, ".chr_", cri, ".annot"))
    target_col <- which(colnames(thick_annot) == paste0("geneset_", cli))
    thin_annot <- thick_annot[,..target_col]
    colnames(thin_annot) <- "ANNOT"
    fwrite(thin_annot, paste0(ldsc_directory, "custom_genesets/cts_thin_annot/geneset_", cli, ".chr_", cri, ".annot"))
    fwrite(thin_annot, paste0(ldsc_directory, "custom_genesets/cts_thin_annot/geneset_", cli, ".chr_", cri, ".annot.gz"))
  }
}

#ID target GWAS
gwas_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/imputed_gwas_hg38_1.1/"
gwas_summary_files <- list.files(gwas_dir)
gwas_summary_files <- gwas_summary_files[-grep(gwas_summary_files, pattern = "README")]

#and make cluster_cts.ldcts file
ldcts <- cbind(paste0("Cluster_", cluster_ids, "  "), paste0("custom_genesets/cts_thin_annot/thingeneset_", geneset_ids, ".chr_"), ",", paste0("custom_genesets/cts_thin_annot/cluster_control.chr_"))
ldcts <- sapply(1:nrow(ldcts), function(ri) paste0(ldcts[ri,], collapse = ""))

sink(paste0(ldsc_directory, "orthologs.ldcts"))
cat(paste0(ldcts, "\n"), sep = "")
sink()

# partition heritability
if(exists("cl")){
  stopCluster(cl)
  remove(cl)
}
cl <- makeCluster(8, outfile="")
registerDoParallel(cl)
getDoParWorkers()


foreach(gi=1:length(gwas_summary_files)) %dopar% {
  
  if(readLines("~/foreach_continue.txt") != "T"){
    next()
  }
  
  for(cli in cluster_ids[2]){ 
    cat(paste0(" (", gi, " / ", length(gwas_summary_files), ")")) 
    command_ldsc <- paste0("python2 ldsc.py ",
                           "--n-blocks 1000 ",
                           "--h2-cts ", "gwas_sumstats/proper_format/", gwas_summary_files[gi], ".sumstats.gz ",
                           "--ref-ld-chr custom_genesets/geneset_", cli, ".chr_ ",
                           #whoops, bc the baseline is already included here, don't need to list it again (results in singular matrix error bc of 100% redundancy)
                           #see /Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/unused_scripts/GxE_annotation_enrichment_motrpac.R for another example
                           # "--ref-ld-chr custom_genesets/geneset_", cli, ".chr_,custom_genesets/cts_thin_annot/baseline.chr_ ",
                           "--w-ld-chr weights_hm3_no_hla/weights. ",
                           "--overlap-annot ",
                           "--frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. ",
                           "--out output/original_command/", 
                           strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1], "_cluster_", cli, " --print-coefficients")
    
    command_ldsc <- paste0("python2 ldsc.py ", 
                           "--h2-cts ", "gwas_sumstats/proper_format/", gwas_summary_files[gi], ".sumstats.gz ", 
                           "--ref-ld-chr custom_genesets/cts_thin_annot/baseline.chr_ ", 
                           "--w-ld-chr weights_hm3_no_hla/weights. ", 
                           "--overlap-annot ", 
                           "--ref-ld-chr-cts orthologs.ldcts ", 
                           "--frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. ", 
                           "--n-blocks 200 ", 
                           "--out output/cts/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1], "_cluster_", cluster_ids[1], "-", cluster_ids[length(cluster_ids)], "")
    
    command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
    system(command)
  }
}

#read in cts command results
gwas_names <- gsub(".txt.gz", "", gwas_summary_files)
ldsc_output_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/output/cts/"
ldsc_results_paths <- list.files(ldsc_output_dir, full.names = T)
hist(sapply(ldsc_results_paths[grepl("orthologs", ldsc_results_paths) & grepl("results", ldsc_results_paths)], function(x){
  read.csv(x, sep = "\t")$Coefficient_P_value
}))

#copied from /Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/unused_scripts/GxE_annotation_enrichment_motrpac.R

#read the results back in
gwas_names <- gsub(".txt.gz", "", gwas_summary_files)
ldsc_output_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/output/original_command/"
ldsc_results_paths <- list.files(ldsc_output_dir)
ldsc_log_paths <- paste0(ldsc_output_dir, ldsc_results_paths[grep(ldsc_results_paths, pattern = "log")])
ldsc_results_paths <- paste0(ldsc_output_dir, ldsc_results_paths[grep(ldsc_results_paths, pattern = "results")])
cluster_ids <- geneset_ids
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
# log_files$gwas <- stringr::str_replace_all(log_files$gwas, "imputed_", "")


ldsc_results <- as.data.frame(matrix(data = NA, ncol = ncol(fread(ldsc_results_paths[[1]][1])), 
                                     nrow = length(gwas_names) * length(cluster_ids)))
colnames(ldsc_results) <- colnames(fread(ldsc_results_paths[[1]][1]))
ldsc_results$cluster <- rep(cluster_ids, length(gwas_names))
ldsc_results$gwas <- c(sapply(gwas_names, function(gwas_name) rep(gwas_name, length(cluster_ids))))
for(i in 1:nrow(ldsc_results)){
  filename <- paste0(ldsc_output_dir, ldsc_results$gwas[i], "_cluster_", ldsc_results$cluster[i], ".results")
  if(!file.exists(filename)){
    print(paste0("<<", filename, ">> not found."))
    next()
  }
  output <- fread(paste0(ldsc_output_dir, ldsc_results$gwas[i], "_cluster_", ldsc_results$cluster[i], ".results"))
  output$Category
  ldsc_results[i,1:ncol(output)] <- output[grep(pattern = ldsc_results$cluster[i], output$Category),]
}


ldsc_results$logPVal_enrichment <- log10(ldsc_results$Enrichment_p)
hist(ldsc_results$Enrichment_p[!grepl("not", ldsc_results$cluster)], breaks = 100, xlim = c(0,1),
     main = "enrichment p-values for tissue-specific DE genes", xlab = "p-value")
abline(v = 0.05 / sum(!grepl("not", ldsc_results$cluster)), col = 2, lwd = 2)
text(x = 0.05 / sum(!grepl("not", ldsc_results$cluster)), y = par('usr')[4], pos = 3, labels = "bonferroni alpha = 0.05", xpd = NA, col = 2)
