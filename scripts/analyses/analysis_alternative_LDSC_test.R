# libraries
library(foreach)
library(doParallel)
library(parallel)
library(EnsDb.Hsapiens.v79)
library(MotrpacRatTraining6mo) # v1.6.0
library(data.table)

# initialize parallelization
if(!exists("cl")){
  cl <- makeCluster(8, outfile="")
  registerDoParallel(cl)
}
getDoParWorkers()

#load tissue-specific gene sets
node_sets = MotrpacRatTraining6moData::GRAPH_COMPONENTS$node_sets
nodes_to_look_at <- c("8w_F1_M1", "8w_F-1_M-1")
node_metadata <- lapply(setNames(nodes_to_look_at, nodes_to_look_at), function(node_to_look_at)
  cbind(do.call(rbind, strsplit(node_sets[node_to_look_at][[node_to_look_at]][grep(
    node_sets[node_to_look_at][[node_to_look_at]], pattern = "TRNSCRPT")], ";")), node_to_look_at))
node_metadata <- as.data.table(do.call(rbind, node_metadata))
colnames(node_metadata) <- c("ome","tissue","ensembl_gene", "node")

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

#process rat orthologs
node_metadata$human_orthologs_ENSEMBL <- setNames(gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID, gene_map$RAT_ENSEMBL_ID)[node_metadata$ensembl_gene]
node_metadata$human_orthologs_ENSEMBL <- gsub("\\..*", "", node_metadata$human_orthologs_ENSEMBL)
node_metadata <- node_metadata[!is.na(node_metadata$human_orthologs_ENSEMBL),]
node_metadata <- node_metadata[node_metadata$human_orthologs_ENSEMBL %in% ENSG_coord$GENE,]

#specify useful directories
ldsc_directory <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/"
geneset_directory <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets_revisions/"
geneset_directory_thin <- paste0(geneset_directory, "cts_thin_annot/")

#specify genesets
genesets <- split(node_metadata$human_orthologs_ENSEMBL, node_metadata$tissue)
n_rand_genesets <- 5
n_genes_per_tissue <- sapply(genesets, length)
sim_genesets <- unlist(lapply(1:n_rand_genesets, function(i){
  lapply(1:length(n_genes_per_tissue), function(j){
    sim_geneset_name <- paste0(names(n_genes_per_tissue)[j], ".random.", i)
    if(file.exists(paste0(geneset_directory, "geneset_", sim_geneset_name, ".GeneSet"))){
      out <- list(readLines(paste0(geneset_directory, "geneset_", sim_geneset_name, ".GeneSet")))
    } else {
      out <- list(sample(human_ensembl_genes_in_DEA, n_genes_per_tissue[j], F))  
    }
    names(out) <- sim_geneset_name
    out
  })
}), recursive = F)
genesets <- c(list(all_rat_human_orthologs = human_ensembl_genes_in_DEA), genesets)
genesets <- c(genesets, unlist(sim_genesets, recursive = F))
geneset_ids <- names(genesets)

#write gene sets to file
if(!dir.exists(geneset_directory)){dir.create(geneset_directory)}
for(cli in geneset_ids){
  sink(paste0(geneset_directory, "/geneset_", cli, ".GeneSet"))
  cat(paste0(unique(genesets[[cli]]), "\n"), sep = "")
  sink()
}

#write control gene set to file
sink(paste0(geneset_directory, "/Control.GeneSet"))
cat(paste0(human_ensembl_genes_in_DEA, "\n"), sep = "")
sink()

#generic commands
command_changedir <- paste0("source ~/.bash_profile; cd ", ldsc_directory, "; ")
command_initiateLDSC <- "source activate ldsc; "

#ID target GWAS
gwas_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/imputed_gwas_hg38_1.1/"
gwas_summary_files <- list.files(gwas_dir)
gwas_summary_files <- gwas_summary_files[-grep(gwas_summary_files, pattern = "README")]

#construct annotation files
windowsize <- "100000"
n_chromosomes <- 22

#make .annot files
foreach(cri=1:n_chromosomes) %dopar% {
  for(cli in geneset_ids){
    
    if(readLines("~/foreach_continue.txt") != "T"){
      next()
    }
    
    cat(paste0("(chr: ", cri, ", clus: ", cli, ") "))
    command_ldsc <- paste0("python2 make_annot.py --gene-set-file custom_genesets_revisions/geneset_", cli, ".GeneSet --gene-coord-file ENSG_coord.txt --windowsize ", windowsize, " --bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri,".bim --annot-file ", geneset_directory, "thin_geneset_", cli, ".chr_", cri,".annot")
    command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
    system(command, ignore.stdout = T)
  }
}

#thicken up those .annots
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
    
    if(file.exists(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets_revisions/cts_thin_annot/thin_geneset_", cli, 
                          ".chr_", cri, ".annot"))){
      next()
    }
    
    geneset_code <- fread(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets_revisions/thin_geneset_", cli, ".chr_", cri, ".annot"))
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
    
    #write to disk
    # fwrite(merged_annot, paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets_revisions/geneset_", cli, ".chr_", cri, ".annot"), sep = "\t")
    
    #also write thin annot with thicker columns to disk
    merged_annot_semithin <- merged_annot[,1:5]
    fwrite(merged_annot_semithin, 
           paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/custom_genesets_revisions/cts_thin_annot/thin_geneset_", cli, 
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

# compute ld scores for all these annots
foreach(cri=1:n_chromosomes) %dopar% {
  
  if(readLines("~/foreach_continue.txt") != "T"){
    next()
  }
  
  for(cli in geneset_ids){
    
    print(paste0("now processing geneset ", cli))
    
    if(file.exists(paste0(geneset_directory_thin, "thin_geneset_", cli, ".chr_", cri, ".l2.ldscore.gz"))){
      next()
    }
    
    command_ldsc <- paste0("python2 ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri,
                           " --ld-wind-cm 1 --annot ", geneset_directory_thin, "thin_geneset_", cli, 
                           ".chr_", cri,".annot --out ", geneset_directory_thin, "thin_geneset_", cli, ".chr_", cri, 
                           " --print-snps hapmap3_snps/hm.", cri,".snp")
    command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
    system(command)
  }
}  


# make correction to geneset annotation file
foreach(cri=1:n_chromosomes, .packages = c("data.table")) %dopar% {
  # for(cri in 1:n_chromosomes){
  for(cli in geneset_ids){
    cat(paste0("(", cri, ", ", cli, ")"))
    #make geneset annotation thin
    thick_annot <- fread(paste0(ldsc_directory, "custom_genesets_revisions/cts_thin_annot/thin_geneset_", cli, ".chr_", cri, ".annot"))
    target_col <- which(colnames(thick_annot) == paste0("geneset_", cli))
    thin_annot <- thick_annot[,..target_col]
    colnames(thin_annot) <- "ANNOT"
    fwrite(thin_annot, paste0(ldsc_directory, "custom_genesets_revisions/cts_thin_annot/extra_thin_geneset_", cli, ".chr_", cri, ".annot"))
    fwrite(thin_annot, paste0(ldsc_directory, "custom_genesets_revisions/cts_thin_annot/extra_thin_geneset_", cli, ".chr_", cri, ".annot.gz"))
  }
}

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
  
  for(cli in geneset_ids){ 
    
    if(readLines("~/foreach_continue.txt") != "T"){
      next()
    }
    
    cat(paste0(" (", gi, " / ", length(gwas_summary_files), ")")) 
    command_ldsc <- paste0("python2 ldsc.py ",
                           "--n-blocks 1000 ",
                           "--h2 ", "gwas_sumstats/proper_format/", gwas_summary_files[gi], ".sumstats.gz ",
                           "--ref-ld-chr custom_genesets_revisions/cts_thin_annot/thin_geneset_", cli, ".chr_,custom_genesets_revisions/cts_thin_annot/baseline.chr_ ",
                           "--w-ld-chr weights_hm3_no_hla/weights. ",
                           "--overlap-annot ",
                           "--frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. ",
                           "--out output_revisions/group_analysis/", 
                           strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1], "_geneset_", cli, 
                           " --print-coefficients")
    
    command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
    system(command)
    
  }
  
  # make geneset_cts.ldcts file
  subset_geneset_ids <- geneset_ids[!grepl("random", geneset_ids) & 
                                      !grepl("all_rat_human", geneset_ids)]
  ldcts <- cbind(paste0("geneset_", subset_geneset_ids, "  "), 
                 paste0("custom_genesets_revisions/cts_thin_annot/thin_geneset_", subset_geneset_ids, ".chr_"), ",", 
                 paste0("custom_genesets_revisions/cts_thin_annot/thin_geneset_all_rat_human_orthologs.chr_"))
  ldcts <- sapply(1:nrow(ldcts), function(ri) paste0(ldcts[ri,], collapse = ""))
  
  sink(paste0(ldsc_directory, "motrpac_tissues.ldcts"))
  cat(paste0(ldcts, "\n"), sep = "")
  sink()
  
  command_ldsc <- paste0("python2 ldsc.py ", 
                         "--h2-cts ", "gwas_sumstats/proper_format/", gwas_summary_files[gi], ".sumstats.gz ", 
                         "--ref-ld-chr custom_genesets_revisions/cts_thin_annot/baseline.chr_ ", 
                         "--w-ld-chr weights_hm3_no_hla/weights. ", 
                         "--overlap-annot ", 
                         "--ref-ld-chr-cts motrpac_tissues.ldcts ", 
                         "--frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. ", 
                         "--n-blocks 200 ", 
                         "--out output_revisions/specific_analysis/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1], "_geneset_", subset_geneset_ids[1], "-", subset_geneset_ids[length(subset_geneset_ids)], "")
  
  command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
  system(command)
  
  #also iterate over random genesets
  for(ri in 1:n_rand_genesets){
    
    # make geneset_cts.ldcts file
    subset_geneset_ids <- geneset_ids[grepl(paste0("random.", ri), geneset_ids)]
    ldcts <- cbind(paste0("geneset_", subset_geneset_ids, "  "), 
                   paste0("custom_genesets_revisions/cts_thin_annot/thin_geneset_", subset_geneset_ids, ".chr_"), ",", 
                   paste0("custom_genesets_revisions/cts_thin_annot/thin_geneset_all_rat_human_orthologs.chr_"))
    ldcts <- sapply(1:nrow(ldcts), function(ri) paste0(ldcts[ri,], collapse = ""))
    
    sink(paste0(ldsc_directory, "motrpac_tissues_random.", ri, ".ldcts"))
    cat(paste0(ldcts, "\n"), sep = "")
    sink()
    
    command_ldsc <- paste0("python2 ldsc.py ", 
                           "--h2-cts ", "gwas_sumstats/proper_format/", gwas_summary_files[gi], ".sumstats.gz ", 
                           "--ref-ld-chr custom_genesets_revisions/cts_thin_annot/baseline.chr_ ", 
                           "--w-ld-chr weights_hm3_no_hla/weights. ", 
                           "--overlap-annot ", 
                           "--ref-ld-chr-cts motrpac_tissues_random.", ri,".ldcts ", 
                           "--frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. ", 
                           "--n-blocks 200 ", 
                           "--out output_revisions/specific_analysis/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1], "_geneset_", subset_geneset_ids[1], "-", subset_geneset_ids[length(subset_geneset_ids)])
    
    command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
    system(command)
    
  }
  
}

gwas_names <- gsub(".txt.gz", "", gwas_summary_files)


par(mfrow = c(2,2))
tissues_to_examine <- names(n_genes_per_tissue)#[1]
#read in and quickly plot cts command results
ldsc_output_revisions_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/output_revisions/specific_analysis/"
ldsc_results_paths <- list.files(ldsc_output_revisions_dir, full.names = T)

hist(sapply(ldsc_results_paths[grepl("random", ldsc_results_paths) & 
                                 grepl("results", ldsc_results_paths)
  ], function(x){
    res <- read.csv(x, sep = "\t")
    res$Coefficient_P_value[apply(do.call(rbind, lapply(tissues_to_examine, function(tiss) 
      grepl(tiss, res$Name))), 2, any)]
  }), breaks = 0:100/100, xlab = "enrichment p-value", main = "\"Cell-type specific analysis\", random gene sets")

hist(sapply(ldsc_results_paths[!grepl("random", ldsc_results_paths) & 
                                 grepl("results", ldsc_results_paths)
  ], function(x){
    res <- read.csv(x, sep = "\t")
    res$Coefficient_P_value[apply(do.call(rbind, lapply(tissues_to_examine, function(tiss) 
        grepl(tiss, res$Name))), 2, any)]
    
  }), breaks = 0:100/100, xlab = "enrichment p-value", main = "\"Cell-type specific analysis\", actual gene sets")

#read in and quickly plot orig command results
ldsc_output_revisions_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/output_revisions/group_analysis/"
ldsc_results_paths <- list.files(ldsc_output_revisions_dir, full.names = T)

hist(sapply(ldsc_results_paths[grepl("random", ldsc_results_paths) & 
                                 grepl("results", ldsc_results_paths) &
                                 apply(do.call(rbind, lapply(tissues_to_examine, function(tiss) 
                                   grepl(tiss, ldsc_results_paths))), 2, any)
  ], function(x){
    res <- read.csv(x, sep = "\t")
    res$Enrichment_p[res$Category == "L2_0"]
  }), breaks = 0:100/100, xlab = "enrichment p-value", main = "\"Cell-type group analysis\", random gene sets")

hist(sapply(ldsc_results_paths[!grepl("random", ldsc_results_paths) & 
                                 grepl("results", ldsc_results_paths) &
                                 apply(do.call(rbind, lapply(tissues_to_examine, function(tiss) 
                                   grepl(tiss, ldsc_results_paths))), 2, any)
  ], function(x){
    res <- read.csv(x, sep = "\t")
    res$Enrichment_p[res$Category == "L2_0"]
  }), breaks = 0:100/100, xlab = "enrichment p-value", main = "\"Cell-type group analysis\", actual gene sets")

#copied from /Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/unused_scripts/GxE_annotation_enrichment_motrpac.R

#read the results back in
gwas_names <- gsub(".txt.gz", "", gwas_summary_files)
ldsc_output_revisions_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/output_revisions/group_analysis/"
ldsc_results_paths <- list.files(ldsc_output_revisions_dir)
ldsc_log_paths <- paste0(ldsc_output_revisions_dir, ldsc_results_paths[grep(ldsc_results_paths, pattern = "log")])
ldsc_results_paths <- paste0(ldsc_output_revisions_dir, ldsc_results_paths[grep(ldsc_results_paths, pattern = "results")])
geneset_ids <- geneset_ids
ldsc_results_paths <- lapply(gwas_names, function(gwas) ldsc_results_paths[grep(pattern = gwas, ldsc_results_paths)])
ldsc_log_paths <- lapply(gwas_names, function(gwas) ldsc_log_paths[grep(pattern = gwas, ldsc_log_paths)])
names(ldsc_results_paths) <- names(ldsc_log_paths) <- gwas_names

log_files <- as.data.frame(matrix(data = NA, ncol = 5, nrow = length(gwas_names) * length(geneset_ids)))
colnames(log_files) <- c("h2", "h2se", "chi2", "gwas", "geneset")
log_files$geneset <- rep(geneset_ids, length(gwas_names))
log_files$gwas <- c(sapply(gwas_names, function(gwas_name) rep(gwas_name, length(geneset_ids))))
for(i in 1:nrow(log_files)){
  if(!file.exists(paste0(ldsc_output_revisions_dir, log_files$gwas[i], "_geneset_", log_files$geneset[i], ".log"))){next()}
  logfile <- readLines(paste0(ldsc_output_revisions_dir, log_files$gwas[i], "_geneset_", log_files$geneset[i], ".log"))
  h2 <- logfile[grep(logfile, pattern = "Total Observed scale h2")]
  log_files$h2[i] <- as.numeric(strsplit(h2, " ")[[1]][5])
  log_files$h2se[i] <- as.numeric(substr(x = strsplit(h2, " ")[[1]][6], start = 2, stop = nchar(strsplit(h2, " ")[[1]][6]) - 1))
  chi2 <- as.numeric(strsplit(logfile[grep(logfile, pattern = "Mean Chi")], split = " ")[[1]][3])
  log_files$chi2[i] <- chi2
}
# log_files$gwas <- stringr::str_replace_all(log_files$gwas, "imputed_", "")


ldsc_results <- as.data.frame(matrix(data = NA, ncol = ncol(fread(ldsc_results_paths[[1]][1])), 
                                     nrow = length(gwas_names) * length(geneset_ids)))
colnames(ldsc_results) <- colnames(fread(ldsc_results_paths[[1]][1]))
ldsc_results$geneset <- rep(geneset_ids, length(gwas_names))
ldsc_results$gwas <- c(sapply(gwas_names, function(gwas_name) rep(gwas_name, length(geneset_ids))))
for(i in 1:nrow(ldsc_results)){
  filename <- paste0(ldsc_output_revisions_dir, ldsc_results$gwas[i], "_geneset_", ldsc_results$geneset[i], ".results")
  if(!file.exists(filename)){
    print(paste0("<<", filename, ">> not found."))
    next()
  }
  output_revisions <- fread(paste0(ldsc_output_revisions_dir, ldsc_results$gwas[i], "_geneset_", ldsc_results$geneset[i], ".results"))
  output_revisions$Category
  ldsc_results[i,1:ncol(output_revisions)] <- output_revisions[grep(pattern = ldsc_results$geneset[i], output_revisions$Category),]
}


ldsc_results$logPVal_enrichment <- log10(ldsc_results$Enrichment_p)
hist(ldsc_results$Enrichment_p[!grepl("not", ldsc_results$geneset)], breaks = 100, xlim = c(0,1),
     main = "enrichment p-values for tissue-specific DE genes", xlab = "p-value")
abline(v = 0.05 / sum(!grepl("not", ldsc_results$geneset)), col = 2, lwd = 2)
text(x = 0.05 / sum(!grepl("not", ldsc_results$geneset)), y = par('usr')[4], pos = 3, labels = "bonferroni alpha = 0.05", xpd = NA, col = 2)
