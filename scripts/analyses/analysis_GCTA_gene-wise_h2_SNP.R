#specify paths
gtex_pipeline_directory <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/gtex-pipeline/"
gtex_covariates_directory <- paste0(gtex_pipeline_directory, "GTEx_Analysis_v8_eQTL_covariates/")
gtex_covariates_files <- list.files(gtex_covariates_directory)
gtex_log2expr_directory <- paste0(gtex_pipeline_directory, "log2-normalized-expression/")
expression_files <- list.files(gtex_log2expr_directory)
expression_files <- setdiff(expression_files, expression_files[grep("tbi", expression_files)])
gcta_directory <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/gcta_1.93.2beta_mac/"
gcta_covariates_directory <- paste0(gcta_directory, "GTEx_Analysis_v8_covariates/")
if(!dir.exists(gcta_covariates_directory)){dir.create(gcta_covariates_directory)}
gcta_phenotypes_directory <- paste0(gcta_directory, "log2_gene_expression/")
if(!dir.exists(gcta_phenotypes_directory)){dir.create(gcta_phenotypes_directory)}
peer_factor_recommendation <- cbind(N_Greater_Than = c(0,150,250,350), n_peer_factors_to_use = c(1:4*15)) #from https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl#2-calculate-peer-factors

#load name map
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

#munge and write covariates files
for(covariate_file in gtex_covariates_files){
  print(covariate_file)
  covariates <- as.data.frame(fread(paste0(gtex_covariates_directory, covariate_file)))
  rownames(covariates) <- covariates$ID
  covariates <- covariates[,-1]
  pf2u <- paste0("InferredCov", 1:peer_factor_recommendation[sum((ncol(covariates)) > peer_factor_recommendation[,"N_Greater_Than"]), "n_peer_factors_to_use"])
  if(length(unique(as.numeric(covariates["platform",]))) > 2){print("more than two platforms"); break}
  if(length(unique(as.numeric(covariates["sex",]))) > 2){print("more than two sexes"); break}
  if(length(unique(as.numeric(covariates["pcr",]))) > 2){print("more than two pcrs"); break}
  quantitative_gcta_covariates <- cbind(0, 
                                       colnames(covariates), 
                                       t(covariates[pf2u,]),
                                       t(covariates[paste0("PC", 1:5),])
                                       )
  
  #check for uniformity
  identical_covariates_1 <- which(sapply(apply(quantitative_gcta_covariates[,-1], 2, unique), function(x) length(x) == 1)) + 1
  if(length(identical_covariates_1) > 0){
    quantitative_gcta_covariates <- quantitative_gcta_covariates[,-identical_covariates_1]
  }
  
  qualitative_gcta_covariates <- cbind(0, 
                                        colnames(covariates), 
                                       as.numeric(covariates["pcr",]),
                                       as.numeric(covariates["platform",]),
                                       as.numeric(covariates["sex",]-1)
  )
  
  #check for uniformity
  identical_covariates_2 <- which(sapply(apply(qualitative_gcta_covariates[,-1], 2, unique), function(x) length(x) == 1)) + 1
  if(length(identical_covariates_2) > 0){
    qualitative_gcta_covariates <- qualitative_gcta_covariates[,-identical_covariates_2]
  }
  
  fwrite(quantitative_gcta_covariates, paste0(gcta_covariates_directory, strsplit(covariate_file, ".v8.")[[1]][1], ".qcovar"), quote = F, col.names = F, sep = "\t")
  fwrite(qualitative_gcta_covariates, paste0(gcta_covariates_directory, strsplit(covariate_file, ".v8.")[[1]][1], ".covar"), quote = F, col.names = F, sep = "\t")
}

#munge and write expression files
for(tissue in names(motrpac_gtex_map)){
  
  print(tissue)
  
  tissue_expression_directory <- paste0(gcta_phenotypes_directory, tissue, "/")
  if(!dir.exists(tissue_expression_directory)){dir.create(tissue_expression_directory)}
  
  expression_file <- expression_files[grep(tissue, expression_files)]
  expression <- fread(paste0(gtex_log2expr_directory, expression_file))
  expression$gene_id <- gsub('\\..*', '', expression$gene_id)
  expression <- split(expression, f = expression$gene_id)
  for(gene in expression){
    phenotype <- data.table(0,
          colnames(gene)[grep("GTEX-", colnames(gene))],
          as.numeric(as.data.frame(gene)[,colnames(gene)[grep("GTEX-", colnames(gene))]])
          )
    fwrite(phenotype, paste0(tissue_expression_directory, gene$gene_id, ".phen"), quote = F, col.names = F, sep = "\t")
  }
  
}

#now calculate heritability for all the genes in all the tissues
# gcta64 --grm GTEx_v8 --qcovar GTEx_Analysis_v8_covariates/Whole_Blood.qcovar --covar GTEx_Analysis_v8_covariates/Whole_Blood.covar --pheno log2_gene_expression/t30-blood-rna/ENSG00000230312.phen --reml --out output/test --thread-num 1
command_changedir <- paste0("source ~/.zshenv; cd ", gcta_directory, "; ")

#make cluster
library(foreach)
library(doParallel)
library(parallel)
if(exists("cl")){
  stopCluster(cl)
  remove(cl)
} 
cl <- makeCluster(9, outfile="")
registerDoParallel(cl)

#iterate over tissues and genes
foreach(tissue=names(motrpac_gtex_map)) %dopar% {
  
  print(tissue)
  tissue_output_directory <- paste0(gcta_directory, "output/", tissue)
  if(!dir.exists(tissue_output_directory)){dir.create(tissue_output_directory)}
  
  genes <- list.files(paste0(gcta_phenotypes_directory, "/", tissue))
  genes <- gsub(".phen", "", genes)
  
  for(gene in genes){ 
    command_gcta <- paste0("gcta64 --grm GTEx_v8 ", 
                           "--qcovar GTEx_Analysis_v8_covariates/", motrpac_gtex_map[tissue], ".qcovar ",
                           "--covar GTEx_Analysis_v8_covariates/", motrpac_gtex_map[tissue], ".covar ",
                           "--pheno log2_gene_expression/", tissue, "/", gene, ".phen ", 
                           "--reml --out output/", tissue, "/", gene, " ", 
                           "--thread-num 1"
                           )
    command <- paste0(command_changedir, command_gcta)
    system(command)
  
  }
}

#read in results and put into single file
library(data.table)
gcta_output <- list()
for(tissue in names(motrpac_gtex_map)){
  
  print(tissue)
  tissue_output_directory <- paste0(gcta_directory, "output/", tissue, "/")
  
  genes <- list.files(tissue_output_directory)
  genes <- genes[grep("hsq", genes)]
  
  gcta_output_tissue <- data.frame(matrix(0, nrow = length(genes), ncol = 5))
  colnames(gcta_output_tissue) <- c("tissue", "ENSG", "h2", "SE", "p")
  for(gi in 1:length(genes)){ 
    if(gi %% round(length(genes) / 100) == 0){cat(paste0(" ", round(gi / length(genes) * 100)))}
    txt <- readLines(paste0(tissue_output_directory, genes[gi]))  
    gcta_output_tissue[gi,] <- c(tissue, gsub(".hsq", "", genes[gi]), strsplit(txt[5], "\t")[[1]][2:3], strsplit(txt[10], "\t")[[1]][2])
  }
  gcta_output_tissue$h2 <- as.numeric(gcta_output_tissue$h2)
  gcta_output_tissue$SE <- as.numeric(gcta_output_tissue$SE)
  gcta_output_tissue$p <- as.numeric(gcta_output_tissue$p)
  
  gcta_output[[tissue]] <- gcta_output_tissue
}
save(gcta_output, file = paste0(gcta_directory, "gcta_output_GTEx_allTissues.RData"))
load(paste0(gcta_directory, "gcta_output_GTEx_allTissues.RData"))

sum(log(unlist(sapply(gcta_output, function(tx) tx$p))) < (log(0.05) - log(length(unlist(sapply(gcta_output, function(tx) tx$p))))))
