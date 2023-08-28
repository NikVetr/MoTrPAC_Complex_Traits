#### revTWMR ####
library(MotrpacRatTraining6mo)
library(foreach)
library(doParallel)
library(parallel)

#initialize parallelization
if(!exists("cl")){
  cl <- makeCluster(8, outfile="")
  registerDoParallel(cl)
}
getDoParWorkers()

library(data.table)

#specify basic commands
revTWMR_directory <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/revTWMR/"
gtex_pipeline_directory <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/gtex-pipeline/"
gwas_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/imputed_gwas_hg38_1.1/"
gwas_summary_files <- list.files(gwas_dir)
gwas_summary_files <- gwas_summary_files[-grep(gwas_summary_files, pattern = "README")]
command_changedir <- paste0("source ~/.bash_profile; cd ", revTWMR_directory, "; ")

#take a look at what these are doing
example_betaGWAS <- fread(paste0(revTWMR_directory, "bmi.matrix.betaGWAS"))
example_genes.N <- fread(paste0(revTWMR_directory, "genes.N"), header = T)
dim(example_betaGWAS)
# each call to revTWMR takes the form of:
"R < revTWMR.R --no-save bmi.matrix.betaGWAS bmi genes.N"
# where example_betaGWAS is a space-separated table of SNP (rs#), ensembl_id_1, ensembl_2, ..., BETA_GWAS, SE, N
#with eqtl effect sizes for SNPs being in the columns and GWAS effect size, standard error, and sample size at the end
#there is also a genes.N file with the columns of bmi.matrix.betaGWAS in the first column (space separated) and sample size in the second
bgw_cols <- c(1:5, (ncol(example_betaGWAS)-5):ncol(example_betaGWAS))
example_betaGWAS[1:3, ..bgw_cols]
example_genes.N
readLines(paste0(revTWMR_directory, "genes.N"), 3)

#retrieve GTEx eQTLs
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
                     't66-lung'='Lung',
                     't67-small-intestine'='Small_Intestine_Terminal_Ileum',
                     't68-liver'='Liver',
                     't70-white-adipose'='Adipose_Subcutaneous')

GTEx_eQTLs = c("Adipose_Subcutaneous.allpairs.txt.gz", "
                Adipose_Visceral_Omentum.allpairs.txt.gz", "
                Adrenal_Gland.allpairs.txt.gz", "
                Artery_Aorta.allpairs.txt.gz", "
                Artery_Coronary.allpairs.txt.gz", "
                Artery_Tibial.allpairs.txt.gz", "
                Brain_Amygdala.allpairs.txt.gz", "
                Brain_Anterior_cingulate_cortex_BA24.allpairs.txt.gz", "
                Brain_Caudate_basal_ganglia.allpairs.txt.gz", "
                Brain_Cerebellar_Hemisphere.allpairs.txt.gz", "
                Brain_Cerebellum.allpairs.txt.gz", "
                Brain_Cortex.allpairs.txt.gz", "
                Brain_Frontal_Cortex_BA9.allpairs.txt.gz", "
                Brain_Hippocampus.allpairs.txt.gz", "
                Brain_Hypothalamus.allpairs.txt.gz", "
                Brain_Nucleus_accumbens_basal_ganglia.allpairs.txt.gz", "
                Brain_Putamen_basal_ganglia.allpairs.txt.gz", "
                Brain_Spinal_cord_cervical_c-1.allpairs.txt.gz", "
                Brain_Substantia_nigra.allpairs.txt.gz", "
                Breast_Mammary_Tissue.allpairs.txt.gz", "
                Cells_Cultured_fibroblasts.allpairs.txt.gz", "
                Cells_EBV-transformed_lymphocytes.allpairs.txt.gz", "
                Colon_Sigmoid.allpairs.txt.gz", "
                Colon_Transverse.allpairs.txt.gz", "
                Esophagus_Gastroesophageal_Junction.allpairs.txt.gz", "
                Esophagus_Mucosa.allpairs.txt.gz", "
                Esophagus_Muscularis.allpairs.txt.gz", "
                Heart_Atrial_Appendage.allpairs.txt.gz", "
                Heart_Left_Ventricle.allpairs.txt.gz", "
                Kidney_Cortex.allpairs.txt.gz", "
                Liver.allpairs.txt.gz", "
                Lung.allpairs.txt.gz", "
                Minor_Salivary_Gland.allpairs.txt.gz", "
                Muscle_Skeletal.allpairs.txt.gz", "
                Nerve_Tibial.allpairs.txt.gz", "
                Ovary.allpairs.txt.gz", "
                Pancreas.allpairs.txt.gz", "
                Pituitary.allpairs.txt.gz", "
                Prostate.allpairs.txt.gz", "
                Skin_Not_Sun_Exposed_Suprapubic.allpairs.txt.gz", "
                Skin_Sun_Exposed_Lower_leg.allpairs.txt.gz", "
                Small_Intestine_Terminal_Ileum.allpairs.txt.gz", "
                Spleen.allpairs.txt.gz", "
                Stomach.allpairs.txt.gz", "
                Testis.allpairs.txt.gz", "
                Thyroid.allpairs.txt.gz", "
                Uterus.allpairs.txt.gz", "
                Vagina.allpairs.txt.gz", "
                Whole_Blood.allpairs.txt.gz")

GTEx_eQTLs <- GTEx_eQTLs[sapply(motrpac_gtex_map, function(x) grep(x, GTEx_eQTLs))]
GTEx_eQTLs <- gsub(x = GTEx_eQTLs, pattern = "\\n                ", replacement = "")

# #run this if you haven't split them up already
# for(tissue in GTEx_eQTLs){
#   print(tissue)
#   if(!file.exists(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/GTEx_Analysis_v8_eQTL_all_associations/", tissue))){
#     system(paste0("scp -r nikgvetr@smsh11dsu-srcf-d15-38.scg.stanford.edu:/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_all_associations/", tissue, " /Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/GTEx_Analysis_v8_eQTL_all_associations/", tissue))  
#   }
# }

#load tissue-specific gene sets
node_sets = MotrpacRatTraining6moData::GRAPH_COMPONENTS$node_sets
nodes_to_look_at <- c("8w_F1_M1", "8w_F-1_M-1")
node_metadata <- lapply(setNames(nodes_to_look_at, nodes_to_look_at), function(node_to_look_at)
  cbind(do.call(rbind, strsplit(node_sets[node_to_look_at][[node_to_look_at]][grep(
    node_sets[node_to_look_at][[node_to_look_at]], pattern = "TRNSCRPT")], ";")), node_to_look_at))
node_metadata <- as.data.table(do.call(rbind, node_metadata))
colnames(node_metadata) <- c("ome","tissue","ensembl_gene", "node")

# load mapping of rat -> human genes
gene_map <- MotrpacRatTraining6moData::RAT_TO_HUMAN_GENE

#process rat orthologs
ENSG_coord <- read.table(file = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/ENSG_coord.txt", header = T)
node_metadata$human_orthologs_ENSEMBL <- setNames(gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID, gene_map$RAT_ENSEMBL_ID)[node_metadata$ensembl_gene]
node_metadata$human_orthologs_ENSEMBL <- gsub("\\..*", "", node_metadata$human_orthologs_ENSEMBL)
node_metadata <- node_metadata[!is.na(node_metadata$human_orthologs_ENSEMBL),]
node_metadata$tissue_code <- TISSUE_ABBREV_TO_CODE[node_metadata$tissue]
DEGs <- split(node_metadata$human_orthologs_ENSEMBL, node_metadata$tissue_code)
# node_metadata <- node_metadata[node_metadata$human_orthologs_ENSEMBL %in% ENSG_coord$GENE,]

#create directories
if(!dir.exists(paste0(revTWMR_directory, "betaGWAS"))){dir.create(paste0(revTWMR_directory, "betaGWAS"))}
if(!dir.exists(paste0(revTWMR_directory, "genesN"))){dir.create(paste0(revTWMR_directory, "genesN"))}

if(!all(file.exists(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/RSID_POS_MAP_",c(1:22, "X"),".txt")))){
  #parse SNP IDs
  RSID_POS_MAP <- fread("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/RSID_POS_MAP.txt")
  colnames(RSID_POS_MAP)[1] <- "CHROM"
  RSID_POS_MAP <- lapply(setNames(unique(RSID_POS_MAP$CHROM), unique(RSID_POS_MAP$CHROM)), function(cri) RSID_POS_MAP[CHROM == cri])
  
  #slice up the positional map for easier processing
  for(cri in names(RSID_POS_MAP)){
    print(cri)
    fwrite(RSID_POS_MAP[[cri]], paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/RSID_POS_MAP_",cri,".txt"))
  }
} else {
  #read it back in
  RSID_POS_MAP <- list()
  for(cri in c(1:22, "X")){
    print(cri)
    RSID_POS_MAP[[cri]] <- fread(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/RSID_POS_MAP_",cri,".txt"))
  }  
}

#slice up the eQTL files by chromosome for easier processsing
for(tissue in names(motrpac_gtex_map)){
  print(tissue)
  if(!file.exists(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/",
                         "data/external/GTEx_Analysis_v8_eQTL_all_associations/", 
                         tissue, "_", cri, ".txt.gz"))){
    
    eQTL_sumstats <- fread(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/GTEx_Analysis_v8_eQTL_all_associations/", 
                                  GTEx_eQTLs[grep(motrpac_gtex_map[match(tissue, names(motrpac_gtex_map))], GTEx_eQTLs)])[1],
                           select = c("gene_id", "variant_id", "slope"))
    CHROMs <- substr(eQTL_sumstats$variant_id, 4, 5)
    CHROMs <- gsub("_", "", CHROMs)
    CHROMs <- lapply(setNames(c(1:22, "X"), c(1:22, "X")), function(cri) which(CHROMs == cri))
    for(cri in c(1:22, "X")){
      fwrite(eQTL_sumstats[CHROMs[[cri]]], 
             file = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/GTEx_Analysis_v8_eQTL_all_associations/", tissue, "_", cri, ".txt.gz"))
    }
    
    rm(eQTL_sumstats)
    gc()
  }
}

#need to retrieve independent genome-wide significant SNPs
# As instrumental variables, we used independent (r2 < 0.01) significant (PGWAS < 5 × 10−08) 
# SNPs chosen among the 10 K preselected trait-associated SNPs included in a trans-eQTL dataset 
# from eQTLGen Consortium (31,684 whole blood samples). As we are using only strongly independent 
# SNPs, we use the identity matrix to approximate C. The SNPs with larger effects on the outcome
# than on the exposure were removed, as these would indicate a violation of MR assumptions 
# (likely reverse causality and/or confounding).

#filter gwas for genome-wide significance (5E-8) and independence
for(gi in 1:length(gwas_summary_files)){
  cat(paste0(" (", gi, " / ", length(gwas_summary_files), ")")) 
  gwas_sumstats <- fread(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/gwas_sumstats/proper_format/signed/", gwas_summary_files[gi], ".sumstats.gz"))
  gwas_sumstats <- gwas_sumstats[!(gwas_sumstats$A1 == "" | is.na(gwas_sumstats$Z)),]
  #A1 = REF
  signif_gwas_sumstats <- gwas_sumstats[abs(gwas_sumstats$Z) > -(qnorm(2.5E-8)),]
  signif_gwas_sumstats$pval <- abs(pnorm(signif_gwas_sumstats$Z) - c(0,1)[as.numeric(sign(signif_gwas_sumstats$Z) == 1)+1]) * 2
  #can get beta and SE from z-score and sample size
  
  
}

#process all the input files needed for revTWMR
for(tissue in names(motrpac_gtex_map)){
  
  print(tissue)
  
  for(cri in c(1:22, "X")){
    
    cat(paste0(" ", cri))
    
    RSID_POS_MAP <- fread(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/RSID_POS_MAP_",cri,".txt"))
    eQTL_sumstats <- fread(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/GTEx_Analysis_v8_eQTL_all_associations/", tissue, "_", cri, ".txt.gz"))
    vids <- eQTL_sumstats$variant_id
    vids <- substr(vids, start = 4, nchar(vids))
    vids <- strsplit(vids, "_", T)
    vids <- as.data.table(data.frame(data.table::transpose(vids)))
    colnames(vids) <- c("CHROM", "POS", "REF", "ALT", "BUILD")
    vids$RSID <- RSID_POS_MAP$ID[match(vids$POS, RSID_POS_MAP$POS)]
    eQTL_sumstats <- cbind(eQTL_sumstats, vids)
    eQTL_sumstats <- eQTL_sumstats[!is.na(eQTL_sumstats$RSID)]
    eQTL_sumstats$ENSG <- gsub('\\..*','',eQTL_sumstats$gene_id)
    
    for(gi in 1:length(gwas_summary_files)){
      cat(paste0(" (", gi, " / ", length(gwas_summary_files), ")")) 
      gwas_sumstats <- fread(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/gwas_sumstats/proper_format/signed/", gwas_summary_files[gi], ".sumstats.gz"))
      gwas_sumstats <- gwas_sumstats[!(gwas_sumstats$A1 == "" | is.na(gwas_sumstats$Z)),]
      #A1 = REF
      
      #need to retrieve independent SNPs
      # As instrumental variables, we used independent (r2 < 0.01) significant (PGWAS < 5 × 10−08) 
      # SNPs chosen among the 10 K preselected trait-associated SNPs included in a trans-eQTL dataset 
      # from eQTLGen Consortium (31,684 whole blood samples). As we are using only strongly independent 
      # SNPs, we use the identity matrix to approximate C. The SNPs with larger effects on the outcome
      # than on the exposure were removed, as these would indicate a violation of MR assumptions 
      # (likely reverse causality and/or confounding).
      signif_gwas_sumstats <- gwas_sumstats[abs(gwas_sumstats$Z) > -(qnorm(2.5E-8)),]
      
      #get LD correlations
      change_dir_command <- paste0("cd /Volumes/SSD500GB/all.EUR.hg38a/; ")
      initialize_profile <- paste0("source ~/.zshenv; ")
      run_plink_command <- paste0("plink --r gz --bfile all.EUR.hg38a --ld-window 1000000 --ld-window-r2 0 --ld-window-kb 10000 --extract ", 
                                  paste0("/Volumes/SSD500GB/all.EUR.hg38a/", gwas_traitname,"/chr", cri, "/", 
                                         gene_name,"_", as.integer(window_size),"BP-Window_SNPs.txt"),
                                  " --out ", paste0("/Volumes/SSD500GB/all.EUR.hg38a/", gwas_traitname,"/chr", cri, "/", 
                                                    gene_name,"_", as.integer(window_size),"BP-Window_LD; "))
      system_command <- paste0(change_dir_command, initialize_profile, run_plink_command)
      system(system_command)
      

      example_betaGWAS[1:3, ..bgw_cols]
      example_genes.N
      #for genes in 
      DEGs
      
      
      
    }
    
    #clean up memory
    rm(vids)
    rm(RSID_POS_MAP)
    rm(eQTL_sumstats)
    gc()
    
  }
  
}


for(gi in 1:length(gwas_summary_files)){
  cat(paste0(" (", gi, " / ", length(gwas_summary_files), ")")) 
  gwas_sumstats <- fread(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/ldsc/gwas_sumstats/proper_format/signed/", gwas_summary_files[gi], ".sumstats.gz"))
  gwas_sumstats <- gwas_sumstats[!(gwas_sumstats$A1 == "" | is.na(gwas_sumstats$Z)),]
  #A1 = REF
}
