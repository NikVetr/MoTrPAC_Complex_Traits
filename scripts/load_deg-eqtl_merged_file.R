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

# get tissue labels 
# load MotrpacBicQC conversion tables 
bic_animal_tissue_code = as.data.table(MotrpacBicQC::bic_animal_tissue_code)
bic_animal_tissue_code = bic_animal_tissue_code[tissue_name_release!='']
bic_animal_tissue_code[,my_tissue := tolower(gsub(" Powder", "", bic_tissue_name))]
bic_animal_tissue_code[,my_tissue := gsub(' ','_',my_tissue)]
bic_animal_tissue_code[my_tissue == 'blood_rna', my_tissue := 'paxgene_rna']
tissue_codes = bic_animal_tissue_code[, tissue_name_release]
tissue_codes = tissue_codes[tissue_codes != '']

## Define data paths

#download data to local drive
# system("scp nikgvetr@smsh11dsu-srcf-d15-38.scg.stanford.edu:/oak/stanford/groups/smontgom/shared/motrpac/mawg_data/pass1b-06/transcript-rna-seq/mapping/pass1b-06_transcript-rna-seq_feature-mapping_20201002.txt /Users/nikolai/data/smontgom/")

GTEx_logo <- readJPEG("~/Documents/Documents - nikolai/GTEx_logo.jpg")
eqtl = '~/data/smontgom/GTEx_Analysis_v8_eQTL'
#deg = 'gs://mawg-data/pass1b-06/transcript-rna-seq/dea/'
deg = '~/data/smontgom/dea/'
#map = dl_read_gcp('gs://mawg-data/pass1b-06/transcript-rna-seq/mapping/pass1b-06_transcript-rna-seq_feature-mapping_20201002.txt', sep='\t')
# map = fread('~/data/smontgom/pass1b-06_transcript-rna-seq_feature-mapping_20201002.txt', sep='\t', header=T)
# map = fread('~/data/smontgom/pass1b-06_transcript-rna-seq_feature-mapping_20210721.txt', sep='\t', header=T)
gene_map <- fread("~/data/smontgom/gencode.v39.RGD.20201001.human.rat.gene.ids.txt")
gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID <- gsub(gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID, pattern = "\\..*", replacement = "")
map <- unique(gene_map[,c("RAT_ENSEMBL_ID", "HUMAN_ORTHOLOG_ENSEMBL_ID", "HUMAN_ORTHOLOG_SYMBOL")])
colnames(map) <- c("feature_ID", "human_ensembl_gene", "human_gene_symbol")

gwas = '~/data/smontgom/imputed_gwas_hg38_1.1'
coloc = '~/data/smontgom/results_enloc_priors'

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

## Download and format RNA-seq differential expression analysis results for easy loading
# timewise_dea_list = list()
# training_dea_list = list()
# # for (file in system(sprintf('gsutil ls %s | grep "transcript-rna-seq_dea_20201028.txt"', deg), intern=T)){
# #   tissue_code = gsub('_.*','',gsub('.*pass1b-06_','',file))
# #   dea_list[[tissue_code]] = dl_read_gcp(file, sep='\t')
# # }
# for (file in list.files(deg, full.names=T, pattern='timewise-dea')){
#   tissue_code = gsub('_.*','',gsub('.*pass1b-06_','',file))
#   timewise_dea_list[[tissue_code]] = fread(file, sep='\t', header=T)
# }
# for (file in list.files(deg, full.names=T, pattern='training-dea')){
#   tissue_code = gsub('_.*','',gsub('.*pass1b-06_','',file))
#   training_dea_list[[tissue_code]] = fread(file, sep='\t', header=T)
# }
# training_dea_list = lapply(training_dea_list, function(x) x[,removed_samples := as.character(removed_samples)])
# timewise_dea_list = lapply(timewise_dea_list, function(x) x[,removed_samples := as.character(removed_samples)])
# training_dea = rbindlist(training_dea_list)
# timewise_dea = rbindlist(timewise_dea_list)
# # remove aorta (BAT contamination in F 1w,2w samples)
# training_dea = training_dea[tissue!='t65-aorta']
# timewise_dea = timewise_dea[tissue!='t65-aorta']
# # adjust tests across all tissues
# training_dea[,adj_p_value := p.adjust(p_value, method='BH')]
# table(training_dea[,adj_p_value] < 0.05, training_dea[,tissue])
# rna_dea = list(training_dea = training_dea,
#                timewise_dea = timewise_dea)
# save(rna_dea, file='/oak/stanford/groups/smontgom/shared/motrpac/shared_rdata/rna_dea_20210114.RData')

if(!file.exists("~/data/smontgom/rna_dea_20210114.RData")){
  system("scp nikgvetr@smsh11dsu-srcf-d15-38.scg.stanford.edu:/oak/stanford/groups/smontgom/shared/motrpac/shared_rdata/rna_dea_20210114.RData /Users/nikolai/data/smontgom/",
         show.output.on.console = TRUE)
}

if(!exists("rna_dea")){
  load('~/data/smontgom/rna_dea_20210114.RData')
  load('~/data/smontgom/dea/transcript_rna_seq_20210804.RData')
  rna_dea$training_dea <- as.data.table(transcript_rna_seq$training_dea)
  rna_dea$timewise_dea <- as.data.table(transcript_rna_seq$timewise_dea)
  rm(transcript_rna_seq)
}
# list with two data.tables:
# training_dea: use `adj_p_value` to select genes that are differential due to training (use 0.1 (10% FDR) threshold for now)
# timewise_dea: cross-reference list of DEGs with this table to ID sex- and time- specific log fold-changes 

## Intersect DA results with GTEx eQTLs 
if(!exists("deg_eqtl_list")){
  deg_eqtl_list = list()
  for(motrpac_tissue in unique(rna_dea$timewise_dea$tissue)){
    
    if(!motrpac_tissue %in% names(motrpac_gtex_map)){next}
    
    cat(paste0(motrpac_tissue, "\n"))
    
    # read in eQTLs
    gtex_tissue = motrpac_gtex_map[[motrpac_tissue]]
    gtex_egene = fread(sprintf('%s/%s.v8.egenes.txt.gz',eqtl,gtex_tissue), sep='\t', header=T)
    gtex_egene[, human_ensembl_gene := gsub('\\..*','',gene_id)]
    # gtex_egene$human_gene_symbol <- gtex_egene$gene_name
    
    # match human genes with rat ensembl genes 
    gtex_egene = merge(gtex_egene, map, by='human_ensembl_gene')
    gtex_motrpac = merge(x = gtex_egene, y = rna_dea$timewise_dea[tissue == motrpac_tissue], by='feature_ID')
    gtex_motrpac$abs_slope <- abs(gtex_motrpac$slope)
    
    cat(paste0("prop of feature_IDs matched: ", round(length(unique(gtex_motrpac$feature_ID)) / 
                                                        length(unique(rna_dea$timewise_dea[tissue == motrpac_tissue]$feature_ID)), 3),
               "\nprop gene symbols in map: ", round(mean(gtex_egene$human_gene_symbol %in% map$human_gene_symbol), 2), "\n"))
    
    deg_eqtl_list[[motrpac_tissue]] = gtex_motrpac
    # gtex_motrpac$gene_symbol
    
    # de-duplicate gene entries, slowly but straightforwardly, while also preserving ensemble IDs
    # dup_gene_IDs <- table(gtex_motrpac$gene_id)
    # dup_gene_IDs <- names(dup_gene_IDs)[dup_gene_IDs > ifelse(any(motrpac_tissue == c("t63-testes", "t64-ovaries")), 4, 8)]
    # for(dgid in dup_gene_IDs){
    #   dgid_inds <- which(gtex_motrpac$gene_id == dgid)
    #   # print(length(unique(gtex_motrpac[dgid_inds,]$slope)))
    #   # print(length(unique(gtex_motrpac[dgid_inds,]$logFC)))
    #   a <- deg_eqtl_list[[tss]][deg_eqtl_list[[tss]]$gene_id == dgid,]
    #   print(paste0(dgid, " has a max unique entry count of ", max(sapply(colnames(a), function(i) length(unique(a[[i]]))))))
    # }
  }
}

# gene_ID_map <- unique.data.frame(as.data.frame(do.call(rbind, deg_eqtl_list)[,c("feature_ID","human_ensembl_gene")]))
# write.table(gene_ID_map, file = "~/data/smontgom/motrpac_geneID_map.txt", col.names = T)
# deg_eqtl_list[["t68-liver"]]
# for(i in 1:length(names(motrpac_gtex_map))){
#   print(nrow(deg_eqtl_list[[names(motrpac_gtex_map)[i]]]))
# }
# 
# length(unique(deg_eqtl_list[[motrpac_tissue]]$gene_id))
# length(deg_eqtl_list[[motrpac_tissue]]$gene_id)
# 
# deg_eqtl_list[[motrpac_tissue]]$logFC
# deg_eqtl_list[[motrpac_tissue]]$slope
# 
# ind1 <- sample(1:nrow(deg_eqtl_list[[motrpac_tissue]]), 1)
# ind2 <- sample(which(deg_eqtl_list[[motrpac_tissue]]$gene_id == deg_eqtl_list[[motrpac_tissue]]$gene_id[ind1]), 1)
# if(length(ind2) > 0){
#   colnames(deg_eqtl_list[[motrpac_tissue]])[which(abs(as.numeric(deg_eqtl_list[[motrpac_tissue]][ind1,]) -
#                                                       as.numeric(deg_eqtl_list[[motrpac_tissue]][ind2,])) > 1E-10)]
# }

#prop of duplicate gene names
# sapply(names(motrpac_gtex_map), function(name) length(deg_eqtl_list[[name]]$gene_id) / length(unique(deg_eqtl_list[[name]]$gene_id)))
# 
# sapply(names(motrpac_gtex_map), function(name) length(deg_eqtl_list[[name]]$gene_id) / length(unique(deg_eqtl_list[[name]]$gene_id)))
# 
# sapply(1:ncol(deg_eqtl_list[[tss]][deg_eqtl_list[[tss]]$gene_id == "ENSG00000070087.13",]), 
#        function(x) length(unique(deg_eqtl_list[[tss]][deg_eqtl_list[[tss]]$gene_id == "ENSG00000070087.13",x])))
# 
#some quick EDA
# counts <- table(deg_eqtl_list[[tss]]$gene_id)
# ind <- sample(which(counts > 8), 1)
# a <- deg_eqtl_list[[tss]][deg_eqtl_list[[tss]]$gene_id == names(ind),]
# max(sapply(colnames(a), function(i) length(unique(a[[i]]))))
# 
# for(i in colnames(a)){
# 
# }
# dev.off()
# par(mfrow = c(4,4))
# a = 0
# for(tss in names(motrpac_gtex_map)){
#   counts <- table(deg_eqtl_list[[tss]]$gene_id)
#   counts <- counts[counts > 8]
#   print(paste0(tss, " = ", length(counts)))
# }
# 
# 
# tss <- sample(names(motrpac_gtex_map), size = 1)
# par(mfrow = c(4,4))
# for(tss in names(motrpac_gtex_map)){
# inds <- sample(x = length(deg_eqtl_list[[tss]]$logFC), 1E4, replace = F)
# plot(deg_eqtl_list[[tss]]$logFC[inds], deg_eqtl_list[[tss]]$slope[inds], main = tss)
# # plot(deg_eqtl_list[[tss]]$shrunk_logFC[inds], deg_eqtl_list[[tss]]$logFC[inds], main = tss)
# }



# this now includes the intersection of GTEx and MoTrPAC expressed genes annotated in both species
# there are some duplicate gene IDs (both human and rat). figure out how to deal with these
# direction of DEG: logFC
# direction of eQTL: "slope" ("variant_id", "alt")
#save(deg_eqtl_list, file='rdata/deg_eqtl_list_20201109.RData')

calc_gene_intersect = F
if(calc_gene_intersect){
  sapply(names(deg_eqtl_list), function(tissue_1) sapply(names(deg_eqtl_list), function(tissue_2) 
    round(length(intersect(deg_eqtl_list[[tissue_1]]$gene_name, deg_eqtl_list[[tissue_2]]$gene_name)) / 
            length(union(deg_eqtl_list[[tissue_1]]$gene_name, deg_eqtl_list[[tissue_2]]$gene_name)) * 100, 2)))
}

# add in merged gonads
# if(!any(names(deg_eqtl_list) == "t1000-gonads")){
#   deg_eqtl_list$`t1000-gonads` = rbind(deg_eqtl_list$`t63-testes`, deg_eqtl_list$`t64-ovaries`)
# }

#compare old and new rna_dea files
# plot(quantile(old_rna_dea$timewise_dea$logFC, 1:99 / 100), quantile(rna_dea$timewise_dea$logFC, 1:99 / 100))



cols = list(Tissue=MotrpacBicQC::tissue_cols[names(deg_eqtl_list)], 
            Time=MotrpacBicQC::group_cols[paste0(c(1,2,4,8), "w")],
            Sex=MotrpacBicQC::sex_cols[c('male','female')])
# cols$Tissue[which(is.na(cols$Tissue))] <- '#C0C0C0'
# names(cols$Tissue)[which(is.na(names(cols$Tissue)))] <- "t1000-gonads"
