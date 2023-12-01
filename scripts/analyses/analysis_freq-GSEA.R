#load in some metadata
trait_categories <- read.csv("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/gwas_metadata.csv", header = T)
trait_categories$Category[trait_categories$new_Phenotype == "Multiple_Sclerosis"] <- "Psychiatric-neurologic" #correct original coding from "Cardiometabolic"

library(fgsea)
trait_name_key <- setNames(trait_categories$new_Phenotype, trait_categories$Tag)

#make directory
gsea_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/gsea_output/"
if(!dir.exists(gsea_dir)){dir.create(gsea_dir)}

if(file.exists(paste0(gsea_dir, "gsea_trait_x_tissue.txt"))){
  gsea_output <- fread(paste0(gsea_dir, "gsea_trait_x_tissue.txt"))
  gsea_output <- gsea_output[,-"leadingEdge"]
} else {
  
  gsea_output <- do.call(rbind, lapply(tissues, function(tiss){
    
    
    gsea_results <- parallel::mclapply(split(twas_by_tissue[[tiss]], twas_by_tissue[[tiss]]$trait), 
                                       function(trait_twas){
                                         
                                         tryCatch({
                                           
                                           trait <- trait_twas$trait[1]
                                           mcprint(paste0(which(tissues == tiss), ".", tiss, ".", trait))
                                           trait_twas$log10_pvalue <- -log10(trait_twas$pvalue)
                                           ranked_genes <- setNames(object = trait_twas$log10_pvalue, 
                                                                    nm = gsub("\\..*", "", trait_twas$gene))
                                           DE_genes <- list(intersect(DE_genes_x_tissue[[tiss]], names(ranked_genes)))
                                           names(DE_genes) <- paste0(tiss, ".", trait)
                                           gsea_out <- fgsea(pathways = DE_genes, stats = ranked_genes, scoreType = "pos")
                                           
                                           if("data.frame" %in% class(gsea_out)){
                                             gsea_out$tissue <- tiss
                                             gsea_out$trait <- trait  
                                           }
                                           
                                           gsea_out
                                           
                                         }, error = function(e){
                                           print(paste("Error msg: ", e$message))
                                           return(NULL)
                                         })
                                       }, mc.cores = 8)
    
    # gsea_results <- gsea_results[lapply(gsea_results, class) != "try-error"]
    gsea_results <- do.call(rbind, gsea_results)
    gsea_results
  }))
  
  gsea_output$trait_name <- trait_name_key[gsea_output$trait]
  gsea_output$pathway_name <- paste0(gsea_output$tissue, ".", gsea_output$trait_name)
  hist(gsea_output$pval, breaks = 0:100/100)
  # gsea_output$adj.pval <- p.adjust(gsea_output$pval, method = "bh")
  gsea_output$adj.pval <- p.adjust(gsea_output$pval, method = "bonf")
  gsea_output <- gsea_output[,-"leadingEdge"]
  
  fwrite(gsea_output, paste0(gsea_dir, "gsea_trait_x_tissue.txt"))
  
}

hist(gsea_output$adj.pval, breaks = 0:100/100)
head(gsea_output[order(gsea_output$pval),c("trait_name", "tissue", "pval", "adj.pval")], 20)

#try meta-analysis (using HMP?)
library(harmonicmeanp)
L <- nrow(gsea_output)
tissue_ps <- sapply(split(gsea_output$pval, gsea_output$tissue), function(tissps){
  (length(tissps)/L)/sum((1/L)/tissps)
})
tissue_thresh <- qharmonicmeanp(0.05, L) * table(gsea_output$tissue)/L
tissue_thresh <- tissue_thresh[names(tissue_ps)]
any(tissue_ps < tissue_thresh)
head(sort(tissue_ps), n = 10)

trait_ps <- sapply(split(gsea_output$pval, gsea_output$trait_name), function(traitps){
  (length(traitps)/L)/sum((1/L)/traitps)
})
trait_thresh <- qharmonicmeanp(0.05, L) * table(gsea_output$trait_name)/L
trait_thresh <- trait_thresh[names(trait_ps)]
any(trait_ps < trait_thresh)
head(sort(trait_ps), n = 10)
sort(trait_ps / (trait_thresh / 0.05))


trait_cat_key <- setNames(trait_categories$Category, 
                          trait_categories$new_Phenotype)[names(trait_ps)]
L_traits <- length(trait_ps)
trait_category_ps <- sapply(split(trait_ps, trait_cat_key), function(traitcatps){
  (length(traitcatps)/L_traits)/sum((1/L_traits)/traitcatps)
})
trait_category_thresh <- qharmonicmeanp(0.05, L_traits) * 
  table(trait_cat_key) / L_traits
trait_category_thresh <- trait_category_thresh[names(trait_category_ps)]
any(trait_category_ps < trait_category_thresh)
sort(trait_category_ps / (trait_category_thresh / 0.05))

sort(trait_category_ps)

# #try just set unioning the tissue DE genes and combining the predixcan pvals? with HMP also?
tissue_x_twas <- lapply(setNames(unique(some.twas$trait),
                                 unique(some.twas$trait)),
                        function(traiti) {
                          print(traiti)
                          some.twas[some.twas$trait == traiti & compatible_twas_genes,]
                        })
# all_DE_genes <- unlist(lapply(seq_along(DE_genes_x_tissue), function(i) paste0(names(DE_genes_x_tissue)[i], ".", DE_genes_x_tissue[[i]])))
# gsea_output_pooled <- parallel::mclapply(tissue_x_twas, 
#                                          function(trait_twas){
#                                            
#                                            tryCatch({
#                                              
#                                              trait <- trait_twas$trait[1]
#                                              mcprint(trait)
#                                              
#                                              trait_twas$log10_pvalue <- -log10(trait_twas$pvalue)
#                                              ranked_genes <- setNames(object = trait_twas$log10_pvalue, 
#                                                                       nm = paste0(tissue_abbr[trait_twas$tissue], ".", gsub("\\..*", "", trait_twas$gene)))
#                                              # DE_genes <- list(intersect(unique(unlist(DE_genes_x_tissue)), names(ranked_genes)))
#                                              # DE_genes <- list(intersect(names(table(unlist(DE_genes_x_tissue))[table(unlist(DE_genes_x_tissue)) >= 3]), 
#                                              #                            names(ranked_genes)))
#                                              
#                                              #only use 1-to-1 genes here?
#                                              DE_genes <- list(intersect(all_DE_genes, names(ranked_genes)))
#                                              names(DE_genes) <- trait
#                                              gsea_out <- fgseaMultilevel(pathways = DE_genes, 
#                                                                          stats = ranked_genes, 
#                                                                          scoreType = "pos",
#                                                                          eps = 1E-10)
#                                              
#                                              if("data.frame" %in% class(gsea_out)){
#                                                gsea_out$trait <- trait  
#                                              }
#                                              
#                                              return(gsea_out)
#                                              
#                                            }, error = function(e){
#                                              print(paste("Error msg: ", e$message))
#                                              return(NULL)
#                                            })
#                                            
#                                          }, mc.cores = 8)
# 
# gsea_output_pooled <- gsea_output_pooled[lapply(gsea_output_pooled, class) != "NULL"]
# gsea_output_pooled <- do.call(rbind, gsea_output_pooled)
# gsea_output_pooled$trait_name <- setNames(trait_categories$new_Phenotype, trait_categories$Tag)[
#   gsea_output_pooled$trait]
# gsea_output_pooled$trait_fullname <- setNames(trait_categories$Phenotype, trait_categories$Tag)[
#   gsea_output_pooled$trait]
# gsea_output_pooled$adj.pval <- p.adjust(gsea_output_pooled$pval, method = "BH")
# hist(gsea_output_pooled$adj.pval)
# head(gsea_output_pooled[order(gsea_output_pooled$pval),
#                         c("trait_name", "adj.pval")], 20)

#or try combining predixcan probabilities?
if(file.exists(paste0(gsea_dir, "gsea_trait_pooled.txt"))){
  gsea_output_pooled_alt <- fread(paste0(gsea_dir, "gsea_trait_pooled.txt"))
} else {
  gsea_output_pooled_alt <- parallel::mclapply(tissue_x_twas, 
                                               function(trait_twas){
                                                 
                                                 tryCatch({
                                                   
                                                   trait <- trait_twas$trait[1]
                                                   mcprint(trait)
                                                   
                                                   trait_twas$log10_pvalue <- -log10(trait_twas$pvalue)
                                                   trait_twas$ensembl <- gsub("\\..*", "", trait_twas$gene)
                                                   tissue_x_genes <- split(trait_twas$pvalue, trait_twas$ensembl)
                                                   L_txg <- length(unlist(tissue_x_genes))
                                                   ranked_genes <- -log10(sapply(tissue_x_genes, function(x){
                                                     (length(x)/L_txg)/sum((1/L_txg)/x)
                                                   }))
                                                   DE_genes <- list(intersect(names(table(unlist(DE_genes_x_tissue))[table(unlist(DE_genes_x_tissue)) >= 3]),
                                                                              names(ranked_genes)))
                                                   names(DE_genes) <- trait
                                                   gsea_out <- fgseaMultilevel(pathways = DE_genes, 
                                                                               stats = ranked_genes, 
                                                                               scoreType = "pos",
                                                                               eps = 1E-10)
                                                   
                                                   if("data.frame" %in% class(gsea_out)){
                                                     gsea_out$trait <- trait  
                                                   }
                                                   
                                                   return(gsea_out)
                                                   
                                                 }, error = function(e){
                                                   print(paste("Error msg: ", e$message))
                                                   return(NULL)
                                                 })
                                                 
                                               }, mc.cores = 8)
  
  gsea_output_pooled_alt <- gsea_output_pooled_alt[lapply(gsea_output_pooled_alt, class) != "NULL"]
  gsea_output_pooled_alt <- do.call(rbind, gsea_output_pooled_alt)
  gsea_output_pooled_alt$trait_name <- setNames(trait_categories$new_Phenotype, trait_categories$Tag)[
    gsea_output_pooled_alt$trait]
  gsea_output_pooled_alt$adj.pval <- p.adjust(gsea_output_pooled_alt$pval, method = "BH")
  
  fwrite(gsea_output_pooled_alt, paste0(gsea_dir, "gsea_trait_pooled.txt"))
  
  
}
hist(gsea_output_pooled_alt$adj.pval)
head(gsea_output_pooled_alt[order(gsea_output_pooled_alt$pval),
                            c("trait_name", "pval", "adj.pval")], 20)


#ANOTHER IDEA: create a supplemental figure that compares signed negative log10(pvalues) to log-odds of positive effect

#### supp fig 3 preprocessing ####
#first load in the MCMC output

base = paste0("deviation_from_expected_logodds_split_the_difference_informative-MVN-matrix-normal_priors") 
load(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/mcmc_output/", 
            base, 
            ifelse(T, "_allCats", ""), 
            ifelse(F, "_randomgenes", ""),
            ifelse(T, "", "_intersecting-genes-marginal-corrmats"),
            ifelse(T, "_pairwise-interactions", ""),
            ".cmdStanR.fit"))
samps <- data.frame(as_draws_df(out$draws()))
prop_greater_than_0 <- function(x) mean(x>0)

celltotalbias <- apply(subset_samps("cell_total_prob_bias", c("raw", "sd", "logodds", "L"), samps = samps), 2, prop_greater_than_0)
celltotalbias <- setNames(celltotalbias, paste0(data_subset$tissue, ".", 
                                                trait_name_key[data_subset$trait]))
celltotalbias <- sort(celltotalbias, T)
#adjust 0s and 1s to min and max non-0 and non-1 values
celltotalbias <- celltotalbias + (celltotalbias == 0) * 1/nrow(samps) + (celltotalbias == 1) * -1/nrow(samps)
head(celltotalbias, 20)

rowbias <- apply(subset_samps("row_bias", c("raw", "sd", "L"), samps = samps) + samps$overall_bias, 2, prop_greater_than_0)
rowbias <- setNames(sort(rowbias), tissues[order(rowbias)])
rowbias

colbias <- apply(subset_samps("col_bias", c("raw", "sd", "L"), samps = samps) + 
                   samps$overall_bias +
                   subset_samps("colcat_bias", c("raw", "sd", "L"), 
                                samps = samps)[,match(traitwise_partitions$Category[
                                  match(traits, traitwise_partitions$Tag)], trait_cats)], 
                 2, prop_greater_than_0)
trait_key <- setNames(trait_categories$new_Phenotype, trait_categories$Tag)
colbias <- setNames(sort(colbias), trait_key[traits[order(colbias)]])
tail(colbias, 20)

# colcatbias <- apply(subset_samps("colcat_bias", c("raw", "sd", "L"), samps = samps), 2, prop_greater_than_0)
colcatbias <- apply(subset_samps("colcat_bias", c("raw", "sd", "L"), samps = samps) + samps$overall_bias, 
                    2, prop_greater_than_0)
colcatbias <- setNames(sort(colcatbias), trait_cats[order(colcatbias)])
colcatbias
