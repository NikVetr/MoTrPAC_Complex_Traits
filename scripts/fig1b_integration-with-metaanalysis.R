#libraries
library(ks)
library(arrow)
library(edgeR)
library(limma)
library(DESeq2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(org.Rn.eg.db)
library(biomaRt)
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
library(data.table)
# library(MotrpacBicQC)
library(MotrpacRatTraining6mo) # v1.6.0
# also attaches MotrpacRatTraining6moData v1.8.0


# TODO: add these scripts to this repo 
#load de data
source(file = "~/scripts/montgomery_lab/load_deg-eqtl_merged_file.R")

#load functions
source(file = "~/scripts/montgomery_lab/deg-trait_functions.R")


#rat <-> human gene map 
# gencode_gene_map <- rdg_mapping <- fread("~/data/smontgom/gencode.v39.RGD.20201001.human.rat.gene.ids.txt")
gencode_gene_map <- rdg_mapping <- MotrpacRatTraining6moData::RAT_TO_HUMAN_GENE

# load("~/data/smontgom/genes_tested_in_transcriptome_DEA.RData")
genes_tested_in_transcriptome_DEA <- unique(unlist(MotrpacRatTraining6moData::GENE_UNIVERSES$ensembl_gene$TRNSCRPT))
gencode_gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID <- gsub(gencode_gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID, pattern = "\\..*", replacement = "")
gene_map <- gencode_gene_map
all_orthologs_tested <- gene_map$HUMAN_ORTHOLOG_SYMBOL[match(genes_tested_in_transcriptome_DEA, gene_map$RAT_ENSEMBL_ID)]
all_orthologs_tested <- all_orthologs_tested[!is.na(all_orthologs_tested)]


tissues <- c("t55-gastrocnemius", "t56-vastus-lateralis", "t30-blood-rna")
tissues <- setNames(tissues,tissues)


#### snag & process meta-analysis results ####
load("~/data/smontgom/meta_analysis_results.RData")

mouse_symb_to_human = tapply(rdg_mapping$HUMAN_ORTHOLOG_SYMBOL,
                             rdg_mapping$RAT_SYMBOL,
                             function(x)x[1])
human_entrez_to_symbol = tapply(rdg_mapping$HUMAN_ORTHOLOG_SYMBOL,
                                as.character(rdg_mapping$HUMAN_ORTHOLOG_NCBI_GENE_ID),
                                function(x)x[1])

meta_d = all_meta_analysis_res$`longterm,muscle`
meta_d_base_models = lapply(
  meta_d,
  function(x)x[[which(grepl("base",names(x)))[1]]]
)
# extract the simple base model, useful for I2 scores
meta_d_base_models_simple = lapply(meta_d,function(x)x[["simple:base_model"]])
# create a summary data frame of the human data
human_muscle_res = sapply(meta_d_base_models,
                          function(x)c(x$coeffs[1:2],x$mod_p,x$het_p,x$I2))
human_muscle_res = t(human_muscle_res)
colnames(human_muscle_res) = c("beta","beta_se","p","het_p","I2")
human_muscle_res = as.data.frame(human_muscle_res)
# extract the I2 scores from the simple models
human_muscle_res$I2 = sapply(meta_d_base_models_simple,
                             function(x){
                               if(is.null(x$I2) || is.na(x$I2)){return(0)}
                               return(x$I2)
                             }
)
# map entrez ids to gene symbols, will be used to merge with the 
# rat results
human_muscle_res$symbol = human_entrez_to_symbol[rownames(human_muscle_res)]
human_muscle_res = human_muscle_res[!is.na(human_muscle_res$symbol),]
rownames(human_muscle_res) = human_muscle_res$symbol
human_muscle_res$adj_p <- p.adjust(human_muscle_res$p, method = "BH")
hist(human_muscle_res$adj_p, breaks = 100)


meta_d = all_meta_analysis_res$`longterm,blood`
meta_d_base_models = lapply(
  meta_d,
  function(x)x[[which(grepl("base",names(x)))[1]]]
)
# extract the simple base model, useful for I2 scores
meta_d_base_models_simple = lapply(meta_d,function(x)x[["simple:base_model"]])
# create a summary data frame of the human data
human_blood_res = sapply(meta_d_base_models,
                         function(x)c(x$coeffs[1:2],x$mod_p,x$het_p,x$I2))
human_blood_res = t(human_blood_res)
colnames(human_blood_res) = c("beta","beta_se","p","het_p","I2")
human_blood_res = as.data.frame(human_blood_res)
# extract the I2 scores from the simple models
human_blood_res$I2 = sapply(meta_d_base_models_simple,
                            function(x){
                              if(is.null(x$I2) || is.na(x$I2)){return(0)}
                              return(x$I2)
                            }
)
# map entrez ids to gene symbols, will be used to merge with the 
# rat results
human_blood_res$symbol = human_entrez_to_symbol[rownames(human_blood_res)]
human_blood_res = human_blood_res[!is.na(human_blood_res$symbol),]
rownames(human_blood_res) = human_blood_res$symbol
human_blood_res$adj_p <-p.adjust(human_blood_res$p, method = "BH")
hist(human_blood_res$adj_p, breaks = 100)

#process to gene symbols
alpha_fdr <- 0.05
i2_filter <- 101
blood_genes <- human_blood_res$symbol[(human_blood_res$adj_p < alpha_fdr) & (human_blood_res$I2 < i2_filter)]
muscle_genes <- human_muscle_res$symbol[(human_muscle_res$adj_p < alpha_fdr) & (human_muscle_res$I2 < i2_filter)]
blood_signs <- setNames(sign(human_blood_res$beta[human_blood_res$adj_p < alpha_fdr & (human_blood_res$I2 > i2_filter)]), blood_genes)
muscle_signs <- setNames(sign(human_muscle_res$beta[human_muscle_res$adj_p < alpha_fdr & (human_muscle_res$I2 > i2_filter)]), muscle_genes)

#alternatively, use the longterm genes IDed in the paper
longterm_genes <- fread("~/data/smontgom/41467_2021_23579_MOESM6_ESM.csv")
longterm_genes <- cbind(longterm_genes, do.call(rbind, strsplit(longterm_genes$`Discovered in`, ",")))
colnames(longterm_genes)[colnames(longterm_genes) == "V1"] <- "timeperiod"
colnames(longterm_genes)[colnames(longterm_genes) == "V2"] <- "tissue"
longterm_genes <- longterm_genes[longterm_genes$timeperiod == "longterm",]
blood_genes <- longterm_genes$Symbol[longterm_genes$tissue == "blood"]
muscle_genes <- longterm_genes$Symbol[longterm_genes$tissue == "muscle"]
blood_signs <- setNames(sign(human_blood_res$beta[match(blood_genes, human_blood_res$symbol)]), blood_genes)
muscle_signs <- setNames(sign(human_muscle_res$beta[match(muscle_genes, human_muscle_res$symbol)]), muscle_genes)


#### get motrpac clustering results ####

# if(!exists("node_sets")){
#   load("~/data/smontgom/graphical_analysis_results_20211220.RData")  
# }
node_sets <- MotrpacRatTraining6moData::GRAPH_COMPONENTS$node_sets

# nodes_to_look_at_list <- list(c("1w_F1_M1", "1w_F-1_M-1"),
#                               c("2w_F1_M1", "2w_F-1_M-1"),
#                               c("4w_F1_M1", "4w_F-1_M-1"),
#                               c("8w_F1_M1", "8w_F-1_M-1"))

node_shorthand <- expand.grid(-1:1, -1:1)
colnames(node_shorthand) <- c("male", "female")
focal_node_sets <- list(both = which((node_shorthand$male == 1 & node_shorthand$female == 1) | 
                                       (node_shorthand$male == -1 & node_shorthand$female == -1)), 
                        male = which(node_shorthand$male != 0), 
                        female = which(node_shorthand$female != 0), 
                        neither = which(node_shorthand$male == 0 & node_shorthand$female == 0))
# focal_node_sets <- list(both = which((node_shorthand$Var1 == 1 & node_shorthand$Var2 == 1) | 
#                                        (node_shorthand$Var1 == -1 & node_shorthand$Var2 == -1)), 
#                         neither = which(node_shorthand$Var1 == 0 & node_shorthand$Var2 == 0))

nodes_to_look_at_list <- lapply(setNames(paste0(2^(0:3), "w"), paste0(2^(0:3), "w")), function(tpt){
  lapply(focal_node_sets, function(node_set){
    paste0(tpt, "_F", node_shorthand[node_set,"female"], "_M", node_shorthand[node_set,"male"])
  })
})

node_metadata_list <- lapply(setNames(paste0(2^(0:3), "w"), paste0(2^(0:3), "w")), function(timepoint){
  lapply(setNames(names(focal_node_sets), names(focal_node_sets)), function(node_type){
    
    node_metadata <- lapply(nodes_to_look_at_list[[timepoint]][[node_type]], function(node_to_look_at){
      cbind(do.call(rbind, strsplit(node_sets[node_to_look_at][[node_to_look_at]][grep(
        node_sets[node_to_look_at][[node_to_look_at]], pattern = "TRNSCRPT")], ";")), node_to_look_at)
    })
    
    node_metadata <- as.data.table(do.call(rbind, node_metadata))
    colnames(node_metadata) <- c("ome","tissue","ensembl_gene", "node")
    
    #snag signs
    node_metadata_sign <- do.call(rbind, strsplit(node_metadata$node, "_"))
    node_metadata_sign[grep("F-1", node_metadata_sign[,2]),2] <- -1
    node_metadata_sign[grep("F1", node_metadata_sign[,2]),2] <- 1
    node_metadata_sign[grep("F0", node_metadata_sign[,2]),2] <- 0
    node_metadata_sign[grep("M-1", node_metadata_sign[,3]),3] <- -1
    node_metadata_sign[grep("M1", node_metadata_sign[,3]),3] <- 1
    node_metadata_sign[grep("M0", node_metadata_sign[,3]),3] <- 0
    
    if(node_type == "male"){
      out <- setNames(object = as.integer(node_metadata_sign[,3]), node_metadata$ensembl_gene)
    }
    if(node_type == "female"){
      out <- setNames(object = as.integer(node_metadata_sign[,2]), node_metadata$ensembl_gene)
    }
    if(node_type == "both"){
      out <- setNames(object = as.integer(node_metadata_sign[,3]), node_metadata$ensembl_gene)
    }
    if(node_type == "neither"){
      out <- setNames(object = as.integer(node_metadata_sign[,3]), node_metadata$ensembl_gene)
    }
    
    names(out) <- gene_map$HUMAN_ORTHOLOG_SYMBOL[match(names(out), gene_map$RAT_ENSEMBL_ID)]
    
    out <- lapply(setNames(tissues, tissues), function(tiss){
      subout <- out[node_metadata$tissue == MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[tiss]]
      subout[!is.na(names(subout))]
    })
    
    out
  })
})


prop_in_motr <- lapply(setNames(paste0(2^(0:3), "w"), paste0(2^(0:3), "w")), function(timepoint){
  sapply(setNames(names(focal_node_sets), names(focal_node_sets)), function(node_type){
    sapply(setNames(tissues, tissues), function(tiss){
      motr_signs <- node_metadata_list[[timepoint]][[node_type]][[tiss]]
      if(tiss %in% c("t55-gastrocnemius", "t56-vastus-lateralis")){
        meta_signs <- muscle_signs
      } else {
        meta_signs <- blood_signs
      }
      # motr_signs <- motr_signs[motr_signs > 0]
      # meta_signs <- meta_signs[meta_signs > 0]
      
      
      paste0(length(intersect(names(motr_signs), names(meta_signs))),  " / ", length(motr_signs), " = ",
             round(length(intersect(names(motr_signs), names(meta_signs))) / length(motr_signs), 3))
      
      length(intersect(names(motr_signs), names(meta_signs))) / length(motr_signs)
    })
  })
})

#now do the sign matching thing with a sliding window for p-value
log10pval_window_size <- 1
min_log10pval <- -4
nwin <- 100
wins_center <- seq(1E-6, min_log10pval-1E-6, length.out = nwin)
wins <- cbind(wins_center + log10pval_window_size / 2, wins_center - log10pval_window_size / 2)
sexes = c("male", "female")
pval_window_output <- lapply(setNames(tissues, tissues), function(tiss){
  lapply(setNames(sexes, sexes), function(sex_i){
    lapply(setNames(paste0(2^(0:3), "w"), paste0(2^(0:3), "w")), function(timepoint){
      
      motr_signs <- deg_eqtl_list[[tiss]]
      motr_signs <- motr_signs[motr_signs$sex == sex_i & motr_signs$comparison_group == timepoint,c("human_gene_symbol", "p_value", "adj_p_value", "logFC")]
      motr_signs$sign <- sign(motr_signs$logFC)
      if(tiss %in% c("t55-gastrocnemius", "t56-vastus-lateralis")){
        meta_signs <- muscle_signs
      } else {
        meta_signs <- blood_signs
      }
      # motr_signs <- motr_signs[motr_signs > 0]
      shared_genes <- intersect(motr_signs$human_gene_symbol, names(meta_signs))
      motr_signs <- motr_signs[match(shared_genes, motr_signs$human_gene_symbol),]
      meta_signs <- meta_signs[shared_genes]
      signs <- cbind(motr_signs, meta_signs)
      signs$agree <- signs$sign * signs$meta_signs
      signs$log10p <- log10(signs$p_value)
      signs$log10p[signs$log10p < min_log10pval] <- min_log10pval
      
      prop_agree <- sapply(1:nwin, function(wi) mean(signs[signs$log10p < wins[wi,1] & signs$log10p > wins[wi,2], "agree"] == 1)) #can also use cut?
      prop_agree <- zoo::na.approx(prop_agree)      #linearly interpolate missing values
      prop_agree <-c(prop_agree, rep(prop_agree[length(prop_agree)], times = nwin-length(prop_agree)))
      
      #or use a laplace smoother?
      # dexp_smooth(x = signs$log10p[order(signs$log10p, decreasing = T)][1:1000], y = as.numeric(signs$agree == 1)[order(signs$log10p, decreasing = T)][1:1000], r = 5, 
      #             reweight_trunc_tail = F, reweight_n_pts = F, fix_endpoints = F, 
      #             interpolate_at = wins_center)
      
    })
  })
})

pval_window_output_both_sexes <- lapply(setNames(tissues, tissues), function(tiss){
  
    lapply(setNames(paste0(2^(0:3), "w"), paste0(2^(0:3), "w")), function(timepoint){
      
      motr_signs <- deg_eqtl_list[[tiss]]
      motr_signs_m <- motr_signs[motr_signs$sex == "male" & motr_signs$comparison_group == timepoint,c("human_gene_symbol", "p_value", "adj_p_value", "logFC")]
      motr_signs_f <- motr_signs[motr_signs$sex == "female" & motr_signs$comparison_group == timepoint,c("human_gene_symbol", "p_value", "adj_p_value", "logFC")]
      motr_signs_m$sign <- sign(motr_signs_m$logFC)
      motr_signs_f$sign <- sign(motr_signs_f$logFC)
      
      if(tiss %in% c("t55-gastrocnemius", "t56-vastus-lateralis")){
        meta_signs <- muscle_signs
      } else {
        meta_signs <- blood_signs
      }
      
      shared_genes_m <- intersect(motr_signs_m$human_gene_symbol, names(meta_signs))
      motr_signs_m <- motr_signs_m[match(shared_genes_m, motr_signs_m$human_gene_symbol),]
      meta_signs_m <- meta_signs[shared_genes_m]
      signs_m <- cbind(motr_signs_m, meta_signs = meta_signs_m)
      signs_m$agree <- signs_m$sign * signs_m$meta_signs
      signs_m$log10p <- log10(signs_m$p_value)
      signs_m$log10p[signs_m$log10p < min_log10pval] <- min_log10pval
      
      shared_genes_f <- intersect(motr_signs_f$human_gene_symbol, names(meta_signs))
      motr_signs_f <- motr_signs_f[match(shared_genes_f, motr_signs_f$human_gene_symbol),]
      meta_signs_f <- meta_signs[shared_genes_f]
      signs_f <- cbind(motr_signs_f, meta_signs = meta_signs_f)
      signs_f$agree <- signs_f$sign * signs_f$meta_signs
      signs_f$log10p <- log10(signs_f$p_value)
      signs_f$log10p[signs_f$log10p < min_log10pval] <- min_log10pval
      
      signs <- rbind(signs_m, signs_f)
      
      prop_agree <- sapply(1:nwin, function(wi) mean(signs[signs$log10p < wins[wi,1] & signs$log10p > wins[wi,2], "agree"] == 1)) #can also use cut?
      prop_agree <- zoo::na.approx(prop_agree)      #linearly interpolate missing values
      prop_agree <-c(prop_agree, rep(prop_agree[length(prop_agree)], times = nwin-length(prop_agree)))
      
      #or use a laplace smoother?
      # dexp_smooth(x = signs$log10p[order(signs$log10p, decreasing = T)][1:1000], y = as.numeric(signs$agree == 1)[order(signs$log10p, decreasing = T)][1:1000], r = 5, 
      #             reweight_trunc_tail = F, reweight_n_pts = F, fix_endpoints = F, 
      #             interpolate_at = wins_center)
      
    })
})







#### try just a scatterplot? ####

tissues <- c("t55-gastrocnemius", "t56-vastus-lateralis", "t30-blood-rna")
tissues <- setNames(tissues,tissues)
sexes = c("male", "female")
timepoints = paste0(2^(0:3), "w")

if(!exists("testout_list")){
  testout_list <- lapply(setNames(sexes, sexes), function(sex_i){
    lapply(setNames(tissues, tissues), function(tissue){
      # tissue <- tissues[3]
      # sex = sexes[1]
      sex = sex_i
      
      timepoint = timepoints[4]
      motrpac_data <- deg_eqtl_list[[tissue]]
      motrpac_data <- motrpac_data[motrpac_data$comparison_group == timepoint & motrpac_data$sex == sex_i,]     
      
      if(tissue %in% c("t55-gastrocnemius", "t56-vastus-lateralis")){
        xdat <- human_muscle_res
        xval <- xdat$beta
      } else {
        xdat <- human_blood_res
        xval <- xdat$beta
      }
      yval <- motrpac_data$logFC[match(xdat$symbol, motrpac_data$human_gene_symbol)]
      xval <- xval[!is.na(yval)]
      yval <- asinh(yval[!is.na(yval)])
      
      testout <- lapply(setNames(c("spearman", "pearson", "kendall"),c("spearman", "pearson", "kendall")), 
                        function(testmeth) cor.test(xval, sinh(yval), method = testmeth))
      testout
    })
  })
}

#### final figure ####
incl_sex_comparison <- F
list_only_pearson <- T

# cairo_pdf("~/Documents/Documents - nikolai/motrpac_companion/figures/pass1b_fig8_rat-man_comparison.pdf", width = 1200 / 72, height = 675 / 72 * 2, family="Arial Unicode MS", pointsize = 18)
cairo_pdf(paste0("~/Documents/Documents - nikolai/motrpac_companion/figures/pass1b_fig8_rat-man_comparison", ifelse(incl_sex_comparison, "", "_no-direct-sex-effect"), ".pdf"), 
          width = ifelse(incl_sex_comparison, 1200 / 72 * 2, 1200 / 72 * 2 * 5 / 6),
          height = 675 / 72 * 2 * 3 / 4,
          family="Arial Unicode MS", pointsize = 18)
# layout(rbind(c(1,2,3), c(4,5,6), c(7,8,9), c(10,10,10)))
if(incl_sex_comparison){
  layout(rbind(c(1, 4, 7, 10, 13, 16),
               c(2, 5, 8, 11, 14, 16),
               c(3, 6, 9, 12, 15, 16))
  )
} else {
  layout(rbind(c(1, 4, 7, 10, 13),
               c(2, 5, 8, 11, 14),
               c(3, 6, 9, 12, 15))
  )
}

par(xpd = NA, mar = c(3,5,2,3))

for(sex_i in sexes){
  for(tissue in tissues){
    # tissue <- tissues[3]
    # sex = sexes[1]
    sex = sex_i
    
    timepoint = timepoints[4]
    motrpac_data <- deg_eqtl_list[[tissue]]
    motrpac_data <- motrpac_data[motrpac_data$comparison_group == timepoint & motrpac_data$sex == sex_i,]     
    
    if(tissue %in% c("t55-gastrocnemius", "t56-vastus-lateralis")){
      xdat <- human_muscle_res
      xval <- xdat$beta
    } else {
      xdat <- human_blood_res
      xval <- xdat$beta
    }
    yval <- motrpac_data$logFC[match(xdat$symbol, motrpac_data$human_gene_symbol)]
    xval <- xval[!is.na(yval)]
    yval <- asinh(yval[!is.na(yval)])
    
    # plot(xval, 
    #      yval,
    #      pch = 19, col = adjustcolor(1, 0.2), main = "male, 82, rat vastus vs. human muscle")
    # bivkde <- ks::kde(cbind(xval, yval))
    # bivkde$estimate <- bivkde$estimate - min(bivkde$estimate) + min(bivkde$estimate[bivkde$estimate > 0])
    # bivkde$estimate <- log(bivkde$estimate) - min(log(bivkde$estimate)) + 1E-6
    # plot(bivkde)
    
    
    ndens <- 15
    xr <- seq_incl_0(min(xval)-1E-6, max(xval)+1E-6, length.out = ndens)
    yr <- seq_incl_0(min(yval)-1E-6, max(yval)+1E-6, length.out = ndens)
    nincell <- sapply(seq_along(xr)[-1]-1, function(xi) sapply(seq_along(yr)[-1]-1, function(yi){
      sum(xval > xr[xi] & xval < xr[xi+1] & yval > yr[yi] & yval < yr[yi+1])
    }))
    nincell <- nincell[(ndens-1):1,]
    maxlog <- max(log(nincell))
    
    cols <- rev(viridis::mako(100+1))
    alpha_pow <- 1/2
    cols <- sapply(1:101/101, function(colalpha) adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissue], colalpha^(alpha_pow)))
    
    plot(NA, xlim = c(1,ndens), ylim = c(1,ndens), xaxt = "n", frame = F, yaxt = "n", xlab = "", ylab = "")
    for(xi in 1:(ndens-1)){
      for(yi in 1:(ndens-1)){
        rect(xi, ndens-yi+1, xi+1, ndens-yi, col = cols[ceiling(log(nincell[yi,xi]) / maxlog * 100 + 1E-6)], border = adjustcolor(1, 0.1))
      }
    }
    
    #overall frame
    rect(1,1,ndens,ndens)
    
    #axes
    #yax 
    segments(x0 = 1, x1 = 1-ndens/50, y0 = 1:ndens, y1 = 1:ndens)
    text(x = 1-ndens/100, y = 1:ndens, pos = 2, labels = round(sinh(yr), 2), xpd = NA, cex = 0.8)
    #xax
    segments(x0 = 1:ndens, x1 = 1:ndens, y0 = 1, y1 =  1-ndens/50)
    text(x =  1:ndens+0.3, y = 1-ndens/30, pos = 2, labels = round(xr, 2), xpd = NA, srt = 45, cex = 0.8)
    
    #label axes
    text(x = ndens/2, y = -ndens/12, labels = latex2exp::TeX("Human Meta-analysis $\\Beta$"), cex = 1.5)
    text(x = -ndens/10, y = ndens/2, labels = latex2exp::TeX("MoTrPAC $log_{2}FC$"), cex = 1.5, srt = 90)
    
    #mark origin on graph
    segments(x0 = 1, x1 = ndens, y0 = which(yr == 0), y1 = which(yr == 0))
    segments(x0 = which(xr == 0), x1 = which(xr == 0), y0 = 1, y1 = ndens)
    
    #title
    # text(x = ndens/2, y = ndens + ndens/20, labels = MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[tissue], cex = 1.5)
    text_indivcolor(xloc = ndens/2 - ndens/4, y = ndens + ndens/20, 
                    labels_mat = t(matrix(c(MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[tissue], ", ", timepoint, ", ", c(male = "\u2642", female = "\u2640")[sex]))),
                    colors_mat = t(matrix(c(MotrpacRatTraining6moData::TISSUE_COLORS[tissue], 1, MotrpacRatTraining6moData::GROUP_COLORS[timepoint], 1, MotrpacRatTraining6moData::SEX_COLORS[sex]))),
                    cex = 1.75, pos = 4)
    
    #legend for colors
    xl = ndens + ndens/30; xr = ndens + ndens/10; yb = ndens / 2; yt = ndens;
    rect(xleft = xl, ybottom = yb, xright = xr, ytop = yt)
    nblocks <- 30
    nlabs <- 10
    labvals <- c(1,floor(exp(seq(1,maxlog, length.out = nlabs))))
    lablocs <- log(labvals) / (maxlog) 
    lablocs <- lablocs * (yt-yb) + yb
    rect(xleft = xl, xright = xr, 
         yb = seq(yb, yt, length.out = nblocks+1)[-1], 
         yt = seq(yb, yt, length.out = nblocks+1)[-(nblocks+1)],
         border = NA, col = cols[round(seq(1, length(cols), length.out = nblocks))])
    text(x = xr - ndens/100, y = lablocs, labels = labvals, pos = 4, cex = 0.5)
    text(x = mean(c(xl,xr)), y = yt, pos = 3, labels = latex2exp::TeX("$n_{genes}$"))
    
    #get correlation coefficients in there too
    testout <- testout_list[[sex]][[tissue]]
    if(!list_only_pearson){
      text(labels = latex2exp::TeX(
        paste0("Spearman's $\\rho$ = ", round(testout$spearman$estimate,2), ", p = ", "$10^{", round(log10(testout$spearman$p.value), 2), "}$\n")
      ), x = 1 - ndens/100, y = ndens - ndens/20, pos = 4, cex = 0.7)  
      text(labels = latex2exp::TeX(
        paste0("Kendall's $\\tau$ = ", round(testout$kendall$estimate,2), ", p = ", "$10^{", round(log10(testout$kendall$p.value), 2), "}$\n")
      ), x = 1 - ndens/100, y = ndens - ndens/20 - 2, pos = 4, cex = 0.7)
    }
    text(labels = latex2exp::TeX(
      paste0("Pearson's $r$ = ", round(testout$pearson$estimate,2), ", p = ", "$10^{", round(log10(testout$pearson$p.value), 2), "}$\n")
    ), x = 1 - ndens/100, y = ndens - ndens/20 - 1, pos = 4, cex = 0.7)
    
    #figure label
    if(tissue == tissues[1] & sex == sexes[1]){
      fig_label("a)", cex = 2, shrinkX = 0.8)
    }
    
  }
}

#now plot the proportion figure
node_type_cols <- c(both = "purple", MotrpacRatTraining6moData::SEX_COLORS["male"], MotrpacRatTraining6moData::SEX_COLORS["female"], neither = "grey50")
# par(mar = c(3,5.5,3,3.5))

for(tissue in tissues){
  
  #lay down basic plot
  ylims = range(sapply(names(focal_node_sets), function(node_type) sapply(paste0(2^(0:3), "w"), function(tpt)  prop_in_motr[[tpt]][tissue, node_type])))
  ylims = c(floor(ylims[1]*100)/100, ceiling(ylims[2]*100)/100)
  # ylims = list(c(0.1,0.3), c(0.18, 0.4), c(0,0.02))[[match(tissue, tissues)]]
  plot(NA, xlim = c(1,4), ylim = ylims, xaxt = "n", frame = F, yaxt = "n", xlab = "", ylab = "")
  #title
  text(x = 2.5, y = ylims[2], labels = MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[tissue], pos = 3, cex = 1.75, col = MotrpacRatTraining6moData::TISSUE_COLORS[tissue], xpd = NA)
  
  #axes
  #yax
  yax <-seq(ylims[1], ylims[2], by = diff(ylims)/10)
  segments(x0 = 1, x1 = 1-1/20, y0 = yax, y1 = yax)
  text(x = 1-1/20, y = yax, pos = 2, labels = yax, xpd = NA, cex = 0.8)
  #xax
  segments(x0 = 1:4, x1 = 1:4, y0 = ylims[1], y1 = ylims[1]-diff(ylims)/50)
  shadowtext(x =  1:4 + 0.1, y = rep(ylims[1]-diff(ylims)/20, 4), pos = 2, labels = timepoints, xpd = NA, srt = 45, cex = 1.2, col = MotrpacRatTraining6moData::GROUP_COLORS[timepoints], theta = 1)
  
  #label axes
  text(x = 2.5, y = ylims[1]-diff(ylims)/6, labels = latex2exp::TeX("Timepoint"), cex = 1.5)
  text(x = 0.3, y = mean(ylims), labels = latex2exp::TeX("Prop. Genes Replicated in Meta-analysis"), cex = 1, srt = 90)
  text(x = 0.45, y = mean(ylims), labels = latex2exp::TeX("(IHW $\\alpha$ = 0.05)"), cex = 1, srt = 90)
  
  #add legend
  legend(x = 4, y = ylims[2], lwd = 4, col = node_type_cols, legend = names(node_type_cols), bty = "n", seg.len = 1)
  
  for(node_type in names(focal_node_sets)){
    vals <- sapply(paste0(2^(0:3), "w"), function(tpt)  prop_in_motr[[tpt]][tissue, node_type])
    lines(1:4, vals, col = node_type_cols[node_type], lwd = 4)
  }  
  
  #overall frane
  rect(xleft = 1,ybottom = ylims[1],xright = 4,ytop = ylims[2], xpd = NA)
  
  #figure label
  if(tissue == tissues[1]){
    fig_label("b)", cex = 2, shrinkX = 0.8, shrinkY = 0.995)
  }
  
}

#lay down basic plot
for(sex_i in sexes){
  for(tissue in tissues){
    
    ylims = c(-0.02,1.02)
    xlims = range(-wins_center) + c(-1E-2,1E-2)
    # ylims = list(c(0.1,0.3), c(0.18, 0.4), c(0,0.02))[[match(tissue, tissues)]]
    plot(NA, xlim = xlims, ylim = ylims, xaxt = "n", frame = F, yaxt = "n", xlab = "", ylab = "")
    #title
    text_indivcolor(xloc = mean(xlims) - diff(xlims)/5, y = ylims[2] + diff(ylims)/20,pos = 4,
                    labels_mat = t(matrix(c(MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[tissue], ", ", c(male = "\u2642", female = "\u2640")[sex_i]))),
                    colors_mat = t(matrix(c(MotrpacRatTraining6moData::TISSUE_COLORS[tissue], 1, MotrpacRatTraining6moData::SEX_COLORS[sex_i]))),
                    cex = 1.75)
    
    #axes
    #yax
    yax <-seq(0, 1, by = 0.1)
    segments(x0 = xlims[1], x1 = xlims[1]-2/20, y0 = yax, y1 = yax)
    text(x = xlims[1]-1/20, y = yax, pos = 2, labels = yax, xpd = NA, cex = 0.8)
    #xax
    xtickloc <- round(seq(xlims[1], xlims[2], by = 1))
    xtickval <- c(xtickloc[-length(xtickloc)], paste0(">", xtickloc[length(xtickloc)]))
    segments(x0 = xtickloc, x1 = xtickloc, y0 = ylims[1], y1 = ylims[1]-diff(ylims)/50)
    text(x =  xtickloc + 0.1, y = rep(ylims[1]-diff(ylims)/20, 4), pos = 2, labels = xtickval, xpd = NA, srt = 45, cex = 1.2, col = 1)
    
    #label axes
    text(x = mean(xlims), y = ylims[1]-diff(ylims)/6.5, labels = latex2exp::TeX("-log$_{10}$(pval)"), cex = 1.5)
    text(x = xlims[1] - diff(xlims)/6, y = mean(ylims), labels = latex2exp::TeX("Prop. Matching Signs"), cex = 1.5, srt = 90)
    
    #add horizontal line at 0.5
    
    #add legend
    legend(x = xlims[2], y = ylims[2], lwd = 4, col = MotrpacRatTraining6moData::GROUP_COLORS[timepoints], legend = timepoints, bty = "n", seg.len = 1)
    
    #add actual lines
    for(timepoint in timepoints){
      lines(-wins_center, pval_window_output[[tissue]][[sex_i]][[timepoint]], col = MotrpacRatTraining6moData::GROUP_COLORS[timepoint], lwd = 4)
    }  
    
    #overall frane
    rect(xleft = xlims[1],ybottom = ylims[1],xright = xlims[2],ytop = ylims[2], xpd = NA)
    
    #figure label
    if(tissue == tissues[1] & sex_i == sexes[1]){
      fig_label("c)", cex = 2, shrinkX = 0.8)
    }
  }
}


if(incl_sex_comparison){
  par(mar = c(4.25,4.5,2.5,2))
  # output from fig1a_ratman-sex-comparison.R
  samps <- read.table("~/data/smontgom/ratman_sex_comparison_samps.txt")
  colnames(samps) <- gsub("\\.", "-", colnames(samps))
  tord <- order(apply(samps, 2, mean))
  qi_95 <- apply(samps, 2, quantile, probs = c(0.05, 0.95))[,tord]
  qi_100 <- apply(samps, 2, quantile, probs = c(0.0, 1.0))[,tord]
  posterior_means <- apply(samps, 2, mean)[tord]
  inside_cols <- rep("black", length(tord))
  inside_cols[overlaps_with_zero(qi_95-0.5)] <- "white"
  tmp <- vioplot::vioplot(x = samps[,tord], T, 
                          col = MotrpacRatTraining6moData::TISSUE_COLORS[colnames(samps)][tord], outline=FALSE, yaxt = "n",
                          names = colnames(samps)[tord], range = 0, ylab = "", lineCol = MotrpacRatTraining6moData::TISSUE_COLORS[colnames(samps)][tord],
                          xlab = "", cex.lab = 2, plotCentre = "point", 
                          colMed = MotrpacRatTraining6moData::TISSUE_COLORS[colnames(samps)][tord],
                          horizontal = T)
  
  #axes
  xtickvals <- round(seq(min(samps), max(samps), length.out = 8), 2)
  segments(x0 = xtickvals, x1 = xtickvals, y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4])/100, xpd = NA)
  text(x = xtickvals, y = par("usr")[3] - diff(par("usr")[3:4])/50, labels = xtickvals)
  text(x = mean(par("usr")[1:2]), y = par("usr")[3]-diff(par("usr")[3:4])/20, labels = "Posterior Proportion Same Direction", cex = 1.5, xpd = NA)
  tick <- seq_along(colnames(samps))
  axis(2, at = tick, labels = F)
  for(i in 1:length(colnames(samps))){
    segments(y0 = i, y1 = i, x0 = qi_95[1,i], x1 = qi_95[2,i], lwd = 4, col = inside_cols[i])
    points(y = i, x = posterior_means[i], col = inside_cols[i], pch = 19, cex = 1.5)
    
    #guiding lines
    segments(y0 = i, y1 = i, x0 = 0, x1 = qi_100[1,i], lwd = 1, col = "black", xpd = T, lty = 3)
    segments(y0 = i, y1 = i, x0 = qi_100[2,i], x1 = max(qi_100) * 0.5 + 1 * 0.5, lwd = 1, col = "black", xpd = T, lty = 3)
    
  }
  text(y = tick, x = rep(par("usr")[1] - 0.01, length(tick)), colnames(samps)[tord], srt = 0, xpd = T, pos = 2,
       col = MotrpacRatTraining6moData::TISSUE_COLORS[colnames(samps)[tord]])
  
  abline(v=0.5,lwd=3,lty=2, col = adjustcolor(1,0.5), xpd = F)
  
  #figure label
  if(T){
    fig_label("d)", cex = 2, shrinkX = 0.8, shrinkY = 0.985)
  }
}

dev.off()

#### get alternative second set of figures ####
motrpac_expressed_genes <- lapply(tissues, function(ti) unique(deg_eqtl_list[[ti]]$feature_ID))
human_expressed_genes <- list("t55-gastrocnemius" = human_muscle_res$symbol, "t56-vastus-lateralis" = human_muscle_res$symbol, "t30-blood-rna" = human_blood_res$symbol)
human_diff_expressed_genes <- list("t55-gastrocnemius" = muscle_genes, "t56-vastus-lateralis" = muscle_genes, "t30-blood-rna" = blood_genes)

possible_intersecting_genes <- lapply(tissues, function(ti){ 
  motrgenes <- map$human_gene_symbol[match(motrpac_expressed_genes[[ti]], map$feature_ID)]
  motrgenes <- motrgenes[!is.na(motrgenes)]
  intersect(motrgenes, human_expressed_genes[[ti]])
})

included_in_graphical_clustering <- unlist(node_sets[!grepl("0w", names(node_sets))])
included_in_graphical_clustering <- included_in_graphical_clustering[grepl("TRNSCRPT", included_in_graphical_clustering)]
included_in_graphical_clustering <- do.call(rbind, strsplit(included_in_graphical_clustering, ";"))
included_in_graphical_clustering <- lapply(tissues, function(ti){
  x <- unique(included_in_graphical_clustering[included_in_graphical_clustering[,2] == MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[ti],3])
  x <- map$human_gene_symbol[match(x, map$feature_ID)]
  x <- x[!is.na(x)]
  x <- intersect(x, possible_intersecting_genes[[ti]])
  x
})

not_included_in_graphical_clustering <- lapply(tissues, function(ti){
  setdiff(possible_intersecting_genes[[ti]], included_in_graphical_clustering[[ti]])
})

focal_node_sets <- focal_node_sets[c("both", "male", "female")]
timepoint == "8w"
motrpac_DE_nodes <- lapply(setNames(names(focal_node_sets), names(focal_node_sets)), function(node_type){
  
  node_metadata <- lapply(nodes_to_look_at_list[[timepoint]][[node_type]], function(node_to_look_at){
    cbind(do.call(rbind, strsplit(node_sets[node_to_look_at][[node_to_look_at]][grep(
      node_sets[node_to_look_at][[node_to_look_at]], pattern = "TRNSCRPT")], ";")), node_to_look_at)
  })
  
  node_metadata <- as.data.table(do.call(rbind, node_metadata))
  colnames(node_metadata) <- c("ome","tissue","ensembl_gene", "node")
  
  #snag signs
  node_metadata_sign <- do.call(rbind, strsplit(node_metadata$node, "_"))
  node_metadata_sign[grep("F-1", node_metadata_sign[,2]),2] <- -1
  node_metadata_sign[grep("F1", node_metadata_sign[,2]),2] <- 1
  node_metadata_sign[grep("F0", node_metadata_sign[,2]),2] <- 0
  node_metadata_sign[grep("M-1", node_metadata_sign[,3]),3] <- -1
  node_metadata_sign[grep("M1", node_metadata_sign[,3]),3] <- 1
  node_metadata_sign[grep("M0", node_metadata_sign[,3]),3] <- 0
  
  if(node_type == "male"){
    out <- setNames(object = as.integer(node_metadata_sign[,3]), node_metadata$ensembl_gene)
  }
  if(node_type == "female"){
    out <- setNames(object = as.integer(node_metadata_sign[,2]), node_metadata$ensembl_gene)
  }
  if(node_type == "both"){
    out <- setNames(object = as.integer(node_metadata_sign[,3]), node_metadata$ensembl_gene)
  }
  if(node_type == "neither"){
    out <- setNames(object = as.integer(node_metadata_sign[,3]), node_metadata$ensembl_gene)
  }
  
  names(out) <- map$human_gene_symbol[match(names(out), map$feature_ID)]
  
  out <- lapply(setNames(tissues, tissues), function(tiss){
    subout <- out[node_metadata$tissue == MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[tiss]]
    subout <- subout[!is.na(names(subout))]
    subout <- subout[intersect(possible_intersecting_genes[[tiss]], names(subout))]
  })
  
  out
})


motrpac_DE_nodes$neither <- lapply(not_included_in_graphical_clustering, function(x) setNames(rep(0, length(x)), x))

proportions_of_ratDE_in_humans <- do.call(rbind, lapply(motrpac_DE_nodes, function(group){
  sapply(setNames(names(group), names(group)), function(ti){
    length(intersect(names(group[[ti]]), human_diff_expressed_genes[[ti]])) / length(group[[ti]])
  })
}))

count_of_rat_DE_in_humans <- do.call(rbind, lapply(motrpac_DE_nodes, function(group){
  sapply(setNames(names(group), names(group)), function(ti){
    length(intersect(names(group[[ti]]), human_diff_expressed_genes[[ti]]))
  })
}))

count_of_rat_notDE_in_humans <- do.call(rbind, lapply(motrpac_DE_nodes, function(group){
  sapply(setNames(names(group), names(group)), function(ti){
    length(intersect(names(group[[ti]]), setdiff(possible_intersecting_genes[[ti]], human_diff_expressed_genes[[ti]])))
  })
}))


FET_out <- sapply(tissues, function(ti){
  sapply(setNames(names(focal_node_sets),names(focal_node_sets)), function(nodeset){
    ct <- rbind(c(DE = count_of_rat_DE_in_humans[nodeset, ti], nDE = count_of_rat_notDE_in_humans[nodeset, ti]),
                c(DE = count_of_rat_DE_in_humans["neither", ti], nDE = count_of_rat_notDE_in_humans["neither", ti]))
    rownames(ct) <- c(nodeset, "neither")
    fisher.test(ct)$p.value
  })
})

#compare to open targets
if(!exists("dismat")){
  direct_associations <- fread("~/data/smontgom/opentargets/associationByOverallDirect.csv")
  genes <- intersect(unique(direct_associations$targetId), map$human_ensembl_gene[match(unique(unlist(possible_intersecting_genes)), map$human_gene_symbol)])
  genes <- genes[!is.na(genes)]
  diseases <- unique(direct_associations$diseaseId)
  disvec <- setNames(rep(0, length(diseases)), diseases)
  associations_to_use <- split(direct_associations, direct_associations$targetId)
  dismat <- do.call(rbind, mclapply(setNames(genes,genes), function(gene){
    locdisvec <- disvec
    locdisvec[associations_to_use[[gene]]$diseaseId] <- associations_to_use[[gene]]$score
    locdisvec
  }, mc.cores = 8))
  nonzero_entries <- which(apply(dismat, 2, sum) > 1E-6)
  dismat <- dismat[,nonzero_entries]
}


genes_in_intersect <- lapply(motrpac_DE_nodes, function(group){
  lapply(tissues, function(ti){
    map$human_ensembl_gene[match(intersect(names(group[[ti]]), human_diff_expressed_genes[[ti]]), map$human_gene_symbol)]
  })
})

evidence_score_threshold <- 0.8
n_genes_above_evidence_score <- sapply(setNames(names(genes_in_intersect), names(genes_in_intersect)), function(gi){
  sapply(tissues, function(ti){
    geneset <- genes_in_intersect[[gi]][[ti]]
    geneset <- intersect(geneset, rownames(dismat))
    sum(apply(dismat[geneset,] > evidence_score_threshold, 1, sum) > 0.1)
  })
})



#####

par(mfrow = c(2,1), xpd = NA)
ylims <- c(0,0.8)
plot(NA, xlim = c(0.5,3.5), ylim = ylims, xaxt = "n", frame = F, yaxt = "n", xlab = "", ylab = "")
rect(0.5, -0.025, 3.5, ylims[2])

for(i in 1:nrow(proportions_of_ratDE_in_humans)){
  
  for(j in 1:ncol(proportions_of_ratDE_in_humans)){
    rect(xleft = j + 0.5 - i / 5 - 1 / 10, 
         xright = j + 0.5 - i / 5 + 1 / 10, 
         ybottom = 0, ytop = proportions_of_ratDE_in_humans[i,j],
         col = node_type_cols[rownames(proportions_of_ratDE_in_humans)[i]])
    if(i == 1){
      text(MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[colnames(proportions_of_ratDE_in_humans)[j]], x = j, y = -0.025, pos = 1)
    }
    if(i != 4){
      segments(x0 = j + 0.5 - i / 5 , x1 = j + 0.5 - 4 / 5, 
               y0 = max(proportions_of_ratDE_in_humans[,j]) + (4-i) / 20, 
               y1 = max(proportions_of_ratDE_in_humans[,j]) + (4-i) / 20)
      segments(x0 = j + 0.5 - i / 5 , x1 = j + 0.5 - i / 5, 
               y0 = max(proportions_of_ratDE_in_humans[,j]) + (4-i) / 20, 
               y1 = max(proportions_of_ratDE_in_humans[,j]) + (4-i) / 20 - 0.01)
      segments(x0 = j + 0.5 - 4 / 5 , x1 = j + 0.5 - 4 / 5, 
               y0 = max(proportions_of_ratDE_in_humans[,j]) + (4-i) / 20, 
               y1 = max(proportions_of_ratDE_in_humans[,j]) + (4-i) / 20 - 0.01)
      text(paste0(ifelse(i==1, "FET p-val = ", ""), round(FET_out[rownames(proportions_of_ratDE_in_humans)[i], colnames(proportions_of_ratDE_in_humans)[j]], 7)),
           x = mean(c(j + 0.5 - i / 5 , j + 0.5 - 4 / 5)), y = max(proportions_of_ratDE_in_humans[,j]) + (4-i) / 20 - 0.015, pos = 3, cex = 0.75)
      
    }
  }
}

#yax
yax <-seq(0, ylims[2], by = 0.2)
segments(x0 = 0.5, x1 = 0.45, y0 = yax, y1 = yax)
text(x = 0.45, y = yax, pos = 2, labels = yax, xpd = NA, cex = 0.8)

#label axes
text(x = 2, y = -diff(ylims) / 10, labels = "Tissue", cex = 1.5, pos = 1)
text(x = 0.2, y = mean(ylims), labels = latex2exp::TeX("Prop. MoTrPAC Genes\nReplicated in Meta-analysis"), cex = 1, srt = 90)

#add legend
legend(x = 3.05, y = ylims[2], col = node_type_cols, legend = names(node_type_cols), pt.cex = 2, pch = 15, bty = "n", seg.len = 1)


#plot counts of intersect figure

ylims <- c(0,ceiling(log10(max(count_of_rat_DE_in_humans))))
plot(NA, xlim = c(0.5,3.5), ylim = ylims, xaxt = "n", frame = F, yaxt = "n", xlab = "", ylab = "")
rect(0.5, -0.025, 3.5, ylims[2])

for(i in 1:nrow(count_of_rat_DE_in_humans)){
  
  for(j in 1:ncol(count_of_rat_DE_in_humans)){
    text(labels = paste0(n_genes_above_evidence_score[colnames(count_of_rat_DE_in_humans)[j], rownames(count_of_rat_DE_in_humans)[i]], " genes"),
         x = j + 0.5 - i / 5, y = log10(count_of_rat_DE_in_humans[i,j]) - 0.05, pos = 3, cex = 0.6)
    rect(xleft = j + 0.5 - i / 5 - 1 / 10, 
         xright = j + 0.5 - i / 5 + 1 / 10, 
         ybottom = 0, ytop = log10(count_of_rat_DE_in_humans[i,j]),
         col = node_type_cols[rownames(count_of_rat_DE_in_humans)[i]])
    if(i == 1){
      text(MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[colnames(count_of_rat_DE_in_humans)[j]], x = j, y = -0.025, pos = 1)
    }
  }
}

top_label <- paste0("# genes corresponds to Open Targets hits with > ", evidence_score_threshold, " evidence scores in at least one trait")
text(labels = top_label,
     x = par("usr")[2] - diff(par("usr")[1:2]) / 30, y = ylims[2] + diff(ylims) / 40, pos = 2, cex = 0.75)

#yax
yax <-seq(0, ylims[2], by = 1)
segments(x0 = 0.5, x1 = 0.45, y0 = yax, y1 = yax)
segments(x0 = 0.5, x1 = 0.475, y0 = log10(c((t(t(1:10)) %*% 10^(ylims[1]:(ylims[2]-1)))[-c(1,10),])), y1 = log10(c((t(t(1:10)) %*% 10^(ylims[1]:(ylims[2]-1)))[-c(1,10),])))
text(x = 0.45, y = yax, pos = 2, labels = latex2exp::TeX(paste0("$10^{", yax, "}$")), xpd = NA, cex = 0.8)

#label axes
text(x = 2, y = -diff(ylims) / 10, labels = "Tissue", cex = 1.5, pos = 1)
text(x = 0.2, y = mean(ylims), labels = latex2exp::TeX("Size of Rat-Human DE Intersect"), cex = 1, srt = 90)

#add legend
legend(x = 3.05, y = ylims[2], col = node_type_cols, legend = names(node_type_cols), pt.cex = 2, pch = 15, bty = "n", seg.len = 1)



#####




#### final figure, redux ####
incl_sex_comparison <- F
list_only_pearson <- T

# cairo_pdf("~/Documents/Documents - nikolai/motrpac_companion/figures/pass1b_fig8_rat-man_comparison.pdf", width = 1200 / 72, height = 675 / 72 * 2, family="Arial Unicode MS", pointsize = 18)
cairo_pdf(paste0("~/Documents/Documents - nikolai/motrpac_companion/figures/pass1b_fig8_rat-man_comparison_redux.pdf"), 
          width = 1200 / 72 * 2 * 5 / 6 * 4 / 5 * 3 / 4,
          height = 675 / 72 * 2 * 3 / 4,
          family="Arial Unicode MS", pointsize = 18)
layout(rbind(c(1, 4, 7),
             c(2, 5, 8),
             c(3, 6, 9)))

par(xpd = NA, mar = c(3,5,2,3))

for(tissue in tissues){
  sex = sex_i
  
  timepoint = timepoints[4]
  motrpac_data <- deg_eqtl_list[[tissue]]
  motrpac_data <- motrpac_data[motrpac_data$comparison_group == timepoint,]     
  
  if(tissue %in% c("t55-gastrocnemius", "t56-vastus-lateralis")){
    xdat <- human_muscle_res
    xval <- xdat$beta
  } else {
    xdat <- human_blood_res
    xval <- xdat$beta
  }
  yval <- motrpac_data$logFC
  xval <- xdat$beta[match(motrpac_data$human_gene_symbol, xdat$symbol)]
  yval <- yval[!is.na(xval)]
  xval <- xval[!is.na(xval)]
  yval <- asinh(yval[!is.na(yval)])
  
  # plot(xval, 
  #      yval,
  #      pch = 19, col = adjustcolor(1, 0.2), main = "male, 82, rat vastus vs. human muscle")
  # bivkde <- ks::kde(cbind(xval, yval))
  # bivkde$estimate <- bivkde$estimate - min(bivkde$estimate) + min(bivkde$estimate[bivkde$estimate > 0])
  # bivkde$estimate <- log(bivkde$estimate) - min(log(bivkde$estimate)) + 1E-6
  # plot(bivkde)
  
  
  ndens <- 15
  xr <- seq_incl_0(min(xval)-1E-6, max(xval)+1E-6, length.out = ndens)
  yr <- seq_incl_0(min(yval)-1E-6, max(yval)+1E-6, length.out = ndens)
  nincell <- sapply(seq_along(xr)[-1]-1, function(xi) sapply(seq_along(yr)[-1]-1, function(yi){
    sum(xval > xr[xi] & xval < xr[xi+1] & yval > yr[yi] & yval < yr[yi+1])
  }))
  nincell <- nincell[(ndens-1):1,]
  maxlog <- max(log(nincell))
  
  cols <- rev(viridis::mako(100+1))
  alpha_pow <- 1/2
  cols <- sapply(1:101/101, function(colalpha) adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissue], colalpha^(alpha_pow)))
  
  plot(NA, xlim = c(1,ndens), ylim = c(1,ndens), xaxt = "n", frame = F, yaxt = "n", xlab = "", ylab = "")
  for(xi in 1:(ndens-1)){
    for(yi in 1:(ndens-1)){
      rect(xi, ndens-yi+1, xi+1, ndens-yi, col = cols[ceiling(log(nincell[yi,xi]) / maxlog * 100 + 1E-6)], border = adjustcolor(1, 0.1))
    }
  }
  
  #overall frame
  rect(1,1,ndens,ndens)
  
  #axes
  #yax 
  segments(x0 = 1, x1 = 1-ndens/50, y0 = 1:ndens, y1 = 1:ndens)
  text(x = 1-ndens/100, y = 1:ndens, pos = 2, labels = round(sinh(yr), 2), xpd = NA, cex = 0.8)
  #xax
  segments(x0 = 1:ndens, x1 = 1:ndens, y0 = 1, y1 =  1-ndens/50)
  text(x =  1:ndens+0.3, y = 1-ndens/30, pos = 2, labels = round(xr, 2), xpd = NA, srt = 45, cex = 0.8)
  
  #label axes
  text(x = ndens/2, y = -ndens/12, labels = latex2exp::TeX("Human Meta-analysis $\\Beta$"), cex = 1.5)
  text(x = -ndens/10, y = ndens/2, labels = latex2exp::TeX("MoTrPAC $log_{2}FC$"), cex = 1.5, srt = 90)
  
  #mark origin on graph
  segments(x0 = 1, x1 = ndens, y0 = which(yr == 0), y1 = which(yr == 0))
  segments(x0 = which(xr == 0), x1 = which(xr == 0), y0 = 1, y1 = ndens)
  
  #title
  # text(x = ndens/2, y = ndens + ndens/20, labels = MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[tissue], cex = 1.5)
  labels_mat <- t(matrix(c(MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[tissue], ", ", timepoint, ", ", "\u2642", " + ", "\u2640")))
  text_indivcolor(xloc = ndens/2 - strwidth(paste0(labels_mat, collapse = ""), "user", cex = 1.75) / 2, y = ndens + ndens/20, 
                  labels_mat = labels_mat,
                  colors_mat = t(matrix(c(MotrpacRatTraining6moData::TISSUE_COLORS[tissue], 1, MotrpacRatTraining6moData::GROUP_COLORS[timepoint], 1, MotrpacRatTraining6moData::SEX_COLORS["male"], 1, MotrpacRatTraining6moData::SEX_COLORS["female"]))),
                  cex = 1.75, pos = 4)
  
  #legend for colors
  xl = ndens + ndens/30; xr = ndens + ndens/10; yb = ndens / 2; yt = ndens;
  rect(xleft = xl, ybottom = yb, xright = xr, ytop = yt)
  nblocks <- 30
  nlabs <- 10
  labvals <- c(1,floor(exp(seq(1,maxlog, length.out = nlabs))))
  lablocs <- log(labvals) / (maxlog) 
  lablocs <- lablocs * (yt-yb) + yb
  rect(xleft = xl, xright = xr, 
       yb = seq(yb, yt, length.out = nblocks+1)[-1], 
       yt = seq(yb, yt, length.out = nblocks+1)[-(nblocks+1)],
       border = NA, col = cols[round(seq(1, length(cols), length.out = nblocks))])
  text(x = xr - ndens/100, y = lablocs, labels = labvals, pos = 4, cex = 0.5)
  text(x = mean(c(xl,xr)), y = yt, pos = 3, labels = latex2exp::TeX("$n_{genes}$"))
  
  #get correlation coefficients in there too
  testout <- testout_list[[sex]][[tissue]]
  if(!list_only_pearson){
    text(labels = latex2exp::TeX(
      paste0("Spearman's $\\rho$ = ", round(testout$spearman$estimate,2), ", p = ", "$10^{", round(log10(testout$spearman$p.value), 2), "}$\n")
    ), x = 1 - ndens/100, y = ndens - ndens/20, pos = 4, cex = 0.7)  
    text(labels = latex2exp::TeX(
      paste0("Kendall's $\\tau$ = ", round(testout$kendall$estimate,2), ", p = ", "$10^{", round(log10(testout$kendall$p.value), 2), "}$\n")
    ), x = 1 - ndens/100, y = ndens - ndens/20 - 2, pos = 4, cex = 0.7)
    text(labels = latex2exp::TeX(
      paste0("Pearson's $r$ = ", round(testout$pearson$estimate,2), ", p = ", "$10^{", round(log10(testout$pearson$p.value), 2), "}$\n")
    ), x = 1 - ndens/100, y = ndens - ndens/20 - 1, pos = 4, cex = 0.7)
  } else {
    text(labels = latex2exp::TeX(
      paste0("$r_{i,j}$ = ", round(testout$pearson$estimate,2))
    ), x = 1 - ndens/100, y = ndens - ndens/20 - 0.25, pos = 4, cex = 1.5)
    text(labels = latex2exp::TeX(
      paste0("p = ", "$10^{", round(log10(testout$pearson$p.value), 2), "}$")
    ), x = 1 - ndens/100, y = ndens - ndens/20 - 1.35, pos = 4, cex = 1.5)
  }
  
  #figure label
  if(tissue == tissues[1]){
    fig_label("a)", cex = 2, shrinkX = 0.8)
  }
  
}


#now plot the proportion figure
node_type_cols <- c(both = "purple", MotrpacRatTraining6moData::SEX_COLORS["male"], MotrpacRatTraining6moData::SEX_COLORS["female"], neither = "grey50")
# par(mar = c(3,5.5,3,3.5))

ylims <- c(0,0.8)
par(xpd = NA, mar = c(3,5,2,2))
plot(NA, xlim = c(0.5,3.5), ylim = ylims, xaxt = "n", frame = F, yaxt = "n", xlab = "", ylab = "")
rect(0.5, -0.025, 3.5, ylims[2])

for(i in 1:nrow(proportions_of_ratDE_in_humans)){
  
  for(j in 1:ncol(proportions_of_ratDE_in_humans)){
    rect(xleft = j + 0.5 - i / 5 - 1 / 10, 
         xright = j + 0.5 - i / 5 + 1 / 10, 
         ybottom = 0, ytop = proportions_of_ratDE_in_humans[i,j],
         col = node_type_cols[rownames(proportions_of_ratDE_in_humans)[i]])
    if(i == 1){
      text(MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[colnames(proportions_of_ratDE_in_humans)[j]], x = j, y = -0.025, pos = 1)
    }
    if(i != 4){
      segments(x0 = j + 0.5 - i / 5 , x1 = j + 0.5 - 4 / 5, 
               y0 = max(proportions_of_ratDE_in_humans[,j]) + (4-i) / 20, 
               y1 = max(proportions_of_ratDE_in_humans[,j]) + (4-i) / 20)
      segments(x0 = j + 0.5 - i / 5 , x1 = j + 0.5 - i / 5, 
               y0 = max(proportions_of_ratDE_in_humans[,j]) + (4-i) / 20, 
               y1 = max(proportions_of_ratDE_in_humans[,j]) + (4-i) / 20 - 0.01)
      segments(x0 = j + 0.5 - 4 / 5 , x1 = j + 0.5 - 4 / 5, 
               y0 = max(proportions_of_ratDE_in_humans[,j]) + (4-i) / 20, 
               y1 = max(proportions_of_ratDE_in_humans[,j]) + (4-i) / 20 - 0.01)
      text(paste0(ifelse(i==1, "FET p-val = ", ""), round(FET_out[rownames(proportions_of_ratDE_in_humans)[i], colnames(proportions_of_ratDE_in_humans)[j]], 3)),
           x = mean(c(j + 0.5 - i / 5 , j + 0.5 - 4 / 5)), y = max(proportions_of_ratDE_in_humans[,j]) + (4-i) / 20 - 0.015, pos = 3, cex = 0.75)
      
    }
  }
}

#yax
yax <-seq(0, ylims[2], by = 0.2)
segments(x0 = 0.5, x1 = 0.45, y0 = yax, y1 = yax)
text(x = 0.45, y = yax, pos = 2, labels = yax, xpd = NA, cex = 0.8)

#label axes
text(x = 2, y = -diff(ylims) / 10, labels = "Tissue", cex = 1.5, pos = 1)
text(x = 0, y = mean(ylims), labels = latex2exp::TeX("Prop. MoTrPAC Genes\nReplicated in Meta-analysis"), cex = 1.5, srt = 90)

#add legend
legend(x = 2.75, y = ylims[2], col = node_type_cols, legend = names(node_type_cols), pt.cex = 2, pch = 15, bty = "n", seg.len = 1)

#add figure label
fig_label("b)", cex = 2, shrinkX = 0.8)


#plot counts of intersect figure

ylims <- c(0,ceiling(log10(max(count_of_rat_DE_in_humans))))
plot(NA, xlim = c(0.5,3.5), ylim = ylims, xaxt = "n", frame = F, yaxt = "n", xlab = "", ylab = "")
rect(0.5, -0.025, 3.5, ylims[2])

for(i in 1:nrow(count_of_rat_DE_in_humans)){
  
  for(j in 1:ncol(count_of_rat_DE_in_humans)){
    text(labels = paste0(n_genes_above_evidence_score[colnames(count_of_rat_DE_in_humans)[j], rownames(count_of_rat_DE_in_humans)[i]], "\ngenes"),
         x = j + 0.5 - i / 5, y = log10(count_of_rat_DE_in_humans[i,j]) - 0.05, pos = 3, cex = 0.5)
    rect(xleft = j + 0.5 - i / 5 - 1 / 10, 
         xright = j + 0.5 - i / 5 + 1 / 10, 
         ybottom = 0, ytop = log10(count_of_rat_DE_in_humans[i,j]),
         col = node_type_cols[rownames(count_of_rat_DE_in_humans)[i]])
    if(i == 1){
      text(MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[colnames(count_of_rat_DE_in_humans)[j]], x = j, y = -0.025, pos = 1)
    }
  }
}

top_label <- paste0("\"# genes\" corresponds to\nOpen Targets hits with > ", evidence_score_threshold, " evidence scores in at least one trait")
text(labels = top_label,
     x = par("usr")[2] - diff(par("usr")[1:2]) / 30, y = ylims[2] + diff(ylims) / 30, pos = 2, cex = 0.75)

#yax
yax <-seq(0, ylims[2], by = 1)
segments(x0 = 0.5, x1 = 0.45, y0 = yax, y1 = yax)
segments(x0 = 0.5, x1 = 0.475, y0 = log10(c((t(t(1:10)) %*% 10^(ylims[1]:(ylims[2]-1)))[-c(1,10),])), y1 = log10(c((t(t(1:10)) %*% 10^(ylims[1]:(ylims[2]-1)))[-c(1,10),])))
text(x = 0.45, y = yax, pos = 2, labels = latex2exp::TeX(paste0("$10^{", yax, "}$")), xpd = NA, cex = 0.8)

#label axes
text(x = 2, y = -diff(ylims) / 10, labels = "Tissue", cex = 1.5, pos = 1)
text(x = 0, y = mean(ylims), labels = latex2exp::TeX("Size of Rat-Human DE Intersect"), cex = 1.5, srt = 90)

#add legend
legend(x = 2.75, y = ylims[2], col = node_type_cols, legend = names(node_type_cols), pt.cex = 2, pch = 15, bty = "n", seg.len = 1)

plot.new()

#trajectories through time figure
#lay down basic plot
for(tissue in tissues){
  
  ylims = c(-0.02,1.02)
  xlims = range(-wins_center) + c(-1E-2,1E-2)
  # ylims = list(c(0.1,0.3), c(0.18, 0.4), c(0,0.02))[[match(tissue, tissues)]]
  plot(NA, xlim = xlims, ylim = ylims, xaxt = "n", frame = F, yaxt = "n", xlab = "", ylab = "")
  #title
  labels_mat <- t(matrix(c(MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[tissue], ", 8w, ", "\u2642", " & ", "\u2640")))
  text_indivcolor(xloc = mean(xlims) - strwidth(paste0(labels_mat, collapse = ""), "user", cex = 1.75) / 2, y = ylims[2] + diff(ylims)/20,pos = 4,
                  labels_mat = labels_mat,
                  colors_mat = t(matrix(c(MotrpacRatTraining6moData::TISSUE_COLORS[tissue], 1, MotrpacRatTraining6moData::SEX_COLORS["male"], 1, MotrpacRatTraining6moData::SEX_COLORS["female"]))),
                  cex = 1.75)
  
  #axes
  #yax
  yax <-seq(0, 1, by = 0.1)
  segments(x0 = xlims[1], x1 = xlims[1]-2/20, y0 = yax, y1 = yax)
  text(x = xlims[1]-1/20, y = yax, pos = 2, labels = yax, xpd = NA, cex = 0.8)
  #xax
  xtickloc <- round(seq(xlims[1], xlims[2], by = 1))
  xtickval <- c(xtickloc[-length(xtickloc)], paste0(">", xtickloc[length(xtickloc)]))
  segments(x0 = xtickloc, x1 = xtickloc, y0 = ylims[1], y1 = ylims[1]-diff(ylims)/50)
  text(x =  xtickloc + 0.1, y = rep(ylims[1]-diff(ylims)/20, 4), pos = 2, labels = xtickval, xpd = NA, srt = 45, cex = 1.2, col = 1)
  
  #label axes
  text(x = mean(xlims), y = ylims[1]-diff(ylims)/6.5, labels = latex2exp::TeX("-log$_{10}$(pval)"), cex = 1.5)
  text(x = xlims[1] - diff(xlims)/6, y = mean(ylims), labels = latex2exp::TeX("Prop. Matching Signs"), cex = 1.5, srt = 90)
  
  #add horizontal line at 0.5
  
  #add legend
  legend(x = xlims[2], y = ylims[2], lwd = 4, col = MotrpacRatTraining6moData::SEX_COLORS[sexes], legend = c("m", "f"), bty = "n", seg.len = 1)
  
  #add actual lines
  for(sex_i in sexes){
    lines(-wins_center, pval_window_output[[tissue]][[sex_i]][["8w"]], col = MotrpacRatTraining6moData::SEX_COLORS[sex_i], lwd = 4)
  }  
  
  #overall frane
  rect(xleft = xlims[1],ybottom = ylims[1],xright = xlims[2],ytop = ylims[2], xpd = NA)
  
  #figure label
  if(tissue == tissues[1]){
    fig_label("c)", cex = 2, shrinkX = 0.8)
  }
}


dev.off()






