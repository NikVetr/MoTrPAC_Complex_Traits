#### run preprocessing script ####
source("/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/helper_scripts/figure_set_3_preprocessing.R")

#### figure plotting ####
#graphical parameters
incl_significance <- T
incl_cell_totals <- F
trait_category_legend_below = T
use_tissue_cols_for_cols <- F
opacity_power_scaler <- 0.25
opacity_white_threshold <- 0.85
use_counts <- T
prop_TWAS <- T
order_by_counts <- F
order_by_posterior_enrichment <- T
group_by_tissue_type <- T
use_range_for_maginal_labels <- F
subset_to_traits <- T
focal_traitcats <- c("Cardiometabolic", "Immune", "Endocrine system", "Anthropometric", "Allergy", "Aging")

#for grouping by tissues
tissue_cats <- list(circulation = c("BLOOD", "HEART", "SPLEEN"),
                    skeletal_muscle = c("SKM-GN", "SKM-VL"),
                    other = rev(c("ADRNL", "KIDNEY", "LUNG", "LIVER")),
                    adipose = c("WAT-SC"),
                    brain = c("CORTEX", "HYPOTH", "HIPPOC"),
                    GI = c("SMLINT", "COLON"))
tissue_cats <- rev(tissue_cats)
disp_amount <- 0.5

#some easy graphical params
cols = list(Tissue=MotrpacRatTraining6moData::TISSUE_COLORS[names(motrpac_gtex_map)], 
            Time=MotrpacRatTraining6moData::GROUP_COLORS[paste0(c(1,2,4,8), "w")],
            Sex=MotrpacRatTraining6moData::SEX_COLORS[c('male','female')])
mars <- list(c(3,
               0 + ifelse(subset_to_traits, 2, 0),
               6 + ifelse(group_by_tissue_type, disp_amount * 4.5, 0),
               6.5 + ifelse(subset_to_traits, 1, 0)),
             c(5,6,7,7.5)+0, 
             c(5,1.5,7,13)+0,
             c(4.25,2.5,2.5,2.5), 
             c(3,2.5,1.25,2.5), 
             c(10,6.5,1.5,2.5),
             c(3,2,3,2))

layout_mat <- rbind(
  c(1,1,1,1,4),
  c(2,2,3,3,5),
  c(6,6,6,6,6)
)
layout_heights = c(1,1,0.6)

layout_mat <- rbind(
  c(1,1,1,1,1,1,1,1,7,7),
  c(1,1,1,1,1,1,1,1,4,4),
  c(2,2,2,2,3,3,3,3,4,4),
  c(2,2,2,2,3,3,3,3,5,5),
  c(2,2,2,2,3,3,3,3,5,5),
  c(6,6,6,6,6,6,6,6,6,6),
  c(6,6,6,6,6,6,6,6,6,6)
)
layout_heights = c(1,1.5,1,1,1,1,0.75)


#preliminary processing for additional graphical params
if(subset_to_traits){
  trait_subset <- colnames(n_deg_sigtwas_intersect)[
    traitwise_partitions$Category[match(colnames(n_deg_sigtwas_intersect), traitwise_partitions$Tag)] %in% 
      focal_traitcats]
  categories_represented <- unique(traitwise_partitions$Category[match(trait_subset, traitwise_partitions$Tag)])
} else {
  trait_subset <- colnames(n_deg_sigtwas_intersect)
}
categories_represented <- unique(traitwise_partitions$Category[match(trait_subset, traitwise_partitions$Tag)])

if(use_counts){
  table_to_use <- n_deg_sigtwas_intersect  
} else {
  if(prop_TWAS){
    table_to_use <- round(prop_twas_are_degs * 1000)  
  } else{
    table_to_use <- round(prop_degs_are_twas * 1000)  
  }
}

if(order_by_counts){
  table_to_use <- table_to_use[,colnames(n_deg_sigtwas_intersect)]
  signif_matrix_to_use <- signif_matrix[,colnames(n_deg_sigtwas_intersect)]
} else if(order_by_posterior_enrichment) {
  trait_order <- traits[order(colbias, decreasing = T)]
  trait_order <- trait_order[trait_order %in% colnames(table_to_use)]
  table_to_use <- table_to_use[,trait_order]
  signif_matrix_to_use <- signif_matrix[,match(trait_order, colnames(signif_matrix))]
} else {
  table_to_use <- table_to_use[,colnames(prop_twas_are_degs)]
  signif_matrix_to_use <- signif_matrix[,colnames(prop_twas_are_degs)]
}

if(subset_to_traits){
  table_to_use <- table_to_use[,trait_subset]
  signif_matrix_to_use <- signif_matrix_to_use[,trait_subset]
}


if(group_by_tissue_type){
  tissue_disps <- unlist(lapply(1:length(tissue_cats), function(tci) rep(disp_amount * (tci), length(tissue_cats[[tci]]))))
  tissue_cats_bars_ylocs <- cbind(start = (c(0, cumsum(unlist(lapply(tissue_cats, function(tc) length(tc))) + disp_amount)) + 2 * disp_amount)[-(length(tissue_cats)+1)], 
                                  end = cumsum(unlist(lapply(tissue_cats, function(tc) length(tc))) + disp_amount) + disp_amount)
  
  #get in correct order
  table_to_use <- table_to_use[unlist(tissue_cats),]
  signif_matrix_to_use <- signif_matrix_to_use[unlist(tissue_cats),]
  
}

if(use_range_for_maginal_labels){
  n_genes_in_nodes_label <- apply(apply(n_genes_in_nodes_matrix, 1, range), 2, paste0, collapse = " - ")
  sig_twas_by_trait_genes_label <- apply(apply(sig_twas_by_trait_genes_matrix, 2, range), 2, paste0, collapse = " - ")
} else {
  n_genes_in_nodes_label <- apply(n_genes_in_nodes_matrix, 1, max)
  sig_twas_by_trait_genes_label <- apply(sig_twas_by_trait_genes_matrix, 2, max)
}

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
cols$category <- category_colors

if(use_tissue_cols_for_cols){
  heatmap_cols <- sapply((1:max(table_to_use) / max(table_to_use))^opacity_power_scaler, function(opcty) 
    adjustcolor("black", opcty))
} else {
  heatmap_cols <- viridisLite::viridis(n = max(table_to_use, na.rm = T)*100+1)
  heatmap_cols <- heatmap_cols[round(log(1:max(table_to_use, na.rm = T)) / log(max(table_to_use, na.rm = T)) * max(table_to_use, na.rm = T) * 100 + 1)]
}

#load appropriate model fit
if(base != "deviation_from_expected_logodds_split_the_difference_informative-MVN-matrix-normal_priors"){
  
  base = "deviation_from_expected_logodds_split_the_difference_informative-MVN-matrix-normal_priors"
  load(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/mcmc_output/", 
              base, 
              ifelse(use_all_cats, "_allCats"), 
              ifelse(use_random_DE_genes, "_randomgenes", ""),
              ifelse(!use_kronecker_interactions, "_pairwise-interactions", ""),
              ".cmdStanR.fit"))
  
  samps <- data.frame(as_draws_df(out$draws()))
  prop_greater_than_0 <- function(x) mean(x>0)
  cellbias <- apply(subset_samps("cell_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)
  celltotalbias <- apply(subset_samps("cell_total_prob_bias", c("raw", "sd", "logodds"), samps = samps), 2, prop_greater_than_0)
  rowbias <- apply(subset_samps("row_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)
  colbias <- apply(subset_samps("col_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)
  trait_key <- setNames(trait_categories$new_Phenotype, trait_categories$Tag)
  colcatbias <- apply(subset_samps("colcat_bias", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)
  
  #generate a "significance" matrix
  signif_threshold <- 0.05
  signif_df <- data_subset[1:(nrow(data_subset)), c("tissue", "trait")]
  signif_df$prob_diff_is_positive <- apply(subset_samps("cell_total_prob_bias", c("raw", "sd", "logodds"), samps = samps), 2, prop_greater_than_0)
  signif_df$signif <- 0
  signif_df$signif[signif_df$prob_diff_is_positive > (1-signif_threshold)] <- 1
  signif_df$signif[signif_df$prob_diff_is_positive < signif_threshold] <- -1
  
  signif_df[signif_df$signif != 0,]
  
  signif_matrix <- reshape(signif_df[,-match("prob_diff_is_positive", colnames(signif_df))], idvar = "tissue", timevar = "trait", direction = "wide")
  rownames(signif_matrix) <- signif_matrix$tissue
  signif_matrix <- signif_matrix[,-match("tissue", colnames(signif_matrix))]
  colnames(signif_matrix) <- gsub(colnames(signif_matrix), pattern = "signif.", replacement = "")
  
  #add back in the HYPOTHALAMUS if we removed it
  if(nrow(signif_matrix) != nrow(n_deg_sigtwas_intersect)){
    signif_matrix <- signif_matrix[rownames(n_deg_sigtwas_intersect),]
    rownames(signif_matrix) <- rownames(n_deg_sigtwas_intersect)
    signif_matrix[is.na(signif_matrix)] <- 0
  }
  
  if(!all(colnames(signif_matrix) == colnames(n_deg_sigtwas_intersect)) & all(rownames(signif_matrix) == rownames(n_deg_sigtwas_intersect))){
    stop("something's wrong with the significance matrix")
  }
}

#start the plotting!
cairo_pdf(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/figures/fig4_intersect-enrichment_redux.pdf"), 
          width = 2000 / 72 * ncol(table_to_use) / 80, 
          height = 1300 / 72 + ifelse(group_by_tissue_type, disp_amount * 0.75, 0), 
          family="Arial Unicode MS", pointsize = 18.5)
par(xpd = NA, 
    mar = mars[[1]])
layout(layout_mat, heights = layout_heights)


#initialize blank plot
plot(1, xaxt="n",yaxt="n",bty="n",pch="",ylab="",xlab="", main="", sub="", 
     xlim= c(-5,ncol(table_to_use)), ylim = c(-5,nrow(table_to_use)))


if(group_by_tissue_type){
  tissue_cats_bars_xlocs <- sapply(tissue_cats, function(tc) max(strwidth(tc, units = "user"))) + 0.2
  tissue_cats_bars_xlocs <- rep(max(tissue_cats_bars_xlocs), length(tissue_cats_bars_xlocs))
}

for(ri in 1:nrow(table_to_use)){
  # text(x = 0.5, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0), pos = 2, labels = rownames(table_to_use)[ri], col = MotrpacRatTraining6moData::TISSUE_COLORS[rownames(table_to_use)[ri]])
  text(x = 0.5, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0), pos = 2, labels = rownames(table_to_use)[ri], col = 1)
  for(ci in 1:ncol(table_to_use)){
    if(ri == 1){
      #trait names
      text(x = ci+0.5, y = -0.9, pos = 2, srt = 45, cex = 0.85,
           labels = gsub("craneal", "cranial", gsub("_Scatter", "", trait_categories$new_Phenotype[match(colnames(table_to_use)[ci], 
                                                                                                         trait_categories$Tag)])))
      trait_category <- trait_categories$Category[match(colnames(table_to_use)[ci], trait_categories$Tag)]
      rect(xleft = ci + 1/2,
           xright = ci - 1/2,
           ybottom = ri - 1/2 - 1,
           ytop =  ri + 1/2 - 1,
           col = category_colors[trait_category])
      text(labels = substr(trait_category, 1, 2), y = ri - 1, x = ci,
           col = "white", cex = 0.8)
    }
    
    #vertical total # options
    if(ri == nrow(table_to_use)){
      text(x = ci-0.5, y = ri+1 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), pos = 4, srt = 45,
           labels = sig_twas_by_trait_genes_label[colnames(table_to_use)[ci]], cex = 0.9)
    }
    
    #horiz total # of options
    if(ci == ncol(table_to_use)){
      text(x = ci+0.45, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0), pos = 4,
           labels = n_genes_in_nodes_label[rownames(table_to_use)[ri]])
    }
    
    #the actual cells
    if(use_tissue_cols_for_cols){
      rect(xleft = ci + 1/2,
           xright = ci - 1/2,
           ybottom = ri - 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0),
           ytop =  ri + 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0),
           col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[rownames(table_to_use)[ri]], (table_to_use[ri, ci] / max(table_to_use, na.rm = T))^opacity_power_scaler ))  
    } else {
      rect(xleft = ci + 1/2,
           xright = ci - 1/2,
           ybottom = ri - 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0),
           ytop =  ri + 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0),
           col = heatmap_cols[table_to_use[ri, ci]])
    }
    
    
    
    # rect(xleft = ci + 1/2,
    #      xright = ci - 1/2,
    #      ybottom = ri - 0.475,
    #      ytop =  ri + 0.475,
    #      col = heatmap_cols[table_to_use[ri, ci]],
    #      border = MotrpacRatTraining6moData::TISSUE_COLORS[rownames(table_to_use)[ri]])
    
    #text inside of cells
    if(table_to_use[ri, ci] != 0){
      text_in_cell <- table_to_use[ri, ci]
      ndigs <- nchar(text_in_cell)
      if(incl_significance){
        signif_dir <- c("₋", "","⁺")[match(signif_matrix_to_use[ri, ci], -1:1)]
        # text_in_cell <- paste0(text_in_cell, signif_dir)
        
        #try using corner triangles
        ctr <- corner_triangle_ratio <- 1/3
        cxr = ci + 1/2
        cxl = ci - 1/2
        cyb = ri - 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0)
        cyt =  ri + 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0)
        if(signif_dir == "⁺"){
          polygon(x = c(cxl * (ctr) + cxr * (1-ctr), cxr, cxr),
                  y = c(cyt, cyt, cyt * (1-ctr) + cyb * ctr),
                  col = "red")  
        }
        if(signif_dir == "₋"){
          polygon(x = c(cxl * (ctr) + cxr * (1-ctr), cxr, cxr),
                  y = c(cyb, cyb, cyb * (1-ctr) + cyt * ctr),
                  col = "lightblue")  
        }
      }
      text(text_in_cell, x = ci, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0), 
           col = ifelse((table_to_use[ri, ci] / max(table_to_use, na.rm = T))^opacity_power_scaler > opacity_white_threshold, "black", "white"), 
           cex = 0.85 + ifelse(subset_to_traits, 0.1, 0) - (ndigs-1)/6)
      if(incl_cell_totals){
        text(sig_twas_by_trait_genes_matrix[,colnames(table_to_use)][ri, ci], 
             x = ci - 0.7375, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0) - 0.325, 
             col = "black", cex = 0.2, pos = 4, srt = 90)
        text(n_genes_in_nodes_matrix[,colnames(table_to_use)][ri, ci], 
             x = ci - 0.7, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0) - 0.4, 
             col = "black", cex = 0.2, pos = 4, srt = 0)
      }
    }
    
  }
}

#tissue category bars
if(group_by_tissue_type){
  for(bi in 1:length(tissue_cats_bars_xlocs)){
    segments(x0 = -tissue_cats_bars_xlocs[bi],
             x1 = -tissue_cats_bars_xlocs[bi],
             y0 = tissue_cats_bars_ylocs[bi,1],
             y1 = tissue_cats_bars_ylocs[bi,2],
             lwd = 3)
    bracket_length <- 0.2
    segments(x0 = -tissue_cats_bars_xlocs[bi],
             x1 = -tissue_cats_bars_xlocs[bi]+bracket_length,
             y0 = tissue_cats_bars_ylocs[bi,1],
             y1 = tissue_cats_bars_ylocs[bi,1],
             lwd = 3)
    segments(x0 = -tissue_cats_bars_xlocs[bi],
             x1 = -tissue_cats_bars_xlocs[bi]+bracket_length,
             y0 = tissue_cats_bars_ylocs[bi,2],
             y1 = tissue_cats_bars_ylocs[bi,2],
             lwd = 3)
    
    text(x = -tissue_cats_bars_xlocs[bi], y = mean(tissue_cats_bars_ylocs[bi,]) + ifelse(grepl("_", names(tissue_cats)[bi]), -0.25, 0), pos = 2, 
         labels = gsub("Gi", "GI", stringr::str_to_title(gsub("_", "\n", names(tissue_cats)[bi]))), cex = 1.25)
  }  
}

#legend for heatmap
x_adj <- 2.25
y_adj <- 0
yb_adjust <- ifelse(incl_significance, 3, 0)
n_legend_rects_to_use <- 30
n_legend_labels_to_use <- 10
legend_yvals <- round(seq(0, max(table_to_use), length.out = n_legend_labels_to_use))
legend_ylocs <- seq(yb_adjust, 1 + nrow(table_to_use) + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), length.out = n_legend_rects_to_use)
legend_ycols <- round(seq(1, max(table_to_use), length.out = n_legend_rects_to_use))
for(i in 1:(n_legend_rects_to_use-1)){
  yb = legend_ylocs[i]
  yt = legend_ylocs[i+1]
  # print(paste(c(yb,yt)))
  rect(xleft = ncol(table_to_use) + x_adj + 2 + 1/2,
       xright = ncol(table_to_use) + x_adj + 2 - 1/2,
       ybottom = yb,
       ytop =  yt,
       col = heatmap_cols[legend_ycols[i]], border = NA)
}
for(i in 1:n_legend_labels_to_use){
  text(labels = legend_yvals[i], x = ncol(table_to_use) + x_adj + 2.4, pos = 4, cex = 0.75,
       y = -0.25 + yb_adjust + legend_yvals[i] / max(table_to_use) * (1+nrow(table_to_use) - yb_adjust + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0)))
}

#overall rect for 0
rect(xleft = ncol(table_to_use) + x_adj + 2 + 1/2,
     xright = ncol(table_to_use) + x_adj + 2 - 1/2,
     ybottom = min(legend_ylocs) - diff(range(legend_ylocs))/50,
     ytop =  max(legend_ylocs))

#legend for significance
if(incl_significance){
  # text(labels = paste0("Pr(Enr.) > ", (1 - signif_threshold), " : X⁺\nPr(Dep.) > ", (1 - signif_threshold), " : X₋"),
  #      x = ncol(table_to_use) + x_adj - 1.75,
  #      y = yb_adjust-3.25, pos = 4, cex = 0.75)
  
  signif_label <- paste0("Pr(Dep.) > ", (1 - signif_threshold), " : ")
  ctr <- 0.5
  cxr = ncol(table_to_use) + x_adj - 2.25 + strwidth(signif_label) + 0.45
  cxl = ncol(table_to_use) + x_adj - 2.25 + strwidth(signif_label) - 0.45
  cyb = yb_adjust - 3.75 - 0.45
  cyt = yb_adjust - 3.75 + 0.45
  
  text(labels = signif_label,
       x = ncol(table_to_use) + x_adj - 1.75,
       y = cyb * 0.325 + cyt * 0.675, pos = 4, cex = 0.75)
  rect(xleft = cxl,
       xright = cxr,
       ybottom = cyb,
       ytop =  cyt)
  polygon(x = c(cxl * (ctr) + cxr * (1-ctr), cxr, cxr),
          y = c(cyb, cyb, cyb * (1-ctr) + cyt * ctr),
          col = "lightblue")  
  
  cyb <- cyb + 1.25
  cyt <- cyt + 1.25
  # cxl <- cxl - 0.125
  # cxr <- cxr - 0.125
  
  signif_label2 <- paste0("Pr(Enr.) > ", (1 - signif_threshold), " : ")
  text(labels = signif_label2,
       x = ncol(table_to_use) + x_adj - 1.75,
       y = cyb * 0.675 + cyt * 0.325, pos = 4, cex = 0.75 * strwidth(signif_label) / strwidth(signif_label2))
  rect(xleft = cxl,
       xright = cxr,
       ybottom = cyb,
       ytop =  cyt)
  polygon(x = c(cxl * (ctr) + cxr * (1-ctr), cxr, cxr),
          y = c(cyt, cyt, cyt * (1-ctr) + cyb * ctr),
          col = "red")  
  
} else {
  text(labels = paste0("* indicate \nIHW-adj. \np-vals < 0.05 \nfrom Fisher's Exact Test"),
       x = ncol(table_to_use) + x_adj - 1.75,
       y = yb_adjust-1.25, pos = 4, cex = 0.75, )
}



#legend for trait categories
if(trait_category_legend_below){
  x_adj2 <- 0
  y_adj2 <- -0.5
  for(i in 1:length(categories_represented)){
    rect(xleft = -1/2 + i + cumsum(c(0,strwidth(categories_represented, units = "user") + 1))[i] + x_adj2,
         xright = 1/2 + i + cumsum(c(0,strwidth(categories_represented, units = "user") + 1))[i] + x_adj2,
         ybottom = -11 + y_adj2,
         ytop = -10 + y_adj2,
         col = category_colors[categories_represented[i]], border = 1)
    text(labels = categories_represented[i], pos = 4, y = -10.5 + y_adj2, 
         x = x_adj2 + 0.35 + i + cumsum(c(0,strwidth(categories_represented, units = "user") + 1))[i])
    text(labels = substr(categories_represented, 1, 2)[i], y = -10.5 + y_adj2, 
         x = i + cumsum(c(0,strwidth(categories_represented, units = "user") + 1))[i] + x_adj2,
         col = "white", cex = 0.8)
    #draw circle?
    arc(t1 = 3*pi/2, t2 = pi/2, r1 = (10.5-y_adj2) / 2, r2 = (10.5-y_adj2) / 2, center = c(0,(-10.5 + y_adj2)/2), 
        lwd = 2, res = 100, adjx = ifelse(order_by_counts, 1, 1.1))
    points(0, 0, pch = -9658, cex = 2)
  }
} else {
  x_adj2 <- x_adj - 2.5
  y_adj2 <- y_adj - 2.5
  for(i in 1:length(category_colors)){
    rect(xleft = 0 + i,
         xright = 1 + i,
         ybottom = -10,
         ytop = -9,
         col = category_colors[i], border = 1)
    text(labels = names(category_colors)[i], pos = 4, y = i + y_adj2 - 7, x = ncol(table_to_use) - 1/2 + x_adj2 + 1/3)
  }
}

#labels
#horiz label for total
segments(x0 = -2, x1 = 0.5, y0 = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0),
         y1 = nrow(table_to_use) + 1 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), lwd = 2)
segments(x0 = -2, x1 = 0.5, y0 = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), 
         y1 = nrow(table_to_use) + 2.5 + ifelse(use_range_for_maginal_labels, 0, -0.5) + 
           ifelse(order_by_counts, 0, -0.8) + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0) + 
           ifelse(subset_to_traits, 1, 1.5), lwd = 2)
text(x = -2, y = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), pos = 2, 
     labels = paste0(ifelse(use_range_for_maginal_labels, "total", "max"), " # of\nPrediXcan hits"))

#vertical label for total
segments(x0 = ncol(table_to_use) + 2.5, x1 = ncol(table_to_use) + 1, 
         y0 = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), 
         y1 = nrow(table_to_use) + 0.75 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), lwd = 2)
segments(x0 = ncol(table_to_use) + 2.5, x1 = ncol(table_to_use) + 2 + 
           ifelse(subset_to_traits, 0, 1.25), 
         y0 = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), 
         y1 = nrow(table_to_use) + 0.75 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), lwd = 2)
text(x = ncol(table_to_use) + 2 + ifelse(order_by_counts, 0, 2), y = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), pos = 3, 
     labels = paste0(ifelse(use_range_for_maginal_labels, "total", "max"), " # of DEGs"))


#legend and title labels
text(labels = ifelse(use_counts, latex2exp::TeX("$n_{intersect}$"), latex2exp::TeX("‰_{ intersect}")), pos = 3, font = 2, cex = 1.25,
     x = ncol(table_to_use) + x_adj + 2, y = nrow(table_to_use) + y_adj + 0.875 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0))
if(use_counts){
  text(latex2exp::TeX(paste0("number of genes in 8w - F↓M↓+ 8w - F↑M↑ with IHW-significant PrediXcan at $\\alpha$ = 0.05")), 
       x = 2 + ifelse(subset_to_traits, 0, 20), 
       y = nrow(table_to_use) + 3.25 + 
         ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0) + 
         ifelse(subset_to_traits, 0, 1), pos = 4, cex = 1.5, font = 2)
} else {
  text(latex2exp::TeX(paste0("proportion (‰) of IHW significant TWAS hits at $\\alpha$ = 0.05 in 8w - F↓M↓ or 8w - F↑M↑")), 
       x = 0, y = nrow(table_to_use) + 4.25 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), pos = 4, cex = 2.35, font = 2)
}


# ~~~~~~~~~~~~ #
# scatterplots #
# ~~~~~~~~~~~~ #

use_focal_vs_compl <- T

data1 <- data.frame(count = as.integer(unlist(c(n_deg_sigtwas_intersect))))
data1$tissue <- rep(rownames(n_deg_sigtwas_intersect), ncol(n_deg_sigtwas_intersect))
data1$trait <- unlist(lapply(colnames(n_deg_sigtwas_intersect), function(tri) rep(tri, nrow(n_deg_sigtwas_intersect))))
data1$total <- n_genes_in_nodes[data1$tissue]
data1$TWAS_Hit <- "YES"

traits <- unique(data1$trait)
tissues <- unique(data1$tissue)
d <- list(cell_count = data1$count,
          total = sapply(1:nrow(data1), function(i) total_number_of_possible_hits_matrix[data1$tissue[i], data1$trait[i]]),
          row_count = sapply(1:nrow(data1), function(i) n_genes_in_nodes_matrix[data1$tissue[i], data1$trait[i]]),
          col_count = sapply(1:nrow(data1), function(i) sig_twas_by_trait_genes_matrix[data1$tissue[i], data1$trait[i]]),
          row_index = match(data1$tissue, tissues),
          col_index = match(data1$trait, traits),
          row_n = length(unique(data1$tissue)),
          col_n = length(unique(data1$trait)))




# category_shapes <- setNames(15:19, salient.categories)
category_shapes <- setNames(15:19, c("Cardiometabolic", "Aging", "Anthropometric", 
                                     "Immune", "Psychiatric-neurologic")
)
category_shapes <- setNames(c(24,19,15,43,23,25)[1:length(focal_traitcats)], focal_traitcats)



par(mar = mars[[2]], pty="s")
# layout(t(c(rep(1,10), rep(2,10), 3)))
# layout(t(c(1,2,3)), widths = c(1,1,0.1))
#first plot


# sapply(1:nrow(signif_df), function(i) adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index[i]]], as.numeric(signif_df$signif[i] != 0) * 0.625 + 0.125))
if(use_focal_vs_compl){
  xvals <- (d$col_count - d$cell_count) / (d$total - d$row_count)
  yvals <- d$cell_count / d$row_count
} else {
  xvals <- d$row_count / d$total * d$col_count / d$total  
  yvals <- d$cell_count / d$total
}

# plot(logit(xvals), logit(yvals), xlim = ifelse2(use_focal_vs_compl, c(-14,-1), c(-15,-4.25)), 
#      ylim = ifelse2(use_focal_vs_compl, c(-6, -1), c(-9.5,-4.25)),
#      pch = category_shapes[traitwise_partitions$Category[match(traits[d$col_index], traitwise_partitions$Tag)]], 
#      col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5), bg = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5),
#      cex = 1.5,
#      xlab = "", ylab = "", yaxt = 'n', xaxt = "n")
plot(logit(xvals), logit(yvals), xlim = ifelse2(use_focal_vs_compl, c(-14,-1), c(-15,-4.25)), 
     ylim = ifelse2(use_focal_vs_compl, c(-6, -1), c(-9.5,-4.25)),
     xlab = "", ylab = "", yaxt = 'n', xaxt = "n", type = "n")
#can try rasterizing points?
# raster_points(x = logit(xvals), y =  logit(yvals),
#               pch = category_shapes[traitwise_partitions$Category[match(traits[d$col_index], traitwise_partitions$Tag)]], 
#               col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5), bg = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5),
#               cex = 1.5, raster_res = 200)
#smaller file with regular points, not that much data fundamentally
points(x = logit(xvals), y =  logit(yvals),
       pch = category_shapes[traitwise_partitions$Category[match(traits[d$col_index], traitwise_partitions$Tag)]], 
       col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5), bg = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5),
       cex = 1.5)
xw <- diff(par("usr")[1:2])
yh <- diff(par("usr")[3:4])
box(which = "plot", lwd = 2)
#axes
text(labels = "focal logit(frequency)", x = par("usr")[2] + diff(par("usr")[1:2])/8, srt = 270, pos = 1, 
     y = mean(par("usr")[3:4]), xpd = NA, cex = 1.25)
text(labels = "complement logit(frequency)", x = mean(par("usr")[1:2]), srt = 0, pos = 1, 
     y = par("usr")[3] - diff(par("usr")[3:4])/12, xpd = NA, cex = 1.25)
yaxlocs <- ifelse2(use_focal_vs_compl, c(-6:-1), c(-5:-9))
segments(x0 = par("usr")[2], x1 = par("usr")[2] + diff(par("usr")[1:2])/100, y0 = yaxlocs, y1 = yaxlocs, lwd = 2, xpd = NA)
text(labels = yaxlocs, x = par("usr")[2] + diff(par("usr")[1:2])/200, srt = 0, pos = 4, 
     y = yaxlocs, xpd = NA, cex = 1.25)
xaxlocs <- ifelse2(use_focal_vs_compl, -7:-1*2, -7:-3*2)
segments(x0 = xaxlocs, x1 = xaxlocs, y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4])/80, lwd = 2, xpd = NA)
text(labels = xaxlocs, x = xaxlocs, srt = 0, pos = 1, 
     y = par("usr")[3] - diff(par("usr")[3:4])/80, xpd = NA, cex = 1.25)
#1 to 1 line
abline(0,1, lty = 2, lwd = 4, col = adjustcolor(1,0.75), xpd = F)

#add arrows for enriched vs depleted
arrow_locs_enrich <- c(-3, -3, -2.25, -1.5)
grad_arrow_curve(arrow_locs_enrich, direc = "v", cols = c("white", "red"), 
                 prop_head_width = 3, prop_shaft_length = 0.2, taper_ratio = 0.7, outline_col = "white",
                 w = 0.6)
text(x = mean(arrow_locs_enrich[1:2]), y = sum(arrow_locs_enrich[3:4] * c(0.35,0.65)), 
     labels = "enriched", col = "white", srt = 90, cex = 0.5)

arrow_locs_deplete <- c(-2, -2, -3.5, -4.25)
grad_arrow_curve(arrow_locs_deplete, direc = "v", cols = c("white", "blue"), 
                 prop_head_width = 3, prop_shaft_length = 0.2, taper_ratio = 0.7, outline_col = "white",
                 w = 0.6)

text(x = mean(arrow_locs_deplete[1:2]), y = sum(arrow_locs_deplete[3:4] * c(0.35,0.65)), 
     labels = "depleted", col = "white", srt = 270, cex = 0.5)



#secondplot
xl <- par("usr")[1]
xr <- par("usr")[1] + diff(par("usr")[1:2]) / 2
yt <- par("usr")[4]
yb <- par("usr")[4] - diff(par("usr")[3:4]) / 2

buffer <- 0.075
bx <- (xr-xl) * buffer
by <- (yt-yb) * buffer
rect(xleft = xl, ybottom = yb, xright = xr, ytop = yt, lwd = 2, xpd = NA)

points(x = ((xvals - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
       y = ((yvals - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, 
       pch = category_shapes[traitwise_partitions$Category[match(traits[d$col_index], traitwise_partitions$Tag)]], 
       col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5), bg = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues[d$row_index]], 0.5), 
       cex = 1)
one2one_coords <- c(x0 = 0 * (xr - xl - 2*bx) + xl + bx,
                    y0 = 0 * (yt - yb - 2*by) + yb + by,
                    x1 = ((max(yvals, na.rm = T) - 
                             min(xvals, na.rm = T)) / 
                            max(xvals, na.rm = T)) * 
                      (xr - xl - 2*bx) + xl + bx,
                    y1 = ((max(yvals, na.rm = T) - 
                             min(yvals, na.rm = T)) / 
                            max(yvals, na.rm = T)) * 
                      (yt - yb - 2*by) + yb + by)
one2one_slope <- (one2one_coords["y1"] - one2one_coords["y0"]) / 
  (one2one_coords["x1"] - one2one_coords["x0"])
one2one_coords["y1"] <- one2one_coords["y1"] - (one2one_coords["x1"] - xr) * one2one_slope
one2one_coords["x1"] <- xr

segments(x0 = one2one_coords["x0"],
         y0 = one2one_coords["y0"],
         x1 = one2one_coords["x1"],
         y1 = one2one_coords["y1"],
         lty = 2, lwd = 2, col = adjustcolor(1,0.75))

#axes
text(labels = "focal frequency", x = par("usr")[1] - diff(par("usr")[1:2])/ifelse2(use_focal_vs_compl, 8, 7), srt = 90, pos = 1, 
     y = par("usr")[4] - diff(par("usr")[3:4])/5, xpd = NA, cex = 1)
text(labels = "complement frequency", x = par("usr")[1] + diff(par("usr")[1:2])/4, srt = 0, pos = 1, 
     y = par("usr")[4] + diff(par("usr")[3:4])/8, xpd = NA, cex = 1)
xaxlocs <- ifelse2(use_focal_vs_compl, c(0:4/20), c(0:6/500))
yaxlocs <- ifelse2(use_focal_vs_compl, c(0:6/20), c(0:6/500))

segments(x0 = par("usr")[1], x1 = par("usr")[1] - diff(par("usr")[1:2])/100, 
         y0 = ((yaxlocs - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, 
         y1 = ((yaxlocs - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, 
         lwd = 2, xpd = NA)
segments(y0 = par("usr")[4], y1 = par("usr")[4] + diff(par("usr")[3:4])/100, 
         x0 = ((xaxlocs - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
         x1 = ((xaxlocs - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
         lwd = 2, xpd = NA)
text(labels = xaxlocs, x = ((xaxlocs - min(xvals, na.rm = T)) / max(xvals, na.rm = T)) * (xr - xl - 2*bx) + xl + bx, 
     srt = 0, pos = 3, 
     y = par("usr")[4], xpd = NA, cex = 0.75)
text(labels = yaxlocs, x = par("usr")[1], 
     srt = 0, pos = 2, 
     y = ((yaxlocs - min(yvals, na.rm = T)) / max(yvals, na.rm = T)) * (yt - yb - 2*by) + yb + by, xpd = NA, cex = 0.75)
# text(labels = 0:6/500, x = -7:-3*2, srt = 0, pos = 1, 
#      y = par("usr")[3] - diff(par("usr")[3:4])/80, xpd = NA, cex = 1.25)
# abline(0,1, lty = 2, lwd = 6, col = adjustcolor(1,0.75))

#legends
legend(x = xl, y = yb, legend = names(category_shapes), pch = category_shapes, col = adjustcolor(1, 0.5), pt.bg = adjustcolor(1, 0.5), cex = 0.75, pt.cex = 1.4)
legend(x = xr, y = yt, legend = tissues, pch = 19, col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues], 0.5), cex = 0.55, pt.cex = 1.2)
segments(x0 = par("usr")[1] + diff(par("usr")[1:2]) / 1.55, 
         y0 = par("usr")[4] - diff(par("usr")[3:4])/50,
         x1 = par("usr")[1] + diff(par("usr")[1:2]) / 1.45, 
         y1 = par("usr")[4] - diff(par("usr")[3:4])/50, 
         lty = 2, lwd = 2, col = adjustcolor(1,0.75))
text(x = par("usr")[1] + diff(par("usr")[1:2]) / 1.45, 
     y = par("usr")[4] - diff(par("usr")[3:4])/45, 
     pos = 4, cex = 0.75, labels = "1-to-1 line")


#make second plot with proportion disease-like effects data
par(mar = mars[[3]])
data <- fread("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/prop_positive_effects_in_DEG-TWAS_intersection.txt")
xvals <- data$count_all / data$total_all
yvals <- data$count_inters / data$total_inters


cex_pow <- 1/2
plot(xvals, yvals, xlim = c(0.48, 0.52), 
     ylim = c(0,1),
     pch = category_shapes[traitwise_partitions$Category[match(data$trait, traitwise_partitions$Tag)]], 
     col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[data$tissue], 0.5), bg = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[data$tissue], 0.5), 
     cex = (data$total_inters / max(data$total_inters))^cex_pow * 5,
     xlab = "", ylab = "", yaxt = 'n', xaxt = "n")
box(which = "plot", lwd = 2)

#axes
text(labels = "focal frequency (+ effects)", x = par("usr")[1] - diff(par("usr")[1:2])/7, srt = 90, pos = 1, 
     y = mean(par("usr")[3:4]), xpd = NA, cex = 1.25)
text(labels = "complement frequency (+ effects)", x = mean(par("usr")[1:2]), srt = 0, pos = 1, 
     y = par("usr")[3] - diff(par("usr")[3:4])/12, xpd = NA, cex = 1.25)
yaxlocs <- seq(0,1,by=0.1)
segments(x0 = par("usr")[1], x1 = par("usr")[1] - diff(par("usr")[1:2])/100, y0 = yaxlocs, y1 = yaxlocs, lwd = 2, xpd = NA)
text(labels = yaxlocs, x = par("usr")[1] - diff(par("usr")[1:2])/200, srt = 0, pos = 2, 
     y = yaxlocs, xpd = NA, cex = 1.25)
xaxlocs <- round(seq(par("usr")[1], par("usr")[2], by = 0.01), 2)
segments(x0 = xaxlocs, x1 = xaxlocs, y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4])/80, lwd = 2, xpd = NA)
text(labels = xaxlocs, x = xaxlocs, srt = 0, pos = 1, 
     y = par("usr")[3] - diff(par("usr")[3:4])/80, xpd = NA, cex = 1.25)
# abline(0,1, lty = 2, lwd = 6, col = adjustcolor(1,0.75))
abline(h=0.5,lwd=2,lty=2,col=adjustcolor(1,0.5),xpd=F)
abline(v=0.5,lwd=2,lty=2,col=adjustcolor(1,0.5),xpd=F)

#legend
xl <- par("usr")[2]
xr <- par("usr")[2] + diff(par("usr")[1:2]) / 2
yt <- par("usr")[4]
yb <- par("usr")[4] - diff(par("usr")[3:4]) / 2.5
legend(x = xl, y = yb + (yt-yb)/20, legend = names(category_shapes), bty = "n",
       pch = category_shapes, col = adjustcolor(1, 0.5), pt.bg = adjustcolor(1, 0.5), cex = 0.75, pt.cex = 1.4, xpd = NA)
legend(x = xl + (xr-xl)/100, y = yt, legend = tissues, pch = 19, bty = "n", xpd = NA,
       col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues], 0.5), cex = 0.55, pt.cex = 1.2)

legcexs <- c(0.5, 1:5)
legvals <- (round((legcexs / 5) ^ (1/cex_pow) * max(data$total_inters)))
legcexs <- (legvals / max(data$total_inters))^cex_pow * 5
# points(x = rep(xl + (xr-xl)/11, 6) + (rev(legcexs)) / 4000, y = yb - (yt-yb) * 1.4 + 
#          cumsum(rep((yt-yb)/60, length(legcexs))) + cumsum(legcexs + c(0,legcexs[-length(legcexs)])) / 120, 
#        cex = legcexs, xpd = NA, pch = 19, col = adjustcolor(1, 0.5))
text(labels = latex2exp::TeX("$n_{intersect}$"), x = xl - (xr-xl)/50, y = yb - (yt-yb) * 0.6, pos = 4, xpd = NA)
fixed_incr <- (yt-yb)/25
y_disp <- (yt-yb) * 1.5
points(x = rep(xl + (xr-xl)/10, 6), y = yb - y_disp + 
         cumsum(rep(fixed_incr, length(legcexs))) + cumsum(legcexs + c(0,legcexs[-length(legcexs)])) / 120, 
       cex = legcexs, xpd = NA, pch = 19, col = adjustcolor(1, 0.5))
text(x = rep(xl + (xr-xl)/10, 6), y = yb - y_disp + 
       cumsum(rep(fixed_incr, length(legcexs))) + cumsum(legcexs + c(0,legcexs[-length(legcexs)])) / 120, 
     cex = legcexs / 5, xpd = NA, labels = legvals, col = "white")
text(labels = legvals[1:2], x = xl + (xr-xl)/15 + cumsum(legcexs)[1:2]/4000, y = yb - y_disp + 
       cumsum(rep(fixed_incr, length(legcexs)))[1:2] + cumsum(legcexs + c(0,legcexs[-length(legcexs)]))[1:2] / 120, 
     pos = 4, xpd = NA, cex = 0.5)

#add arrows for enriched vs depleted
arrow_locs_enrich <- c(0.519, 0.519, 0.525, 0.525 + 0.75 * diff(par("usr")[3:4]) / yh)
grad_arrow_curve(arrow_locs_enrich, direc = "v", cols = c("white", "red"), 
                 prop_head_width = 3, prop_shaft_length = 0.2, taper_ratio = 0.7, outline_col = "white",
                 w = 0.6 * diff(par("usr")[1:2]) / xw)
text(x = mean(arrow_locs_enrich[1:2]), y = sum(arrow_locs_enrich[3:4] * c(0.35,0.65)), 
     labels = "enriched", col = "white", srt = 90, cex = 0.5)

arrow_locs_deplete <- c(0.519, 0.519, 0.475, 0.475 - 0.75 * diff(par("usr")[3:4]) / yh)
grad_arrow_curve(arrow_locs_deplete, direc = "v", cols = c("white", "blue"), 
                 prop_head_width = 3, prop_shaft_length = 0.2, taper_ratio = 0.7, outline_col = "white",
                 w = 0.6 * diff(par("usr")[1:2]) / xw)

text(x = mean(arrow_locs_deplete[1:2]), y = sum(arrow_locs_deplete[3:4] * c(0.35,0.65)), 
     labels = "depleted", col = "white", srt = 270, cex = 0.5)



# ~~~~~~~~~~~~ #
# violin plots #
# ~~~~~~~~~~~~ #

par(mar = mars[[4]], pty="m")
focal_samps <- subset_samps("row_bias", c("raw", "sd"), samps = samps) + samps$overall_bias
n_to_trim <- round(nrow(focal_samps) * 0.002)
focal_samps <- apply(focal_samps, 2, trim_n, n_to_trim)

tissues <- tissues_intersect.model
colnames(focal_samps) <- tissues
tord <- order(apply(focal_samps, 2, mean))
qi_95 <- apply(focal_samps, 2, quantile, probs = c(0.05, 0.95))[,tord]
qi_100 <- apply(focal_samps, 2, quantile, probs = c(0.0, 1.0))[,tord]
posterior_means <- apply(focal_samps, 2, mean)[tord]
inside_cols <- rep("black", length(tord))
inside_cols[overlaps_with_zero(qi_95)] <- "white"
tmp <- vioplot::vioplot(x = focal_samps[,tord], T, 
                        col = MotrpacRatTraining6moData::TISSUE_COLORS[colnames(focal_samps)][tord], outline=FALSE, yaxt = "n",
                        names = rep("", length(tord)), range = 0, ylab = "", lineCol = MotrpacRatTraining6moData::TISSUE_COLORS[colnames(focal_samps)][tord],
                        xlab = "", cex.lab = 2, plotCentre = "point", 
                        colMed = MotrpacRatTraining6moData::TISSUE_COLORS[colnames(focal_samps)][tord],
                        horizontal = T) 


#axes
xtickvals <- seq3(range(qi_100), 5, 0)
segments(x0 = xtickvals, x1 = xtickvals, y0 = par("usr")[3], 
         y1 = par("usr")[3] - diff(par("usr")[3:4])/50, xpd = NA)
text(x = xtickvals, y = par("usr")[3] - diff(par("usr")[3:4])/25, labels = xtickvals)
text(x = mean(par("usr")[1:2]), y = par("usr")[3]-diff(par("usr")[3:4])/9, 
     labels = latex2exp::TeX("multilevel \\textbf{tissue} enrichment"), cex = 1.25, xpd = NA)
tick <- seq_along(colnames(focal_samps))
axis(2, at = tick, labels = F)
for(i in 1:length(colnames(focal_samps))){
  segments(y0 = i, y1 = i, x0 = qi_95[1,i], x1 = qi_95[2,i], lwd = 4, col = inside_cols[i])
  points(y = i, x = posterior_means[i], col = inside_cols[i], pch = 19, cex = 1.5)
  
  #guiding lines
  segments(y0 = i, y1 = i, x0 = min(qi_100), x1 = qi_100[1,i], lwd = 1, col = "black", xpd = T, lty = 3)
  segments(y0 = i, y1 = i, x0 = qi_100[2,i], x1 = max(qi_100), lwd = 1, col = "black", xpd = T, lty = 3)
  
}
#group labels
text(y = tick, x = rep(par("usr")[1] - 0.01, length(tick)), colnames(focal_samps)[tord], srt = 0, xpd = T, pos = 2,
     col = MotrpacRatTraining6moData::TISSUE_COLORS[colnames(focal_samps)[tord]], xpd = NA, cex = 1.25)

abline(v=0,lwd=3,lty=2, col = adjustcolor(1,0.5), xpd = F)

#legend for violin plots
xval <- -1.5
yval <- 5.75 #2
yh <- 4
yt <- yval + yh/2
yb <- yval - yh/2
xr <- -30:30/10
xvals_disp <- (dnorm(xr) / 6) * 2
xvals <- cbind(xval + xvals_disp, xval - xvals_disp)
yvals <- seq(yb,yt,length.out = nrow(xvals))
# lines(xvals[,1], yvals, xpd = NA)
# lines(xvals[,2], yvals, xpd = NA)
polygon(x = c(xvals[,1], xvals[,2]), y = c(yvals, rev(yvals)), xpd = NA, col = adjustcolor(1,0.15))
ybloc <- yb + (qnorm(p = c(0.05,0.95)) - min(xr)) / diff(range(xr)) * yh
segments(x0 = xval, x1 = xval, y0 = ybloc[1], y1 = ybloc[2], xpd = NA, lwd = 2)
points(xval, yval, pch = 19, xpd = NA)
xlab_loc <- max(xvals) + diff(range(xvals)) / 5
segments(x0 = xval, x1 = xlab_loc, 
         y0 = c(max(yvals), min(yvals), yval, ybloc),
         y1 = c(max(yvals), min(yvals), yval, ybloc),
         xpd = NA, lty = 2)
violabs <- paste0("$P_{", c("99.9", "0.01", "\\mu", "5", "95"), "}$")
violabs[[3]] <- "$\\mu$"
text(latex2exp::TeX(violabs), x = xlab_loc - diff(range(xvals)) / 10, y = c(max(yvals), min(yvals), yval, ybloc), pos = 4, xpd = NA, cex = 0.875)

segments(x0 = min(xvals) - diff(range(xvals))/100, x1 = max(xvals) + diff(range(xvals))/100, 
         y0 = min(yvals) - diff(range(yvals)) / 8, y1 = min(yvals) - diff(range(yvals)) / 8, xpd = NA, lwd = 3.25)
points(xval, min(yvals) - diff(range(yvals)) / 8, pch = 19, xpd = NA, cex = 1.2)
segments(x0 = min(xvals), x1 = max(xvals), y0 = min(yvals) - diff(range(yvals)) / 8, y1 = min(yvals) - diff(range(yvals)) / 8, xpd = NA, lwd = 2, col = "white")
points(xval, min(yvals) - diff(range(yvals)) / 8, pch = 19, xpd = NA, cex = 1, col = "white")
text(latex2exp::TeX("$CI_{90}$ includes 0$"), x = max(xvals) + diff(range(xvals))/100, 
     y = min(yvals) - diff(range(yvals)) / 7, pos = 4, xpd = NA, cex = 0.875)

segments(x0 = min(xvals) - diff(range(xvals))/100, x1 = max(xvals) + diff(range(xvals))/100, 
         y0 = min(yvals) - diff(range(yvals)) / 8 * 2, y1 = min(yvals) - diff(range(yvals)) / 8 * 2, xpd = NA, lwd = 3.25)
points(xval, min(yvals) - diff(range(yvals)) / 8 * 2, pch = 19, xpd = NA, cex = 1.2)
text(latex2exp::TeX("$CI_{90}$ excludes 0$"), x = max(xvals) + diff(range(xvals))/100, 
     y = min(yvals) - diff(range(yvals)) / 7 - diff(range(yvals)) / 8, 
     pos = 4, xpd = NA, cex = 0.875)


###
#colcat bias violin plots
###

salient.categories <- trait_cats
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
cols$category <- category_colors

par(mar = mars[[5]])
focal_samps <- subset_samps("colcat_bias", c("raw", "sd"), samps = samps) + samps$overall_bias
focal_samps <- apply(focal_samps, 2, trim_n, n_to_trim)
colnames(focal_samps) <- trait_cats
tord <- order(apply(focal_samps, 2, mean))
qi_95 <- apply(focal_samps, 2, quantile, probs = c(0.05, 0.95))[,tord]
qi_100 <- apply(focal_samps, 2, quantile, probs = c(0.0, 1.0))[,tord]
posterior_means <- apply(focal_samps, 2, mean)[tord]
inside_cols <- rep("black", length(tord))
inside_cols[overlaps_with_zero(qi_95)] <- "white"

par(xpd=F) #vioplot keeps plotting extra stuff
tmp <- vioplot::vioplot(x = focal_samps[,tord], T, 
                        col = cols$category[colnames(focal_samps)][tord], outline=FALSE, yaxt = "n",
                        names = rep("", length(tord)), range = 0, ylab = "", 
                        lineCol = cols$category[colnames(focal_samps)][tord],
                        xlab = "", cex.lab = 2, plotCentre = "point", 
                        colMed = cols$category[colnames(focal_samps)][tord],
                        horizontal = T)

#axes
xtickvals <- seq3(range(qi_100), 5, 0)
segments(x0 = xtickvals, x1 = xtickvals, y0 = par("usr")[3], 
         y1 = par("usr")[3] - diff(par("usr")[3:4])/40, xpd = NA)
text(x = xtickvals, y = par("usr")[3] - diff(par("usr")[3:4])/20, labels = xtickvals, xpd = NA)
text(x = mean(par("usr")[1:2]), y = par("usr")[3]-diff(par("usr")[3:4])/7.5, 
     labels = latex2exp::TeX("multilevel \\textbf{category} enrichment"), cex = 1.25, xpd = NA)
tick <- seq_along(colnames(focal_samps))
axis(2, at = tick, labels = F)
for(i in 1:length(colnames(focal_samps))){
  segments(y0 = i, y1 = i, x0 = qi_95[1,i], x1 = qi_95[2,i], lwd = 4, col = inside_cols[i])
  points(y = i, x = posterior_means[i], col = inside_cols[i], pch = 19, cex = 1.5)
  
  #guiding lines
  segments(y0 = i, y1 = i, x0 = min(qi_100), x1 = qi_100[1,i], lwd = 1, col = "black", xpd = T, lty = 3)
  segments(y0 = i, y1 = i, x0 = qi_100[2,i], x1 = max(qi_100), lwd = 1, col = "black", xpd = T, lty = 3)
  
}
catlabs <- colnames(focal_samps)[tord]
# catlabs <- gsub("system ", "system\n", gsub("-", "-\n", catlabs))
catlabs <- gsub("system ", "", gsub("-neurologic", "", catlabs))
text(y = tick, x = rep(par("usr")[1] - 0.01, length(tick)), 
     labels = catlabs, 
     srt = 0, xpd = T, pos = 2, cex = 1.25,
     col = cols$category[colnames(focal_samps)[tord]], xpd = NA)


gsub("_", "\n", colnames(focal_samps))
gsub("Gi", "GI", stringr::str_to_title(gsub("_", "\n", colnames(focal_samps))))

abline(v=0,lwd=3,lty=2, col = adjustcolor(1,0.5), xpd = F)

##
# traits
##

par(mar = mars[[6]])
trait_subset_bool <- traitwise_partitions$Category[match(traits, traitwise_partitions$Tag)] %in% salient.categories
focal_samps <- subset_samps("col_bias", c("raw", "sd"), samps = samps) + samps$overall_bias +
  subset_samps("colcat_bias", c("raw", "sd"), samps = samps)[,match(trait_categories$Category[match(traits, trait_categories$Tag)], trait_cats)]
focal_samps <- apply(focal_samps, 2, trim_n, n_to_trim)
colnames(focal_samps) <- traits
tord <- order(apply(focal_samps, 2, mean))
qi_95 <- apply(focal_samps, 2, quantile, probs = c(0.05, 0.95))[,tord]
qi_100 <- apply(focal_samps, 2, quantile, probs = c(0.0, 1.0))[,tord]
posterior_means <- apply(focal_samps, 2, mean)[tord]
inside_cols <- rep("black", length(tord))
inside_cols[overlaps_with_zero(qi_95)] <- "white"

tmp <- vioplot::vioplot(x = focal_samps[,tord], T,
                        col = cols$category[trait_categories$Category[match(colnames(focal_samps)[trait_subset_bool][tord], trait_categories$Tag)]],
                        outline=FALSE, yaxt = "n",
                        names = rep("", length(tord)), range = 0, ylab = "",
                        lineCol = cols$category[trait_categories$Category[match(colnames(focal_samps)[trait_subset_bool][tord], trait_categories$Tag)]],
                        xlab = "", cex.lab = 2, plotCentre = "point",
                        colMed = cols$category[colnames(focal_samps)][tord],
                        horizontal = F,
                        xlim = c(4,ncol(focal_samps)-3))

ytickvals <- seq3(range(qi_100), 5, 0)
segments(y0 = ytickvals, y1 = ytickvals, x0 = par("usr")[1], x1 = par("usr")[1] - diff(par("usr")[1:2])/200, xpd = NA)
text(y = ytickvals, x = par("usr")[1] - diff(par("usr")[1:2])/75, labels = ytickvals, xpd = NA)
text(y = mean(par("usr")[3:4]), x = par("usr")[1]-diff(par("usr")[1:2])/30,
     labels = latex2exp::TeX("multilevel \\textbf{trait} enrichment"), cex = 1.25, xpd = NA, srt = 90)
tick <- seq_along(colnames(focal_samps))


tick <- seq_along(traits[trait_subset_bool])
axis(1, at = tick, labels = F)
for(i in 1:length(traits[trait_subset_bool])){
  segments(x0 = i, x1 = i, y0 = qi_95[1,i], y1 = qi_95[2,i], lwd = 2, col = inside_cols[i])
  points(x = i, y = posterior_means[i], col = inside_cols[i], pch = 19)
  #guiding lines
  segments(x0 = i, x1 = i, y0 = min(qi_100) - diff(range(qi_100)) / 30, y1 = qi_100[1,i], lwd = 1, col = "black", xpd = T, lty = 3)
  segments(x0 = i, x1 = i, y0 = qi_100[2,i], y1 = max(qi_100) * 0.75 + 1 * 0.25, lwd = 1, col = "black", xpd = T, lty = 3)
}
trait_labs <- trait_categories$new_Phenotype[match(colnames(focal_samps)[trait_subset_bool][tord], trait_categories$Tag)]
trait_labs[which.max(nchar(trait_labs))]
trait_labs[trait_labs == "Attention_Deficit_Hyperactivity_Disorder"] <- "ADHD"
trait_labs <- gsub("craneal", "cranial", trait_labs)
text(tick + 0.55, par("usr")[3] - diff(par("usr")[3:4])/12, cex = 0.8,
     labels = trait_labs, 
     srt = 45, xpd = T, pos = 2,
     col = cols$category[trait_categories$Category[match(colnames(focal_samps)[trait_subset_bool][tord], trait_categories$Tag)]], xpd = NA)
abline(h=0,lwd=3,lty=2, col = adjustcolor(1,0.5), xpd = F)

#histogram for trait x gene effects
par(mar = mars[[7]])
par(xpd = NA)
hist(celltotalbias, breaks = 0:20/20, main = "", ylab = "",
     xlab = "", freq = T, yaxt = "n")
axis(2, at = 0:3/20*length(celltotalbias), labels = 0:3/20, xpd = NA)
text(x = par("usr")[1]-diff(par("usr")[1:2])/5, y = mean(par("usr")[3:4]), 
     "Density", cex = 1.25, srt = 90)
text(y = par("usr")[3]-diff(par("usr")[3:4])/2.25, x = mean(par("usr")[1:2]), 
     "Proportion Posterior Mass > 0" , cex = 1.25)
hist(celltotalbias[celltotalbias > 0.95], breaks = 0:20/20, main = "",
     xlab = "", freq = T, add = T, col = "red")
hist(celltotalbias[celltotalbias < 0.05], breaks = 0:20/20, main = "",
     xlab = "", freq = T, add = T, col = "lightblue")
text(latex2exp::TeX("\\textbf{Tissue} × \\textbf{Trait} Enrichment"), x = 0.5, y = par("usr")[4],
     pos = 3, xpd = NA, cex = 1.25)

#figure labels
fig_lab(label = "a)", xp = 2, yp = 99, draw_grid = F)
fig_lab(label = "b)", xp = 79, yp = 99)
fig_lab(label = "c)", xp = 4, yp = 59.5)
fig_lab(label = "d)", xp = 39.5, yp = 59.5)
fig_lab(label = "e)", xp = 81, yp = 84)
fig_lab(label = "f)", xp = 81, yp = 51)
fig_lab(label = "g)", xp = 5, yp = 24)
# fig_label("d)", addX = 0, addY = 0, region = "device", cex = 2, xpd = NA)
# fig_label("e)", addX = 0, addY = 0, region = "device", cex = 2, xpd = NA)
# fig_label("f)", addX = 0, addY = 0, region = "device", cex = 2, xpd = NA)
# fig_label("g)", addX = 0, addY = 0, region = "device", cex = 2, xpd = NA)


dev.off()


