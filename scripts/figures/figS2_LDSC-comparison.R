#### run preprocessing script ####
run_preprocessing_scripts <- F #or load the data directly
if(run_preprocessing_scripts){
  figure_id <- "S2"
  source("/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/helper_scripts/figure_set_2_preprocessing.R")
} else {
  library(Cairo)
  library(data.table)
  library(MotrpacBicQC)
  source(file = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/helper_scripts/deg-trait_functions.R")
  load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/figures/figS2_LDSC-comparison.RData")
}

#### figure plotting ####
change_names = F
change_names_in_plot = T
nicole_mods = T
arnold_mods = F

grDevices::cairo_pdf(filename = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/figures/fig_S2_LDSC_output.pdf"), 
                     width = 1500 / 72, height = 2000 / 72 * 18 / 17 / 114 * length(coloc_phenotypes_sub) * 2, family="Arial Unicode MS", pointsize = 16)

par(mar = c(3,3,4,3), xpd = NA)
layout(matrix(c(1,2,1,3), 2, 2), heights = c(1,0.75))

for(type_of_plot in 2:2){
  
  use_enrichment <- c(T, F)[type_of_plot]
  
  max_horiz_axis <- 2.05
  
  
  #first get metainfo
  if(reorder_vertical){
    
    if(use_enrichment){
      order_traits <- order(sapply(coloc_phenotypes_sub, function(coloc_phenotype) mean(ldsc_results_sub[ldsc_results_sub$gwas == coloc_phenotype,]$Enrichment)), decreasing = T)
    } else {
      order_traits <- order(sapply(coloc_phenotypes_sub, function(coloc_phenotype) mean(ldsc_results_sub[ldsc_results_sub$gwas == coloc_phenotype,]$Prop._h2)), decreasing = T)
    }
    
  } else {
    order_traits <- 1:length(coloc_phenotypes_sub)
  }
  
  if(partition_by_category){
    order_traits <- order_traits[order(traitwise_partitions[,2][match(coloc_phenotypes_sub[order_traits], traitwise_partitions[,1])], sort(order_traits))]
    trait_category_counts <- as.data.table(table(traitwise_partitions[,2][match(coloc_phenotypes_sub[order_traits], traitwise_partitions[,1])] ))
    trait_category_counts$cN <- cumsum(trait_category_counts$N)
  }
  
  if(use_enrichment){
    if(fix_axes){
      logprob_range <- fix_axes_bounds_enrichment
    } else {
      logprob_range <- range(ldsc_results_sub$Enrichment)
    }
  } else {
    if(fix_axes){
      logprob_range <- fix_axes_bounds_heritability
    } else {
      logprob_range <- range(ldsc_results_sub$Prop._h2)
    }
  }
  
  logprob_range <- c(logprob_range[1] - diff(logprob_range) * buffer_min_and_max / 2, logprob_range[2] + diff(logprob_range) * buffer_min_and_max / 2)
  logprob_range <- c(0,ifelse(use_enrichment, 
                              max(ldsc_results_sub$Enrichment[ldsc_results_sub$adj_p < 0.05]) * 1.025, 
                              max(ldsc_results_sub$Prop._h2[ldsc_results_sub$adj_p < 0.05])) * 1.025)
  max_prob <- diff(logprob_range)
  
  nameloc <- 0.3 + ifelse(change_names | change_names_in_plot, 0, 0.3)
  
  plot(1,1,xlim = c(0,2.1), ylim = c(0,10), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
  
  ylocs <- seq(10, 0, length.out = length(coloc_phenotypes_sub))
  
  #put boxes around categories, if we're doing that  
  if(partition_by_category){
    segments(x0 = -0.1, x1 = -0.1, y0 = 0+diff(ylocs)[1]/2, y1 = 10-diff(ylocs)[1]/2, col = "grey35")
    for(traitcat in 1:(nrow(trait_category_counts)-1)){
      segments(y0 = ylocs[trait_category_counts$cN[traitcat]] + diff(ylocs)[1]/2, y1 = ylocs[trait_category_counts$cN[traitcat]] + diff(ylocs)[1]/2, x0 = -0.1, x1 = nameloc, col = "grey35")
    }
    
    catylocs <- (ylocs[trait_category_counts$cN] + c(ylocs[1], ylocs[trait_category_counts$cN][-length(trait_category_counts$cN)])) / 2 - diff(ylocs)[1]/1.5
    catnames <- as.vector(trait_category_counts$V1)
    catnames <- sapply(catnames, function(x) strsplit(x, split = c("-"))[[1]][1])
    catnames <- sapply(catnames, function(x) strsplit(x, split = c(" "))[[1]][1])
    
    #manually fudge category locations
    catylocs[catnames == "Hair"] <- catylocs[catnames == "Hair"] + diff(ylocs)[1]/1.5
    catylocs[catnames == "Aging"] <- catylocs[catnames == "Aging"] - diff(ylocs)[1]/3
    catylocs[catnames == "Endocrine"] <- catylocs[catnames == "Endocrine"] - diff(ylocs)[1]/10
    catylocs[catnames == "Allergy"] <- catylocs[catnames == "Allergy"] - diff(ylocs)[1]*1.2
    catylocs[catnames == "Anthropometric"] <- catylocs[catnames == "Anthropometric"] - diff(ylocs)[1]*1.2
    catylocs[catnames == "Cardiometabolic"] <- catylocs[catnames == "Cardiometabolic"] - diff(ylocs)[1]*2
    for(catname in 1:length(catnames)){
      text(x = -0.1, y = catylocs[catname] + ifelse(catnames[catname] == "Endocrine", 0.1, 0), labels = catnames[catname], pos = 2, srt = 90, col = "grey35")
    }
    
  }
  
  #plot names of traits
  if(change_names_in_plot){
    coloc_phenotypes_sub_newnames <- traitname_map[match(coloc_phenotypes_sub, traitname_map[,1]),2]
  } else {
    coloc_phenotypes_sub_newnames <- coloc_phenotypes_sub
  }
  
  if(nicole_mods){
    # text(paste0(1:length(coloc_phenotypes_sub_newnames), ": ", 
    #             coloc_phenotypes_sub_newnames)[order_traits], 
    #      col = bad_boys_cols[order_traits], x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
    text(paste0(coloc_phenotypes_sub_newnames)[order_traits], 
         col = bad_boys_cols[order_traits], x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
  } else {
    text(paste0(1:length(coloc_phenotypes_sub_newnames), ": ", 
                coloc_phenotypes_sub_newnames)[order_traits], 
         col = pheno_cols, x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
  }
  
  #vert axis label
  # text(labels = "Phenotypes", x = nameloc, y = 10.45, pos = 2, cex = 2.5, xpd = NA, family="Courier", col = "grey25")
  text(labels = "Phenotypes", x = nameloc, y = 10.45, pos = 2, cex = 2.5, xpd = NA, col = "grey25")
  #horiz axis label
  text(labels = paste0("Heritability ", ifelse(use_enrichment, "Enrichment", "Proportion"), " Across Tissues"), 
       x = 1.2, y = -0.675, cex = 2.25, pos = 1, xpd = NA)
  
  
  #guiding lines for traits
  if(nicole_mods){
    segments(x0 = nameloc, x1 = max_horiz_axis, y0 = ylocs, y1 = ylocs, col = "grey80", lty = 3, lwd = 0.5)
  }
  
  #axes
  segments(x0 = nameloc, x1 = nameloc, y0 = 10.675, y1 = -0.1, lwd = 2)
  segments(x0 = nameloc, x1 = max_horiz_axis, y0 = -0.1, y1 = -0.1, lwd = 2)
  
  #horiz axis ticks and nums
  # if(use_geometric_mean){
  #   
  #   probs <- ceiling(logprob_range[1]):ceiling(logprob_range[2])
  #   segments(x0 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
  #            x1 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
  #            y0 = -0.1, y1 = -0.15, lwd = 2)
  #   text(x = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
  #        y = -0.135, pos = 1, cex = 1.25,
  #        labels =  c(sapply(probs, function(num) as.expression(bquote(10^.(num))))))
  #   
  # } else {
  #   
  #   segments(x0 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
  #            x1 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
  #            y0 = -0.1, y1 = -0.15, lwd = 2)
  #   text(x = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
  #        y = -0.135, pos = 1, cex = 1.25,
  #        labels =  c(sapply(seq(ceiling(logprob_range[1]), ceiling(logprob_range[2]), by = 1), function(num) as.expression(bquote(10^.(num))))))
  # }
  
  if(use_enrichment){
    probs <- seq(ceiling(logprob_range[1]), floor(logprob_range[2]), by = ifelse(fix_axes, 0.5, 1))  
  } else {
    probs <- seq(ceiling(logprob_range[1]), (logprob_range[2]), by = ifelse(fix_axes, 0.1, 0.1))
  }
  
  segments(x0 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc,
           x1 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc,
           y0 = -0.1, y1 = -0.25, lwd = 2)
  text(x = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc,
       y = -0.25, pos = 1, cex = 1.5,
       labels =  probs)
  
  
  # #let's get minor tick marks in there too
  # if(use_geometric_mean){
  #   minticks <- log10(as.vector(sapply(1:(length(probs)-1), function(p) seq(10^probs[p], 10^probs[p+1], length.out = 10))))
  #   minticks <- (minticks - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc
  #   minticks <- minticks[minticks < max_horiz_axis]
  #   segments(x0 = minticks, x1 =  minticks, y0 = -0.1, y1 = -0.125, lwd = 1)
  #   #hmm being weird
  # } else {
  #   minticks <- log10(as.vector(t(t(1:9)) %*% t(10^seq(0, max_prob, by = 10^round(log(max_prob))))))
  #   minticks <- minticks / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025
  #   minticks <- minticks[minticks < max_horiz_axis]
  #   segments(x0 = minticks, x1 =  minticks, y0 = -0.1, y1 = -0.125, lwd = 1)
  # }
  
  
  #dashed lines for ticks
  segments(x0 = (probs[-1] - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
           x1 = (probs[-1] - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y0 = 10, y1 = -0.075, lwd = 2, lty = 2, col = "grey75")
  segments(x0 = (1 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
           x1 = (1 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
           y0 = 10, y1 = -0.075, lwd = 3, lty = 1, col = "grey75")
  
  
  #helpful hint line about direction
  
  # if(use_geometric_mean){
  #   segments(x0 = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
  #            x1 = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y0 = 10, y1 = 10.4, lwd = 2, lty = 2, col = "grey75")
  #   text(labels = "non-DE Genes", srt = 90, x = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = 10.4, pos = 2, col = "grey75", cex = 0.75)
  #   text(labels = "DE Genes", srt = 270, x = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = 10.4, pos = 4, col = "grey75", cex = 0.75)
  #   text(labels = "Stronger Coloc\nSignal In", srt = 0, x = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = 10.3875, pos = 3, col = "grey75", cex = 0.75)
  # } else {
  #   segments(x0 = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
  #            x1 = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y0 = 10, y1 = 10.4, lwd = 2, lty = 2, col = "grey75")
  #   text(labels = "eQTL > DEG", srt = 90, x = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = 10.4, pos = 2, col = "grey75", cex = 0.75)
  #   text(labels = "DEG > eQTL", srt = 270, x = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = 10.4, pos = 4, col = "grey75", cex = 0.75)
  #   text(labels = "Stronger Coloc\nSignal In", srt = 0, x = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = 10.3875, pos = 3, col = "grey75", cex = 0.75)
  # }
  
  #plot legend for colors etc.
  lower_legend_by <- 1.6
  rect(xleft = 1.97, ybottom = 8.15 - lower_legend_by - 0.05*n_points_for_legend, 
       ytop = 10.1 - lower_legend_by, xright = 2.2, border = NA, col = "white")
  # text(labels = stringr::str_replace_all(cluster_names, "_", " "),
  #      y = seq(9.95,8.75,length.out = length(cluster_names)), x = 2, pos = 4)
  text(labels = sapply(cluster_names, function(cli) strsplit(cli, "-")[[1]][1]),
       y = seq(9.95,8.75 - lower_legend_by,length.out = length(cluster_names)) - 0.02, x = 2, pos = 4)
  points(x = rep(1.9925, length(cluster_names)), y = seq(9.955,8.755 - lower_legend_by,length.out = length(cluster_names)), col = cols$cluster, pch = 19, cex = 1.75)
  
  #plot legend for points
  lower_legend_by <- lower_legend_by + 0.3
  pt_loc_expander <- 8
  point_legend_cexes <- seq(from = 0.4, to = max_point_cex, length.out = n_points_for_legend)
  points_legend_pchs <- rep(19, n_points_for_legend)
  point_legend_cexes_log10_pvals <- round((point_legend_cexes / max_point_cex)^(1/point_cex_power) * minimum_enrichment_logPval, 2)
  # points_legend_pchs[point_legend_cexes_log10_pvals < (log10(0.05) - log10(nrow(ldsc_results_sub)))] <- 18
  points_legend_pchs[point_legend_cexes_log10_pvals < (log10(0.05))] <- 18
  
  points(y = 8.55 - lower_legend_by - cumsum(point_legend_cexes + pt_loc_expander) / sum(point_legend_cexes) * 0.05*n_points_for_legend, 
         x = rep(1.9925, n_points_for_legend), cex = point_legend_cexes, col = "grey50", pch = points_legend_pchs)
  text(labels = latex2exp::TeX("$log_{10}$(enrichment p-val)"), y = 8.6 - lower_legend_by, x = 1.975, cex = 1.1, pos = 4,  font = 2)
  # text(labels = paste0("0.05 FDR @ ", round(log10(0.05) - log10(nrow(ldsc_results_sub)),2), "\n           (Bonferroni)"), y = 8.45 , x = 2.1, cex = 0.8, pos = 4,  font = 2)
  text(labels = latex2exp::TeX(paste0("$\\alpha = 0.05$, IHW")), 
       y = 8.35 - lower_legend_by, x = 2.1, cex = 0.8, pos = 4,  font = 2)
  # text(labels = latex2exp::TeX(paste0("(Benjamini-Hochberg)")), 
  #      y = 8.35 - lower_legend_by, x = 2.07, cex = 0.8, pos = 4,  font = 2)
  # text(labels = paste0("< ", round(log10(0.05) - log10(nrow(ldsc_results_sub)),2), "\n> ", round(log10(0.05) - log10(nrow(ldsc_results_sub)),2)), 
  #      y = 8.25 , x = 2.175, cex = 1, pos = 4,  font = 2)
  text(labels = paste0("< ", round(log10(0.05),2), "\n> ", round(log10(0.05),2)), 
       y = 8.25 - lower_legend_by - 0.1 - c(0.15), x = 2.145, cex = 1, pos = 4,  font = 2)
  points(pch = c(19,18), y = c(8.41, 8.2) - lower_legend_by - 0.1 - 0.15, x = rep(2.145,2), cex = c(1.25, 1.75), 
         col = sapply(c(opacity_insig_points, opacity_sig_points), function(opcty) adjustcolor("grey50", alpha.f = opcty)))
  
  
  text(labels = point_legend_cexes_log10_pvals,y = 8.545 - lower_legend_by - cumsum(point_legend_cexes + pt_loc_expander) / sum(point_legend_cexes) * 0.05*n_points_for_legend, 
       x = rep(1.99, n_points_for_legend) + point_legend_cexes / 200, cex = 1, pos = 4, pch = 19)
  
  
  #figure out cex params
  for(cluster in cluster_names){
    
    # horizontal lines for colocalizing traits
    trait_locs <- ylocs
    if(use_enrichment){
      trait_probs <- ldsc_results_sub$Enrichment[ldsc_results_sub$cluster == cluster]
      #for bonferroni
      # trait_logPvals <- (ldsc_results_sub$logPVal_enrichment[ldsc_results_sub$cluster == cluster])[order_traits]
      #for BH adjustment
      trait_logPvals <- log10((ldsc_results_sub$adj_p[ldsc_results_sub$cluster == cluster])[order_traits])
    } else {
      trait_probs <- ldsc_results_sub$Prop._h2[ldsc_results_sub$cluster == cluster]
      trait_logPvals <- (ldsc_results_sub$logPVal_enrichment[ldsc_results_sub$cluster == cluster])[order_traits]
    }
    trait_probs <- trait_probs[order_traits]
    
    point_cex <- (-trait_logPvals / -minimum_enrichment_logPval)^point_cex_power
    point_cex <- point_cex * max_point_cex
    
    good_points <- trait_probs >= logprob_range[1] & trait_probs <= logprob_range[2]
    
    pchs <- rep(19, length(trait_probs))
    opacities <- rep(opacity_insig_points, length(trait_probs))
    
    #mark "significant" points
    # pchs[trait_logPvals < (log10(0.05) - log10(nrow(ldsc_results_sub)))] <- 18
    # opacities[trait_logPvals < (log10(0.05) - log10(nrow(ldsc_results_sub)))] <- opacity_sig_points
    pchs[trait_logPvals < log10(0.05)] <- 23
    opacities[trait_logPvals < log10(0.05)] <- opacity_sig_points
    
    point_cols <- sapply(opacities, function(opcty) grDevices::adjustcolor(cols$Tissue[match(cluster, cluster_names)], alpha.f = opcty))
    point_bgs <- point_cols
    point_cols[trait_logPvals < log10(0.05)] <- 1
    
    
    points(x = ((trait_probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc)[good_points], y = (trait_locs)[good_points], 
           pch = pchs[good_points], bg = point_bgs[good_points], col = point_cols[good_points], cex = point_cex[good_points])
    
    
  }
  
  
  
  if(arnold_mods){
    addImg(obj = arnold_bayes, x = 1.9, y = 1.6775, width = 0.5)
  }
  
}


sub_1 <- ldsc_results_cts_alt_fullyconditional[ldsc_results_cts_alt_fullyconditional$gwas %in% traits_to_keep,]
sub_2 <- ldsc_results_cts_alt[ldsc_results_cts_alt$gwas %in% traits_to_keep,]
tissues <- tissue_order[tissue_order %in% unique(gsub(x = sub_1$cluster, "-sex_homogeneous_changing", ""))]


par(mar = c(4.25,4.25,4,2) * 1.5, xpd = F)
xlims <- range(sub_1$Enrichment)
ylims <- range(sub_2$Enrichment)
plot(sub_1$Enrichment,  sub_2$Enrichment, main = latex2exp::TeX("$h^2_{SNP}$ Enrichment"),
     xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "", pch = 19, cex.main = 2.5,
     col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[gsub(x = sub_1$cluster, "-sex_homogeneous_changing", "")], 0.5),
     xlim = xlims, ylim = ylims, cex = 2)
rect(xleft = par("usr")[1], ybottom = par("usr")[3], xright = par("usr")[2], ytop = par("usr")[4], lwd = 2, xpd = NA)
text(x = mean(xlims), y = ylims[1] - diff(ylims) / 7, labels = "Enrichment (Conditional on Tissue Annotation)", cex = 1.75, pos = 1, xpd = NA)
text(x = xlims[1] - diff(xlims) / 6, y = mean(ylims), labels = "Enrichment (Unconditional on Tissue Annotation)", cex = 1.75, srt = 90, xpd = NA)

#axes
xvals <- round(seq(ceiling(xlims[1]), floor(xlims[2]), length.out = 5))
xvals <- seq(ceiling(xlims[1]), floor(xlims[2]), by = diff(xvals)[1])
segments(x0 = xvals, x1 = xvals, y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4]) / 50, lwd = 2, xpd = NA)
text(xvals, y = par("usr")[3] - diff(par("usr")[3:4]) / 40, labels = xvals, cex = 1.75, pos = 1, xpd = NA)

yvals <- round(seq(ceiling(ylims[1]), floor(ylims[2]), length.out = 5))
yvals <- seq(ceiling(ylims[1]), floor(ylims[2]), by = diff(yvals)[1])
segments(y0 = yvals, y1 = yvals, x0 = par("usr")[1], x1 = par("usr")[1] - diff(par("usr")[1:2]) / 50, lwd = 2, xpd = NA)
text(yvals, x = par("usr")[1] - diff(par("usr")[1:2]) / 50, labels = yvals, cex = 1.75, pos = 2, xpd = NA)


legend("topleft", pch = 19, col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues], 0.5), legend = tissues, cex = 1.25, pt.cex = 2.5, bty="n")
legend("bottomright", lty = 2, col = "red", legend = "1-to-1 line", bty="n", lwd = 3, cex = 1.5)
abline(0,1, col = "red", lty = 2, lwd = 4)


#second figure, prop h2

xlims <- range(sub_1$Prop._h2)
ylims <- range(sub_2$Prop._h2)
plot(sub_1$Prop._h2,  sub_2$Prop._h2, main = latex2exp::TeX("Proportion $h^2_{SNP}$"),
     xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "", pch = 19, cex.main = 2.5,
     col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[gsub(x = sub_1$cluster, "-sex_homogeneous_changing", "")], 0.5),
     xlim = xlims, ylim = ylims, cex = 2)
rect(xleft = par("usr")[1], ybottom = par("usr")[3], xright = par("usr")[2], ytop = par("usr")[4], lwd = 2, xpd = NA)
text(x = mean(xlims), y = ylims[1] - diff(ylims) / 7, labels = latex2exp::TeX("Proportion $h^2_{SNP}$ (Conditional on Tissue Annotation)"), cex = 1.75, pos = 1, xpd = NA)
text(x = xlims[1] - diff(xlims) / 4.5, y = mean(ylims), labels = latex2exp::TeX("Proportion $h^2_{SNP}$ (Unconditional on Tissue Annotation)"), cex = 1.75, srt = 90, xpd = NA)

#axes
xvals <- round(seq((xlims[1]), (xlims[2]), length.out = 5), 2)
segments(x0 = xvals, x1 = xvals, y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4]) / 50, lwd = 2, xpd = NA)
text(xvals, y = par("usr")[3] - diff(par("usr")[3:4]) / 40, labels = xvals, cex = 1.75, pos = 1, xpd = NA)

yvals <- round(seq((ylims[1]), (ylims[2]), length.out = 5), 2)
segments(y0 = yvals, y1 = yvals, x0 = par("usr")[1], x1 = par("usr")[1] - diff(par("usr")[1:2]) / 50, lwd = 2, xpd = NA)
text(yvals, x = par("usr")[1] - diff(par("usr")[1:2]) / 50, labels = yvals, cex = 1.75, pos = 2, xpd = NA)

legend("topleft", pch = 19, col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues], 0.5), legend = tissues, cex = 1.25, pt.cex = 2.5, bty="n")
legend("bottomright", lty = 2, col = "red", legend = "1-to-1 line", bty="n", lwd = 3, cex = 1.5)
abline(0,1, col = "red", lty = 2, lwd = 4)

fig_label("c)", xpd = NA, cex = 3, shrinkX = 0.75, shrinkY = 0.95)
fig_label("b)", xpd = NA, cex = 3, shrinkX = 7.5, shrinkY = 0.95)
fig_label("a)", xpd = NA, cex = 3, shrinkX = 7.75, shrinkY = 2.65)

dev.off()
