#### run preprocessing script ####
source("/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/helper_scripts/figure_set_3_preprocessing.R")
source("/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/helper_scripts/figure_set_4_preprocessing.R")

#### figure plotting ####

load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/deg_sigtwas_proportion.txt")
subset_to_traits <- T
cairo_pdf(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/figures/fig5_prop_pos_hits_bayesian_scatterplot", 
                 ifelse(subset_to_traits, "_sub.pdf", "")), 
          width = 2000 / 72 * length(twas_with_hits) / 2 / 73, height = 700 / 72 * 2, family="Arial Unicode MS")
par(mfrow = c(2,1), xpd = NA, 
    mar = c(19,
            9,
            5,
            7 + ifelse(subset_to_traits, 2, 0)))
lwd <- 3
tissue_name_cex <- 1.2


for(subset_i in 1:2){
  if(subset_to_traits){
    trait_subset <- colnames(n_deg_sigtwas_intersect)[
      traitwise_partitions$Category[match(colnames(n_deg_sigtwas_intersect), traitwise_partitions$Tag)] %in% 
        ifelse2(subset_to_traits, c("Cardiometabolic", "Psychiatric-neurologic"), salient.categories)]
    
    #alternatively
    if(subset_i == 1){
      trait_subset <- names(trait_means)[trait_means > 0.5]  
    } else {
      trait_subset <- names(trait_means)[trait_means < 0.5]
    }
    
    traits_to_plot <- trait_subset
    
  } else {
    traits_to_plot <- twas_with_hits
  }
  categories_represented <- unique(traitwise_partitions$Category[match(traits_to_plot, traitwise_partitions$Tag)])
  categories_represented <- gsub("\\ .*", "", gsub("\\-.*", "", categories_represented))
  category_colors <- c(category_colors, setNames(category_colors[sapply(categories_represented, function(cati) grep(cati, names(category_colors))[1])], 
                                                 categories_represented)[!(categories_represented %in% names(category_colors))])
  
  cols = list(Tissue=MotrpacRatTraining6moData::TISSUE_COLORS[names(motrpac_gtex_map)], 
              Time=MotrpacRatTraining6moData::GROUP_COLORS[paste0(c(1,2,4,8), "w")],
              Sex=MotrpacRatTraining6moData::SEX_COLORS[c('male','female')])
  cols$Tissue<- cols$Tissue[order(match(MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[names(cols$Tissue)], MotrpacRatTraining6moData::TISSUE_ORDER))]
  
  tissue_names <- sapply(strsplit(names(cols$Tissue), "-"), function(x) paste0(x[ifelse(length(x) == 2, c(2), list(2:3))[[1]]], collapse = " "))
  
  #other plotting params
  trait_category_legend_below <- T
  
  #do the plotting
  plot(100,100, xlim = c(3,length(traits_to_plot)), ylim = range(deg_sigtwas_proportion_posterior_mean[,traits_to_plot]), xpd=NA, 
       ylab = "", xlab = "", xaxt = "n", 
       yaxt = "n", bty="n", cex.lab = 2, cex.axis = 1.25)
  ry <- diff(par("usr")[3:4])
  
  text(label = "Posterior Mean Difference (logit-scale)\nin Positive Effects on GWAS Trait", x = par("usr")[1] - diff(par("usr")[1:2])/15, 
       y = mean(par("usr")[3:4]), srt = 90, pos = 3,
       cex = 1.75)
  
  for(sex in c("male", "female")[1]){
    
    # mean_freq <- apply(sapply(traits_to_plot, function(trait) as.numeric(deg_sigtwas_proportion[, trait,"8w","p",sex])), 2, weighted.mean, na.rm = T)
    
    mean_freq <- trait_means[traits_to_plot]
    
    order_traits_to_plot <- order(mean_freq, decreasing = T)
    
    #plot faded positive and negative regions
    yseq <- round(seq(min(deg_sigtwas_proportion_posterior_mean[,traits_to_plot]) - 0.0025, max(deg_sigtwas_proportion_posterior_mean[,traits_to_plot]) + 0.0025, length.out = 9), 2)
    rect(xl = 1, xr = length(traits_to_plot) + 2, yb = 0.5, ytop = max(yseq),
         col = grDevices::adjustcolor("red", 0.1), border = NA)
    rect(xl = 1, xr = length(traits_to_plot) + 2, yb = min(yseq), ytop = 0.5,
         col = grDevices::adjustcolor("blue", 0.1), border = NA)
    
    #vertical axis
    segments(x0 = 1 - length(traits_to_plot) * 0.005, x1 = 1 - 3 * 0.04, y0 = yseq, y1 = yseq, lwd = 2)
    segments(x0 = 1 - 3 * 0.02, x1 = 1 - 3 * 0.02, y0 = min(yseq), y1 = max(yseq), lwd = 2)
    text(x = 1 - 3 * 0.05, y = yseq, labels = round(logit(yseq), 2), pos = 2, cex = 1.5)
    
    #horizontal axis
    text(1:length(traits_to_plot) + 1.5, min(yseq) - ry / 10, srt = 45, pos = 2,
         labels = gsub("raneal", "ranial", paste0(ifelse(trait_sigs[traits_to_plot][order_traits_to_plot], "♦ ", ""), 
                                                  trait_categories$new_Phenotype[match(traits_to_plot[order_traits_to_plot], trait_categories$Tag)])), 
         cex = 1.125, col = 1)
    segments(1:length(traits_to_plot) + 1, min(yseq) - ry / 15, 1:length(traits_to_plot) + 1, max(yseq), col = adjustcolor(1, 0.2), lty = 2)
    
    #plot points
    points(1:length(traits_to_plot)+1, y = sort(mean_freq[traits_to_plot], decreasing = T), pch = "*", cex = 3)
    
    for(trait_i in (1:length(traits_to_plot))){
      
      #plot category blocks
      rect(xleft = which(order_traits_to_plot == trait_i) + 1/2, 
           xright = which(order_traits_to_plot == trait_i) + 3/2,
           ybottom = min(yseq)- ry / 30, 
           ytop = min(yseq)-ry / 14,
           col = category_colors[trait_categories$Category[match(traits_to_plot[trait_i], trait_categories$Tag)]])
      text(labels = substr(trait_categories$Category[match(traits_to_plot[trait_i], trait_categories$Tag)], 1, 2), 
           x = which(order_traits_to_plot == trait_i) + 1,
           y = min(yseq)-ry / 19.5, 
           col = "white", cex = 1)
      
      
      d <- unlist(deg_sigtwas_proportion_posterior_mean[,traits_to_plot[trait_i]])
      d_sig <- unlist(cell_sigs[,traits_to_plot[trait_i]])
      tissue_abbr_rev <- MotrpacRatTraining6moData::TISSUE_ABBREV_TO_CODE
      names(d) <- tissue_abbr_rev[rownames(deg_sigtwas_proportion_posterior_mean)]
      
      dn <- as.numeric(deg_sigtwas_proportion[, traits_to_plot[trait_i],"8w","n",sex])
      names(dn) <- rownames(deg_sigtwas_proportion)
      dn <- dn[names(d)]
      
      #plot white points first
      points(x = rep(which(order_traits_to_plot == trait_i) + 1, length(d)), y = d, pch = 19 - d_sig, col = "white", cex = dn^0.25)
      #and then the actual points
      points(x = rep(which(order_traits_to_plot == trait_i) + 1, length(d)), y = d, pch = 19 - d_sig, 
             col = adjustcolor(cols$Tissue[names(d)], 0.5), cex = dn^0.25)
      #and then the significant points
      points(x = rep(which(order_traits_to_plot == trait_i) + 1, length(d)), y = d, pch = 5, 
             col = 1, cex = dn^0.25*d_sig*0.8)
      
    }
    
    
    #legend for tissues
    tisscols <- cols$Tissue[!(names(cols$Tissue) %in% c("t64-ovaries", "t63-testes"))]
    points(x = rep(length(traits_to_plot) + 3, length(tisscols)), 
           y = seq(mean(yseq)-ry / 15, max(yseq), length.out = length(tisscols)), 
           col = adjustcolor(tisscols, 0.5),
           pch = 19, cex = 2)
    text(x = rep(length(traits_to_plot) + 3.25, length(tisscols)), 
         y = seq(mean(yseq)-ry / 15, max(yseq), length.out = length(tisscols)), 
         pos = 4, pch = 19, cex = 1,
         labels = MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[names(tisscols)])
    
    #legend for tissue size
    n_pts <- 5
    pt_size_logs <- seq(1, log(max(as.numeric(deg_sigtwas_proportion[,,"8w","n",sex]), 
                                   na.rm = T)) / log(2), length.out = n_pts)
    pt_size_legend <- round(2^pt_size_logs)
    pty <-  par("usr")[3] + cumsum(pt_size_legend^0.25/800) * ry * 15 + 
      cumsum(rep(ry / 200, n_pts)) - ry / 6
    text(x = length(traits_to_plot) + 2.25, 
         y = max(pty) + diff(pty)[4], 
         pos = 4, pch = 19, cex = 1.1,
         labels = "Sample Size")
    points(x = rep(length(traits_to_plot) + 3.25, n_pts), 
           y = pty, 
           col = adjustcolor(1, 0.5),
           pch = 19, cex = pt_size_legend^0.25)
    text(x = rep(length(traits_to_plot) + 3.25, n_pts) + pt_size_legend^0.25/10, 
         y = pty, 
         pos = 4, pch = 19, cex = 1,
         labels = pt_size_legend)
    
    #legend for stars
    points(x = length(traits_to_plot) + 3, 
           y = par("usr")[4]-ry / 1.525, 
           col = "black", pch = "*", cex = 3)
    text(x = length(traits_to_plot) + 3.25, cex = 1.1,
         y = par("usr")[4]-ry / 1.5, 
         pos = 4, labels = "Posterior\nMean")
    
    
    #legend for diamonds
    text(x = length(traits_to_plot) * 1.055, cex = 1,
         y = par("usr")[3] + ry / 4.9, 
         pos = 4, labels = paste0("♦ indicates\nPr(enr. or dep.)\n> ", 1-alpha_post))
    
    
    #legend for categories
    if(trait_category_legend_below){
      yadj <- 0
      yloc_factor <- ifelse(subset_i == 1, 1.8, 2)
      for(i in 1:length(categories_represented)){
        rect(xleft = 0 + i + cumsum(c(0,strwidth(categories_represented, cex = 1.25, units = "user") + 1))[i],
             xright = 1 + i + cumsum(c(0,strwidth(categories_represented, cex = 1.25, units = "user") + 1))[i],
             ybottom = min(yseq) - ry / yloc_factor - ry / 14,
             ytop = min(yseq) - ry / yloc_factor - ry / 30,
             col = category_colors[categories_represented[i]], border = 1)
        text(labels = categories_represented[i], pos = 4, 
             y = min(yseq) - ry / yloc_factor - ry / 18, 
             cex = 1.25, 
             x = 0.85 + i + cumsum(c(0,strwidth(categories_represented, cex = 1.25, units = "user") + 1))[i])
        text(labels = substr(categories_represented[i], 1, 2), 
             x = 0.5 + i + cumsum(c(0,strwidth(categories_represented, cex = 1.25, units = "user") + 1))[i],
             y = min(yseq) - ry / yloc_factor - ry / 19.5, 
             col = "white", cex = 1)
        #draw circle?
        ends_y <- c(min(yseq) - ry / yloc_factor - ry / 18, min(yseq) - ry / 18)
        arc(t1 = 3*pi/2, t2 = pi/2+1E-6, r1 = diff(range(ends_y))/2, r2 = diff(range(ends_y))/2, 
            center = c(0.5, mean(ends_y)), lwd = 3, res = 50, adjx = 210 + ifelse(subset_i == 1, 50, 70))
        points(0.5, ends_y[2], pch = -9658, cex = 3)
      }
    } else {
      x_adj2 <- x_adj - 2.5
      y_adj2 <- y_adj - 2.5
      for(i in 1:length(categories_represented)){
        rect(xleft = 0 + i,
             xright = 1 + i,
             ybottom = -10,
             ytop = -9,
             col = category_colors[categories_represented[i]], border = 1)
        text(labels = categories_represented[i], pos = 4, y = i + y_adj2 - 7, x = ncol(table_to_use) - 1/2 + x_adj2 + 1/3)
      }
    }
    
  }
  
  title_text <- paste0("Traits \\textbf{", ifelse(subset_i == 1, "More", "Less"), 
                       "} Similar to Transcriptomic Adaptation to Exercise")
  text(latex2exp::TeX(title_text), 
       x = (par("usr")[1]), y = par("usr")[4] + ry / 15, 
       cex = 2.25, pos = 4)
  rect(xleft = mean(par("usr")[1]) + strwidth(latex2exp::TeX("Traits  "), cex = 2.25),
       xright = mean(par("usr")[1]) + strwidth(latex2exp::TeX(paste0("Traits  \\textbf{", 
                                                                     ifelse(subset_i == 1, "More", "Less"), "}")), 
                                               cex = 2.25),
       ybottom = par("usr")[4] + ry / 20,
       ytop = par("usr")[4] + ry / 20 + 
         strheight(latex2exp::TeX(paste0("\\textbf{", ifelse(subset_i == 1, "More", "Less"), "}")), cex = 2.8),
       col = "white", border = "white")
  text(latex2exp::TeX(paste0("\\textbf{", ifelse(subset_i == 1, "More", "Less"), "}")), 
       x = (par("usr")[1]) + strwidth("Traits ", cex = 2.25), 
       y = par("usr")[4] + ry / 15 + (
         strheight(latex2exp::TeX(title_text), cex = 2.25) - 
           strheight(latex2exp::TeX(paste0("\\textbf{", ifelse(subset_i == 1, "More", "Less"), "}")), cex = 2.25)
       ) / 2, 
       cex = 2.25, pos = 4, col = adjustcolor(c("red", "blue")[subset_i], "0.75"))
  
  
}
dev.off()
