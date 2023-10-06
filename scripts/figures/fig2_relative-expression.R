#### run preprocessing script ####
source("/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/helper_scripts/figure_set_1_preprocessing.R")
source(file = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/analyses/analysis_GREx_RelativeEffectSize.R")

#### figure plotting ####

# actual plotting
grDevices::cairo_pdf(filename = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/figures/fig2_relative-expression.pdf"), 
                     width = 1350 / 72, height = 1125 / 72 * 1.5 * 2 / 6, family="Arial Unicode MS", pointsize = 20.5)

layout_mat <- kronecker(rbind(
  c(1,1,2,3),
  c(1,1,4,5),
  c(6,13,9,11),
  c(7,8,10,12)
), matrix(1,2,2))
layout_mat[5:6, 1:3] <- layout_mat[5:6, 3:1]
layout_mat <- layout_mat + 1
layout_mat <- rbind(matrix(1, nrow = 4, ncol = 8), layout_mat)
layout_mat <- layout_mat[9:12,]
layout_mat <- layout_mat - min(layout_mat) + 1

layout(layout_mat, heights = c(1,1,1,1))


# ~~~~~~~~~~~~~~~~~~~~~~~~
#relative expression plots
# ~~~~~~~~~~~~~~~~~~~~~~~~

par(mar = c(3,3,1,3) + 1)
#logFCs
xr <- c(-2,2)
x_logfc <- seq(from = xr[1], to = xr[2], length.out = 256)
logfc_dens <- sapply(deg_eqtl_list, function(del)
  density((del$logFC), na.rm = T, from = xr[1], to = xr[2], n = 256)$y)
plot(NA,NA, xlim = xr, ylim = c(0, quantile(logfc_dens, 1)), xlab = "", ylab = "",
     xaxt = "n", xpd = NA)
xvals <- seq(xr[1], xr[2], length.out = 5)
segments(y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4])/50, x0 = xvals, x1 = xvals, xpd = NA)
text(y = par("usr")[3] - diff(par("usr")[3:4])/50, pos = 1, x = xvals, labels = round((xvals), 2), xpd = NA)
text(x = par("usr")[1] - diff(par("usr")[1:2])/4.5, y = mean(par("usr")[3:4]), labels = "", srt = 90, xpd = NA)
text(x = par("usr")[1] - diff(par("usr")[1:2])/4.5, y = mean(par("usr")[3:4]), labels = "Density", srt = 90, xpd = NA)
text(x = mean(par("usr")[1:2]), y = par("usr")[3] - diff(par("usr")[3:4])/4.5, labels = latex2exp::TeX("Exercise-induced $log_2FC$"), srt = 0, xpd = NA)

for(tissue in colnames(logfc_dens)){
  polygon(c(x_logfc, rev(x_logfc)), c(logfc_dens[,tissue], rep(0,length(x_logfc))), col = adjustcolor(cols$Tissue[tissue], 0.5))  
}

#surrounding text
# text(x = -6, y = -6.5, labels = "(", cex = 16, xpd = NA, family = "Helvetica Narrow", srt = 0, col = 1)
par(xpd = NA)
arc(t1 = pi/2, t2 = 3*pi/2, r1 = 3, r2 = 3, center = c(-6.1, -6.875), adjx = 0.2, lwd = 3)
arc(t1 = pi/2-1E-5, t2 = 3*pi/2, r1 = 3, r2 = 3, center = c(5.6, -6.875), adjx = 0.2, lwd = 3)
text(x = 6.25, -3.75, labels = "1/2", cex = 2)
segments(x0 = -5.5, x1 = 5, y0 = -2.75, y1 = -2.75, lwd = 5, xpd = NA)
segments(x0 = 6, x1 = 7, y0 = -2.5, y1 = -2.5, lwd = 10, xpd = NA, col = 2)
points(7.25, -2.5, pch = -9658, cex = 5, xpd = NA, col = 2)
par(xpd = T)

fig_label(text = "a)", region = "plot", cex = 3, shrinkX = 2, shrinkY = 1.065)
# fig_label(text = "b)", region = "plot", cex = 3, shrinkX = 3.55, shrinkY = -0.5)
# fig_label(text = "c)", region = "plot", cex = 3, shrinkX = -4.75, shrinkY = -0.45)
fig_label(text = "b)", region = "plot", cex = 3, shrinkX = -3.75, shrinkY = 1.065)
# fig_label(text = "e)", region = "plot", cex = 3, shrinkX = -4.75, shrinkY = -0.45)


#snp heritability
par(mar = c(3,3.5,1,3.5)+1.5)
if(!exists("gcta_output")){
  load(file = paste0(gcta_directory, "gcta_output_GTEx_allTissues_list_IHW.RData"))
  load(paste0(gcta_directory, "gcta_output_GTEx_allTissues.RData"))
}
n_h2s <- sum(sapply(gcta_output, function(i) nrow(i)))
gcta_output_h2s <- lapply(setNames((names(gcta_output)), (names(gcta_output))), function(tissue) gcta_output[[tissue]]$h2)
if(!exists("h2_freqs")){
  h2_freqs <- lapply(setNames(rev(names(gcta_output)), rev(names(gcta_output))), function(tissue){
    out <- hist(do.call(c, gcta_output_h2s[1:match(tissue, names(gcta_output_h2s))]), breaks = 0:5/5, plot = F)
    out$counts <- out$counts / n_h2s
    out
  })
}
tissues <- rev(names(gcta_output))
plot(h2_freqs[[tissues[1]]], col = adjustcolor(cols$Tissue[tissues[1]], 0.5), xlab = "", ylab = "",
     main = "")
for(tissue in tissues[-1]){
  plot(h2_freqs[[tissue]], col = adjustcolor(cols$Tissue[tissue], 0.5), add = T)
}
text(x = par("usr")[1] - diff(par("usr")[1:2])/4, y = mean(par("usr")[3:4]), labels = "Density", srt = 90, xpd = NA)
text(x = mean(par("usr")[1:2]), y = par("usr")[3] - diff(par("usr")[3:4])/4.5, labels = latex2exp::TeX("Estimated $h^2_{SNP}$"), srt = 0, xpd = NA)
box("plot")


#invgamma hyperpriors
x_var <- 1:500/1000
load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/GREx_sds_expression")
load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/GREx_invgamma_estimate")
var_dens <- sapply(invgamma_estimates, function(invgamma_estimate)
  dinvgamma(x_var, shape = invgamma_estimate[1], scale = invgamma_estimate[2]))
colnames(var_dens) <- names(invgamma_estimates)
plot(NA,NA, xlim = range(x_var), ylim = c(0, quantile(var_dens, 1)), xlab = "", ylab = "")
for(tissue in colnames(var_dens)){
  polygon(c(x_var, rev(x_var)), c(var_dens[,tissue], rep(0,length(x_var))), col = adjustcolor(cols$Tissue[tissue], 0.5))  
}
text(x = par("usr")[1] - diff(par("usr")[1:2])/4, y = mean(par("usr")[3:4]), labels = "Density", srt = 90, xpd = NA)
text(x = par("usr")[1] - diff(par("usr")[1:2])/2, y = mean(par("usr")[3:4]), labels = "•", srt = 90, xpd = NA, cex = 4)
text(x = mean(par("usr")[1:2]), y = par("usr")[3] - diff(par("usr")[3:4])/4.5, labels = latex2exp::TeX("Expression Variance"), srt = 0, xpd = NA)

#quantile plots
par(mar = c(4,5,2,0), xpd = NA)
nnmap <- as.data.frame(MotrpacBicQC::bic_animal_tissue_code[,4:5])
qs2use <- round(invlogit(seq(-9,9, length.out = 75)), 4)


for(sex_i in c("male", "female")){
  for(type in 3:4){ #the two sexhomo categories, which subsets to the 8w_M1/-1_F1/-1 nodes
    
    EZ_PZ <- relative_expression_data[[sex_i]][[type]]
    EZ_PZ <- lapply(EZ_PZ, function(x) {
      out <- x[match(paste0(sprintf("%.2f", qs2use*100), "%"), paste0(sprintf("%.2f", as.numeric(gsub("%", "", rownames(x)))), "%")),]
      out[is.na(out)] <- 0
      out
    })
    
    ti = "8w"
    EZ_PZ[[ti]] <- EZ_PZ[[ti]][,which(!apply(apply(EZ_PZ[[ti]], 2, is.na), 2, any))]
    f_p <- 0.4
    f_x <- 2
    ylims = c(min(sort(EZ_PZ[[ti]])[sort(EZ_PZ[[ti]]) != -Inf]),max(sort(EZ_PZ[[ti]],T)[sort(EZ_PZ[[ti]],T) != Inf]))
    ylims <- squish_middle_x(ylims, f_x)
    
    plot(100,100,xlim = c(0,1.3), ylim = ylims, xpd=NA, ylab = "", xlab = "", xaxt = "n", yaxt = "n", bty="n", cex.lab = 1.25, cex.axis = 1.25)
    text(x = par("usr")[1] - diff(par("usr")[1:2])/6, y = mean(par("usr")[3:4]), labels = "Standardized Effect Size (SD)", srt = 90, xpd = NA)
    
    text("Quantile", x = 0.5, y = ylims[1] - diff(ylims)/5, pos = 1, cex = 1)
    if(ti == "2w"){
      text(latex2exp::TeX(paste0("Ratio of Exercise DE to \\sqrt{\\textbf{", ifelse(type == 1, "Phenotypic", "Genetic"), "} Variance in $log_2$(Gene Expression)}")), 
           x = 1.35, y = ylims[2] + diff(ylims) / 15, pos = 3, cex = 3)
    }
    if(ti == "1w"){
      fig_label(ifelse(type == 1, "a)", "b)"), region = "figure", pos = "topleft", cex = 3.5)
    }
    
    xylocs_tissue_names <- cbind(rep(1.05, ncol(EZ_PZ[[ti]])), redistribute(as.numeric(squish_middle_x(tail(EZ_PZ[[ti]], 1), f_x)), diff(ylims) / 20))
    rownames(xylocs_tissue_names) <- colnames(EZ_PZ[[ti]]); colnames(xylocs_tissue_names) <- c("x", "y")
    
    #horiz axis
    segments(x0 = 0, y0 = ylims[1] - diff(ylims) / 100, x1 = 1, y1 = ylims[1] - diff(ylims) / 100, lwd = 2)
    segments(x0 = 0:10/10, y0 = ylims[1] - diff(ylims) / 100, x1 = 0:10/10, y1 = ylims[1] - diff(ylims) / 50, lwd = 2, xpd = NA)
    horiz_axis_labels <- round(unsquish_middle_p(0:10/10, f_p), 3)
    horiz_axis_labels[1] <- 0; horiz_axis_labels[length(horiz_axis_labels)] <- 1;
    text(labels = horiz_axis_labels, x = 0:10/10 + 0.035, y = rep(ylims[1] - diff(ylims) / 25, 10), pos = 2, xpd = NA, srt = 90, cex = 0.75)
    segments(x0 = 0.5, y0 = ylims[1], x1 = 0.5, y1 = ylims[2], lwd = 2, lty = 2, col = "grey50")
    segments(x0 = 0:10/10, y0 = ylims[1], x1 = 0:10/10, y1 = ylims[2], lwd = 1, lty = 3, col = "grey75")
    
    #vert axis
    segments(x0 = -1/100, y0 = ylims[1] + diff(ylims)/100, x1 = -1/100, y1 = ylims[2], lwd = 2)
    segments(x0 = -1/100, y0 = seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), x1 = -1/50,
             y1 = seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), lwd = 2)
    text(labels = round(unsquish_middle_x(seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), f_x), 2),
         x = -1/50, y = seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), pos = 2, xpd = NA, cex = 0.75)
    segments(x0 = 0, y0 = 0, x1 = 1, y1 = 0, lwd = 2, lty = 2, col = "grey50")
    segments(x0 = 0, y0 = seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), x1 = 1, 
             y1 = seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), lwd = 1, lty = 3, col = "grey75")
    
    
    for(tissue in colnames(EZ_PZ[[ti]])){
      if(all(abs(EZ_PZ[[ti]][,tissue]) < 1E-5)){next()}
      lines(squish_middle_p(qs2use, f_p), squish_middle_x(EZ_PZ[[ti]][,tissue], f_x), col = cols$Tissue[tissue])
      text(nnmap$abbreviation[match(tissue, nnmap$tissue_name_release)], x = xylocs_tissue_names[tissue,"x"], y = xylocs_tissue_names[tissue,"y"], pos = 4, xpd = NA,
           col = cols$Tissue[tissue], cex = 0.75)
      segments(x0 = squish_middle_p(tail(qs2use, 1), f_p), x1 = xylocs_tissue_names[tissue,"x"]+0.01,
               y0 = squish_middle_x(tail(EZ_PZ[[ti]][,tissue], 1), f_x), y1 = xylocs_tissue_names[tissue,"y"],
               lty = 3, col = cols$Tissue[tissue])
    }
    text(x = -0.01, y = ylims[2] - diff(ylims) / 8, labels = ti, 
         cex = 3, col = cols$Time[ti], pos = 4)
    text(x = 0.275, y = ylims[2] - diff(ylims) / 8, labels = ",", 
         cex = 3, col = cols$Time[ti], pos = 4)
    text(labels = c(female = "\u2640", male = "\u2642")[sex_i], x = 0.35, y = ylims[2] - diff(ylims) / 8, pos = 4, cex = 3, 
         col = MotrpacRatTraining6moData::SEX_COLORS[sex_i], family = "Arial Unicode MS")
    
    if(sex_i == "male"){
      text(latex2exp::TeX(paste0("Ratio of Exercise DE to \\sqrt{", 
                                 ifelse(type == 3, "Phenotypic", "Genetic"), " Variance in $log_2$(Gene Expression)}")), 
           x = 1.5, y = ylims[2] + diff(ylims) / 50, pos = 3, cex = 1)
    }
  }
  
  
}

dev.off()
