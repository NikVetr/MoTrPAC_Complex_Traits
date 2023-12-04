#### run preprocessing script ####
run_preprocessing_scripts <- F #or load the data directly
if(run_preprocessing_scripts){
  figure_id <- "S1"
  source("/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/helper_scripts/figure_set_1_preprocessing.R")
} else {
  library(Cairo)
  library(pracma)
  source(file = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/helper_scripts/deg-trait_functions.R")
  load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/figures/figS1_trait-x-tissue-correlations.RData")
}

#### figure plotting ####
vertically_oriented <- F
grDevices::cairo_pdf(filename = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/figures/figS1_high-level-overview",
                                       ifelse(vertically_oriented, "_vertical", "_horizontal"),".pdf"), 
                     width = ifelse(vertically_oriented, 775 / 72, 1300 / 72), height = ifelse(vertically_oriented, 1300/72, 600 / 72), family="Arial Unicode MS")
if(vertically_oriented){
  layout(t(t(c(rep(1,1),rep(2,1)))))
} else {
  layout(t(c(rep(1,1),rep(2,1))))  
}

if(vertically_oriented){
  par(mar = c(6,8,3,10), xpd = NA)
} else {
  par(mar = c(6,7,3,6), xpd = NA)
}

# par(mar = c(6,5,3,4), xpd = NA)
# layout(mat = t(as.matrix(c(rep(1,5), rep(2,6)))))

#plotting params
traits <- rownames(gcor_mat)
ncols <- 101
rate = 0.001
incl_h2s <- T
exp_dist_cols <- round(cumsum(c(1, dexp(1:(ncols-1), rate = rate) / min(dexp(1:(ncols-1), rate = rate)))))
heatcols <- viridis::magma(max(exp_dist_cols))[exp_dist_cols]
heatcols <- RColorBrewer::brewer.pal(11, "RdBu")[-c(1,11)]
heatcols <- rev(colorRampPalette(heatcols)(ncols))

pow <- 1

h2_cols <- rev(viridis::mako(n = ncols))
h2_cols <- sapply(1:ncols, function(coli) adjustcolor(h2_cols[coli], alpha.f = ((10+coli)/(ncols+10))^pow))
minrangeh2 <- c(min(log10(c(estimated_h2_ldsc, estimated_h2_mesc))), diff(range(log10(c(estimated_h2_ldsc, estimated_h2_mesc)))))
h2med_cols <- rev(viridis::rocket(n = ncols))
h2med_cols <- sapply(1:ncols, function(coli) adjustcolor(h2med_cols[coli], alpha.f = ((10+coli)/(ncols+10))^pow))
minrangeh2med <- c(min(log10(estimated_h2med_mesc)), diff(range(log10(estimated_h2med_mesc))))
h2medoh2_cols <- rev(viridis::rocket(n = ncols))
h2medoh2_cols <- sapply(1:ncols, function(coli) adjustcolor(h2medoh2_cols[coli], alpha.f = ((10+coli)/(ncols+10))^pow))
minrangeh2medoh2 <- c(min(log10(estimated_h2med_over_h2_mesc)), diff(range(log10(estimated_h2med_over_h2_mesc))))

# plot(1:length(heatcols), 1:length(heatcols), cex = 2, pch = 19, col = heatcols)

plot(1,1,xlim = c(0,length(traits)), ylim = c(0,length(traits)), xpd = NA,
     col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
sorted_row_inds <- order(cmdscale(1-gcor_mat, k = 1))
for(rowi in 1:length(traits)){
  text(labels = traits[sorted_row_inds[rowi]], x = length(traits) + ifelse(incl_h2s, 5, 0), y = length(traits) - rowi + 1 - 0.1, 
       col = cols$category[traitwise_partitions$Category[match(traitname_map$Tag[
         match(traits[sorted_row_inds[rowi]], traitname_map$new_Phenotype)],traitwise_partitions$Tag)]]
       , pos = 4, xpd = NA, cex = 0.45, font = 1)
  for(colj in 1:length(traits)){
    rect(yb = length(traits) - rowi + 1 - 1/2, yt = length(traits) - rowi + 1 + 1/2, 
         xl = colj - 1/2, xr = colj + 1/2, pch = 15, cex = 1, border = NA,
         col = heatcols[round((gcor_mat[sorted_row_inds[rowi], sorted_row_inds[colj]] + 1) / 2 * 100) + 1])
    
    if(gcor_mat_sig[sorted_row_inds[rowi], sorted_row_inds[colj]]){
      points(y = length(traits) - rowi + 1, x = colj, pch = 19, col = "white", cex = 0.2)
    }
    
    if(rowi == 1){
      text(labels = traits[sorted_row_inds[colj]], x = colj + 1.5, y = -0.2, 
           col = cols$category[traitwise_partitions$Category[match(traitname_map$Tag[
             match(traits[sorted_row_inds[colj]], traitname_map$new_Phenotype)],traitwise_partitions$Tag)]]
           , pos = 2, srt = 90, xpd = NA, cex = 0.45, font = 1)
    }
    
  }
  
  if(incl_h2s){
    trait_tag <- trait_categories$Tag[match(traits[sorted_row_inds[rowi]], trait_categories$new_Phenotype)]
    h2lab_cex <- 0.65
    h2lab_height <- strheight("$\\textit{h}^{2}_{SNP}-LDSC$", cex = h2lab_cex)
    h2lab_dispint <- 1.25
    h2lab_xloc <- length(traits) + 2.25 - (h2lab_height + h2lab_dispint) * 1.5 + (0:3) * (h2lab_height + h2lab_dispint)
    
    rect(yb = length(traits) - rowi + 1 - 1/2, yt = length(traits) - rowi + 1 + 1/2, 
         xl = length(traits) + 2 - 1/2, xr = length(traits) + 2 + 1/2, pch = 15, cex = 1, border = NA,
         col = h2_cols[ceiling((log10(estimated_h2_ldsc[trait_tag]) - minrangeh2[1]) / minrangeh2[2] * 100)])
    if(estimated_h2_ldsc_pval[trait_tag] < bonf_p_h2){
      points(y = length(traits) - rowi + 1, x = length(traits) + 2, pch = 19, col = "white", cex = 0.2)
    }
    if(rowi == 1){
      text(y = length(traits) + 2, x = h2lab_xloc[1], pos = 4, srt = 90,
           labels = latex2exp::TeX(paste0("$\\textit{h}^{2}_{SNP}-LDSC$")), cex = h2lab_cex)
      segments(y1 = length(traits) + 1.75, x1 = h2lab_xloc[1] + 3/4,
               y0 = length(traits) + 0.75, x0 = length(traits) + 2,
               lwd = 1, lty = 3, col = "grey50")
    }
    
    
    rect(yb = length(traits) - rowi + 1 - 1/2, yt = length(traits) - rowi + 1 + 1/2, 
         xl = length(traits) + 3 - 1/2, xr = length(traits) + 3 + 1/2, pch = 15, cex = 1, border = NA,
         col = h2_cols[ceiling((log10(estimated_h2_mesc[trait_tag]) - minrangeh2[1]) / minrangeh2[2] * 100)])
    if(estimated_h2_mesc_pval[trait_tag] < bonf_p_h2){
      points(y = length(traits) - rowi + 1, x = length(traits) + 3, pch = 19, col = "white", cex = 0.2)
    }
    if(rowi == 1){
      text(y = length(traits) + 2, x = h2lab_xloc[2], pos = 4, srt = 90,
           labels = latex2exp::TeX(paste0("$\\textit{h}^{2}_{SNP}-MESC$")), cex = h2lab_cex)
      segments(y1 = length(traits) + 1.75, x1 = h2lab_xloc[2] + 3/4,
               y0 = length(traits) + 0.75, x0 = length(traits) + 3,
               lwd = 1, lty = 3, col = "grey50")
    }
    
    rect(yb = length(traits) - rowi + 1 - 1/2, yt = length(traits) - rowi + 1 + 1/2, 
         xl = length(traits) + 4 - 1/2, xr = length(traits) + 4 + 1/2, pch = 15, cex = 1, border = NA,
         col = h2med_cols[ceiling((log10(estimated_h2med_mesc[trait_tag]) - minrangeh2med[1]) / minrangeh2med[2] * 100)])
    if(estimated_h2med_mesc_pval[trait_tag] < bonf_p_h2){
      points(y = length(traits) - rowi + 1, x = length(traits) + 4, pch = 19, col = "white", cex = 0.2)
    }
    if(rowi == 1){
      text(y = length(traits) + 2, x = h2lab_xloc[3], pos = 4, srt = 90,
           labels = latex2exp::TeX(paste0("$\\textit{h}^{2}_{mediated}$")), cex = h2lab_cex)
      segments(y1 = length(traits) + 1.75, x1 = h2lab_xloc[3] + 3/4,
               y0 = length(traits) + 0.75, x0 = length(traits) + 4,
               lwd = 1, lty = 3, col = "grey50")
    }
    
    rect(yb = length(traits) - rowi + 1 - 1/2, yt = length(traits) - rowi + 1 + 1/2, 
         xl = length(traits) + 5 - 1/2, xr = length(traits) + 5 + 1/2, pch = 15, cex = 1, border = NA,
         col = h2medoh2_cols[ceiling((log10(estimated_h2med_over_h2_mesc[trait_tag]) - minrangeh2medoh2[1]) / minrangeh2medoh2[2] * 100)])
    if(estimated_h2med_over_h2_mesc_pval[trait_tag] < bonf_p_h2){
      points(y = length(traits) - rowi + 1, x = length(traits) + 6, pch = 19, col = "white", cex = 0.2)
    }
    if(rowi == 1){
      text(y = length(traits) + 2, x = h2lab_xloc[4], pos = 4, srt = 90,
           labels = latex2exp::TeX(paste0("$\\textit{h}^{2}_{mediated}$ $\\textit{h}^{-2}_{SNP}$")), cex = h2lab_cex)
      segments(y1 = length(traits) + 1.75, x1 = h2lab_xloc[4] + 3/4,
               y0 = length(traits) + 0.75, x0 = length(traits) + 5,
               lwd = 1, lty = 3, col = "grey50")
    }
    
  }
}

#legend for correlation matrix
xl = -2.5; xr = -0.5; yb = length(traitnames) / 3; yt = length(traitnames) / 1.5
ydisp <- -3
rect(xleft = xl, xright = xr, col = heatcols, border = NA,
     yb = seq(yb, yt, length.out = length(heatcols)+1)[-(length(heatcols)+1)] + ydisp, yt= seq(yb, yt, length.out = length(heatcols)+1)[-1] + ydisp)
rect(xleft = xl, xright = xr, ybottom = yb + ydisp, ytop = yt + ydisp)
text(labels = -5:5/5, x = xl, pos = 2, y = seq(yb, yt, length.out = 11) + ydisp, cex = 1)
text(x = mean(c(xl, xr)), y = yt + ydisp - 0.25, labels = latex2exp::TeX("$\\textit{r}_g$"), pos = 3, cex = 2, font = 2)

#title
text("Genetic Correlation Matrix", x = length(traits) / 2 + 0.5, y = length(traits) + 0.5, pos = 3, cex = 2.5, font = 2)

#legend for trait categories
points(rep((xl+xr)/2 + 0.5, length(salient.categories)), 1:length(salient.categories)*1.5, pch = 15, col = cols$category[salient.categories], cex = 1.75)
text(rep((xl+xr)/2 + 0.5, length(salient.categories)), 1:length(salient.categories)*1.5, pos = 2, col = 1, 
     labels = gsub("-.*", "", salient.categories), cex = 0.75)

#legend for bonferroni correction
text(labels = latex2exp::TeX(paste0("• mark p-val < $10^{", round(log10(0.025  / choose(dim(gcor_mat)[1], 2)), 2), "}$")), 
     x = -0.5, y = yt + ydisp + nrow(gcor_mat)/5, cex = 0.75, srt = 90)
text(labels = "bonferroni", x = -1.5, y = yt + ydisp + nrow(gcor_mat)/5, cex = 0.75, srt = 90, font = 3, family = "Arial")

if(!vertically_oriented){
  fig_label(text = "a)", region = "plot", cex = 3, shrinkX = -20, shrinkY = 1.065)
}

#legend for h2 estimates
xl = length(traitnames) - 4; xr = length(traitnames) - 3; yb = -12.75; yt = -0.55
ydisp <- 0
rect(xleft = xl, xright = xr, col = h2_cols, border = NA,
     yb = seq(yb, yt, length.out = length(heatcols)+1)[-(length(heatcols)+1)] + ydisp, yt= seq(yb, yt, length.out = length(heatcols)+1)[-1] + ydisp)
rect(xleft = xl, xright = xr, ybottom = yb + ydisp, ytop = yt + ydisp)
text(labels = latex2exp::TeX(paste0("$10^{", round(seq(minrangeh2[1], minrangeh2[1] + minrangeh2[2], length.out = 5), 2), "}$")), 
     x = xl, pos = 4, y = seq(yb, yt, length.out = 5) + ydisp, cex = 0.5)
text(x = mean(c(xl, xr)), y = yb + ydisp + 0.25, labels = latex2exp::TeX("$\\textit{h}^{2}_{SNP}$"), pos = 1, cex = 0.751, font = 2)

#legend for h2med estimates
xl = length(traitnames) + 3; xr = length(traitnames) + 4; yb = -12.75; yt = -0.55
ydisp <- 0
rect(xleft = xl, xright = xr, col = h2med_cols, border = NA,
     yb = seq(yb, yt, length.out = length(heatcols)+1)[-(length(heatcols)+1)] + ydisp, yt= seq(yb, yt, length.out = length(heatcols)+1)[-1] + ydisp)
rect(xleft = xl, xright = xr, ybottom = yb + ydisp, ytop = yt + ydisp)
text(labels = latex2exp::TeX(paste0("$10^{", round(seq(minrangeh2med[1], minrangeh2med[1] + minrangeh2med[2], length.out = 5), 2), "}$")), 
     x = xl, pos = 4, y = seq(yb, yt, length.out = 5) + ydisp, cex = 0.5)
text(x = mean(c(xl, xr)), y = yb + ydisp + 0.25, labels = latex2exp::TeX("$\\textit{h}^{2}_{mediated}$"), pos = 1, cex = 0.751, font = 2)

#legend for h2med estimates
xl = length(traitnames) + 10; xr = length(traitnames) + 11; yb = -12.75; yt = -0.55
ydisp <- 0
rect(xleft = xl, xright = xr, col = h2medoh2_cols, border = NA,
     yb = seq(yb, yt, length.out = length(heatcols)+1)[-(length(heatcols)+1)] + ydisp, yt= seq(yb, yt, length.out = length(heatcols)+1)[-1] + ydisp)
rect(xleft = xl, xright = xr, ybottom = yb + ydisp, ytop = yt + ydisp)
text(labels = latex2exp::TeX(paste0("$10^{", round(seq(minrangeh2medoh2[1], minrangeh2medoh2[1] + minrangeh2medoh2[2], length.out = 5), 2), "}$")), 
     x = xl, pos = 4, y = seq(yb, yt, length.out = 5) + ydisp, cex = 0.5)
text(x = mean(c(xl, xr)) + 2.5, y = yb + ydisp + 0.25, labels = latex2exp::TeX("$\\textit{h}^{2}_{mediated}$ $\\textit{h}^{-2}_{SNP}$"), pos = 1, cex = 0.751, font = 2)

#now do the tissue transcriptome correlation matrix
col_df = data.frame(row.names=rownames(zcor))
col_df$Tissue = gsub('_.*','',rownames(col_df))
col_df$Time = sapply(rownames(col_df), function(x) unname(unlist(strsplit(x, '_')))[3])
col_df$Sex = sapply(rownames(col_df), function(x) unname(unlist(strsplit(x, '_')))[2])

colours = list(Tissue=MotrpacRatTraining6moData::TISSUE_COLORS[unique(col_df$Tissue)], 
               Time=MotrpacRatTraining6moData::GROUP_COLORS[unique(col_df$Time)],
               Sex=MotrpacRatTraining6moData::SEX_COLORS[c('male','female')])
colours$Tissue["t56-vastus-lateralis"] <- colorRampPalette(c(as.character(colours$Tissue["t56-vastus-lateralis"]), "black"))(3)[2]
nnmap <- as.data.frame(MotrpacBicQC::bic_animal_tissue_code[,4:5])
nnmap <- nnmap[nnmap$tissue_name_release != "",]
tiss_ord <- nnmap$tissue_name_release[match(rev(MotrpacRatTraining6moData::TISSUE_ORDER), nnmap$abbreviation)]
colours$Tissue <- colours$Tissue[match(tiss_ord, names(colours$Tissue))]
colours$Tissue <- colours$Tissue[!is.na(colours$Tissue)]

#figure parameters
nnmap <- as.data.frame(MotrpacBicQC::bic_animal_tissue_code[,4:5])
nice_names <- sapply(names(colours$Tissue), function(tissue) nnmap$abbreviation[match(tissue, nnmap$tissue_name_release)])

corr_thresh <- 0.3

axis.length <- 1.5
center_rescaler <- 1.25
inner_shifter <- 0.96
no_concentric_arcs <- F
opacity_concentric_arcs <- 1
opacity_nonconcentric_arcs <- 0.5
outer_shifter <- 1.04
line_weight_power <- 2
line_weight_multiplier <- 15
numbers_in_squares <- F
tissue_names_not_colors <- T
tissue_name_cex <- 0.29
across_relationships <- T
adjacent_relationships <- T
opacity_multiplier_within_tissue <- 0.35

tissues <- as.list(names(colours$Tissue))
tissues <- unlist(tissues)


tissues_to_include <- tissues
print(tissues_to_include)


if(vertically_oriented){
  par(mar = c(7,12.5,5.5,8), xpd = NA)
} else {
  par(mar = c(7,11.5,5.5,4.5), xpd = NA)
}



if(tissue_names_not_colors){
  nice_names[nice_names == "LUNG"] <- "LUNGS"
  nice_names[nice_names == "BAT"] <- "BR-AT"
  inner_shifter <- 0.9
  outer_shifter <- 1.105
}

et <- as.data.frame(do.call(rbind, sapply(1:(nrow(zcor)-1), function(ri) t(sapply((ri+1):(ncol(zcor)), function(ci) 
  c(tiss1 = col_df$Tissue[ri], tiss2 = col_df$Tissue[ci], corr = zcor[ri,ci], color = (sign(zcor[ri,ci]) + 1) / 2 + 1,
    sex1 = which(col_df$Sex[ri] == c("male", "female")), sex2 = which(col_df$Sex[ci] == c("male", "female")), 
    time1 = which(col_df$Time[ri] == paste0(c(1,2,4,8), "w")), time2 = which(col_df$Time[ci] == paste0(c(1,2,4,8), "w")))
)))))

et$weight <- abs(as.numeric(et$corr))
et <- et[et$weight > corr_thresh,]
n_tiss <- length(unique(c(et$tiss1, et$tiss2)))
aas <- list(c(pi/4 + 1E-3, 5*pi/4 - 1E-3), 
            c(3*pi/4 - 1E-3, 7*pi/4 + 1E-3), 
            rev(c(pi/4 - 1E-3, 5*pi/4 + 1E-3)), 
            rev(c(3*pi/4 + 1E-3, 7*pi/4 - 1E-3)))
rs <- 1:(n_tiss)*axis.length/(n_tiss)

# seq(axis.length/(n_tiss)/2, axis.length, by = axis.length/n_tiss)
# rs <- rs - rs[1] + diff(rs)[1]/2

et$theta1 <- match(et$sex1, c(1,2))
et$theta2 <- match(et$sex2, c(1,2))
et$r1 <- rs[match(et$tiss1, names(colours$Tissue))]
et$r2 <- rs[match(et$tiss2, names(colours$Tissue))]
et$color <- c("#2096be", "#e28a4a")[as.numeric(et$color)]

et$concentric <- et$tiss1 != et$tiss2
et$opacity <- opacity_nonconcentric_arcs
et$opacity[et$tiss1 == et$tiss2] <- et$opacity[et$tiss1 == et$tiss2] * opacity_multiplier_within_tissue
et$opacity[et$concentric] <- opacity_concentric_arcs

et_wt <- lapply(1:4, function(time) et[et$time1 == time & et$time2 == time & (et$tiss1 %in% tissues_to_include | et$tiss2 %in% tissues_to_include),])
# et_wt <- lapply(et_wt, function(et_sub) et_sub[1:3,])
if(no_concentric_arcs){
  et_wt <- lapply(1:4, function(time) et_wt[[time]][et_wt[[time]]$tiss1 != et_wt[[time]]$tiss2,])
}

#get within-tissue opacities multiplied
for(i in 1:length(et_wt)){
  et_wt[[i]]$opacity[et_wt[[i]]$tiss1 == et_wt[[i]]$tiss2] <- et_wt[[i]]$opacity[et_wt[[i]]$tiss1 == et_wt[[i]]$tiss2] * opacity_multiplier_within_tissue
}


plot(1,1,xlim = c(-2,2), ylim = c(-2,2), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
centers <- list(c(-1,1), c(1,1), c(1,-1), c(-1,-1))
centers <- lapply(centers, function(x) x * center_rescaler)
for(time in 1:4){
  if(nrow(et_wt[[time]]) == 0){next()}
  for(ri in 1:nrow(et_wt[[time]])){
    arc(t1 = aas[[time]][et_wt[[time]][ri,"theta1"]], t2 = aas[[time]][et_wt[[time]][ri,"theta2"]], r1 = et_wt[[time]][ri,"r1"], center = centers[[time]] * outer_shifter,
        r2 = et_wt[[time]][ri,"r2"], lwd = et_wt[[time]]$weight[ri]^line_weight_power * line_weight_multiplier, col = adjustcolor(et_wt[[time]]$color[ri], et_wt[[time]]$opacity[ri]), 
        random_selfing = F, clockwise_selfing = (et_wt[[time]][ri,"sex1"] == "1") * c(1,-1,1,-1)[time]*-1 + rev(c(0,1,0,1))[time], self_adjust = 0.3)
  }
  for(sex in 1:2){
    # line(t = aas[[time]][sex], 0, axis.length, col = colours$Sex[sex], lwd = 4, center = centers[[time]])
    # for(tissue in 1:n_tiss){
    #   polarp(t = aas[[time]][sex], r = rs[tissue], col = colours$Tissue[tissue], pch = 19, cex = 0.85, center = centers[[time]])
    # }
  }
}

max_outer_radii <- sapply(et_wt, function(etwt) max(etwt$r1 + etwt$r2) / 2)
week_label_locs <- sapply(1:4, function(w) polar2cart(r = max_outer_radii[w], t = (-(0:3)*pi/2 + 3*pi/4)[w]) + centers[[w]])
week_label_locs[week_label_locs == Inf] <- -center_rescaler * 1.2
week_label_locs[week_label_locs == -Inf] <- center_rescaler * 1.2
timelab_nudge <- 0.25

for(i in 1:4){
  text(c("\u2642", "\u2640", "\u2642", "\u2640")[i], col = c(colours$Sex, colours$Sex)[i], 
       x = c(0,1,0,-1)[i]*(2 - sqrt(2) + axis.length)*1.15 - c(2,0,-2.15,0)[i] * 0.06, 
       y = c(1,0,-1,0)[i]*(2 - sqrt(2) + axis.length)*1.275 - c(-1,-1,-1,-1)[i] * 0.06, 
       cex = 4.5, srt = c(45,90,225,270)[i], pos = 1)
  # shadowtext(paste0(c(1,2,4,8), "w")[i], col = colours$Time[i], x = c(-1,1,1,-1)[i]*(1 + axis.length/sqrt(2)) * 0.99 * center_rescaler, 
  #            y = c(1,1,-1,-1)[i]*(1 + axis.length/sqrt(2)) * 0.99 * center_rescaler, cex = 4,
  #            srt = c(45,-45,235,135)[i], pos = 1)
  
  #plot weeks
  shadowtext(paste0(c(1,2,4,8), "w")[i], col = colours$Time[i], x = week_label_locs[1,i] + c(-1,1,1,-1)[i]*timelab_nudge, 
             y = week_label_locs[2,i] + c(1,1,-1,-1)[i]*timelab_nudge, cex = 4,
             srt = c(45,-45,235,135)[i], pos = 1)
  
}

#two / to is where, 1: counter-clockwise, 2: clockwise, 3: opposite
et$tiw_diff <- (as.numeric(et$time2) - as.numeric(et$time1))
et$tiw <- NA
et$tiw[et$tiw_diff == c(-1, 3)[1] | et$tiw_diff == c(-1, 3)[2]] <- 1
et$tiw[et$tiw_diff == c(1, -3)[1] | et$tiw_diff == c(1, -3)[2]] <- 2
et$tiw[et$tiw_diff == c(-2, 2)[1] | et$tiw_diff == c(-2, 2)[2]] <- 3
et$tiw[et$tiw_diff == 0] <- 0

centers_bt <- list(list(c(-2,0), c(0,2)),
                   list(c(0,2), c(2,0)),
                   list(c(2,0), c(0,-2)),
                   list(c(0,-2), c(-2,0)))
centers_bt <- lapply(centers_bt, function(x1) lapply(x1, function(x2) x2 * center_rescaler * inner_shifter))

#specify angles between new axes
aas_bt <- list(list(c(1*pi/4, 7*pi/4), c(5*pi/4, 7*pi/4)),
               list(c(5*pi/4, 7*pi/4), c(3*pi/4,5*pi/4)),
               list(c(3*pi/4,5*pi/4), c(1*pi/4,3*pi/4)),
               list(c(1*pi/4,3*pi/4), c(1*pi/4, 7*pi/4)))

et$theta1 <- ifelse(et$tiw == 1, 2, 1)
et$theta2 <- ifelse(et$theta1 == 1, 2, 1)


rs_bt <- seq(sqrt(2) * center_rescaler * inner_shifter - axis.length, 
             axis.length + sqrt(2) * center_rescaler * inner_shifter, 
             length.out = n_tiss * 2 + 1)[-(n_tiss+1)]
# rs_bt <- c(rs, rs + max(rs) + diff(rs)[1]) - diff(rs)[1] + sqrt(2) * center_rescaler - axis_length
# polarp(r = rs_bt, t = rep(pi/4, length(rs_bt)), center = c(-2.5,-0), cex = , pch = 1, col = "black")



# rs_bt <- seq(sqrt(2) * center_rescaler - max(rs), 
#              sqrt(2) * center_rescaler, 
#              length.out = n_tiss)
# rs_bt <- c(rs_bt, rs_bt + max(rs_bt) - min(rs_bt) + diff(rs_bt)[1])

# rs_bt <- rs_bt - rs_bt[1] + diff(rs_bt)[1]

et$tpairs <- sapply(1:nrow(et), function(ri) paste0(sort(c(et$time1[ri], et$time2[ri])), collapse = ""))
et$close_sex <- NA
et$close_sex[et$tpairs == "12"] <- 1
et$close_sex[et$tpairs == "23"] <- 2
et$close_sex[et$tpairs == "34"] <- 1
et$close_sex[et$tpairs == "14"] <- 2

# et$r1 <- rs_bt[match(et$tiss1, names(colours$Tissue)) + (as.numeric(et$sex1) - 1) * n_tiss]
#gets index in rs_bt for close / far sex
t1_along <- cbind((n_tiss + 1) - match(et$tiss1, names(colours$Tissue)), n_tiss + match(et$tiss1, names(colours$Tissue))) 
et$r1 <- rs_bt[sapply(1:nrow(t1_along), function(ri) t1_along[ri, 2 - as.numeric(et$sex1[ri] == et$close_sex[ri])])]
# et$r2 <- rs_bt[match(et$tiss2, names(colours$Tissue)) + (as.numeric(et$sex2) - 1) * n_tiss]
t2_along <- cbind((n_tiss + 1) - match(et$tiss2, names(colours$Tissue)), n_tiss + match(et$tiss2, names(colours$Tissue)))
et$r2 <- rs_bt[sapply(1:nrow(t2_along), function(ri) t2_along[ri, 2 - as.numeric(et$sex2[ri] == et$close_sex[ri])])]
et$concentric <- et$tiss1 != et$tiss2 | et$sex1 != et$sex2
et$opacity <- opacity_nonconcentric_arcs
et$opacity[et$concentric] <- opacity_concentric_arcs


et_bt <- lapply(1:4, function(time) et[et$time1 == time & et$time2 != time & (et$tiss1 %in% tissues_to_include | et$tiss2 %in% tissues_to_include),])

if(no_concentric_arcs){
  et_bt <- lapply(1:4, function(time) et_bt[[time]][et_bt[[time]]$tiss1 != et_bt[[time]]$tiss2 | et_bt[[time]]$sex1 != et_bt[[time]]$sex2,])
}

#get between-tissue opacities multiplied
for(i in 1:length(et_bt)){
  et_bt[[i]]$opacity[et_bt[[i]]$tiss1 == et_bt[[i]]$tiss2] <- et_bt[[i]]$opacity[et_bt[[i]]$tiss1 == et_bt[[i]]$tiss2] * opacity_multiplier_within_tissue
}

#hack to fix sex mixup -- TODO find originl bug
# for(time in 1:4){
#   temp <- et_bt[[time]]$sex1
#   et_bt[[time]]$sex1 <- et_bt[[time]]$sex2
#   et_bt[[time]]$sex2 <- temp
# }

# et_bt <- lapply(et_bt, function(et_sub) et_sub[1:10,])
# et_bt <- lapply(1:4, function(et_sub) et_bt[[et_sub]])

if(adjacent_relationships){
  
  for(time in 1:4){
    if(nrow(et_bt[[time]]) == 0){next()}
    for(ri in 1:nrow(et_bt[[time]])){
      if(any(et_bt[[time]][ri,"tiw"] == c(1,2))){
        t1 = aas_bt[[time]][[et_bt[[time]][ri,"tiw"]]][et_bt[[time]][ri,"theta1"]]
        t2 = aas_bt[[time]][[et_bt[[time]][ri,"tiw"]]][et_bt[[time]][ri,"theta2"]]
        r1 = et_bt[[time]][ri,"r1"]
        r2 = et_bt[[time]][ri,"r2"]
        
        #hack to get around 1w and 8w axes switching -- should probs find more principled solution sometime
        if(all(sort(c(et_bt[[time]]$time1[ri], et_bt[[time]]$time2[ri])) == c(1,4))){
          tt <- t1; t1 <- t2; t2 <- tt
          # rt <- r1; r1 <- r2; r2 <- rt
        }
        
        arc(t1 = t1,
            t2 = t2,
            r1 = r1,
            r2 = r2,
            center = centers_bt[[time]][et_bt[[time]][ri,"tiw"]][[1]],
            lwd = et_bt[[time]]$weight[ri]^line_weight_power * line_weight_multiplier,
            col = adjustcolor(et_bt[[time]]$color[ri], et_bt[[time]]$opacity[ri])
        )
      }
    }
  }
  
}




# time = 1
# ri = 1
# arc(t1 = aas_bt[[time]][[et_bt[[time]][ri,"tiw"]]][et_bt[[time]][ri,"theta1"]], 
#     t2 = aas_bt[[time]][[et_bt[[time]][ri,"tiw"]]][et_bt[[time]][ri,"theta2"]], 
#     r1 = et_bt[[time]][ri,"r1"], 
#     r2 = et_bt[[time]][ri,"r2"],
#     center = centers_bt[[time]][et_bt[[time]][ri,"tiw"]][[1]],
#     lwd = et_bt[[time]]$weight[ri]^1.5 * 5, 
#     col = adjustcolor(et_bt[[time]]$color[ri], 0.85)
# )

if(!tissue_names_not_colors){
  for(time in 1:4){
    for(sex in 1:2){
      # line(t = aas[[time]][sex], 0, axis.length, col = colours$Sex[sex], lwd = 4, center = centers[[time]])
      for(tissue in 1:n_tiss){
        polarp(t = aas[[time]][sex], r = rs[tissue], col = colours$Sex[sex], pch = 18, cex = 2.65*axis.length, center = centers[[time]])
      }
      for(tissue in 1:n_tiss){
        polarp(t = aas[[time]][sex], r = rs[tissue], col = "white", pch = 18, cex = 1.775*axis.length, center = centers[[time]])
      }
    }
    for(sex in 1:2){
      for(tissue in 1:n_tiss){
        polarp(t = aas[[time]][sex], r = rs[tissue], col = colours$Tissue[tissue], pch = 18, cex = 1.775*axis.length, center = centers[[time]])
        
      }
    }
    
    if(numbers_in_squares){
      for(sex in 1:2){
        for(tissue in 1:n_tiss){
          cartesian_coordinates <- polar2cart(aas[[time]][sex], rs[tissue]) + centers[[time]]
          text(tissue, x = cartesian_coordinates[1], y = cartesian_coordinates[2], col = "white", 
               srt = c(45,-45,-45,45)[time], adj = 0.5, cex = 0.55, font = 2)   
        }
      }
    }
    
    
  }
} else {
  for(time in 1:4){
    for(sex in 1:2){
      for(tissue in 1:n_tiss){
        cartesian_coordinates <- polar2cart(aas[[time]][sex], rs[tissue]) + centers[[time]]
        if((nice_names[tissue] == "TESTES" & sex == 2) | (nice_names[tissue] == "OVARY" & sex == 1)){
          name_color = colours$Sex[sex]
        } else {
          name_color = colours$Sex[sex]
        }
        text(nice_names[tissue], x = cartesian_coordinates[1], y = cartesian_coordinates[2], 
             col = adjustcolor(name_color, ifelse(tissues[tissue] %in% tissues_to_include & 
                                                    !((nice_names[tissue] == "TESTES" & sex == 2) | (nice_names[tissue] == "OVARY" & sex == 1)), 1, 0.5)), 
             srt = c(-45,45,-45,45)[time], adj = 0.5, cex = 1 / strwidth(nice_names[tissue]) * tissue_name_cex, font = 2)   
      }
    }
  }
}

tissues <- names(nice_names)
if(across_relationships){
  for(time in 1:4){
    et_bt_across <- et_bt[[time]][et_bt[[time]][,"tiw"] == 3,]
    if(nrow(et_bt_across) == 0){
      next()
    }
    for(ri in 1:nrow(et_bt_across)){
      cartesian_coordinates_1 <- polar2cart(aas[[as.numeric(et_bt_across$time1[ri])]][as.numeric(et_bt_across$sex1[ri])], 
                                            rs[which(tissues == et_bt_across$tiss1[ri])]) + centers[[as.numeric(et_bt_across$time1[ri])]] * inner_shifter
      cartesian_coordinates_2 <- polar2cart(aas[[as.numeric(et_bt_across$time2[ri])]][as.numeric(et_bt_across$sex2[ri])], 
                                            rs[which(tissues == et_bt_across$tiss2[ri])]) + centers[[as.numeric(et_bt_across$time2[ri])]] * inner_shifter
      # cartesian_coordinates_1 <- (cartesian_coordinates_1 - centers[as.numeric(et_bt_across$time1[ri])][[1]]) * inner_shifter + centers[as.numeric(et_bt_across$time1[ri])][[1]]
      # cartesian_coordinates_2 <- (cartesian_coordinates_2 - centers[as.numeric(et_bt_across$time2[ri])][[1]]) * inner_shifter + centers[as.numeric(et_bt_across$time2[ri])][[1]]
      # cartesian_coordinates_1 <- cartesian_coordinates_1 * inner_shifter
      # cartesian_coordinates_2 <- cartesian_coordinates_2 * inner_shifter
      segments(x0 = cartesian_coordinates_1[1], y0 = cartesian_coordinates_1[2], x1 = cartesian_coordinates_2[1], y1 = cartesian_coordinates_2[2],
               lwd = et_bt_across$weight[ri]^line_weight_power * line_weight_multiplier,
               col = adjustcolor(et_bt_across$color[ri], et_bt_across$opacity[ri]))
    }
  }
}

xl <- -3.2; yb <- -0.425; xr <- -2.9; yt <- 0.575;
# et_wt[[time]]$weight[ri]^1 * 4
line_weights <- c(0:10/10)^line_weight_power * line_weight_multiplier
lws <- line_weights / 96 / (par("pin")[1]  / 4)
# segments(xl,yb,xl,yt,lwd = 10)
corresponding_heights <- seq(0,abs(yt - yb)/2,length.out = 11)


# line_weight_power <- 2
# line_weight_multiplier <- 10

hoffset <- 0 + c(5:0/5, 1:5/5)^line_weight_power * line_weight_multiplier / 300
voffset_rhos <- -0.1
text(labels = round(seq(-1, 1, length.out = 11), 2), y = seq(yb, yt, length.out = 11) + voffset_rhos, 
     x = xl - 0.0075 + hoffset, pos = 4, las=2, cex=0.9)
text(labels = latex2exp::TeX(paste0("|$\\textit{\\rho}$|")), y = yt - 0.03 + voffset_rhos, x = (xl) - (xr-xl)*0.25, pos = 3, cex = 2)
text(labels = paste0("≥ ", corr_thresh), y = yt + 0.0325 + voffset_rhos, x = xl + 0.225, pos = 3, cex = 1.1)
polygon(x = c(xl - lws/2, xl + rev(lws)/2),
        y = c((yb + yt) / 2 + corresponding_heights, yt - corresponding_heights) + voffset_rhos, col = "#e28a4a")
polygon(x = c(xl - lws/2, xl + rev(lws)/2),
        y = c((yb + yt) / 2 - corresponding_heights, yb + corresponding_heights) + voffset_rhos, col = "#2096be")

if(!tissue_names_not_colors){
  points(x = rep(xl, n_tiss), y = seq(yb - 0.25, yb - 0.25 - (yt - yb) / 11 * n_tiss, length.out = n_tiss), col = colours$Tissue, cex = 2, pch = 15)
  if(numbers_in_squares){
    text(1:19, x = rep(xl, n_tiss), y = seq(yb - 0.25, yb - 0.25 - (yt - yb) / 11 * n_tiss, length.out = n_tiss), col = "white", 
         adj = 0.5, cex = 0.55, font = 2)   
  }
  text(x = rep(xl, n_tiss), y = seq(yb - 0.25, yb - 0.25 - (yt - yb) / 11 * n_tiss, length.out = n_tiss), labels = nice_names, pos = 4, cex = 0.7)
}
# addImg(png::readPNG("~/Pictures/deathrats1.png"), 0, 0, width = 1)

fig_label(text = "a)", region = "plot", cex = 3, shrinkX = 1.7, shrinkY = 4.25)
if(vertically_oriented){
  fig_label(text = "b)", region = "plot", cex = 3, shrinkX = 1.7, shrinkY = 1.32)
}

dev.off()
