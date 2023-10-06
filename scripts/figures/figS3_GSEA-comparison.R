#### run preprocessing script ####
source("/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/helper_scripts/figure_set_3_preprocessing.R")
source("/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/analyses/analysis_freq-GSEA.R")

#### figure plotting ####
#snag category colors
{
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

cairo_pdf(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/figures/figS3_intersect-enrichment_Bayes_vs_freq.pdf"), 
          width = 2200 / 72, 
          height = 2000 / 72, 
          family="Arial Unicode MS", pointsize = 28.5)
par(xpd = NA, mar = c(5,2,4,1), pty="s")
layout(rbind(
  c(1,1,2,2),
  c(1,1,2,2),
  c(3,3,4,6),
  c(3,3,5,6)
), heights = c(1,1,1))

#pairs
trait_x_tissue_ps <- setNames(gsea_output$pval, gsea_output$pathway_name)
# trait_x_tissue_adj_ps <- -log10(setNames(gsea_output$adj.pval, gsea_output$pathway_name))
trait_x_tissue_adj_ps <- -log10(setNames(gsea_output$pval * length(gsea_output$pval), gsea_output$pathway_name))

#from locator()... tried a programmatic solution but not enough time
# temp <- locator(); cat(paste0("piks <- cbind(\n\tx = c(", paste0(temp$x, collapse = ", "), "),\n", "\ty = c(", paste0(temp$y, collapse = ", "), ")\n)"))
#for BH-adjusted values
piks <- cbind(
  x = c(-7.04596421988288, -6.0123057215065, -2.19572049673214, -2.43425707328054, -1.40059857490415, 1.26305986321962, 1.81964520849922, 9.25403517759093, 9.21427908149953),
  y = c(-3.18977144364532, -3.24332748037729, -1.70805442739391, -1.40457021924603, -0.0121132642145909, -0.137077349922541, -0.529821619290384, -2.86843522325345)
)
# pt.pos = c(3,1,3,1,1,3,1), pt.cex = 1.75,
# prop_label_push = c(4,4,4,4,4,9,4)/100

xy <- cbind(logit(celltotalbias), trait_x_tissue_adj_ps[names(celltotalbias)])
dists <- as.matrix(dist(rbind(piks, xy)))[1:nrow(piks), nrow(piks)+1:nrow(xy)]
pts.to.label <- c(apply(dists, 1, which.min), 16, 14)

#plot trait x tissue
plot.scatter(x = logit(celltotalbias), 
             y = trait_x_tissue_adj_ps[names(celltotalbias)],
             cex = 1.5, pch = 19, 
             xlabel = latex2exp::TeX("logit(difference term posterior mass > 0)"),
             ylabel = latex2exp::TeX("-log$_{10}$(adjusted p-value)"),
             axis.lab.dev = 0.0875, lab.cex = 1.5,
             main.title = latex2exp::TeX("pairwise (\\textit{\\textbf{trait}} x \\textit{\\textbf{tissue}}) enrichment comparison"),
             pt.labels = names(celltotalbias),
             # pt.labels = NULL,
             pts.to.label = pts.to.label,
             pt.pos = c(3,1,2,3,2,
                        3,2,1,3,3,
                        4,1), 
             pt.cex = 1.75,
             prop_label_push = c(6,4,4,6,4,
                                 6,4,8,4,13,
                                 6,8,10)/100,
             pt.col = adjustcolor(TISSUE_COLORS[gsub("\\..*", "", names(celltotalbias))], 0.5))
# text(xy[pts.to.label,1], xy[pts.to.label,2], labels = 1:length(pts.to.label), col = 4, cex = 2)




#trait-wise aggregation, HMP
trait_adj_ps <- -log10(trait_ps / (trait_thresh / 0.05))
# temp <- locator(); cat(paste0("piks <- cbind(\n\tx = c(", paste0(temp$x, collapse = ", "), "),\n", "\ty = c(", paste0(temp$y, collapse = ", "), ")\n)"))
piks <- cbind(
  x = c(-7.35450884554205, -7.32513250990435, -6.26758442694711, -3.47683254136552, -2.8893058286115, 2.10467122979766, 2.42781092181237, 4.36664907390063, 4.83667044410385, 2.95658496329099, 0.0777040707962936),
  y = c(-1.97611823905796, -1.69248666124853, -1.45818231436247, -0.188006118085423, 0.86019753903641, -0.311324195393874, -0.668946619588382, -0.594955773203312, -1.17455073655303, -1.72948208444106, -2.27208162459825)
)
xy <- cbind(logit(colbias), trait_adj_ps[names(colbias)])
dists <- as.matrix(dist(rbind(piks, xy)))[1:nrow(piks), nrow(piks)+1:nrow(xy)]
pts.to.label <- apply(dists, 1, which.min)

plot.scatter(x = logit(colbias), 
             y = trait_adj_ps[names(colbias)],
             cex = 1.5, pch = 19, 
             xlabel = latex2exp::TeX("logit(difference term posterior mass > 0)"),
             ylabel = latex2exp::TeX("-log$_{10}$(adjusted p-value)"),
             axis.lab.dev = 0.0875, lab.cex = 1.5,
             main.title = latex2exp::TeX("\\textit{\\textbf{trait}} enrichment comparison (HMP)"),
             pt.labels = names(colbias),
             # pt.labels = NULL,
             pts.to.label = pts.to.label,
             pt.pos = c(1,3,3,2,1,2,3,3,1,1,1), pt.cex = 1.75,
             prop_label_push = c(4,8,10,4,4,3,2,8,4,9,4)/100,
             pt.col = adjustcolor(category_colors[trait_cat_key[names(colbias)]], 0.5))
# text(xy[pts.to.label,1], xy[pts.to.label,2], labels = 1:length(pts.to.label), col = 4, cex = 3)

#genes in 3 tissues, HMP for TWAS
# trait_adj_ps <- -log10(setNames(gsea_output_pooled_alt$adj.pval, gsea_output_pooled_alt$trait_name))
trait_adj_ps <- -log10(setNames(gsea_output_pooled_alt$pval * length(gsea_output_pooled_alt$pval), gsea_output_pooled_alt$trait_name))
# temp <- locator(); cat(paste0("piks <- cbind(\n\tx = c(", paste0(temp$x, collapse = ", "), "),\n", "\ty = c(", paste0(temp$y, collapse = ", "), ")\n)"))
piks <- cbind(
  x = c(-7.3179032742577, -7.20766574694757, -6.25962301208041, -5.64229285914366, -2.53359458899785, -1.80602690875097, 2.5814266781924, 3.00032928197091, 4.23498958784442, 4.69798720254698, 3.2428518420532, 2.80190173281267),
  y = c(-2.05367997178105, -1.82772236457882, -1.61158900116798, -1.43475261292275, -0.776528278898839, -0.452328233782585, 1.30621140487831, 0.431853707443563, 0.0192354682046946, -0.206722138997544, -0.648813109610618, -1.52317080704536)
)
xy <- cbind(logit(colbias), trait_adj_ps[names(colbias)])
dists <- as.matrix(dist(rbind(piks, xy)))[1:nrow(piks), nrow(piks)+1:nrow(xy)]
pts.to.label <- apply(dists, 1, which.min)

plot.scatter(x = logit(colbias), 
             y = trait_adj_ps[names(colbias)],
             cex = 1.5, pch = 19, 
             xlabel = latex2exp::TeX("logit(difference term posterior mass > 0)"),
             ylabel = latex2exp::TeX("-log$_{10}$(adjusted p-value)"),
             axis.lab.dev = 0.0875, lab.cex = 1.5,
             main.title = latex2exp::TeX("\\textit{\\textbf{trait}} enrichment comparison (multi-tissue)"),
             pt.labels = names(colbias),
             # pt.labels = NULL,
             pts.to.label = pts.to.label,
             pt.pos = c(1,3,3,3,3,
                        3,2,2,3,3,
                        2,4), pt.cex = 1.75,
             prop_label_push = c(4,11,11,13,8,
                                 8,4,4,10,6,
                                 2,4)/100,
             pt.col = adjustcolor(category_colors[trait_cat_key[names(colbias)]], 0.5))
# text(xy[pts.to.label,1], xy[pts.to.label,2], labels = 1:length(pts.to.label), col = 4, cex = 3)

#change size for smaller plots
par(xpd = NA, mar = c(4,1,3,0), pty="s")

#tissue aggregation
tissue_adj_ps <- -log10(tissue_ps / (tissue_thresh / 0.05))
pts.to.label <- 1:length(tissue_adj_ps)
xy <- cbind(x = logit(rowbias), 
            y = tissue_adj_ps[names(rowbias)])
plot.scatter(x = logit(rowbias), 
             y = tissue_adj_ps[names(rowbias)],
             cex = 1.5, pch = 19, 
             xlabel = latex2exp::TeX("logit(difference term posterior mass > 0)"),
             ylabel = latex2exp::TeX("-log$_{10}$(adjusted p-value)"),
             axis.lab.dev = 0.225, lab.cex = 1.5,
             main.title = latex2exp::TeX("\\textit{\\textbf{tissue}} enrichment comparison (HMP)"),
             pt.labels = tolower(names(rowbias)),
             # pt.labels = NULL,
             pts.to.label = pts.to.label,
             pt.pos = c(3,2,1,1,2,
                        1,2,3,1,4,
                        3,1,3,1,1), 
             pt.cex = 1.75,
             prop_label_push = c(2.5,12,4,14,4,
                                 4,4,25,8,4,
                                 4,17,4,8,4)/100,
             pt.col = adjustcolor(TISSUE_COLORS[names(rowbias)], 0.5),
             prop_rect_expand = 0.1)
# text(xy[pts.to.label,1], xy[pts.to.label,2], labels = 1:length(pts.to.label), col = 4, cex = 1)

#trait category aggregation
trait_category_adj_ps <- -log10(trait_category_ps / (trait_category_thresh / 0.05))
pts.to.label <- 1:length(trait_category_adj_ps)
xy <- cbind(x = logit(colcatbias),
            y = trait_category_adj_ps[names(colcatbias)])
plot.scatter(x = logit(colcatbias),
             y = trait_category_adj_ps[names(colcatbias)],
             cex = 1.5, pch = 19,
             xlabel = latex2exp::TeX("logit(difference term posterior mass > 0)"),
             ylabel = latex2exp::TeX("-log$_{10}$(adjusted p-value)"),
             axis.lab.dev = 0.225, lab.cex = 1.5,
             main.title = latex2exp::TeX("\\textit{\\textbf{trait category}} enrichment comparison (HMP)"),
             pt.labels = gsub("( |\\-).*", "", names(colcatbias)),
             # pt.labels = NULL,
             pts.to.label = pts.to.label,
             pt.pos = c(1,2,2,4,3,
                        2,2,1,1,4,
                        1,3), 
             pt.cex = 1.75,
             prop_label_push = c(6,6,15,6,6,
                                 6,6,3,6,6,
                                 3,6)/100,
             pt.col = adjustcolor(category_colors[names(colcatbias)], 0.5),
             prop_rect_expand = 0.1)
# text(xy[pts.to.label,1], xy[pts.to.label,2], labels = 1:length(pts.to.label), col = 4, cex = 1)

#legends
legend(x = par("usr")[2] + diff(par("usr")[1:2]) / 2, 
       y = par("usr")[3] + diff(par("usr")[3:4]) * 1.1, 
       legend = names(category_colors), pch = 21, 
       col = 1, pt.bg = adjustcolor(category_colors, 0.5), 
       cex = 1.5, pt.cex = 1.75, bty = "n", title = latex2exp::TeX("\\textbf{Trait Categories}"))
legend(x = par("usr")[2] + diff(par("usr")[1:2]) / 2, 
       y = par("usr")[3] + diff(par("usr")[3:4]) * 2.55, 
       legend = tissues, pch = 21, ncol = 2,
       col = 1, pt.bg = adjustcolor(TISSUE_COLORS[tissues], 0.5), 
       cex = 1.5, pt.cex = 1.75, bty = "n", title = latex2exp::TeX("\\textbf{Tissues}"))

#figure subpanel labels
fig_label("a)", cex = 2.5, shrinkX = 30.5, shrinkY = 8.5)
fig_label("b)", cex = 2.5, shrinkX = 30.5, shrinkY = 3.5)
fig_label("c)", cex = 2.5, shrinkX = 0.5, shrinkY = 8.5)
fig_label("d)", cex = 2.5, shrinkX = 3.5, shrinkY = 3.75)
fig_label("e)", cex = 2.5, shrinkX = 3.5, shrinkY = 1.25)

# legend(x = xr, y = yt, legend = tissues, pch = 19, col = adjustcolor(MotrpacRatTraining6moData::TISSUE_COLORS[tissues], 0.5), cex = 0.55, pt.cex = 1.2)
# segments(x0 = par("usr")[1] + diff(par("usr")[1:2]) / 1.55, 
#          y0 = par("usr")[4] - diff(par("usr")[3:4])/50,
#          x1 = par("usr")[1] + diff(par("usr")[1:2]) / 1.45, 
#          y1 = par("usr")[4] - diff(par("usr")[3:4])/50, 
#          lty = 2, lwd = 2, col = adjustcolor(1,0.75))
# text(x = par("usr")[1] + diff(par("usr")[1:2]) / 1.45, 
#      y = par("usr")[4] - diff(par("usr")[3:4])/45, 
#      pos = 4, cex = 0.75, labels = "1-to-1 line")



dev.off()

