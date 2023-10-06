#### run preprocessing script ####
source("/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/helper_scripts/figure_set_1_preprocessing.R")

#### figure plotting ####

# if(any(!(unlist(.Devices) %in% c("", "null device")))){
#   for(i in 1:(sum(!(unlist(.Devices) %in% c("", "null device"))))){
#     dev.off()
#   }
# }

# flowchart node parameters
cols <- c(rgb(255,255,255,m=255),
          rgb(255,255,255,m=255),
          rgb(240,172,1,m=255), 
          rgb(32,48,80,m=255), 
          rgb(139,10,10,m=255), 
          rgb(16,94,38,m=255), 
          rgb(204,85,0,m=255), 
          rgb(54,1,63,m=255)
)
cols <- c(cols, rep(cols[8], 2), 1)
ws <- c(1, 0.65, rep(0.75, 8), 2.25) * 1.5
hs <- c(0.75,
        0.75,
        rep(0.75, 8), 1.25)
hat_prop <- 0
raster <- T
nslices <- 2E3
raster_res <- 81
# raster <- F
nslices <- 1E3
arrow_alpha <- 0.75

locs <- do.call(rbind, list(
  c(0.25, 1.2), 
  c(0.25, 2.4),
  c(2,1),
  c(2,2),
  c(2,3),
  c(4.25,1),
  c(4.25,2),
  c(4.25,3),
  c(6.5,3),
  c(8.75,3),
  c(7.5,1.375))
)

node_text <- c(
  "\\textbf{F344 Rats}\n(MoTrPAC EET Study)",
  "\\textbf{Humans}\n(GTEx, Public GWAS)",
  "\\textbf{Exercise}\nChronic exercise\nexperiment, 19\ntissues collected\nat 1, 2, 4, and 8\nweek timepoints",
  "\\textbf{GTEx}\nGene$\\times$Tissue\nExpression data\n& eQTL mapping\nacross 54 tissues",
  "\\textbf{GWAS}\n114 Genome-\nWide Association\nStudies across 12\ntrait categories",
  "\\textbf{Relative Effects}\nContextualize DE\nrelative to natural\nhuman variation\nin GTEx ($log_2$\nscale)",
  "\\textbf{$h^2_{SNP}$ Estimation}\nQuantify both\nSNP heritability\nenrichment and\nexpression $h^2_{SNP}$",
  "\\textbf{TWAS}\nGene x Tissue\nAssociations from\nGWAS & eQTL\noutput",
  "\\textbf{DEG \u2229 TWAS}\nCross-reference\ntissue-specific DE\ngenes and TWAS\nhits against prior\nexpectation",
  "\\textbf{Directionality}\nEvaluate gene\nintersects for\nenrichment in +/-\nassociations",
  "\\textbf{Gene x Tissue x Trait Prioritization}\nIdentify shared etiology of impactful exercise\nresponse across transcriptomic architecture\nof human phenotypic variation, implicating\nmultiple organs and organ systems"
)

# node_text <- rep("", length(node_text))
text_vdisp <- c(-0.9, 0.475, rep(-0.1125, length(node_text)-3), -0.2)
text_hdisp <- c(0.075, -0.1, rep(0.035, length(node_text)-3), 0.075)
text_cex <- c(rep(0.75, 2), rep(0.7, length(node_text)-3), 0.9) * 2

CairoFonts(
  regular="Arial Unicode MS:style=Regular",
  bold="Arial Unicode MS:style=Bold",
  italic="Arial Unicode MS:style=Italic",
  bolditalic="Arial Unicode MS:style=Bold Italic,BoldItalic",
  symbol="Symbol", usePUA=TRUE
)

pes <- c(NA,
         NA,
         rep(0.25, 8), 0.125) / 2
imgsrc <- c("/Volumes/2TB_External/MoTrPAC_Complex_Traits/figures/treadmill_rat.png",
            "/Volumes/2TB_External/MoTrPAC_Complex_Traits/figures/motrpac_human.png")

# flowchart edge parameters
edges <- do.call(rbind, list(c(1,2),
                             c(1,3),
                             c(2,4),
                             c(2,5),
                             c(3,6),
                             c(3,7),
                             c(4,6),
                             c(4,7),
                             c(4,8),
                             c(5,7),
                             c(5,8),
                             c(8,9),
                             c(9,10),
                             c(6,11),
                             c(7,11),
                             c(9,11),
                             c(10,11))
)

edge_dir <- c("v", rep("h", 14), rep("v", 2))
edge_labels <- c("         RGD Orthologs", rep("", length(edge_dir)-1))
edge_outline_cols <- c(rep("white", 4), rep("black", nrow(edges) - 4))
edge_disp_x <- cbind(0.04 + as.numeric(edges[,1] == 2) * -0.35 + as.numeric(edges[,1] == 1) * -0.4,
                     rep(-0.04, nrow(edges)))
edge_disp_y <- cbind(as.numeric(edge_dir == "v") * -0.025,
                     as.numeric(edge_dir == "v") * 0.025)
prop_shaft_length <- 0.25 + as.numeric(edges[,1] == 1) * 0.1
shaft_width <- 0.25 + as.numeric(edges[,1]==1) * as.numeric(edges[,2]==2) * 0.15
prop_head_width <- 2.1 + as.numeric(edges[,1]==1) * 0.1 + as.numeric(edges[,1] %in% c(9,10)) * as.numeric(edges[,2] %in% c(11)) * 0.7

#check if grex script has been called
if(!exists("called_grex_script")){called_grex_script <- F}
if(!called_grex_script){
  called_grex_script <- T
  source("/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/analyses/analysis_GREx_RelativeEffectSize.R")
}
if(!exists("relative_expression_data")){
  load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/relative_expression_motrpac_gtex")
}

# actual plotting
grDevices::cairo_pdf(filename = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/figures/fig1_high-level-overview_redux_take-2.pdf"), 
                     width = 1350 / 72, height = 1125 / 72 * 1.5 * 4 / 6, family="Arial Unicode MS", pointsize = 18.5)

layout_mat <- kronecker(rbind(
  c(1,1,2,3),
  c(1,1,4,5),
  c(6,13,9,11),
  c(7,8,10,12)
), matrix(1,2,2))
layout_mat[5:6, 1:3] <- layout_mat[5:6, 3:1]
layout_mat <- layout_mat + 1
layout_mat <- rbind(matrix(1, nrow = 4, ncol = 8), layout_mat)
layout_mat <- layout_mat[1:8,]

layout(layout_mat, heights = c(1,1,1,1))

#high level flowchart
CairoFonts(
  regular="Arial Unicode MS:style=Regular",
  bold="Arial Unicode MS:style=Bold",
  italic="Arial Unicode MS:style=Italic",
  bolditalic="Arial Unicode MS:style=Bold Italic,BoldItalic",
  symbol="Symbol", usePUA=TRUE
)

par(xpd = T, mar = c(0,0,0,0))
plot(NULL, xaxt="n",yaxt="n",bty="n",pch="",ylab="",xlab="", main="", sub="", 
     xlim = c(-0.125,9), ylim = c(0.5,3.5))

#first draw the edges
for(i in 1:nrow(edges)){
  
  if(edge_dir[i] == "h"){
    
    grad_arrow_curve(c(locs[edges[i,1],1] + ws[edges[i,1]]/2 + edge_disp_x[i,1], 
                       locs[edges[i,2],1] - ws[edges[i,2]]/2 + edge_disp_x[i,2], 
                       locs[edges[i,1],2], 
                       locs[edges[i,2],2]), w = shaft_width[i],
                     prop_shaft_length = prop_shaft_length[i], prop_head_width = prop_head_width[i],
                     cols = cols[edges[i,]], col_alpha = arrow_alpha, direc = edge_dir[i], outline_col = edge_outline_cols[i], raster = raster, 
                     nslices = nslices, raster_res = raster_res)  
    
  } else if(edge_dir[i] == "v"){
    
    if(edges[i,1]==1){
      edge_cols <- c("white", "black")
    } else {
      edge_cols <- cols[edges[i,]]  
    }
    
    grad_arrow_curve(c(locs[edges[i,1],1], 
                       locs[edges[i,2],1], 
                       locs[edges[i,1],2] - hs[edges[i,1]]/2 * ifelse(edges[i,1]==1, -1, 1) + edge_disp_y[i,1], 
                       locs[edges[i,2],2] + hs[edges[i,2]]/2 * ifelse(edges[i,1]==1, -1, 1) + edge_disp_y[i,2]), 
                     cols = edge_cols, 
                     col_alpha = arrow_alpha, direc = edge_dir[i],
                     prop_shaft_length = prop_shaft_length[i], prop_head_width = prop_head_width[i],
                     w = shaft_width[i], outline_col = edge_outline_cols[i], raster = raster, 
                     nslices = nslices, raster_res = raster_res)
  }
  
  #edge labels
  text(labels = edge_labels[i], x = locs[edges[i,2],1], y = locs[edges[i,2],2] - hs[edges[i,2]]/1.15, 
       srt = ifelse(edge_dir[i] == "h", 0, 90), cex = 0.8, col = "white")
  
}

#then the nodes
for(i in 1:nrow(locs)){
  if(i %in% c(1,2)){
    addImg(png::readPNG(imgsrc[i]), x = locs[i,1], y = locs[i,2], width = ws[i])
  } else {
    rrect(loc = locs[i,], w = ws[[i]], h = hs[[i]], border = cols[i], col = adjustcolor(cols[i], 0.1),
          lwd = 2, pe = pes[i], hat_prop = hat_prop, bold_border = ifelse(i == nrow(locs), 0.05, 0))
    # text(locs[i,1] - ws[i]/1.9, locs[i,2] + hs[i]/2, labels = i, cex = 1.25)  
  }
  
  #write the node text
  # text(locs[i,1] - ws[i]/2, locs[i,2] + hs[i]/2.25 - strheight(node_text[i], units = "u", cex = text_cex[i]) - hs[i] * hat_prop + text_vdisp[i],
  #      labels = latex2exp::TeX(node_text[i]), cex = text_cex[i], pos = 4, font = 3) #font does not work with expressions :/ gotta use \textbf{}, or else it works w/ normal text()
  first_line_text_adj <- strheight(node_text[i], cex = text_cex[i])
  latext(labs = node_text[i], 
         x = locs[i,1] - ws[i]/2 + text_hdisp[i], 
         y = locs[i,2] + hs[i]/2 + text_vdisp[i], 
         pos = 4, cex = text_cex[i], 
         boxh = ifelse(cols[i] == "#FFFFFF", NA, hs[i]), 
         first_line_col = ifelse(cols[i] == "#FFFFFF", 1, cols[i]),
         first_line_hadj = NA, first_line_center_loc = ifelse(i %in% c(1,2), NA, locs[i,1]))
  
}


#node intersects plot
cols = list(Tissue=MotrpacRatTraining6moData::TISSUE_COLORS[names(motrpac_gtex_map)], 
            Time=MotrpacRatTraining6moData::GROUP_COLORS[paste0(c(1,2,4,8), "w")],
            Sex=MotrpacRatTraining6moData::SEX_COLORS[c('male','female')])
cols$Tissue[which(is.na(cols$Tissue))] <- '#C0C0C0'
par(mar = c(4,5,4,1), xpd = NA)
load("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/node_metadata_list.RData")
ensembl_genes <- orig_ensembl_genes <- lapply(split(node_metadata_list$`8w`$human_ensembl_gene[!is.na(node_metadata_list$`8w`$human_ensembl_gene)], 
                                                    node_metadata_list$`8w`$tissue[!is.na(node_metadata_list$`8w`$human_ensembl_gene)]), unique)
symbol_map <- unique(node_metadata_list$`8w`[,c("human_gene_symbol", "human_ensembl_gene")])
all_genes <- unlist(orig_ensembl_genes)
n_tissues_per_gene <- table(all_genes)
ensembl_genes$THREE <- names(n_tissues_per_gene)[n_tissues_per_gene > 2]


tissue_colors <- c(MotrpacRatTraining6moData::TISSUE_COLORS, THREE = "black")
jaccard <- function(x, y) length(intersect(x,y)) / length(union(x,y))
jacmat <- sapply(orig_ensembl_genes, function(x) sapply(orig_ensembl_genes, function(y) jaccard(x,y)))
jacmat_inds <- order(cmdscale(1-jacmat, k = 1))
jacmat <- jacmat[jacmat_inds, jacmat_inds]
n_tissues_per_gene <- table(unlist(orig_ensembl_genes))
nt2pg <- table(n_tissues_per_gene)
ensembl_genes_df <- data.frame(gene = unlist(orig_ensembl_genes),
                               tissue = rep(names(orig_ensembl_genes), sapply(orig_ensembl_genes, length), each = T))
ensembl_genes_df$n_tissue <- n_tissues_per_gene[match(ensembl_genes_df$gene, names(n_tissues_per_gene))]
n_per_cat <- lapply(setNames(sort(unique(n_tissues_per_gene)), sort(unique(n_tissues_per_gene))), 
                    function(nt){x <- table(ensembl_genes_df$tissue[ensembl_genes_df$n_tissue == nt]); x})
prop_per_cat <- lapply(n_per_cat, function(nt) nt / sum(nt))


plot(NA, xlim = c(0.5, length(prop_per_cat)+0.5), ylim = c(0,2), frame = F, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
for(nt in 1:length(prop_per_cat)){
  ntcs <- cumsum(n_per_cat[[nt]])
  ptcs <- ntcs / max(ntcs)
  rect(xleft = nt - 0.5, xright = nt + 0.5, yb = c(0, ptcs[-length(ptcs)]), ytop = ptcs,
       col = tissue_colors[names(ptcs)])
  rect(xleft = nt - 0.5, xright = nt + 0.5, yb = 1.05, ytop = 1.05 + log10(nt2pg[nt]) / log10(max(nt2pg)) * 0.9 + 0.05,
       col = "lightgrey")
  text(nt, 1.05 + log10(nt2pg[nt]) / log10(max(nt2pg)) * 0.9 + 0.05, labels = nt2pg[nt], pos = 3, cex = 1.5, xpd = NA)
}

#axes
#bottom plot
text(x = 1:6, labels = 1:6, y = 0, pos = 1, cex = 1.5, font = 2)
nnum <- 6
segments(x0 = 0.425, x1 = 0.425, y = 0, y1 = 1, lwd = 2)
segments(x0 = 0.425, x1 = 0.35, y = seq(0, 1, length.out = nnum), y1 = seq(0, 1, length.out = nnum), lwd = 2)
text(x = 0.35, y = seq(0, 1, length.out = nnum), labels = seq(0,1,length.out=nnum), pos = 2, cex = 1.25, xpd = NA)
text(labels = "Number of Tissues Differentially Expressing the Same Gene", x = (length(nt2pg) + 1) / 2, y = -0.125, pos = 1, cex = 1.5, xpd = NA)
text(labels = "Tissue Composition of Bin", x = -0.15, y = 0.5, pos = 3, cex = 1.5, xpd = NA, srt = 90)

#top plot
log10_count_ticks <- (0:floor(max(log10(nt2pg))))
log10_count_ylocs <- 1.1 + 0.9 * log10_count_ticks / log10(max(nt2pg))
segments(x0 = 0.425, x1 = 0.425, y = 1.1, y1 = 2, lwd = 2)
segments(x0 = 0.425, x1 = 0.35, y = log10_count_ylocs, y1 = log10_count_ylocs, lwd = 2)
text(x = 0.375, y = log10_count_ylocs, labels = latex2exp::TeX(paste0("$10^{", log10_count_ticks, "}$")), pos = 2, cex = 1.25, xpd = NA)
text(labels = "Number of Genes in Bin", x = -0.15, y = 1.525, pos = 3, cex = 1.5, xpd = NA, srt = 90)


#jaccard matrix
xl = 2; xr = 6; yb = 1.5; yt = 2
ncols <- 101
rate = 0.04
exp_dist_cols <- round(cumsum(c(1, dexp(1:(ncols-1), rate = rate) / min(dexp(1:(ncols-1), rate = rate)))))
colgrad <- viridis::viridis(max(exp_dist_cols))[exp_dist_cols]

xyrat <- diff(par("usr")[1:2]) / diff(par("usr")[3:4])

for(j in 2:ncol(jacmat)){
  rect(xleft = xl + (j-1) / nrow(jacmat) * (xr - xl), xright = xl + j / nrow(jacmat) * (xr - xl), 
       ybottom = yt + 1 / ncol(jacmat) * (yt - yb), ytop =  yt + 2 / ncol(jacmat) * (yt - yb),
       col = MotrpacRatTraining6moData::TISSUE_COLORS[colnames(jacmat)[j]])
  text(x = xl + (j-0.75) / nrow(jacmat) * (xr - xl), y = yt + 2.5 / ncol(jacmat) * (yt - yb),
       labels = colnames(jacmat)[j], pos = 4, cex = 0.75, xpd = NA, srt = 45)
  
  for(i in 1:(j-1)){
    rect(xleft = xl + (j-1) / nrow(jacmat) * (xr - xl), xright = xl + j / nrow(jacmat) * (xr - xl), 
         ybottom = yb + (nrow(jacmat)-i) / ncol(jacmat) * (yt - yb), ytop =  yb + (nrow(jacmat)-i+1) / ncol(jacmat) * (yt - yb),
         col = colgrad[floor(jacmat[i,j] * 100) + 1])
    if(j == ncol(jacmat)){
      rect(xleft = xr + 1 / xyrat / nrow(jacmat) * (xr - xl), xright = xr + (1 + 1 / xyrat) / nrow(jacmat) * (xr - xl), 
           ybottom = yb + (nrow(jacmat)-i) / ncol(jacmat) * (yt - yb), ytop =  yb + (nrow(jacmat)-i+1) / ncol(jacmat) * (yt - yb),
           col = MotrpacRatTraining6moData::TISSUE_COLORS[rownames(jacmat)[i]])
      text(x = xr + (0.9 + 1 / xyrat) / nrow(jacmat) * (xr - xl), y = yb + (nrow(jacmat)-i+0.5) / ncol(jacmat) * (yt - yb),
           labels = rownames(jacmat)[i], pos = 4, cex = 0.75)
    }
  }
}

#legend for heatmap
rect(xleft = xr + (0.5 + 1 / xyrat) / nrow(jacmat) * (xr - xl), xright = xr + (1 + 1 / xyrat) / nrow(jacmat) * (xr - xl),
     ybottom = seq(yb - (yt-yb)/1.5, yb, length.out = ncols+1)[-ncols], ytop =  seq(yb - (yt-yb)/1.5, yb, length.out = ncols+1)[-1],
     col = colgrad, border = colgrad)
rect(xleft = xr + (0.5 + 1 / xyrat) / nrow(jacmat) * (xr - xl), xright = xr + (1 + 1 / xyrat) / nrow(jacmat) * (xr - xl),
     ybottom = yb, ytop =  yb - (yt-yb)/1.5, border = 1)
nnum <- 6
text(x = xr + (0.875 + 1 / xyrat) / nrow(jacmat) * (xr - xl), y = seq(yb - (yt-yb)/1.5, yb, length.out = nnum), labels = seq(0,1,length.out=nnum), pos = 4, cex = 0.75)
text(x = xr + (0.375 + 1 / xyrat) / nrow(jacmat) * (xr - xl), y = yb - (yt-yb)/3, labels = "Jaccard Index", pos = 3, cex = 0.75, srt = 90)


#now do the Open Targets curves
if(!exists("n_traits_above_at_least_1")){
  use_indirect <- F
  load(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/open-targets_tissue-x-disease_", 
              ifelse(use_indirect, "indirect", "direct"), 
              "-associations"))
  breakpoints <- 0:100/100
  n_above <- sapply(names(tissue_x_disease), function(tissue) 
    log10((data.frame(cats = cut(tissue_x_disease[[tissue]][tissue_x_disease[[tissue]] > 1E-6], breaks=c(breakpoints, Inf)), ordered_result=TRUE) %>% 
             dplyr::count(cats, .drop=FALSE) %>% arrange(desc(cats)) %>% mutate(cumfreq = cumsum(n)) %>% arrange(cats))$cumfreq))
  prop_above <- log10(10^n_above %*% (diag(1 / sapply(ensembl_genes, length))))
  colnames(prop_above) <- colnames(n_above)
  
  n_above_at_least_1 <- sapply(names(tissue_x_disease), function(tissue) 
    log10((data.frame(cats = cut(apply(tissue_x_disease[[tissue]], 1, max)[apply(tissue_x_disease[[tissue]], 1, max) > 1E-6], breaks=c(breakpoints, Inf)), ordered_result=TRUE) %>% 
             dplyr::count(cats, .drop=FALSE) %>% arrange(desc(cats)) %>% mutate(cumfreq = cumsum(n)) %>% arrange(cats))$cumfreq))
  
  
  n_traits_above_at_least_1 <- sapply(names(tissue_x_disease), function(tissue) 
    log10((data.frame(cats = cut(apply(tissue_x_disease[[tissue]], 2, max)[apply(tissue_x_disease[[tissue]], 2, max) > 1E-6], breaks=c(breakpoints, Inf)), ordered_result=TRUE) %>% 
             dplyr::count(cats, .drop=FALSE) %>% arrange(desc(cats)) %>% mutate(cumfreq = cumsum(n)) %>% arrange(cats))$cumfreq))
}

#now plot Open Targets curves
main_title <- paste0("Open Targets ", ifelse(use_indirect, "Overall / Indirect", "Direct"), " Associations")

par(xpd = F)


curves_list <- list(n_above, prop_above, n_above_at_least_1, n_traits_above_at_least_1)
vlabs <- c("# trait x gene pairs\n≥ given evidence score",
           "average # traits per gene\n≥ given evidence score",
           "# genes with ≥ one association\n≥ given evidence score",
           "# traits with ≥ one association\n≥ given evidence score")

for(i in 1:4){
  
  if(i < 3){
    par(mar = c(3,4,2,4))
  } else {
    par(mar = c(4,4,1,4))
  }
  
  curves_to_use <- curves_list[[i]]
  vlab <- vlabs[[i]]
  ylims <- range(curves_to_use[is.finite(curves_to_use)])
  
  plot(NA, NA, xlim = rev(c(0,1)), ylim = ylims, xlab = "", ylab = "", 
       yaxt = "n", xaxt = "n", main = "", xpd = NA, frame = F)
  text(x = 1.1675, y = mean(par("usr")[3:4]), label = vlab, srt = 90, pos = 3, xpd = T)
  text(x = mean(par("usr")[1:2]), y = par("usr")[3] - diff(par("usr")[3:4])/10, pos = 1, xpd = T, label = "evidence score")
  
  if(i == 1){
    text(main_title, x = -0.425, y = ylims[2] + diff(ylims) / 20, pos = 3, xpd = NA, cex = 2)
  }
  segments(y0 = ceiling(ylims[1]):floor(ylims[2]), y1 = ceiling(ylims[1]):floor(ylims[2]), x0 = 0, x1 = 1, lty = 3)
  segments(y0 =ceiling(ylims[1]):floor(ylims[2]), y1 = ceiling(ylims[1]):floor(ylims[2]), x0 = 0, x1 = 1, lty = 3)
  minticks <- log10(2:9) + rep(floor(ylims[1]):ceiling(ylims[2]), each = 8)
  minticks <- minticks[minticks < par("usr")[4] & minticks > par("usr")[3]]
  segments(y0 = minticks, y1 = minticks, x0 = 0, x1 = 1, lty = 3, lwd = 0.5)
  segments(y0 = minticks, x0 = 1, y1 = minticks, x1 = 1.0125, lwd = 0.75, xpd = NA)
  segments(y0 = -1E9, x0 = 0:5/5, y1 = 1E9, x1 = 0:5/5, lwd = 0.75, lty = 3)
  
  #yax
  segments(x0 = 1, x1 = 1, y0 = par("usr")[3], y1 = par("usr")[4], lwd = 2)
  segments(x0 = 1.025, x1 = 1, y0 = ceiling(ylims[1]):floor(ylims[2]), y1 = ceiling(ylims[1]):floor(ylims[2]), lwd = 2)
  text(x = 0.99, y = ceiling(ylims[1]):floor(ylims[2]), pos = 2, xpd = T,
       labels = latex2exp::TeX(paste0("$10^{", ceiling(ylims[1]):floor(ylims[2]), "}$")))
  
  #hax
  segments(x0 = 0, x1 = 1, y0 = par("usr")[3], y1 = par("usr")[3], lwd = 5)
  segments(x0 = 0:5/5, x1 = 0:5/5, y0 = par("usr")[3], y1 =  par("usr")[3] - diff(par("usr")[3:4]) / 40, xpd = T, lwd = 2)
  text(x = 0:5/5, y = par("usr")[3] - diff(par("usr")[3:4]) / 100, pos = 1, xpd = T,
       labels = 0:5/5)
  
  #the actual curves
  for(tissue in names(tissue_x_disease)){
    lines(breakpoints, curves_to_use[,tissue], lwd = 2, col = tissue_colors[tissue])
  }
  
  #tissue labels
  tiss_order <- colnames(curves_to_use)[order(curves_to_use[1,], decreasing = T)]
  text(labels = tiss_order, col = tissue_colors[tiss_order], cex = 0.75, font = 2, xpd = NA,
       pos = 4, x = -0.05, y = ylims[2] - cumsum(c(0,rep(strheight(tiss_order[1], units = "user"), length(tissue_x_disease) - 1))))
  segments(x0 = -0.08, x1 = 0, y0 = ylims[2] - cumsum(c(0,rep(strheight(tiss_order[1], units = "user"), length(tissue_x_disease) - 1))),
           y1 = curves_to_use[1,tiss_order], col = tissue_colors[tiss_order], xpd = T)
  
}

fig_label(text = "a)", region = "plot", cex = 3, shrinkX = 6.45, shrinkY = 5.55)
fig_label(text = "b)", region = "plot", cex = 3, shrinkX = 6.45, shrinkY = 2.5)
fig_label(text = "c)", region = "plot", cex = 3, shrinkX = 3, shrinkY = 2.5)

# fig_label(text = "b)", region = "plot", cex = 3, shrinkX = 1.5, shrinkY = 1.32)
# fig_label(text = "c)", region = "plot", cex = 3, shrinkX = -1.5, shrinkY = 1.32, xpd = NA)

dev.off()
