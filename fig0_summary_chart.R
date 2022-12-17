#libraries
library(Cairo)

#### functions ####
source(file = "~/scripts/montgomery_lab/deg-trait_functions.R")

#### node parameters ####
cols <- c(rgb(255,255,255,m=255),
          rgb(255,255,255,m=255),
          rgb(240,172,1,m=255), 
          rgb(32,48,80,m=255), 
          rgb(139,10,10,m=255), 
          rgb(16,94,38,m=255), 
          rgb(204,85,0,m=255), 
          rgb(54,1,63,m=255)
        )
cols <- c(cols, rep(cols[8], 2), "grey50")
ws <- c(0.75, 0.5, rep(0.75, 8), 2) * 1.5
hs <- c(0.75,
        0.75,
        rep(0.75, 8), 1.25)
hat_prop <- 0.15

locs <- do.call(rbind, list(
  c(0.125, 1.25), 
  c(0.125, 2.5),
  c(1,1),
  c(1,2),
  c(1,3),
  c(2,1),
  c(2,2),
  c(2,3),
  c(3,3),
  c(4,3),
  c(3.5,1.25))
)

node_titles <- c(
  "Rat Model\n(Chronic Exercise Training)",
  "Human\nPopulations",
  "Exercise Training Response",
  "GTEx Gene x Tissue Expression",
  "Genome Wide Association Studies",
  "Relative Effect Size Estimates",
  "SNP-Heritability Estimates",
  "Gene x Tissue Associations",
  "DEG \u2229 PrediXcan Enrichment",
  "Protective Enrichment",
  "Gene x Tissue x Trait Prioritization"
)
title_cex <- c(rep(0.75, 2), rep(0.525, length(node_titles)-3), 1.2)
title_disp_mod <- c(5.5,-2, rep(0.9, length(node_titles)-3), 0.25)

node_text <- c(
  "",
  "",
  paste0("\u2022 Male and Female F344 Rats\n  trained across 1, 2, 4, and 8 weeks\n\n",
         "\u2022 874 samples across 19 Tissues\n  subjected to RNA-Seq and DEA\n\n",
         "\u2022 repfdr used to cluster genes"),
  paste0("\u2022 GTEx v8 downloaded representing\n  54 tissues across 948 donors\n\n",
         "\u2022 15 human tissues able to be\n  matched to comparable tissues in rats\n\n"),
  paste0("\u2022 Genome-wide summary statistics\n  across 114 publicly available\n  GWAS downloaded\n\n",
         "\u2022 Represent a variety of plausibly\n  exercise-responsive human traits\n  across many diverse categories"),
  paste0("\u2022 Non-exercise rat transcriptome lacks\n  variability\n\n",
         "\u2022 To contextualize observed DE\n  we residualized out sex, ancestry,\n  and other variation and estimated\n  pseudolog-2 SD, comparing to rat DE"),
  paste0("\u2022 Using annotations based on our\n  DE genesets, we estimated tissue-\n  specific h2_SNP enrichment in LDSC\n\n",
         "\u2022 Additionally estimated expression-\n  mediated h2_SNP enrichment in\n  MESC"),
  paste0("\u2022 Downloaded and filtered S-PrediXcan\n  output, x-referencing eQTL and\n  GWAS association statistics\n\n",
         "\u2022 Causal relations between tissue-\n  specific gene expression and trait\n  expression revealed using SNPs"),
  paste0("\u2022 Evaluate whether S-PrediXcan hits\n  more common in DE gene sets using\n  a multilevel Bayesian model\n\n",
         "\u2022 Allows for partial pooling of signal\n  across tissues, traits, and trait\n  categories"),
  paste0("\u2022 Intersecting gene sets may not be\n  not reflect anticipated (e.g. protective)\n  effects of exercise\n\n",
         "\u2022 Examine whether gene intersects\n  enriched in protective or anti-protective\n  directions"),
  paste0("\u2022 Identify tissues with confident or strong enrichment for directional\n  effects in specific human phenotypes\n\n",
         "\u2022 Reference these against genes whose DE from exercise is large\n  relative to standing variation in human populations\n\n",
         "\u2022 Further integrate biological significance of specific tissues by\n  examing heritability enrichment results\n\n")
  
)
text_cex <- c(rep(0.75, 2), rep(0.525, length(node_titles)-3), 0.85)
CairoFonts(
  regular="Arial Unicode MS:style=Regular",
  bold="Arial Unicode MS:style=Bold",
  italic="Arial Unicode MS:style=Italic",
  bolditalic="Arial Unicode MS:style=Bold Italic,BoldItalic",
  symbol="Symbol", usePUA=TRUE
)

locs[,1] <- 2 * locs[,1]
pes <- c(NA,
         NA,
         rep(0.25, 8), 0.125)
imgsrc <- c("~/repos/MoTrPAC_Complex_Traits/figures/treadmill_rat.png",
            "~/repos/MoTrPAC_Complex_Traits/figures/motrpac_human.png")

#### edge parameters ####
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



#### wordy figure plotting ####

#set up the plot
cairo_pdf(paste0("~/repos/MoTrPAC_Complex_Traits/figures/fig0_summary_chart.pdf"),
          width = 1500/72,
          height = 500/72,
          family="Arial Unicode MS", pointsize = 18.5)
par(xpd = T, mar = c(0,0,0,0))
plot(NULL, xaxt="n",yaxt="n",bty="n",pch="",ylab="",xlab="", main="", sub="", 
     xlim = c(0,8.5), ylim = c(0.5,3.5))

#first draw the edges
for(i in 1:nrow(edges)){
  
  if(edge_dir[i] == "h"){
    
    grad_arrow_curve(c(locs[edges[i,1],1] + ws[edges[i,1]]/2 + 0.02 + ifelse(edges[i,1]==2, -0.25, 0) + ifelse(edges[i,1]==1, -0.2, 0), 
                       locs[edges[i,2],1] - ws[edges[i,2]]/2 - 0.02, 
                       locs[edges[i,1],2], 
                       locs[edges[i,2],2]), 
                     cols = cols[edges[i,]], col_alpha = 0.75, direc = edge_dir[i], outline_col = edge_outline_cols[i])  
  
  } else if(edge_dir[i] == "v"){
    
    if(edges[i,1]==1){
      edge_cols <- c("white", "black")
    } else {
      edge_cols <- cols[edges[i,]]  
    }
    
    grad_arrow_curve(c(locs[edges[i,1],1], 
                       locs[edges[i,2],1], 
                       locs[edges[i,1],2] - hs[edges[i,1]]/2 * ifelse(edges[i,1]==1, -1, 1) - 0.05, 
                       locs[edges[i,2],2] + hs[edges[i,2]]/2 * ifelse(edges[i,1]==1, -1, 1) + 0.05), 
                     cols = edge_cols, 
                     col_alpha = 0.75, direc = edge_dir[i],
                     prop_shaft_length = ifelse(edges[i,1]==1, 0.25, 0.15),
                     w = ifelse(edges[i,1]==1, 0.15, 0.2), outline_col = edge_outline_cols[i])
  }
  
  #edge labels
  text(labels = edge_labels[i], x = locs[edges[i,2],1], y = locs[edges[i,2],2] - hs[edges[i,2]]/1.15, 
       srt = ifelse(edge_dir[i] == "h", 0, 90), cex = 0.45, col = "white")
  
}

#then the nodes
for(i in 1:nrow(locs)){
  if(i %in% c(1,2)){
    addImg(png::readPNG(imgsrc[i]), x = locs[i,1], y = locs[i,2], width = ws[i])
  } else {
    rrect(loc = locs[i,], w = ws[[i]], h = hs[[i]], border = cols[i], col = adjustcolor(cols[i], 0.1),
          lwd = 2, pe = pes[i], hat_prop = hat_prop)
    text(locs[i,1] - ws[i]/1.9, locs[i,2] + hs[i]/2, labels = i, cex = 1.25)  
  }
  
  #write the node titles
  # par(family = "Arial Unicode MS Bold")
  text(locs[i,1], locs[i,2] + hs[i]/2 - hs[i] * hat_prop - strheight(node_titles[i], units = "u", cex = title_cex[i]) * title_disp_mod[i], 
         labels = latex2exp::TeX(node_titles[i]), cex = title_cex[i], pos = 3, font = 2)
  # par(family = "Arial Unicode MS")
  
  #write the node text
  text(locs[i,1] - ws[i]/2, locs[i,2] + hs[i]/2.25 - strheight(node_text[i], units = "u", cex = text_cex[i]) - hs[i] * hat_prop,
       labels = latex2exp::TeX(node_text[i]), cex = text_cex[i], pos = 4)
}


dev.off()


