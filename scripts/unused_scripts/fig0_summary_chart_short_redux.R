#libraries
library(Cairo)

#### functions ####
source(file = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/helper_scripts/deg-trait_functions.R")

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
cols <- c(cols, rep(cols[8], 2), 1)
ws <- c(1, 0.75, rep(0.75, 8), 2.25) * 1.5
hs <- c(0.75,
        0.75,
        rep(0.75, 8), 1.25)
hat_prop <- 0
raster <- F
nslices <- 4E2
raster_res <- 101
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
  "\textbf{F344 Rats}\n$(MoTrPAC EET Study)$",
  "\textbf{Humans}\n(GTEx, Public GWAS)",
  "\textbf{Exercise}\nChronic exercise experiment, 19 tissues collected at 1, 2, 4, and 8 week timepoints",
  "\textbf{GTEx}\nGene$\times$Tissue Expression data & eQTL mapping across 54 tissues",
  "\textbf{GWAS}\n114 Genome- Wide Association Studies across 12 trait categories",
  "\textbf{Relative Effects}\nContextualize DE relative to natural human variation in GTEx ($log_2$ scale)",
  "\textbf{$h^2_{SNP}$ Estimation}\nQuantify both SNP heritability enrichment and expression $h^2_{SNP}$",
  "\textbf{TWAS}\nGene x Tissue Associations from GWAS & eQTL output",
  "\textbf{DEG \u2229 TWAS}\nCross-reference tissue-specific DE genes and TWAS hits against prior expectation",
  "\textbf{Directionality}\nEvaluate gene intersects for enrichment in +/- associations",
  "\textbf{Gene x Tissue x Trait Prioritization}\nIdentify shared etiology of impactful exercise response across transcriptomic architecture of human phenotypic variation, implicating multiple organs and organ systems"
)

# node_text <- rep("", length(node_text))
text_vdisp <- c(-0.85, 0.55, rep(-0.1125, length(node_text)-3), -0.15)
text_hdisp <- c(0, -0.15, rep(-0.04, length(node_text)-3), 0)
text_cex <- c(rep(0.75, 2), rep(0.7, length(node_text)-3), 0.9)

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
edge_disp_x <- cbind(0.04 + as.numeric(edges[,1] == 2) * -0.35 + as.numeric(edges[,1] == 1) * -0.4,
                     rep(-0.04, nrow(edges)))
edge_disp_y <- cbind(as.numeric(edge_dir == "v") * -0.025,
                     as.numeric(edge_dir == "v") * 0.025)
prop_shaft_length <- 0.25 + as.numeric(edges[,1] == 1) * 0.1
shaft_width <- 0.25 + as.numeric(edges[,1]==1) * as.numeric(edges[,2]==2) * 0.15
prop_head_width <- 2.1 + as.numeric(edges[,1]==1) * 0.1 + as.numeric(edges[,1] %in% c(9,10)) * as.numeric(edges[,2] %in% c(11)) * 0.7


#First wrap all the text in boxes to find maximum containing cex
all_rect_coords <- lapply(1:nrow(locs), function(i){
  list(x0 = locs[i,1] - ws[i]/2,
       x1 = locs[i,1] + ws[i]/2,
       y0 = locs[i,2] - hs[i]/2,
       y1 = locs[i,2] + hs[i]/2
  )
})

#iterate through all values
i = 1
tokens <- string_to_tokens(txt_string = node_text[i])
find_optimal_cex_and_lines(txt = tokens$tokens, rect_coords = all_rect_coords[[i]], fixed_cex = NA,
                                                         newlines = tokens$has_newline)$cex
#parentheses are a special character :/ also other problems, whatever, will return to this later

#then put all the words on lines corresponding to the boxes



#### abridged figure plotting ####

#set up the plot
cairo_pdf(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/figures/fig0_summary_chart_short.pdf"),
          width = 1000/72,
          height = 450/72,
          family="Arial Unicode MS", pointsize = 18.5)


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
       srt = ifelse(edge_dir[i] == "h", 0, 90), cex = 0.45, col = "white")
  
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


dev.off()


