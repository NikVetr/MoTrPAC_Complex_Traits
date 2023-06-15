#libraries
library(Cairo)

#### functions ####
source(file = "~/scripts/montgomery_lab/deg-trait_functions.R")

latext <- function(labs, x, y, cex = 1, boxh = NA, first_line_col = 1, first_line_hadj = NA, col = 1, pos = NULL, first_line_center_loc = NA, ...){
  new_labs <- strsplit(labs, split = "\n")[[1]]
  new_labs_no_uflow <- lapply(new_labs, function(l) latex2exp::TeX(gsub("g|j|p|q|y|,|_|\\(|\\)|Q|u", "a", l)))
  new_labs_no_oflow <- lapply(new_labs, function(l) latex2exp::TeX(gsub("\\^", "a", l)))
  new_labs <- lapply(new_labs, function(l) latex2exp::TeX(l))
  
  wsh <- strheight("M\nM", cex = cex) - strheight("M", cex = cex) * 2
  lineh <- sapply(new_labs, strheight, cex = cex)  
  lineh_no_uflow <- sapply(new_labs_no_uflow, strheight, cex = cex)
  lineh_no_oflow <- sapply(new_labs_no_oflow, strheight, cex = cex)
  ebot <- lineh_no_uflow - lineh
  etop <- lineh_no_oflow - lineh
  uflow_adj <- (lineh_no_uflow - lineh) / 2 - (lineh - lineh_no_oflow) / 2
  uflow_adj <- (lineh_no_uflow - lineh)
  flow_adj <- ebot
  flow_adj[etop < -1E-6] <- flow_adj[etop < -1E-6] / 2 + etop[etop < -1E-6] / 2
  
  charh <- rep(strheight("A", cex = cex), length(new_labs)-1)
  yadj <- -cumsum(c(0, charh + wsh)) + lineh/2
  
  if(!is.na(boxh)){
    yadj[-1] <- yadj[-1] - (boxh + tail(yadj,1) + yadj[2]) / 5
  }
  
  for(i in 1:length(new_labs)){
    # print(ifelse(i==1 & !is.na(first_line_hadj), paste0(new_labs[[i]], ": ", (first_line_hadj - strwidth(new_labs[[i]], cex = cex)) / 2), 0))
    if(!is.na(first_line_center_loc) & i == 1){
      text(x = first_line_center_loc, 
           y = y + yadj[i] + flow_adj[i], labels = new_labs[[i]], cex = cex, col = first_line_col,
           pos = NULL, ...)
    } else {
      text(x = x + ifelse(i==1 & !is.na(first_line_hadj), (first_line_hadj - strwidth(new_labs[[i]], cex = cex)) / 2, 0), 
           y = y + yadj[i] + flow_adj[i], labels = new_labs[[i]], cex = cex, col = ifelse(i==1, first_line_col, col),
           pos = pos, ...)  
    }
    
    # points(x = x, y = y+yadj[i]-lineh[i]/2)
    # segments(x0 = x, x1 = x + strwidth(new_labs[[i]], cex = cex) * 1.5, y0 = y+yadj[i]-lineh[i]/2+0.0, y1 = y+yadj[i]-lineh[i]/2+0.0)
    # rect(xleft = x, xright = x + strwidth(new_labs[[i]], cex = cex), ybottom = y+yadj[i]-lineh[i]/2, ytop = y+yadj[i]+lineh[i]/2)
  }
}

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


