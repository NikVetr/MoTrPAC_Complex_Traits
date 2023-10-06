#### run preprocessing scripts ####
source("/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/helper_scripts/figure_set_3_preprocessing.R")
source("/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/helper_scripts/figure_set_4_preprocessing.R")

#### figure plotting ####
#specify which traits are good and which are bad
cols = list(Tissue=MotrpacRatTraining6moData::TISSUE_COLORS[names(motrpac_gtex_map)], 
            Time=MotrpacRatTraining6moData::GROUP_COLORS[paste0(c(1,2,4,8), "w")],
            Sex=MotrpacRatTraining6moData::SEX_COLORS[c('male','female')])
cols$Tissue<- cols$Tissue[order(match(MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV[names(cols$Tissue)], MotrpacRatTraining6moData::TISSUE_ORDER))]
tissue_names <- sapply(strsplit(names(cols$Tissue), "-"), function(x) paste0(x[ifelse(length(x) == 2, c(2), list(2:3))[[1]]], collapse = " "))
trait_goodness <- list("Standing_Height_UKB" = 0,
                       "BMI_UKB" = -1,
                       "Height" = 0,
                       "Body_Fat_Percentage_UKB" = -1,
                       "Hypertension_UKBS" = -1,
                       "Hayfever_UKB" = -1,
                       "Birth_Weight_UKB" = 1,
                       "Neuroticism_UKB" = -1,
                       "Asthma_UKB" = -1,
                       "Asthma_UKBS" = -1,
                       "Schizophrenia" = -1,
                       "Chronotype_UKB" = 0,
                       "Fluid_Intelligence_Score_UKB" = 1,
                       "Rheumatoid_Arthritis" = -1,
                       "Education_Years" = 1,
                       "Inflammatory_Bowel_Disease" = -1,
                       "HDL_Cholesterol" = 1,
                       "High_Cholesterol_UKBS" = -1,
                       "Crohns_Disease" = -1,
                       "Sleep_Duration_UKB" = 0,
                       "Triglycerides" = -1,
                       "Insomnia_UKB" = -1,
                       "Birth_Weight" = 1,
                       "Hip_Circumference_EUR" = -1,
                       "LDL_Cholesterol" = -1,
                       "Waist_Circumference_EUR" = -1,
                       "Ulcerative_Colitis" = -1,
                       "Mothers_Age_At_Death_UKB" = -1,
                       "Heart_Attack_UKB" = -1,
                       "Multiple_Sclerosis" = -1,
                       "BMI_Active_Inds" = 0,
                       "Fathers_Age_At_Death_UKB" = -1,
                       "Systemic_Lupus_Erythematosus" = -1,
                       "BMI_EUR" = 1,
                       "Psoriasis_UKBS" = -1,
                       "Deep_Venous_Thrombosis_UKBS" = -1,
                       "Coronary_Artery_Disease" = -1,
                       "Waist-to-Hip_Ratio_EUR" = -1,
                       "Alzheimers_Disease" = -1,
                       "Heart_Rate" = -1,
                       "Multiple_Sclerosis_UKBS" = -1,
                       "Rheumatoid_Arthritis_UKBS" = -1,
                       "HDL_Cholesterol_NMR" = 1,
                       "Depressive_Symptoms" = -1,
                       "Ankylosing_Spondylitis_UKBS" = -1,
                       "Type_1_Diabetes_UKBS" = -1,
                       "Adiponectin" = 1,
                       "BMI_Childhood" = -1,
                       "Triglycerides_NMR" = -1,
                       "CH2DB_NMR" = 0,
                       "LDL_Cholesterol_NMR" = -1,
                       "Diastolic_Blood_Pressure" = -1,
                       "Chronotype" = 0,
                       "Crohns_Disease_UKBS" = -1,
                       "Pubertal_Height_Male" = 0,
                       "Fasting_Glucose" = -1,
                       "Attention_Deficit_Hyperactivity_Disorder" = -1,
                       "Birth_Length" = 0,
                       "Bone_Mineral_Density" = 1,
                       "Systolic_Blood_Pressure" = -1,
                       "Celiac_Disease" = -1,
                       "Sleep_Duration" = 0,
                       "Smoker" = -1,
                       "Ulcerative_Colitis_UKBS" = -1,
                       "Deep_Venous_Thrombosis_UKB" = -1,
                       "Chronic_Kidney_Disease" = -1,
                       "Insomnia_In_Both_Sexes" = -1,
                       "Pubertal_Height_Female" = 0,
                       "Intracraneal_Volume" = 0,
                       "Epilepsy" = -1,
                       "Stroke" = -1,
                       "Migraine_UKBS" = -1,
                       "Insomnia_UKBS" = -1,
                       "Pulmonary_Embolism_UKB" = -1)
trait_goodness <- setNames(as.numeric(trait_goodness), names(trait_goodness))
mean(trait_goodness == 0)

#plotting params
plot_gene_names <- T
add_in_alltiss <- T
trait_patts <- c("triglycerides", "rheumatoid_arthritis", "high_cholesterol", "asthma_ukbs")[c(3,4)]
n_deg_sigtwas_intersect
# trait_patts <- c("body", "standing")
# trait_patts <- c("triglycerides")


cairo_pdf(paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/figures/fig6_protective_effect_opposing.pdf"), 
          width = 1400 / 72, height = 475 / 72 * length(trait_patts), family="Arial Unicode MS", pointsize = 12.5)
par(mfrow = c(length(trait_patts),1), xpd = NA, mar = c(4,3,4,17.75))
lwd <- 3
tissue_name_cex <- 1.2
gene_cex <- 0.75
gene_location <- 5.5
f_offset <- 7.15

my_tissue_abbr <- MotrpacRatTraining6moData::TISSUE_CODE_TO_ABBREV
my_tissue_abbr["all_tissues"] <- "ALL"
for(trait_patt in trait_patts){
  trait <- intersect(twas_with_hits, 
                     trait_categories$Tag[grep(trait_categories$new_Phenotype, pattern = trait_patt, ignore.case = T)])[1]
  
  for(sex in c("male", "female")){
    
    #retrieve data
    d <- apply(deg_sigtwas_proportion[, trait,,"p",sex], 2, as.numeric)
    dn <-  apply(deg_sigtwas_proportion[, trait,,"n",sex], 2, as.numeric)
    colnames(d) <- colnames(dn) <- colnames(deg_sigtwas_proportion[, trait,,"p",sex])
    rownames(d) <- rownames(dn) <- rownames(deg_sigtwas_proportion[, trait,,"p",sex])
    
    #add in summary of all tissues
    if(add_in_alltiss){
      at_prop <- apply(d * dn, 2, sum, na.rm = T)
      dn <- rbind(dn, all_tissues = apply(dn, 2, sum, na.rm = T))
      d <- rbind(d, all_tissues = at_prop / dn["all_tissues",])
      # d <- rbind(d, all_tissues = apply(d, 2, mean, na.rm = T))
    }
    
    dn <- dn[!apply(apply(d, 1, is.na), 2, all),]
    d <- d[!apply(apply(d, 1, is.na), 2, all),]
    d <- d[!is.nan(d[,"8w"]),]
    d[is.nan(d)] <- 0.5
    dg <- deg_sigtwas_proportion[, trait,,"genes",sex]
    dg <- lapply(dg[,"8w"], function(x) strsplit(x, " ~ ")[[1]])
    
    #munge gene names into composite strings
    dgm <- deg_sigtwas_proportion[, trait,,"genes","male"]
    dgm <- lapply(dgm[,"8w"], function(x) strsplit(x, " ~ ")[[1]])
    dgf <- deg_sigtwas_proportion[, trait,,"genes","female"]
    dgf <- lapply(dgf[,"8w"], function(x) strsplit(x, " ~ ")[[1]])
    dg <- lapply(setNames(names(dgm), names(dgm)), function(tiss){
      if(is.na(dgm[[tiss]][1])){return(NA)}
      mvs <- do.call(rbind, strsplit(dgm[[tiss]], " "))
      fvs <- do.call(rbind, strsplit(dgf[[tiss]], " "))
      if(nrow(mvs) > 1){
        apply(cbind(mvs[,2:3], mvs[,1], fvs[match(mvs[,1], fvs[,1]), 2:3]), 1, paste0, collapse = " ")  
      } else {
        mvs <- c(mvs)
        fvs <- c(fvs)
        paste0(c(mvs[2:3], mvs[1], fvs[2:3]), collapse = " ")
      }
    })
    
    
    #iterate through d to make identical lines parallel
    line_thickness <- lwd / 96 / par("pin")[2] * (par("usr")[4] - par("usr")[3])
    need_to_increment <- matrix(T, nrow = nrow(d), ncol = ncol(d), dimnames = list(rownames(d), colnames(d)))
    while(any(need_to_increment)){
      need_to_increment <- matrix(F, nrow = nrow(d), ncol = ncol(d), dimnames = list(rownames(d), colnames(d)))
      for(coi in 1:(ncol(d)-1)){
        unchanging_tissues <- rownames(unique(d[,c(coi,coi+1)]))
        need_to_increment[setdiff(rownames(d), unchanging_tissues),c(coi,coi+1)] <- T
      }
      d[need_to_increment] <- d[need_to_increment] + line_thickness * 1.05
    }
    
    #find coordinates to plot tissue names
    if(sex == "male"){
      ylocs_scale <- 100
      xylocs_tissue_names <- FField::FFieldPtRep(coords = cbind(rep(4 + 0.25, nrow(d)),
                                                                c((d[,"8w"] + 1:length(d[,"8w"])/1000 - length(d[,"8w"])/2000) * ylocs_scale)),
                                                 rep.fact = 20, adj.max = 15)
      xylocs_tissue_names$y <- xylocs_tissue_names$y / ylocs_scale
      # xylocs_tissue_names <- as.data.frame(cbind(x = rep(4 + 0.25, nrow(d)), y = c(d[,"8w"])))
    }
    
    #find coordinates to plot gene names
    if(plot_gene_names){
      xylocs_tissues_genes <- data.frame(gene = unlist(dg[rownames(xylocs_tissue_names)[order(xylocs_tissue_names$y, decreasing = T)]]))
      xylocs_tissues_genes$tissue <- rownames(xylocs_tissue_names)[sapply(rownames(xylocs_tissues_genes), function(x) 
        grep(strsplit(x, "-")[[1]][1], rownames(xylocs_tissue_names)))]
      xylocs_tissues_genes$x <- gene_location
      xylocs_tissues_genes$y <- xylocs_tissue_names$y[match(xylocs_tissues_genes$tissue, rownames(xylocs_tissue_names))]
      xylocs_tissues_genes$y <- xylocs_tissues_genes$y + 1:length(xylocs_tissues_genes$y)/1E4
      xylocs_tissues_genes$y <- seq(1.15, -0.15, length.out = nrow(xylocs_tissues_genes))
      # xylocs_tissues_genes_locs <- FField::FFieldPtRep(coords = cbind(xylocs_tissues_genes$x,
      #                                                                 xylocs_tissues_genes$y * 500),
      #                                                  rep.fact = 30, adj.max = 5, iter.max = 5E3)
      # xylocs_tissues_genes$y <- xylocs_tissues_genes_locs$y / 500
    }
    
    if(sex == "male"){
      
      plot(100,100,xlim = c(1,9.5), ylim = c(0,1), xpd=NA, 
           ylab = "", xlab = "", xaxt = "n", 
           yaxt = "n", bty="n", cex.lab = 1.5, cex.axis = 1.25)
      
      
      text("Timepoint", x = 2.5 + c(0,f_offset), y = -0.125, pos = 1, cex = 1.5)
      if(trait_goodness[trait_categories$new_Phenotype[trait_categories$Tag == trait]] == -1){
        vlab <- "Proportion Risk Enhancing Effects"
      } else if(trait_goodness[trait_categories$new_Phenotype[trait_categories$Tag == trait]] == 0){
        vlab <- "Proportion Positive Effects on GWAS Trait"
      } else if(trait_goodness[trait_categories$new_Phenotype[trait_categories$Tag == trait]] == 1){
        vlab <- "Proportion Risk Reducing Effects"
      }
      
      text(vlab, x = 0.5, y = 0.5, pos = 3, cex = 1.5, srt = 90, xpd = NA)
      text(vlab, x = 4 + f_offset + 0.6, y = 0.5, pos = 1, cex = 1.5, srt = 270, xpd = NA)
      
      
      #plot faded positive and negative regions
      rect(xl = 1, xr = 4, yb = 0.5, ytop = 1,
           col = grDevices::adjustcolor("red", 0.1), border = NA)
      rect(xl = 1, xr = 4, yb = 0, ytop = 0.5,
           col = grDevices::adjustcolor("blue", 0.1), border = NA)
      
      #trait and sex name
      text(paste0(trait_categories$new_Phenotype[trait_categories$Tag == trait], ","), col = 1, cex = 2, font = 1, pos = 3,
           x = par("usr")[2] * 0.185 + par("usr")[1] * 0.825, y = par("usr")[3] * 0 + par("usr")[4] * 1)
      text(c("\u2642", "\u2640")[c("male", "female") == sex], col = cols$Sex[sex], cex = 3.5, font = 2, pos = 3,
           x = par("usr")[2] * 0.185 + par("usr")[1] * 0.825 + strwidth(paste0(trait_categories$new_Phenotype[trait_categories$Tag == trait], ",     "), units = "use")*2/2, 
           y = par("usr")[3] * 0 + par("usr")[4] * 1)
      
      
      #horizontal axis
      segments(x0 = 1:4, x1 = 1:4, y0 = - 0.02, y1 = - 0.04, lwd = 2)
      segments(x0 = 1, x1 = 4, y0 = - 0.02, y1 = - 0.02, lwd = 2)
      text(x = 1:4, y = - 0.07, labels = paste0(2^(0:3), "w"), pos = , cex = 1.251)
      
      #vertical axis
      segments(x0 = 1 - 3 * 0.02, x1 = 1 - 3 * 0.04, y0 = 0:5/5, y1 = 0:5/5, lwd = 2)
      segments(x0 = 1 - 3 * 0.02, x1 = 1 - 3 * 0.02, y0 = 0, y1 = 1, lwd = 2)
      text(x = 1 - 3 * 0.035, y = 0:5/5, labels = 0:5/5, pos = 2, cex = 1.25)
      
      
      for(tissue in rownames(d)){
        
        if(all(is.na(d[tissue,]))){
          next()
        }
        
        lines(1:4, d[tissue,], lwd = lwd, col = ifelse(is.na(cols$Tissue[tissue]), 1, cols$Tissue[tissue]))
        
        #plot gene names
        maxwidth_genename <- max(strwidth(xylocs_tissues_genes$gene, units = "user", cex = gene_cex))
        if(plot_gene_names & tissue != "all_tissues"){
          tissue_genes <- xylocs_tissues_genes[xylocs_tissues_genes$tissue == tissue,]
          
          #break up gene names to allow for different colors
          genes_df <- do.call(rbind, strsplit(tissue_genes$gene, "(?=[\\(,)])", perl = TRUE))
          tissue_genes$gene <- tissue_genes$gene[order(genes_df[,3], decreasing = T)]
          genes_df <- genes_df[order(genes_df[,3], decreasing = T),]
          if(nrow(tissue_genes) == 1){
            genes_df <- t(as.matrix(genes_df))
          }
          max_ef <- 5
          colors_mat <- cbind(1, 
                              ifelse3(genes_df[,2] == "+", "red", "blue"), 
                              1, 
                              rev(viridisLite::cividis(n = 100))[min2(ceiling(as.numeric(genes_df[,4]) / max_ef * 100), 100)], 
                              1,
                              cols$Tissue[tissue], 
                              1, 
                              ifelse3(genes_df[,8] == "+", "red", "blue"), 
                              1, 
                              rev(viridisLite::cividis(n = 100))[min2(ceiling(as.numeric(genes_df[,10]) / max_ef * 100), 100)], 
                              1)
          text_indivcolor(labels_mat = genes_df, xloc = tissue_genes$x + (maxwidth_genename - strwidth(tissue_genes$gene, units = "user", cex = gene_cex)) / 2, y = tissue_genes$y, 
                          colors_mat = colors_mat, pos = 4, cex = gene_cex)
          
          
          # #all same color text
          # text(labels = tissue_genes$gene,
          #      x = tissue_genes$x,
          #      y = tissue_genes$y,
          #      cex = gene_cex, col = cols$Tissue[tissue], pos = 4)
          
          #plot connecting lines
          for(gene_i in 1:nrow(tissue_genes)){
            segments(x0 = xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue] + 
                       strwidth(paste0(my_tissue_abbr[tissue], " (", dn[tissue,"8w"], ")   "), cex = tissue_name_cex, units = "user"), 
                     y0 = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue],
                     x1 = tissue_genes$x[gene_i] + gene_cex*strwidth("  ", cex = tissue_name_cex, units = "user"),
                     y1 = tissue_genes$y[gene_i],
                     col = cols$Tissue[tissue], lty = 3)
            segments(x0 = tissue_genes$x[gene_i] + gene_cex*strwidth("  ", cex = tissue_name_cex, units = "user"), 
                     y0 = tissue_genes$y[gene_i],
                     x1 = tissue_genes$x[gene_i] + gene_cex*strwidth("  ", cex = tissue_name_cex, units = "user") + 
                       (maxwidth_genename - strwidth(tissue_genes$gene[gene_i], units = "user", cex = gene_cex)) / 2,
                     y1 = tissue_genes$y[gene_i],
                     col = cols$Tissue[tissue], lty = 3)
          }
        }
        
      }
      
      #plot tissue names
      for(tissue in rownames(d)){
        text(x = xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue], cex = tissue_name_cex,
             y = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue], 
             labels = paste0(my_tissue_abbr[tissue], " (", dn[tissue,"8w"], ")"), col = cols$Tissue[tissue], pos = 4)
        segments(x0 = 4, y0 = d[tissue,"8w"], 
                 x1 = xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue] + 0.075, 
                 y1 = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue],
                 col = ifelse(is.na(cols$Tissue[tissue]), 1, cols$Tissue[tissue]), lty = 3)
      }
      
      # legend(x = 1, y = 1.5, legend = tissue_names,
      #        col = cols$Tissue, lwd = 3, ncol = 4, cex = 1, border = NA, seg.len = 1, bg = NA, bty = "n", x.intersp = 0.25, text.width = 0.65)
      # segments(x0 = 0.9, y0 = 0.5, x1 = 4.1, y1 = 0.5, lty = 3, lwd = 3, col = "lightgrey")
    } else if(sex == "female") {
      
      #find desired offset for female plot
      # f_offset <- (gene_location - 4) * 2 + maxwidth_genename
      
      #plot faded positive and negative regions
      rect(xl = 1+f_offset, xr = 4+f_offset, yb = 0.5, ytop = 1,
           col = grDevices::adjustcolor("red", 0.1), border = NA)
      rect(xl = 1+f_offset, xr = 4+f_offset, yb = 0, ytop = 0.5,
           col = grDevices::adjustcolor("blue", 0.1), border = NA)
      
      #trait and sex name
      text(paste0(trait_categories$new_Phenotype[trait_categories$Tag == trait], ","), col = 1, cex = 2, font = 1, pos = 3,
           x = par("usr")[2] * 0.185 + par("usr")[1] * 0.825+f_offset, y = par("usr")[3] * 0 + par("usr")[4] * 1)
      text(c("\u2642", "\u2640")[c("male", "female") == sex], col = cols$Sex[sex], cex = 3.5, font = 2, pos = 3,
           x = par("usr")[2] * 0.185 + par("usr")[1] * 0.825 + strwidth(paste0(trait_categories$new_Phenotype[trait_categories$Tag == trait], ",     "), units = "use")*2/2+f_offset, 
           y = par("usr")[3] * 0 + par("usr")[4] * 1)
      
      #horizontal axis
      segments(x0 = 1:4+f_offset, x1 = 1:4+f_offset, y0 = - 0.02, y1 = - 0.04, lwd = 2)
      segments(x0 = 1+f_offset, x1 = 4+f_offset, y0 = - 0.02, y1 = - 0.02, lwd = 2)
      text(x = 1:4+f_offset, y = - 0.07, labels = rev(paste0(2^(0:3), "w")), pos = , cex = 1.251)
      
      #vertical axis
      segments(x0 = 1 + 3 * 0.02+f_offset+3, x1 = 1 + 3 * 0.04+f_offset+3, y0 = 0:5/5, y1 = 0:5/5, lwd = 2)
      segments(x0 = 1 + 3 * 0.02+f_offset+3, x1 = 1 + 3 * 0.02+f_offset+3, y0 = 0, y1 = 1, lwd = 2)
      text(x = 1 + 3 * 0.035+f_offset+3, y = 0:5/5, labels = 0:5/5, pos = 4, cex = 1.25)
      
      
      for(tissue in rownames(d)){
        
        if(all(is.na(d[tissue,]))){
          next()
        }
        
        lines(rev(1:4)+f_offset, d[tissue,], lwd = lwd, col = ifelse(is.na(cols$Tissue[tissue]), 1, cols$Tissue[tissue]))
        
      }
      
      #plot tissue names
      for(tissue in rownames(d)){
        text(x = -(xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue] - 4) + f_offset + 1 - 0.025, cex = tissue_name_cex,
             y = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue], 
             labels = paste0("(", dn[tissue,"8w"], ") ", my_tissue_abbr[tissue]), col = cols$Tissue[tissue], pos = 2)
        segments(x0 = 1+f_offset, y0 = d[tissue,"8w"], 
                 x1 = -(xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue] - 4) + f_offset + 1 - 0.075, 
                 y1 = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue],
                 col = ifelse(is.na(cols$Tissue[tissue]), 1, cols$Tissue[tissue]), lty = 3)
      }
      
      #plot connecting lines
      
      for(tissue in rownames(d)){
        if(tissue == "all_tissues"){next()}
        tissue_genes <- xylocs_tissues_genes[xylocs_tissues_genes$tissue == tissue,]
        
        for(gene_i in 1:nrow(tissue_genes)){
          segments(x0 = -(xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue] - 4) + f_offset + 1 - 0.075 - 
                     strwidth(paste0(my_tissue_abbr[tissue], " (", dn[tissue,"8w"], ")  "), cex = tissue_name_cex, units = "user"), 
                   y0 = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue],
                   x1 = tissue_genes$x[gene_i] + gene_cex*strwidth("  ", cex = tissue_name_cex, units = "user") + 
                     (maxwidth_genename - strwidth(tissue_genes$gene[gene_i], units = "user", cex = gene_cex)) + 
                     strwidth(tissue_genes$gene[gene_i], units = "user", cex = gene_cex),
                   y1 = tissue_genes$y[gene_i],
                   col = cols$Tissue[tissue], lty = 3)
          
          segments(x0 = tissue_genes$x[gene_i] + gene_cex*strwidth("  ", cex = tissue_name_cex, units = "user") + 
                     (maxwidth_genename - strwidth(tissue_genes$gene[gene_i], units = "user", cex = gene_cex)) + 
                     strwidth(tissue_genes$gene[gene_i], units = "user", cex = gene_cex),
                   y0 = tissue_genes$y[gene_i],
                   x1 = tissue_genes$x[gene_i] + gene_cex*strwidth("  ", cex = tissue_name_cex, units = "user") +
                     (maxwidth_genename - strwidth(tissue_genes$gene[gene_i], units = "user", cex = gene_cex)) / 2 +
                     strwidth(tissue_genes$gene[gene_i], units = "user", cex = gene_cex),
                   y1 = tissue_genes$y[gene_i],
                   col = cols$Tissue[tissue], lty = 3)
          
        }
      }
      
    }
  }
  
}

dev.off()

