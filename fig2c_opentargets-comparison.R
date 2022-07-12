#libraries
library(jsonlite)
library(data.table)
library(dplyr)
library(parallel)
# library(ghql)
library(sparklyr)
library(sparklyr.nested)
library(MotrpacBicQC)

#functions


#dload data
ot_dir <- "~/data/smontgom/opentargets/"
desired_folders <- c("associationByOverallIndirect", "associationByOverallDirect", "evidence", "targets", "diseases")
desired_types <- c("parquet", "json")
for(dt in desired_types){
  for(df in desired_folders){
    if(!dir.exists(paste0(ot_dir, dt, "/", df, "/"))){
      system(paste0("cd ~/data/smontgom/opentargets/", dt, "; ",
                    "wget --recursive --no-parent --no-host-directories --cut-dirs 8 ", 
                    "ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.04/output/etl/", dt, "/", df))
    }  
  }
}

#snag DE genesets
load("~/data/smontgom/node_metadata_list.RData")
ensembl_genes <- orig_ensembl_genes <- lapply(split(node_metadata_list$`8w`$human_ensembl_gene[!is.na(node_metadata_list$`8w`$human_ensembl_gene)], 
                       node_metadata_list$`8w`$tissue[!is.na(node_metadata_list$`8w`$human_ensembl_gene)]), unique)
symbol_map <- unique(node_metadata_list$`8w`[,c("human_gene_symbol", "human_ensembl_gene")])
all_genes <- unlist(ensembl_genes)
n_tissues_per_gene <- table(all_genes)
ensembl_genes$THREE <- names(n_tissues_per_gene)[n_tissues_per_gene > 2]

#paths
ot_dirs <- setNames(paste0(ot_dir, "json/", list.files(paste0(ot_dir, "json/")), "/"), list.files(paste0(ot_dir, "json/")))[1:2]
files <- lapply(ot_dirs, function(i){x <- paste0(i, list.files(i)); x[grep("json", x)]})

#read in OT files, then write them to an easier format
for(fileset in names(files)){
  print(fileset)
  if(!file.exists(paste0(ot_dir, fileset, ".csv"))){
    data <- lapply(files[[fileset]], function(indiv_file){
      stream_in(file(indiv_file))
    })
    data <- do.call(rbind, data)
    fwrite(data, file = paste0(ot_dir, fileset, ".csv"), sep = ",")
  }
}

overall_associations <- fread("~/data/smontgom/opentargets/associationByOverallIndirect.csv")
direct_associations <- fread("~/data/smontgom/opentargets/associationByOverallDirect.csv")

use_indirect <- F
if(use_indirect){
  associations_to_use <- overall_associations   
} else{
  associations_to_use <- direct_associations 
}

associations_to_use <- split(associations_to_use, associations_to_use$targetId)


## path to ClinVar (EVA) evidence dataset 
## directory stored on your local machine
#from tutorials https://platform-docs.opentargets.org/data-access/datasets#accessing-and-querying-datasets and https://therinspark.com/starting.html
evidencePath <- "~/data/smontgom/opentargets/parquet/evidence/sourceId=eva"
indAssPath <- "~/data/smontgom/opentargets/parquet/associationByOverallIndirect/"
diseasePath <- "~/data/smontgom/opentargets/parquet/diseases/"

## establish connection
Sys.setenv(JAVA_HOME="/Library/Java/JavaVirtualMachines/adoptopenjdk-8.jdk/Contents/Home")
conf <- spark_config()
conf$spark.executor.memory <- "32G"
conf$spark.executor.cores <- 2

sc <- spark_connect(master = "local", version = "2.0.0")
spark_web(sc)
## read evidence dataset
# evd <- spark_read_json(sc, path = evidencePath)
evd <- spark_read_parquet(sc, path = evidencePath)
ind_ass <- spark_read_parquet(sc, path = indAssPath)
dis <- spark_read_parquet(sc, path = diseasePath)


## Browse the disease schema
discol <- dis %>%
  sdf_schema() %>%
  lapply(function(x) do.call(tibble, x)) %>%
  bind_rows()

diseaseSelect <- dis %>%
  select(id,
         description,
         name)
# ,
#          ontology,
#          synonyms,
#          therapeuticAreas)

diseases <- diseaseSelect %>% collect()

diseases$name[match(associations_to_use[[ensembl_genes$ADRNL[1]]]$diseaseId, diseases$id)]

tissue_x_disease <- lapply(setNames(names(ensembl_genes),names(ensembl_genes)), function(tissue){
  print(tissue)
  # tissue <- "ADRNL"
  genes <- ensembl_genes[[tissue]]
  genes <- genes[!is.na(genes)]
  disvec <- setNames(rep(0, nrow(diseases)), diseases$id)
  dismat <- do.call(rbind, mclapply(setNames(genes,genes), function(gene){
    locdisvec <- disvec
    locdisvec[associations_to_use[[gene]]$diseaseId] <- associations_to_use[[gene]]$score
    locdisvec
  }, mc.cores = 8))
  nonzero_entries <- which(apply(dismat, 2, sum) > 1E-6)
  dismat <- dismat[,nonzero_entries]
  dismat
})

save(tissue_x_disease, file = paste0("~/data/smontgom/open-targets_tissue-x-disease_", 
                                       ifelse(use_indirect, "indirect", "direct"), 
                                       "-associations"))

#find sumstats for this table
thresh <- 0.8
tissues <- names(tissue_x_disease[1:15])
passing_gene_x_trait <- as.data.frame(do.call(rbind, sapply(tissues, function(tiss){
  x <- tissue_x_disease[[tiss]]
  inds <- which(x > thresh, arr.ind = T)
  return(cbind(tissue = tiss, gene = rownames(x)[inds[,1]], trait = colnames(x)[inds[,2]]))
})))
nrow(passing_gene_x_trait)
length(unique(passing_gene_x_trait$trait))
mean(sapply(split(passing_gene_x_trait, passing_gene_x_trait$tissue), function(x) length(unique(x$gene))))
easy_tissues <- c("BAT", "BLOOD", "SKM-GN", 'SKM-VL', 'WAT-SC')

sub_passing_gene_x_trait <- passing_gene_x_trait[passing_gene_x_trait$tissue %in% setdiff(tissues, easy_tissues),]
nrow(sub_passing_gene_x_trait)
mean(sapply(split(sub_passing_gene_x_trait, sub_passing_gene_x_trait$tissue), function(x) length(unique(x$gene))))




# max(do.call(rbind, associations_to_use)$score)
breakpoints <- 0:100/100
n_above <- sapply(names(tissue_x_disease), function(tissue) 
  log10((data.frame(cats = cut(tissue_x_disease[[tissue]][tissue_x_disease[[tissue]] > 1E-6], breaks=c(breakpoints, Inf)), ordered_result=TRUE) %>% 
           count(cats, .drop=FALSE) %>% arrange(desc(cats)) %>% mutate(cumfreq = cumsum(n)) %>% arrange(cats))$cumfreq))


n_above_at_least_1 <- sapply(names(tissue_x_disease), function(tissue) 
  log10((data.frame(cats = cut(apply(tissue_x_disease[[tissue]], 1, max)[apply(tissue_x_disease[[tissue]], 1, max) > 1E-6], breaks=c(breakpoints, Inf)), ordered_result=TRUE) %>% 
           count(cats, .drop=FALSE) %>% arrange(desc(cats)) %>% mutate(cumfreq = cumsum(n)) %>% arrange(cats))$cumfreq))


n_traits_above_at_least_1 <- sapply(names(tissue_x_disease), function(tissue) 
  log10((data.frame(cats = cut(apply(tissue_x_disease[[tissue]], 2, max)[apply(tissue_x_disease[[tissue]], 2, max) > 1E-6], breaks=c(breakpoints, Inf)), ordered_result=TRUE) %>% 
           count(cats, .drop=FALSE) %>% arrange(desc(cats)) %>% mutate(cumfreq = cumsum(n)) %>% arrange(cats))$cumfreq))


###### do the plotting #####
tissue_colors <- c(tissue_cols, THREE = "black")
main_title <- paste0("OpenTargets ", ifelse(use_indirect, "Overall / Indirect", "Direct"), " Associations")

par(mfrow = c(2,2), mar = c(4,6,4,2))

ylims <- range(n_above[is.finite(prop_above)])
plot(NA, NA, xlim = c(-0.15,1), ylim = ylims, xlab = "evidence score", 
     ylab = "# trait x gene pairs\nat or above given evidence score", yaxt = "n", main = main_title)
segments(y0 = 0:floor(max(n_above)), y1 = 0:floor(max(n_above)), x0 = -1, x1 = 2, lty = 3)
segments(y0 = 0:floor(max(n_above)), y1 = 0:floor(max(n_above)), x0 = -1, x1 = 2, lty = 3)
minticks <- log10(2:9) + rep(0:ceiling(par("usr")[4]), each = 8)
minticks <- minticks[minticks < par("usr")[4]]
segments(y0 = minticks, y1 = minticks, x0 = -1, x1 = 2, lty = 3, lwd = 0.5)
segments(y0 = minticks, x0 = par("usr")[1], y1 = minticks, x1 = par("usr")[1] - diff(par("usr")[1:2])/100, lwd = 0.75, xpd = NA)
segments(y0 = -1, x0 = 0:5/5, y1 = 1E9, x1 = 0:5/5, lwd = 0.75, lty = 3)

axis(2, at = 0:floor(max(n_above)), labels = latex2exp::TeX(paste0("$10^{", 0:floor(max(n_above)), "}$")), las = 1)
for(tissue in names(tissue_x_disease)){
  lines(breakpoints, n_above[,tissue], lwd = 2, col = tissue_colors[tissue])
}

#tissue labels
tiss_order <- colnames(n_above)[order(n_above[1,], decreasing = T)]
text(labels = tiss_order, col = tissue_colors[tiss_order], cex = 0.75, font = 2,
     pos = 2, x = -0.05, y = ylims[2] - cumsum(c(0,rep(strheight(tiss_order[1], units = "user"), length(tissue_x_disease) - 1))))
segments(x0 = -0.065, x1 = 0, y0 = ylims[2] - cumsum(c(0,rep(strheight(tiss_order[1], units = "user"), length(tissue_x_disease) - 1))),
         y1 = n_above[1,tiss_order], col = tissue_colors[tiss_order])

#now the prop plot
prop_above <- log10(10^n_above %*% (diag(1 / sapply(ensembl_genes, length))))
colnames(prop_above) <- colnames(n_above)

ylims <- range(prop_above[is.finite(prop_above)])

plot(NA, NA, xlim = c(-0.15,1), ylim = ylims, xlab = "evidence score", ylab = "average # traits per gene\nat or above given evidence score", 
     yaxt = "n", main = main_title)
segments(y0 = ceiling(ylims[1]):floor(ylims[2]), y1 = ceiling(ylims[1]):floor(ylims[2]), x0 = -1, x1 = 2, lty = 3)
segments(y0 =ceiling(ylims[1]):floor(ylims[2]), y1 = ceiling(ylims[1]):floor(ylims[2]), x0 = -1, x1 = 2, lty = 3)
minticks <- log10(2:9) + rep(floor(ylims[1]):ceiling(ylims[2]), each = 8)
minticks <- minticks[minticks < par("usr")[4] & minticks > par("usr")[3]]
segments(y0 = minticks, y1 = minticks, x0 = -1, x1 = 2, lty = 3, lwd = 0.5)
segments(y0 = minticks, x0 = par("usr")[1], y1 = minticks, x1 = par("usr")[1] - diff(par("usr")[1:2])/100, lwd = 0.75, xpd = NA)
segments(y0 = -1E9, x0 = 0:5/5, y1 = 1E9, x1 = 0:5/5, lwd = 0.75, lty = 3)

axis(2, at = ceiling(ylims[1]):floor(ylims[2]), labels = latex2exp::TeX(paste0("$10^{", ceiling(ylims[1]):floor(ylims[2]), "}$")), las = 1)
for(tissue in names(tissue_x_disease)){
  lines(breakpoints, prop_above[,tissue], lwd = 2, col = tissue_colors[tissue])
}

#tissue labels
tiss_order <- colnames(prop_above)[order(prop_above[1,], decreasing = T)]
text(labels = tiss_order, col = tissue_colors[tiss_order], cex = 0.75, font = 2,
     pos = 2, x = -0.05, y = ylims[2] - cumsum(c(0,rep(strheight(tiss_order[1], units = "user"), length(tissue_x_disease) - 1))))
segments(x0 = -0.065, x1 = 0, y0 = ylims[2] - cumsum(c(0,rep(strheight(tiss_order[1], units = "user"), length(tissue_x_disease) - 1))),
         y1 = prop_above[1,tiss_order], col = tissue_colors[tiss_order])

#now the ngenes with at least one plot
ylims <- range(n_above_at_least_1[is.finite(n_above_at_least_1)])

plot(NA, NA, xlim = c(-0.15,1), ylim = ylims, xlab = "evidence score", ylab = "# genes with at least one association\nat or above given evidence score", 
     yaxt = "n", main = main_title, xpd = NA)
segments(y0 = ceiling(ylims[1]):floor(ylims[2]), y1 = ceiling(ylims[1]):floor(ylims[2]), x0 = -1, x1 = 2, lty = 3)
segments(y0 =ceiling(ylims[1]):floor(ylims[2]), y1 = ceiling(ylims[1]):floor(ylims[2]), x0 = -1, x1 = 2, lty = 3)
minticks <- log10(2:9) + rep(floor(ylims[1]):ceiling(ylims[2]), each = 8)
minticks <- minticks[minticks < par("usr")[4] & minticks > par("usr")[3]]
segments(y0 = minticks, y1 = minticks, x0 = -1, x1 = 2, lty = 3, lwd = 0.5)
segments(y0 = minticks, x0 = par("usr")[1], y1 = minticks, x1 = par("usr")[1] - diff(par("usr")[1:2])/100, lwd = 0.75, xpd = NA)
segments(y0 = -1E9, x0 = 0:5/5, y1 = 1E9, x1 = 0:5/5, lwd = 0.75, lty = 3)

axis(2, at = ceiling(ylims[1]):floor(ylims[2]), labels = latex2exp::TeX(paste0("$10^{", ceiling(ylims[1]):floor(ylims[2]), "}$")), las = 1)
for(tissue in names(tissue_x_disease)){
  lines(breakpoints, n_above_at_least_1[,tissue], lwd = 2, col = tissue_colors[tissue])
}

#tissue labels
tiss_order <- colnames(n_above_at_least_1)[order(n_above_at_least_1[1,], decreasing = T)]
text(labels = tiss_order, col = tissue_colors[tiss_order], cex = 0.75, font = 2,
     pos = 2, x = -0.05, y = ylims[2] - cumsum(c(0,rep(strheight(tiss_order[1], units = "user"), length(tissue_x_disease) - 1))))
segments(x0 = -0.065, x1 = 0, y0 = ylims[2] - cumsum(c(0,rep(strheight(tiss_order[1], units = "user"), length(tissue_x_disease) - 1))),
         y1 = n_above_at_least_1[1,tiss_order], col = tissue_colors[tiss_order])



#now the ntraits with at least one plot
ylims <- range(n_traits_above_at_least_1[is.finite(n_traits_above_at_least_1)])

plot(NA, NA, xlim = c(-0.15,1), ylim = ylims, xlab = "evidence score", ylab = "# traits with at least one association\nat or above given evidence score", 
     yaxt = "n", main = main_title, xpd = NA)
segments(y0 = ceiling(ylims[1]):floor(ylims[2]), y1 = ceiling(ylims[1]):floor(ylims[2]), x0 = -1, x1 = 2, lty = 3)
segments(y0 =ceiling(ylims[1]):floor(ylims[2]), y1 = ceiling(ylims[1]):floor(ylims[2]), x0 = -1, x1 = 2, lty = 3)
minticks <- log10(2:9) + rep(floor(ylims[1]):ceiling(ylims[2]), each = 8)
minticks <- minticks[minticks < par("usr")[4] & minticks > par("usr")[3]]
segments(y0 = minticks, y1 = minticks, x0 = -1, x1 = 2, lty = 3, lwd = 0.5)
segments(y0 = minticks, x0 = par("usr")[1], y1 = minticks, x1 = par("usr")[1] - diff(par("usr")[1:2])/100, lwd = 0.75, xpd = NA)
segments(y0 = -1E9, x0 = 0:5/5, y1 = 1E9, x1 = 0:5/5, lwd = 0.75, lty = 3)

axis(2, at = ceiling(ylims[1]):floor(ylims[2]), labels = latex2exp::TeX(paste0("$10^{", ceiling(ylims[1]):floor(ylims[2]), "}$")), las = 1)
for(tissue in names(tissue_x_disease)){
  lines(breakpoints, n_traits_above_at_least_1[,tissue], lwd = 2, col = tissue_colors[tissue])
}

#tissue labels
tiss_order <- colnames(n_traits_above_at_least_1)[order(n_traits_above_at_least_1[1,], decreasing = T)]
text(labels = tiss_order, col = tissue_colors[tiss_order], cex = 0.75, font = 2,
     pos = 2, x = -0.05, y = ylims[2] - cumsum(c(0,rep(strheight(tiss_order[1], units = "user"), length(tissue_x_disease) - 1))))
segments(x0 = -0.065, x1 = 0, y0 = ylims[2] - cumsum(c(0,rep(strheight(tiss_order[1], units = "user"), length(tissue_x_disease) - 1))),
         y1 = n_traits_above_at_least_1[1,tiss_order], col = tissue_colors[tiss_order])



#### image of genes across tissues ####
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

#order the jaccard matrix
inds <- which(upper.tri(jacmat), T)
jacmat_df <- data.frame(tissue_1 = rownames(jacmat)[inds[,1]], tissue_2 = colnames(jacmat)[inds[,2]], jaccard_index = jacmat[inds])
jacmat_df[order(jacmat_df$jaccard_index, decreasing = T),]

#do some plotting
# dev.off()
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
text(labels = "Tissue Composition of Bin", x = 0, y = 0.5, pos = 3, cex = 1.5, xpd = NA, srt = 90)

#top plot
log10_count_ticks <- (0:floor(max(log10(nt2pg))))
log10_count_ylocs <- 1.1 + 0.9 * log10_count_ticks / log10(max(nt2pg))
segments(x0 = 0.425, x1 = 0.425, y = 1.1, y1 = 2, lwd = 2)
segments(x0 = 0.425, x1 = 0.35, y = log10_count_ylocs, y1 = log10_count_ylocs, lwd = 2)
text(x = 0.375, y = log10_count_ylocs, labels = latex2exp::TeX(paste0("$10^{", log10_count_ticks, "}$")), pos = 2, cex = 1.25, xpd = NA)
text(labels = "Number of Genes in Bin", x = 0, y = 1.525, pos = 3, cex = 1.5, xpd = NA, srt = 90)


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
       col = tissue_cols[colnames(jacmat)[j]])
  text(x = xl + (j-0.75) / nrow(jacmat) * (xr - xl), y = yt + 2.5 / ncol(jacmat) * (yt - yb),
       labels = colnames(jacmat)[j], pos = 4, cex = 0.75, xpd = NA, srt = 45)
  
  for(i in 1:(j-1)){
    rect(xleft = xl + (j-1) / nrow(jacmat) * (xr - xl), xright = xl + j / nrow(jacmat) * (xr - xl), 
         ybottom = yb + (nrow(jacmat)-i) / ncol(jacmat) * (yt - yb), ytop =  yb + (nrow(jacmat)-i+1) / ncol(jacmat) * (yt - yb),
         col = colgrad[floor(jacmat[i,j] * 100) + 1])
    if(j == ncol(jacmat)){
      rect(xleft = xr + 1 / xyrat / nrow(jacmat) * (xr - xl), xright = xr + (1 + 1 / xyrat) / nrow(jacmat) * (xr - xl), 
           ybottom = yb + (nrow(jacmat)-i) / ncol(jacmat) * (yt - yb), ytop =  yb + (nrow(jacmat)-i+1) / ncol(jacmat) * (yt - yb),
           col = tissue_cols[rownames(jacmat)[i]])
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



#### create summary data table ####
sumdf <- unique(ensembl_genes_df[,c("gene", "n_tissue")])
sumdf <- cbind(sumdf, do.call(rbind, lapply(sumdf$gene, function(gene){
  tissues <- ensembl_genes_df$tissue[ensembl_genes_df$gene == gene]
  tissues <- c(tissues, rep("", max(as.numeric(names(nt2pg))) - length(tissues)))
  names(tissues) <- paste0("tissue_", 1:length(nt2pg))
  tissues
  })
))
rownames(sumdf) <- NULL

sumdf <- cbind(sumdf, do.call(rbind, mclapply(sumdf$gene, function(gi){
  ass <- overall_associations[overall_associations$targetId == gi,]
  if(nrow(ass) == 0){return(c(NA, NA))}
  max_ass <- ass[which.max(ass$score), c("diseaseId", "score")]
  names(max_ass) <- paste0("indirect_association.", names(max_ass))
  unlist(max_ass)
}, mc.cores = 12)))

sumdf <- cbind(sumdf, do.call(rbind, mclapply(sumdf$gene, function(gi){
  ass <- direct_associations[direct_associations$targetId == gi,]
  if(nrow(ass) == 0){return(c(NA, NA))}
  max_ass <- ass[which.max(ass$score), c("diseaseId", "score")]
  names(max_ass) <- paste0("direct_association.", names(max_ass))
  unlist(max_ass)
}, mc.cores = 12)))

sumdf$indirect_association.disease <- paste0(diseases$name[match(sumdf$indirect_association.diseaseId, diseases$id)], ": ", diseases$description[match(sumdf$indirect_association.diseaseId, diseases$id)])
sumdf$direct_association.disease <- paste0(diseases$name[match(sumdf$direct_association.diseaseId, diseases$id)], ": ", diseases$description[match(sumdf$direct_association.diseaseId, diseases$id)])

sumdf$ensembl_ID <- sumdf$gene
sumdf$gene <- symbol_map$human_gene_symbol[match(sumdf$ensembl_ID, symbol_map$human_ensembl_gene)]

sumdf <- sumdf[with(sumdf, order(n_tissue, indirect_association.score, decreasing = T)),]
fwrite(sumdf, file = "~/data/smontgom/OpenTargets_MoTrPAC_Top_Associations.csv", sep = ",", col.names = T)
