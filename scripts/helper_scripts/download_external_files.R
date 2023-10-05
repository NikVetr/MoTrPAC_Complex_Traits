# files_map <- rbind(read.csv(file = "~/data/smontgom/pass1b_companion_paper_external_filenames_map.csv")[,c("old_path", "new_dir")],
#                    read.csv(file = "~/data/smontgom/pass1b_companion_paper_external_extra_filenames_map.csv")[,c("old_path", "new_dir")])
# files_map <- rbind(files_map, 
#                    cbind(old_path = "~/repos/MoTrPAC_Complex_Traits/", 
#                          new_dir = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/"),
#                    cbind(old_path = "/Volumes/SSD500GB/gtex-pipeline/expression_genotypes/GTEx_v8*", 
#                          new_dir = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/"),
#                    cbind(old_path = "/Volumes/SSD500GB/gtex-pipeline/expression_data/GTEx_Analysis_v8_*_expected_count.gct.gz", 
#                          new_dir = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/"),
#                    cbind(old_path = "/Volumes/SSD500GB/gtex-pipeline/expression_data/GTEx_Analysis_v8_*_tpm.gct.gz", 
#                          new_dir = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/"))
getDirSizeSystem <- function(dirPath) {
  cmd <- paste("du -sk", dirPath)  # Use '-sk' for size in kilobytes
  output <- system(cmd, intern = TRUE)
  total_size_kb <- as.numeric(strsplit(output, "\t")[[1]][1])
  return(total_size_kb / 1024^2)
}

files <- list.files("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/", recursive = F, full.names = T)
files <- cbind(files, file.info(files)[,c("size", "isdir")])
files$size <- files$size / 1024^3
files$size[files$isdir] <- sapply(files$files[files$isdir], getDirSizeSystem)
files[files$isdir,]

files <- files[order(files$isdir),]
files$files <- gsub("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external//", "", files$files)
# write.csv(x = files[,c(1,3,2)], "~/external_files_motrpac_companion.csv", row.names = F)

read.csv("~/external_files_motrpac_companion.csv", row.names = F)
dest_file <- 
download.file("https://drive.google.com/uc?export=download&id=1COG3UXpdMtfDgF9QQ3LCE4UrgLQQRCNK", destfile = "~/RSID_POS_MAP.txt", mode = "wb")
