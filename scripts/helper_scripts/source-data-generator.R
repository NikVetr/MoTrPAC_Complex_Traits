rm(list = ls())

#initialize overall figure parameters
figure_script_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/figures/"
figure_scripts <- list.files(figure_script_dir, full.names = T)

#retrieve figure script metadata and make source data directory
figure_index <- 9
figure_script <- figure_scripts[figure_index]
figure_name <- gsub(pattern = ".R", replacement = "", tail(strsplit(figure_script, "/")[[1]], 1))
figures_data_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/figures/"
figure_data_dir <- paste0(figures_data_dir, figure_name, "/")
if(!dir.exists(figure_data_dir)){dir.create(figure_data_dir)}

#id last call to preprocessing scripts and add in a check to find all objects in memory
figure_script_text <- readLines(figure_script)
source_calls <- which(startsWith(figure_script_text, "source"))
last_source_call <- max(source_calls)
figure_script_text_preprocessing <- c(figure_script_text[1:last_source_call], "preprocessed_memory_objects <- ls()")

#run this and retrieve memory object metadata
eval(parse(text = paste0(figure_script_text_preprocessing, collapse = "\n")))
object_info <- do.call(rbind, lapply(preprocessed_memory_objects, function(objname){
  obj <- get(objname)
  data.frame(name = objname,
             size = as.numeric(object.size(obj)), #in bytes
             class = class(obj),
             len = length(obj))[1,]
}))
object_uses <- lapply(setNames(preprocessed_memory_objects,preprocessed_memory_objects), function(objname){
  figure_script_text[grepl(objname, figure_script_text)]})
object_info <- object_info[!(object_info$class %in% "function"),]
object_info$name_subset <- sapply(object_info$name, function(ni) sum(grepl(ni, object_info$name)) > 1)
object_info$appearance_count <- sapply(object_info$name, function(ni) length(object_uses[[ni]]))

#id which of these objects are actually used in the plotting script
required_objects <- object_info[sapply(object_info$name, function(objname) any(grepl(objname, figure_script_text))),]
required_objects[order(required_objects$size, decreasing = T),]

#check up on specific objects
object_uses["d"]
# setdiff(object_uses["ldsc_results"], object_uses["ldsc_results_sub"])

#save to disk
if(sum(required_objects$size) < 1E6){
  source_data_path <- paste0(figures_data_dir, figure_name, ".RData")
  save(list = required_objects$name, 
       file = source_data_path,
       envir = .GlobalEnv)
  for(i in required_objects$name){
    print(i)
    save(list = i,
         file = paste0(figure_data_dir, i, ".RData"),
         envir = .GlobalEnv)
  }
  load(source_data_path)
}


