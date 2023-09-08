file_lines <- readLines("~/Downloads/motrpac_study_group.txt", warn = F)

institution <- NULL
author_institution_list <- list()

for (i in seq_along(file_lines)) {
  line <- file_lines[i]
  if (i %% 3 == 1) {
    institution <- line
  } else if (i %% 3 == 2) {
    authors <- strsplit(line, ", ")[[1]]
    author_institution_list[[institution]] <- trimws(authors)
  }
}

names(author_institution_list) <- paste0(names(author_institution_list), ", USA")
d <- data.frame(authors = setNames(unlist(author_institution_list), NULL),
                institution = unlist(sapply(names(author_institution_list), 
                                            function(i) rep(i, length(author_institution_list[[i]]))))
)
rownames(d) <- NULL
d$surname <- sapply(d$authors, function(i) tail(strsplit(i, " ")[[1]], 1))
d$index <- as.integer(rank(d$surname))

institutions_ordered <- unique(d$institution[order(d$surname)])
institutions_ordered <- c("Stanford University, Stanford, CA, USA", 
                          institutions_ordered[institutions_ordered != 
                                                 "Stanford University, Stanford, CA, USA"])

d$institution_index <- match(d$institution, institutions_ordered)
d <- d[order(d$surname),]


#latex
latex_author_list <- "\\textbf{\\large MoTrPAC Study Group}\n\n"
latex_affiliation_list <- ""

for (author in d$authors) {
  affiliations <- d$institution_index[author == d$authors]
  latex_author_list <- paste0(latex_author_list, author, "$^{", affiliations, "}$, ")
}

latex_author_list <- substr(latex_author_list, 1, nchar(latex_author_list) - 2)  # Remove the last comma and space

# Generate Affiliations
for (i in 2:length(institutions_ordered)) {
  latex_affiliation_list <- paste0(latex_affiliation_list, "$^{", i, "}$", institutions_ordered[i], ", ")
}

latex_affiliation_list <- substr(latex_affiliation_list, 1, nchar(latex_affiliation_list) - 2)  # Remove the last comma and space

# Complete LaTeX sections
latex_complete <- paste0(latex_author_list, "\n\n\\textbf{\\small ", latex_affiliation_list, "}\n")

# Print the LaTeX code
cat(latex_complete)
