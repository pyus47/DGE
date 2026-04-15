packages <- c("BiocManager", "tibble")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

bio_packages <- c("DESeq2", "edgeR", "limma")
missing_bio <- bio_packages[!vapply(bio_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_bio) > 0) {
  BiocManager::install(missing_bio, ask = FALSE, update = FALSE)
}

message("R packages are ready.")
