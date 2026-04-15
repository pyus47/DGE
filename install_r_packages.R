# 1. Define a local directory for R packages
local_lib <- file.path(getwd(), "R_libs")
if (!dir.exists(local_lib)) dir.create(local_lib)

# 2. Tell R to use this folder for all installations and loading
.libPaths(c(local_lib, .libPaths()))

# 3. Standard installation logic
packages <- c("BiocManager", "tibble")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg,
      repos = "https://cloud.r-project.org",
      lib = local_lib
    ) # Force local install
  }
}

bio_packages <- c("DESeq2", "edgeR", "limma")
# Check if they exist in our local library
missing_bio <- bio_packages[!vapply(bio_packages, function(x) requireNamespace(x, lib.loc = local_lib, quietly = TRUE), logical(1))]

if (length(missing_bio) > 0) {
  BiocManager::install(missing_bio,
    lib = local_lib,
    ask = FALSE,
    update = FALSE
  )
}

message("R packages are ready in: ", local_lib)
