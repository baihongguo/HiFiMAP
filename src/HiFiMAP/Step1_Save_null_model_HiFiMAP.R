

library(data.table)
library(GMMAT)
library(Matrix)
library(fastmatch)

# Source the random vector conversion helper
# (Ensure this path is correct for your environment)
source("src/HiFiMAP/glmmkin2randomvec.R")

# -------------------------------------------------------------------------
# 1. PARSE ARGUMENTS
# -------------------------------------------------------------------------
# Note: IBD file arguments are REMOVED. 
# This script only needs Phenotype and Kinship data.

cmd <- commandArgs(trailingOnly = TRUE)

if (length(cmd) < 4) {
  stop("Usage: Rscript Step1_Preprocess_HiFiMAP.R <pheno_file> <kinship_file> <outfile_prefix> <seed>")
}

pheno.file <- cmd[1] 
global.kinship.file <- cmd[2] 
outfile_prefix <- cmd[3] 
seed <- as.numeric(cmd[4]) 

cat("------------------------------------------------\n")
cat("Step 1: GLMM Null Model Fitting\n")
cat("Phenotype File:", pheno.file, "\n")
cat("Kinship File:  ", global.kinship.file, "\n")
cat("Output Prefix: ", outfile_prefix, "\n")
cat("Seed:          ", seed, "\n")
cat("------------------------------------------------\n")

# -------------------------------------------------------------------------
# 2. LOAD AND ALIGN DATA
# -------------------------------------------------------------------------

# Load Kinship (Expects an R matrix/GRM)
# Adjust 'get(load(...))' if your file format differs
kin_obj <- load(global.kinship.file)
if(is.character(kin_obj)) {
  kin <- get(kin_obj)
} else {
  stop("Could not load kinship object")
}

# Load Phenotypes
pheno <- fread(pheno.file, data.table = FALSE)

# Filter Phenotypes to match Kinship samples
pheno <- subset(pheno, id %in% rownames(kin))
pheno <- pheno[!apply(is.na(pheno), 1, any), ]

cat("Effective sample size:", nrow(pheno), "\n")

# -------------------------------------------------------------------------
# 3. FIT NULL GLMM (GMMAT)
# -------------------------------------------------------------------------
# Adjust the formula "Phenotype ~ Age + Sex" as needed for your real data

obj <- try(glmmkin(as.formula("pheno ~ covar"), 
                   data = pheno, 
                   id = "id", 
                   kins = kin, 
                   family = gaussian(link = "identity"), 
                   verbose = TRUE))

if(inherits(obj, "try-error")) stop("Error: cannot fit the LMM")

cat("Number of individuals in glmmkin:", length(obj$id_include), "\n")
cat("Null model: theta =", obj$theta, "\n")

if(obj$theta[1] == 0) stop("Error: residual variance component estimate is 0")

# -------------------------------------------------------------------------
# 4. PREPARE RANDOM VECTORS
# -------------------------------------------------------------------------
set.seed(seed)

# Ensure Sigma_i (Inverse Covariance) exists
if (is.null(obj$Sigma_i)) {
  if (!is.null(obj$P)) {
    cat("Calculating Sigma_i from P...\n")
    obj$Sigma_i <- tryCatch(
      solve(obj$P),
      error = function(e) {
        message("obj$P is singular, applying nearPD()...")
        as.matrix(nearPD(obj$P)$mat)
      }
    )
  } else {
    stop("Error: obj$Sigma_i is NULL, and obj$P is also missing.")
  }
}

# Prepare Cholesky of Kinship for Random Vector generation
match.idx1 <- match(pheno$id, rownames(kin))
match.idx2 <- match(pheno$id, colnames(kin))
tmpibd <- kin[match.idx1, match.idx2]
IBD.chol <- chol(tmpibd)

# Generate Random Vectors
N.randomvec <- 100
cat("Generating", N.randomvec, "random vectors...\n")

obj1 <- glmmkin2randomvec(obj, 
                          Z = list(t(IBD.chol)), 
                          N.randomvec = N.randomvec, 
                          robust = FALSE)

# Clean up large objects before saving
rm(obj, tmpibd, kin, IBD.chol)

# -------------------------------------------------------------------------
# 5. SAVE OUTPUT
# -------------------------------------------------------------------------
output_file <- paste0(outfile_prefix, "_glmmkin2randomvec.rds")
cat("Saving processed object to:", output_file, "\n")
saveRDS(obj1, file = output_file)

cat("Step 1 Complete.\n")












