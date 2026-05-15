library(data.table)
library(GMMAT)
library(Matrix)
library(rhdf5)
library(fastmatch)

source("glmmkin2randomvec.R")
source("get_changing_positions_and_h5File.R")


cmd <- commandArgs(T)

pheno.file <- cmd[1] #"toy_pheno.txt"
global.kinship.file <- cmd[2] #"toy_global_IBD.RData"
IBD_file <- cmd[3] #"toy_chr20_IBD.txt.gz"
outfile_prefix <- cmd[4] #"toy"
seed = as.numeric(cmd[5]) #12345
id1.col = as.numeric(cmd[6]) #1
id2.col = as.numeric(cmd[7]) #2
chr.col = as.numeric(cmd[8]) #3
start.col = as.numeric(cmd[9]) #4
end.col = as.numeric(cmd[10]) #5
IBD_sorted = as.numeric(cmd[11]) #0


N.randomvec <- 100
kin <- get(load(global.kinship.file))
pheno <- fread(pheno.file, data.table = FALSE)
pheno <- subset(pheno, id %in% rownames(kin))
pheno <- pheno[!apply(is.na(pheno), 1, any), ]
cat("Effective sample size:", nrow(pheno), "\n")

obj <- try(glmmkin(as.formula(paste0("pheno ~ covar")), data = pheno, id = "id", kins = kin, family = gaussian(link = "identity"), verbose = T))
if(class(obj) == "try-error") stop("Error: cannot fit the LMM")
cat("Number of individuals in glmmkin:", length(obj$id_include), "\n")
cat("Null model: theta =", obj$theta, "\n")
if(obj$theta[1]==0) stop("Error: residual variance component estimate is 0")
set.seed(seed)
match.idx1 <- match(pheno$id, rownames(kin))
match.idx2 <- match(pheno$id, colnames(kin))
tmpibd <- kin[match.idx1, match.idx2]
IBD.chol <- chol(tmpibd)
#Sigma <- obj$theta[1]*Diagonal(nrow(tmpibd))+obj$theta[2]*tmpibd
obj1 <- glmmkin2randomvec(obj, Z = list(t(IBD.chol)), N.randomvec = N.randomvec, robust = FALSE)
rm(obj, tmpibd)
saveRDS(obj1, file = paste0(outfile_prefix, "_glmmkin2randomvec.rds"))


get_changing_positions_and_h5File(obj1, IBD_file, outfile_prefix, id1.col, id2.col, chr.col, start.col, end.col, IBD_sorted)























