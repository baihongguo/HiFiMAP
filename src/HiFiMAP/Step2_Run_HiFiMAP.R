library(data.table)
library(Matrix)
library(methods)
library(Rcpp)

source("src/HiFiMAP/pval_caculation_functions.R")

# --- Helper: Load dK ---
read_diff_as_dK <- function(diff_file, n_hap) {
  if (!file.exists(diff_file)) return(NULL)
  dat <- scan(diff_file, what = list(integer(), integer(), integer()), quiet = TRUE)
  ops <- dat[[1]]; rows <- dat[[2]] + 1; cols <- dat[[3]] + 1
  if (length(ops) == 0) return(NULL)
  vals <- ifelse(ops == 1, 1, -1)
  dK <- sparseMatrix(i = rows, j = cols, x = vals, dims = c(n_hap, n_hap))
  return(forceSymmetric(dK))
}

# --- MAIN ---
cmd <- commandArgs(trailingOnly = TRUE)
chr            <- as.numeric(cmd[1])
outfile_prefix <- cmd[2]
threads        <- as.numeric(cmd[3])
glmm_rds_path  <- cmd[4]
ibd_dir        <- cmd[5] 
start_idx      <- as.numeric(cmd[6])
end_idx        <- as.numeric(cmd[7])

Sys.setenv(MKL_NUM_THREADS = threads)
set.seed(seed)

t1 <- proc.time()

# 1. Load Data
obj1 <- readRDS(glmm_rds_path)
residuals <- obj1$scaled.residuals
random_vecs <- obj1$random.vectors
N <- length(obj1$id_include)
ncovars <- obj1$ncovar
N_randomvec <- ncol(random_vecs)

Q_offset <- 2 * sum(residuals^2)
eigenmat_offset <- 2 * crossprod(random_vecs) - Q_offset/(N-ncovars) * obj1$RSigma_R

# 2. Load Maps
sites <- fread(file.path(ibd_dir, "sites.txt"))
samples_vcf <- fread(file.path(ibd_dir, "samples.txt"), header=F)$V1
N_vcf <- length(samples_vcf)

match_idx <- match(obj1$id_include, samples_vcf)
rows <- rep(1:N, each=2)
cols <- c(rbind(2*match_idx - 1, 2*match_idx))
P_proj <- sparseMatrix(i = rows, j = cols, x = 1, dims = c(N, N_vcf * 2))

# ==========================================
# 3. C++ STREAMING FALLBACK ENGINE
# ==========================================
Rcpp::sourceCpp(code = "
#include <Rcpp.h>
#include <fstream>
#include <string>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
List stream_init_stats(std::string mtx_file, IntegerVector match_idx, int N_vcf,
                       NumericVector residuals, NumericMatrix random_vecs) {
                       
    int N = residuals.size();
    int N_randomvec = random_vecs.ncol();

    std::vector<int> hap_to_subj(2 * N_vcf + 1, 0);
    for (int i = 0; i < N; ++i) {
        int vcf_idx = match_idx[i];
        if (IntegerVector::is_na(vcf_idx)) continue;
        hap_to_subj[2 * vcf_idx - 1] = i + 1; 
        hap_to_subj[2 * vcf_idx]     = i + 1;
    }

    double Q_curr = 0.0;
    NumericMatrix Eigen_curr(N_randomvec, N_randomvec);
    double sum_K = 0.0;

    FILE* fp = fopen(mtx_file.c_str(), \"r\");
    if (!fp) stop(\"Cannot open MTX file.\");

    char line[256];
    bool is_symmetric = false;
    
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '%') {
            std::string s(line);
            if (s.find(\"symmetric\") != std::string::npos) {
                is_symmetric = true;
            }
        } else {
            break; 
        }
    }

    int h1, h2, v_int;
    
    while (fscanf(fp, \"%d %d %d\", &h1, &h2, &v_int) == 3) {
        
        int s1 = hap_to_subj[h1];
        int s2 = hap_to_subj[h2];
        
        if (s1 != 0 && s2 != 0) {
            double val = (s1 == s2) ? (v_int / 2.0) : (double)v_int;
            sum_K += val;
            Q_curr += residuals[s1 - 1] * residuals[s2 - 1] * val;
            for (int a = 0; a < N_randomvec; ++a) {
                double v1 = random_vecs(s1 - 1, a);
                for (int b = 0; b < N_randomvec; ++b) {
                    Eigen_curr(a, b) += v1 * random_vecs(s2 - 1, b) * val;
                }
            }
        }
        
        if (is_symmetric && h1 != h2) {
            int s1_mir = hap_to_subj[h2];
            int s2_mir = hap_to_subj[h1];
            
            if (s1_mir != 0 && s2_mir != 0) {
                double val = (s1_mir == s2_mir) ? (v_int / 2.0) : (double)v_int;
                sum_K += val;
                Q_curr += residuals[s1_mir - 1] * residuals[s2_mir - 1] * val;
                for (int a = 0; a < N_randomvec; ++a) {
                    double v1 = random_vecs(s1_mir - 1, a);
                    for (int b = 0; b < N_randomvec; ++b) {
                        Eigen_curr(a, b) += v1 * random_vecs(s2_mir - 1, b) * val;
                    }
                }
            }
        }
    }
    fclose(fp);

    double n_ibd_curr = std::round(sum_K / 2.0);

    return List::create(
        Named(\"Q_curr\") = Q_curr,
        Named(\"Eigen_curr\") = Eigen_curr,
        Named(\"n_ibd_curr\") = n_ibd_curr
    );
}

// [[Rcpp::export]]
List stream_diff_stats(std::string diff_file, IntegerVector match_idx, int N_vcf,
                       NumericVector residuals, NumericMatrix random_vecs,
                       double Q_curr, NumericMatrix Eigen_curr, double n_ibd_curr) {
                       
    int N = residuals.size();
    int N_randomvec = random_vecs.ncol();
    NumericMatrix Eigen_new = clone(Eigen_curr);

    std::vector<int> hap_to_subj(2 * N_vcf + 1, 0);
    for (int i = 0; i < N; ++i) {
        int vcf_idx = match_idx[i];
        if (IntegerVector::is_na(vcf_idx)) continue;
        hap_to_subj[2 * vcf_idx - 1] = i + 1; 
        hap_to_subj[2 * vcf_idx]     = i + 1;
    }

    FILE* fp = fopen(diff_file.c_str(), \"r\");
    if (!fp) stop(\"Cannot open DIFF file.\");

    int op, h1_0, h2_0;
    double sum_dK = 0.0;
    
    while (fscanf(fp, \"%d %d %d\", &op, &h1_0, &h2_0) == 3) {
        
        int h1 = h1_0 + 1;
        int h2 = h2_0 + 1;
        double v_int = (op == 1) ? 1.0 : -1.0;
        
        int s1 = hap_to_subj[h1];
        int s2 = hap_to_subj[h2];
        
        if (s1 != 0 && s2 != 0) {
            double val = (s1 == s2) ? (v_int / 2.0) : v_int;
            sum_dK += val;
            Q_curr += residuals[s1 - 1] * residuals[s2 - 1] * val;
            for (int a = 0; a < N_randomvec; ++a) {
                double v1 = random_vecs(s1 - 1, a);
                for (int b = 0; b < N_randomvec; ++b) {
                    Eigen_new(a, b) += v1 * random_vecs(s2 - 1, b) * val;
                }
            }
        }
        
        if (h1 != h2) {
            int s1_mir = hap_to_subj[h2];
            int s2_mir = hap_to_subj[h1];
            
            if (s1_mir != 0 && s2_mir != 0) {
                double val = (s1_mir == s2_mir) ? (v_int / 2.0) : v_int;
                sum_dK += val;
                Q_curr += residuals[s1_mir - 1] * residuals[s2_mir - 1] * val;
                for (int a = 0; a < N_randomvec; ++a) {
                    double v1 = random_vecs(s1_mir - 1, a);
                    for (int b = 0; b < N_randomvec; ++b) {
                        Eigen_new(a, b) += v1 * random_vecs(s2_mir - 1, b) * val;
                    }
                }
            }
        }
    }
    fclose(fp);

    n_ibd_curr += std::round(sum_dK / 2.0);

    return List::create(
        Named(\"Q_curr\") = Q_curr,
        Named(\"Eigen_curr\") = Eigen_new,
        Named(\"n_ibd_curr\") = n_ibd_curr
    );
}
")

# ==========================================
# 4. INITIALIZATION (Hybrid)
# ==========================================
mat_file <- file.path(ibd_dir, paste0("ibd_mat_", start_idx, ".mtx"))
if (!file.exists(mat_file)) stop(paste("Missing:", mat_file))

# Attempt fast R matrix logic first
K_subj <- tryCatch({
  K_hap_start <- readMM(mat_file)
  tmp <- P_proj %*% K_hap_start %*% t(P_proj)
  diag(tmp) <- diag(tmp) / 2
  tmp
}, error = function(e) {
  NULL # Return NULL if matrix exceeds 2^31-1 limit
})

if (!is.null(K_subj)) {
  # FAST PATH
  Q_curr <- as.numeric(crossprod(residuals, crossprod(K_subj, residuals)))
  Eigen_curr <- crossprod(random_vecs, crossprod(K_subj, random_vecs))
  n_ibd_curr <- round(sum(K_subj) / 2)
  rm(K_subj); gc(verbose = FALSE)
} else {
  # FALLBACK PATH
  cat("  [Fallback] Init matrix too large. Streaming via C++...\n")
  init_stats <- stream_init_stats(mat_file, match_idx, N_vcf, 
                                  as.numeric(residuals), as.matrix(random_vecs))
  Q_curr <- init_stats$Q_curr
  Eigen_curr <- init_stats$Eigen_curr
  n_ibd_curr <- init_stats$n_ibd_curr
}

# ==========================================
# 5. SCAN LOOP (Hybrid)
# ==========================================
res_file <- paste0(outfile_prefix, ".txt")
cat("chr\tpos\tn.ibd.segs\tp.value\n", file = res_file)

calc_pval <- function(E, Q) {
  lambda <- zapsmall(eigen(E + eigenmat_offset - Q/(N-ncovars)*obj1$RSigma_R, 
                           only.values=T, symmetric=T)$values / N_randomvec)
  pv <- try(.Q_pval(0, lambda), silent=T)
  if (inherits(pv, "try-error")) pv <- NA
  if (all(lambda < 0)) pv <- 0
  return(pv)
}

pval_curr <- calc_pval(Eigen_curr, Q_curr)
cat(chr, sites$BP[start_idx+1], n_ibd_curr, pval_curr, "\n", file=res_file, sep="\t", append=T)

loop_indices <- (start_idx + 1):end_idx

for (i in loop_indices) {
  diff_file <- file.path(ibd_dir, paste0("delta_", i, ".diff"))
  
  if (file.exists(diff_file) && file.info(diff_file)$size > 0) {
    
    dK_hap <- read_diff_as_dK(diff_file, N_vcf * 2)
    
    if (!is.null(dK_hap)) {
      
      # Attempt fast R matrix logic
      dK_subj <- tryCatch({
        tmp <- P_proj %*% dK_hap %*% t(P_proj)
        diag(tmp) <- diag(tmp) / 2
        tmp
      }, error = function(e) {
        NULL # Matrix limit exceeded
      })
      
      if (!is.null(dK_subj)) {
        # FAST PATH (Intel MKL)
        if (nnzero(dK_subj) > 0) {
          Q_curr <- Q_curr + as.numeric(crossprod(residuals, dK_subj %*% residuals))
          Eigen_curr <- Eigen_curr + crossprod(random_vecs, dK_subj %*% random_vecs)
          n_ibd_curr <- n_ibd_curr + round(sum(dK_subj) / 2)
          pval_curr <- calc_pval(Eigen_curr, Q_curr)
        }
      } else {
        # FALLBACK PATH (C++ Streaming)
        cat("  [Fallback] Site", i, "diff too large. Streaming via C++...\n")
        updated_stats <- stream_diff_stats(diff_file, match_idx, N_vcf, 
                                           as.numeric(residuals), as.matrix(random_vecs), 
                                           as.numeric(Q_curr), as.matrix(Eigen_curr), as.numeric(n_ibd_curr))
        
        Q_curr <- updated_stats$Q_curr
        Eigen_curr <- updated_stats$Eigen_curr
        n_ibd_curr <- updated_stats$n_ibd_curr
        pval_curr <- calc_pval(Eigen_curr, Q_curr)
      }
    }
  }
  
  cat(chr, sites$BP[i+1], n_ibd_curr, pval_curr, "\n", file=res_file, sep="\t", append=T)
}

t2 <- proc.time()
time_file <- paste0(outfile_prefix, ".time")
write.table(t(c(t2 - t1)), file = time_file, quote = FALSE, row.names = FALSE, col.names = TRUE)
