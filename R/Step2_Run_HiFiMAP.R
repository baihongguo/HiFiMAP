library(data.table)
library(Matrix)
library(rhdf5)
source("/HGCNT95FS/ADDLIE/work/Test/FiMAP/simulation/toy/FiMAP_process_chunk.R")
source("/HGCNT95FS/ADDLIE/work/Test/FiMAP/simulation/toy/pval_caculation_functions.R")

cmd <- commandArgs(T)
chr <- as.numeric(cmd[1]) # 1-22
cat("Running chromosome", chr, "...\n")
seed <- as.numeric(cmd[2]) # 12345
cat("Using seed", seed, "...\n")
outfile.prefix <- cmd[3]
n_chunks = as.numeric(cmd[4])
chunk_inx = as.numeric(cmd[5])
threads = as.numeric(cmd[6])

Sys.setenv(MKL_NUM_THREADS=threads)
IBD_h5file = paste0("/HGCNT95FS/ADDLIE/work/Test/FiMAP/simulation/toy/toy_ibd_chr20.h5")
pos_file = paste0("/HGCNT95FS/ADDLIE/work/Test/FiMAP/simulation/toy/toy_ibd_chr20_CPos_index.txt")
obj1 <- readRDS(paste0("/HGCNT95FS/ADDLIE/work/Test/FiMAP/simulation/toy/toy_glmmkin2randomvec.rds"))

set.seed(seed)


if (n_chunks == 1) {
changing_pos = fread(pos_file)
chunk <- na.omit(changing_pos[[3]])
chunk_first_index = 1
start <- na.omit(changing_pos[[1]])
start.idx.end <- na.omit(changing_pos[[2]])
start.idx.start <- c(1, start.idx.end[-length(start.idx.end)] + 1)
} else {
#For parallel processing
changing_pos = fread(pos_file)
positions <- na.omit(changing_pos[[3]])
chunks <- split(positions, cut(seq_along(positions), n_chunks, labels = FALSE))

starts <- na.omit(changing_pos[[1]])
start.idx.end <- na.omit(changing_pos[[2]])
start.idx.start <- c(1, start.idx.end[-length(start.idx.end)] + 1)


chunk = chunks[[chunk_inx]]
chunk_first_index = which(positions == chunk[1])

start.idx.start = start.idx.start[which(starts %in% chunk)]
start.idx.end = start.idx.end[which(starts %in% chunk)]
start = starts[starts %in% chunk]
} 

process_chunk <- function(pos, start, start.idx.start, start.idx.end, glmmkin.randomvec.obj, IBD_h5file, chr, id1.col = 1, id2.col = 2, chr.col = 3, start.col = 4, end.col = 5) {
  N <- length(glmmkin.randomvec.obj$id_include)
  ncovars <- glmmkin.randomvec.obj$ncovar
  N.randomvec <- ncol(glmmkin.randomvec.obj$random.vectors)
  Q.offset <- 2*sum(glmmkin.randomvec.obj$scaled.residuals^2)
  eigenmat.offset <- 2*crossprod(glmmkin.randomvec.obj$random.vectors)-Q.offset/(N-ncovars)*glmmkin.randomvec.obj$RSigma_R
  
  out <- data.frame(chr = chr, pos = pos, n.ibd.segs = NA, p.value = NA)
  
  
  for(ii in 1:nrow(out)) {
    if(ii %% 10000 == 0) {print(ii); system("date")}
    if(ii == 1) {
      if (start[1] == pos[1]){
        start.idx <- 1
      }
      else start.idx <- 0
      row_indices = as.numeric(changing_pos[chunk_first_index,4]):as.numeric(changing_pos[chunk_first_index,5])
      tmp <- h5read(IBD_h5file, "data", index = list(row_indices, NULL))
      tmp = tmp[which(tmp[,start.col] <= pos[ii] & tmp[,end.col]>=pos[ii]),]
      rm(row_indices)
      
      ibd <- sparseMatrix(tmp[,id1.col], tmp[,id2.col], x = 1, dims = c(N, N), symmetric = TRUE)
      n.ibd.segs <- sum(ibd@x)
      Q <- as.numeric(crossprod(glmmkin.randomvec.obj$scaled.residuals, crossprod(ibd, glmmkin.randomvec.obj$scaled.residuals)))
      eigenmat <- crossprod(glmmkin.randomvec.obj$random.vectors, crossprod(ibd, glmmkin.randomvec.obj$random.vectors))
    } else {
      if(start.idx < length(start)) {
        if(pos[ii] == start[start.idx+1]) {
          start.idx <- start.idx + 1
          suppressWarnings({tmp.add = h5read(IBD_h5file, "data", index = list(start.idx.start[start.idx]:start.idx.end[start.idx], NULL))})
          ibd.add <- sparseMatrix(tmp.add[,id1.col], tmp.add[,id2.col], x = 1, dims = c(N, N), symmetric = TRUE)
          n.ibd.segs <- n.ibd.segs + sum(ibd.add@x)
          Q <- Q + as.numeric(crossprod(glmmkin.randomvec.obj$scaled.residuals, crossprod(ibd.add, glmmkin.randomvec.obj$scaled.residuals)))
          eigenmat <- eigenmat + crossprod(glmmkin.randomvec.obj$random.vectors, crossprod(ibd.add, glmmkin.randomvec.obj$random.vectors))
          tmp <- rbind(tmp, tmp.add)
        }
      }
      if(any(tmp[[end.col]] < pos[ii])) {
        tmp.subtract <- tmp[tmp[,end.col] < pos[ii], , drop = FALSE]
        tmp <- tmp[tmp[,end.col] >= pos[ii], , drop = FALSE] 
        ibd.subtract <- sparseMatrix(tmp.subtract[,id1.col, drop = FALSE], tmp.subtract[,id2.col, drop = FALSE], x = 1, dims = c(N, N), symmetric = TRUE)
        n.ibd.segs <- n.ibd.segs - sum(ibd.subtract@x)
        Q <- Q - as.numeric(crossprod(glmmkin.randomvec.obj$scaled.residuals, crossprod(ibd.subtract, glmmkin.randomvec.obj$scaled.residuals)))
        eigenmat <- eigenmat - crossprod(glmmkin.randomvec.obj$random.vectors, crossprod(ibd.subtract, glmmkin.randomvec.obj$random.vectors))
      }
    }
    lambda <- zapsmall(eigen(eigenmat+eigenmat.offset-Q/(N-ncovars)*glmmkin.randomvec.obj$RSigma_R, only.values = TRUE, symmetric = TRUE)$values/N.randomvec)
    p.value <- try(.Q_pval(0, lambda))
    if(class(p.value) == "try-error") p.value <- NA
    if(all(lambda < 0)) p.value <- 0
    out$n.ibd.segs[ii] <- n.ibd.segs
    out$p.value[ii] <- p.value
  }
  out
}

out = FiMAP_process_chunk(chunk, start, start.idx.start, start.idx.end, obj1, IBD_h5file, chr, id1.col = 1, id2.col = 2, chr.col = 3, start.col = 4, end.col = 5)
write.table(out, paste0(outfile.prefix, "_chr", chr, "_seed", seed, "_chunk",chunk_inx,".txt"), quote=F, row.names=F, col.names=T)

























