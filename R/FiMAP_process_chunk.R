FiMAP_process_chunk <- function(pos, start, start.idx.start, start.idx.end, glmmkin.randomvec.obj, IBD_h5file, chr, id1.col = 1, id2.col = 2, chr.col = 3, start.col = 4, end.col = 5) {
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