get_changing_positions_and_h5File = function(glmmkin.randomvec.obj, IBD_file, outfile_prefix, chr, id1.col = 1, id2.col = 2, chr.col = 3, start.col = 4, end.col = 5, IBD_sorted = 1){
  find_first_overlap_idx <- function(val, end) {
    index <- which.max(val <= end)
  }
  
  segs = fread(IBD_file, data.table = FALSE)
  ids <- glmmkin.randomvec.obj$id_include
  #N <- length(ids)
  chr <- unique(segs[[chr.col]])
  if(length(chr)>1) stop("Error: more than one chromosomes observed in \"infile\".")
  idx1 <- match(segs[[id1.col]], ids)
  idx2 <- match(segs[[id2.col]], ids)
  idx.keep <- !(is.na(idx1) | is.na(idx2))
  segs <- segs[idx.keep, ]
  idx1 <- idx1[idx.keep]
  idx2 <- idx2[idx.keep]
  rm(idx.keep)
  segs[[id1.col]] <- pmin(idx1, idx2)
  segs[[id2.col]] <- pmax(idx1, idx2)
  rm(idx1, idx2)
  
  if (IBD_sorted == 0) {
    segs <- segs[order(segs[[start.col]], segs[[ end.col]]), ]
  }
  
  # get the start positions index and unique start positions and, the last end position index of IBD segments with the same start position
  start <- unique(segs[[start.col]])
  start.idx.start <- fmatch(start, segs[[start.col]])
  start.idx.end = c(start.idx.start[2:length(start.idx.start)]-1,dim(segs)[1])
  
  end <- unique(segs[[end.col]])
  end.idx.start = fmatch(end, segs[[end.col]])
  
  startAend = intersect(start,end) 
  startAend.start.idx = end.idx.start[vapply(startAend, find_first_overlap_idx, integer(1), end = end)]
  startAend.end.idx = fmatch(startAend, rev(segs[[start.col]]))
  startAend.end.idx = length(segs[[start.col]]) - startAend.end.idx+1
  
  start.only = start[!(start %in% startAend)]
  start.only.noEndBefore = start.only[which(start.only<min(end))]
  start.only.noEndBefore.start.idx = rep(1,length(start.only.noEndBefore))
  
  start.only.WithEndBefore = start.only[which(start.only>min(end))]
  #start.only.WithEndBefore.start.idx = end.idx.start[findInterval(start.only.WithEndBefore, end)+1]
  start.only.WithEndBefore.start.idx = end.idx.start[vapply(start.only.WithEndBefore, find_first_overlap_idx, integer(1), end = end)]
  
  start.only.start.idx = c(start.only.noEndBefore.start.idx,start.only.WithEndBefore.start.idx)
  start.only.end.idx = fmatch(start.only, rev(segs[[start.col]]))
  start.only.end.idx = length(segs[[start.col]]) - start.only.end.idx+1
  
  end.only = end[!(end %in% startAend)]
  end.only.end.idx = start.idx.end[findInterval(end.only, start)]
  end.only.start.idx = end.idx.start[vapply(end.only, find_first_overlap_idx, integer(1), end = end)]
  
  pos <- unique(sort(c(start, end)))
  match.idx <- match(pos, c(startAend,start.only,end.only))
  
  pos.overlap.idx.start = c(startAend.start.idx,start.only.start.idx,end.only.start.idx)[match.idx]
  pos.overlap.idx.end = c(startAend.end.idx,start.only.end.idx,end.only.end.idx)[match.idx]
  
  #many points are both add/substract points 
  out <- data.frame(chr = chr, pos = pos, n.ibd.segs = NA, p.value = NA)
  
  max_length <- length(pos)
  start <- c(start, rep(NA, max_length - length(start)))
  start.idx.end <- c(start.idx.end, rep(NA, max_length - length(start.idx.end)))
  
  write.table(cbind(start,start.idx.end,pos, pos.overlap.idx.start,pos.overlap.idx.end),paste0(outfile_prefix, "_ibd_chr",chr,"_CPos_index.txt"),quote=F, row.names=F, col.names=T, sep="\t")
  
  segs = segs[,c(id1.col, id2.col, chr.col, start.col, end.col)]
  dim_segs <- dim(segs)
  chunk_size <- c(1000, dim_segs[2])
  
  h5file <- paste0(outfile_prefix,"_ibd_chr",chr,".h5")
  h5createFile(h5file)
  h5createDataset(h5file, "data", dims = dim_segs, chunk = chunk_size, storage.mode = "integer64")
  h5write(as.matrix(segs), h5file, "data")
  
}

