library(data.table)
library(Matrix)
library(rhdf5)
source("FiMAP_process_chunk.R")
source("pval_caculation_functions.R")

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
IBD_h5file = paste0("toy_ibd_chr20.h5")
pos_file = paste0("toy_ibd_chr20_CPos_index.txt")
obj1 <- readRDS(paste0("toy_glmmkin2randomvec.rds"))

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

out = FiMAP_process_chunk(chunk, start, start.idx.start, start.idx.end, obj1, IBD_h5file, chr, id1.col = 1, id2.col = 2, chr.col = 3, start.col = 4, end.col = 5)
write.table(out, paste0(outfile.prefix, "_chr", chr, "_seed", seed, "_chunk",chunk_inx,".txt"), quote=F, row.names=F, col.names=T)

























