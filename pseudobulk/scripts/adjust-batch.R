# Masashi Fujita

library(sva)

args <- commandArgs(trailingOnly=T)
batch_file <- args[1]
input_tsv <- args[2]
outdir <- args[3]

# load log2CPM
logcpm  <- read.table(input_tsv, sep="\t", header=T, check.names=F)

# convert log2CPM to matrix
logcpm_id <- logcpm$id
logcpm <- as.matrix(logcpm[, colnames(logcpm) != "id"])
rownames(logcpm) <- logcpm_id

# read batch info
b <- read.table(batch_file, sep="\t", header=T)
batches <- b$batch
names(batches) <- b$donor

# reorder batches to have the same order with log2cpm 
batches <- batches[colnames(logcpm)]
stopifnot(names(batches) == colnames(logcpm))

# adjust for batch effect
if(all(table(batches) == 1)){
  res <- logcpm  # if there is no meaningful batch
} else {
  res <- ComBat(dat = logcpm, batch = batches)
}

# export
df <- data.frame(id = rownames(res), res, check.names = F)
filename <- basename(input_tsv)
outfile <- file.path(outdir, filename)
write.table(df, file=outfile, sep="\t", quote = F, row.names = F)
