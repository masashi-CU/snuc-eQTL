# Masashi Fujita

args <- commandArgs(trailingOnly=T)
input_tsv <- args[1]
outdir <- args[2]

# load log2CPM
logcpm  <- read.table(input_tsv, sep="\t", header=T, check.names=F)

# convert log2CPM to matrix
logcpm_id <- logcpm$id
logcpm <- as.matrix(logcpm[, colnames(logcpm) != "id"])
rownames(logcpm) <- logcpm_id

# quantile normalizarion
logcpm <- t(apply(logcpm, 1, rank, ties.method = "average"))
logcpm <- qnorm(logcpm / (ncol(logcpm) + 1))

# export
df <- data.frame(id = rownames(logcpm), logcpm, check.names = F)
filename <- basename(input_tsv)
outfile <- file.path(outdir, filename)
write.table(df, file=outfile, sep="\t", quote = F, row.names = F)
