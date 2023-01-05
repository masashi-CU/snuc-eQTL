# Masashi Fujita

library(edgeR)
library(limma)

args <- commandArgs(trailingOnly=T)
input_rds <- args[1]
outdir <- args[2]

print(input_rds)
filename <- basename(input_rds)
prefix <- sub("\\.rds", "", filename)

# load a count matrix for the celltype
mat  <- readRDS(input_rds)

# filter out low-expression genes
y <- DGEList(counts = mat)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=F]

# TMM normalization
y <- calcNormFactors(y, method = "TMM")

# voom transform
pngfile <- sprintf("%s/%s.voom.png", outdir, prefix)
png(pngfile)
v <- voom(y, plot=T)
dev.off()
logcpm <- v$E

# draw histogram of mean log2CPM
pngfile <- sprintf("%s/%s.hist.png", outdir, prefix)
png(pngfile)
hist(apply(logcpm, 1, mean), 100)
dev.off()

# remove genes if mean log2CPM < 2.0
mean_logcpm <- apply(logcpm, 1, mean)
logcpm <- logcpm[mean_logcpm > 2.0,]

# save log2CPM
outfile <- sprintf("%s/%s.log2cpm.tsv", outdir, prefix)
cat("id", file=outfile)
write.table(logcpm, file=outfile, sep="\t", quote=F, col.names=NA, append=T)
