# Masashi Fujita

# parse args
args <- commandArgs(trailingOnly = T)
input_rds <- args[1]
input_tsv <- args[2]
outdir <- args[3]

print(input_rds)
print(input_tsv)
filename <- basename(input_rds)
prefix <- sub("\\.rds", "", filename)

# load a count matrix for the celltype
mat <- readRDS(input_rds)

# read the number of cells per donor
cell <- read.table(input_tsv, sep = "\t", header = T)

# find donors who have 10 or more cells.
cnt <- cell[cell$num_cells >= 10, , drop = F]

# export matrix if 10 or more donors passed the filter of cell numbers
if(nrow(cnt) >= 10) {
  submat <- mat[, cnt$donor, drop = F]
  output_rds <- file.path(outdir, filename)
  saveRDS(submat, file = output_rds)
} else {
  cat("[ERROR] <10 donors passed the filter\n")
}
