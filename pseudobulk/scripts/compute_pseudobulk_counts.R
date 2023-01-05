# Masashi Fujita

library(Matrix)

# parse args
args <- commandArgs(trailingOnly = T)
annotation_file <- args[1]
celltype_file <- args[2]
donor_matrix_dir <- args[3]
outdir <- args[4]
prefix <- args[5]

# fetch long name of cell type
abb <- read.delim(celltype_file)
cell.type <- abb[abb$prefix == prefix, "cell.type"]

# load cell type annotation
anno <- read.delim(annotation_file)

# extract barcodes of the query cell type
celltype_barcodes <- anno[anno$cell.type == cell.type, "barcode"]

# list UMI count matrices of donors
mat_files <- Sys.glob(file.path(donor_matrix_dir, "*.rds"))

# get the gene names in a matrix
rn1 <- rownames(readRDS(mat_files[1]))

# extract counts of specified cell type and compute pseudobulk.
# Also count cells per donor.
res <- lapply(X = mat_files, FUN = function(mat_file){
  print(mat_file)
  mat <- readRDS(mat_file)
  mat <- mat[, colnames(mat) %in% celltype_barcodes, drop=F]
  stopifnot(rownames(mat) == rn1)   # ensure the consistency of rownames
  list(sum = rowSums(mat), n_cell = ncol(mat))
})
pb <- sapply(res, "[[", "sum")  # pseudobulk count matrix
nc <- sapply(res, "[[", "n_cell")  # number of cells of the cell type

stopifnot(ncol(pb) == length(mat_files))

# extract donor names from matrix file names
donors <- sub("\\.rds", "", basename(mat_files))

# set donor names
colnames(pb) <- donors
nc.df <- data.frame(donor = donors, num_cells = nc)

# export pseudobulk matrix
output_rds <- sprintf("%s/%s.rds", outdir, prefix)
saveRDS(pb, file = output_rds)

# export number of cells
output_tsv <- sprintf("%s/%s.tsv", outdir, prefix)
write.table(nc.df, file = output_tsv, sep = "\t", quote = F, row.names = F)
