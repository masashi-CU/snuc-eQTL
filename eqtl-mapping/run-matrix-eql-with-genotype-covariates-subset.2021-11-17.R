#!/mnt/mfs/cluster/bin/R-4.0.0/bin/Rscript
#$ -cwd
#$ -j y
#$ -l h_vmem=20G
#$ -S /mnt/mfs/cluster/bin/R-4.0.0/bin/Rscript

library(MatrixEQTL)
library(getopt)

## Parse command line arguments

spec <- matrix(c(
  'input', 'i', 1, "character", "input directory that contains log2cpm, dosage, and snppos.",
  'output', 'o', 1, "character", "output directory",
  'celltype', 't', 1, "character", "cell type",
  'cov', 'c', 1, "character", "covariate file",
  'help', 'h', 0, "logical", "show help message"
), byrow=TRUE, ncol=5)
opt = getopt(spec)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

if ( is.null(opt$input) ) { stop("input directory must be specified") }
if ( is.null(opt$output) ) { stop("output directory must be specified") }
if ( is.null(opt$celltype) ) { stop("celltype must be specified") }

if(! dir.exists(opt$output)){
  dir.create(opt$output, showWarnings = T, recursive = T, mode = "0775")
}
################################################################################

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
# SNP_file_name = sprintf("%s/%s/get-dosage.%s.dosage", opt$input, opt$celltype, opt$celltype);
# snps_location_file_name = sprintf("%s/%s/get-dosage.%s.snppos", opt$input, opt$celltype, opt$celltype);
SNP_file_name = "/mnt/mfs/ctcn/team/masashi/snuc-eqtl/genotype/get-dosage.ALL.dosage"
snps_location_file_name = "/mnt/mfs/ctcn/team/masashi/snuc-eqtl/genotype/get-dosage.ALL.snppos"

# Gene expression file name
expression_file_name = sprintf("%s/%s/%s.log2cpm.tsv", opt$input, opt$celltype, opt$celltype);
# gene_location_file_name = file.path(opt$input, "get-tss-pos.tsv")
gene_location_file_name = "/mnt/mfs/ctcn/team/masashi/snuc-eqtl/transcriptome/get-tss-pos.tsv"

# Covariates file name
# Set to character() for no covariates
covariates_file_name = opt$cov;

# Output file name
output_file_name_cis = NULL;
output_file_name_tra = NULL;

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 0;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 4;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Subset genotype data to have the same donors with gene expression data

matching_cols <- match(gene$columnNames, sub("^\\[[0-9]+\\]", "", snps$columnNames))
stopifnot(!is.na(matching_cols))
snps$ColumnSubsample(matching_cols)
stopifnot(gene$columnNames == sub("^\\[[0-9]+\\]", "", snps$columnNames))

## Subset covariates to have the same donors with gene expression data
if(length(covariates_file_name)>0) {
  matching_cols <- match(gene$columnNames, cvrt$columnNames)
  stopifnot(!is.na(matching_cols))
  cvrt$ColumnSubsample(matching_cols)
  stopifnot(gene$columnNames == cvrt$columnNames)
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

# Copy SNP names
stopifnot(nrow(snps) == nrow(snpspos));
rownames(snps) = snpspos[,1];

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = FALSE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)

## Plot the Q-Q plot of local and distant p-values
pngfile = file.path(opt$output, "qqplot.png");
png(pngfile);
plot(me);
dev.off();

## Export result
rdsfile = file.path(opt$output, "matrix-eqtl.rds")
saveRDS(me, rdsfile)
