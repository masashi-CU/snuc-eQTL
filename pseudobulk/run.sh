#!/bin/bash
#
# Written by Masashi Fujita
#

set -eux

################################################################################
# CONFIGURATION

# Command to run R script
run="qsub -sync y -cwd -l h_vmem=40G -S /mnt/mfs/cluster/bin/R-4.0.0/bin/Rscript"

# Annotation file of cell types. It should a tab-delimited file with a header 
# line and two columns. The "barcode" column has cell barcodes, and the 
# "cell.type" column has the cell type of the cell.
annotation_file="example/annotation.tsv"

# List of cell types to be analyzed. It should a tab-delimited file with a 
# header line and two columns. The first "prefix" column will be used as a file
# name prefix of output matrix. The second "cell.type" column should match the
# same column in the above annotation file.
celltype_file="example/celltypes.tsv"

# Input folder that stores UMI count matrix per donor. File names of the matrix
# must be "DONOR_ID.rds".
donor_matrix_dir="example/mat-per-donor"

# List of experimental batches for donors. This information will be used for 
# batch effect correction using the ComBat software. The "donor" colum is donor
# IDs. The "batch" column is experimental batches.
batch_file="example/batch.tsv"

################################################################################
# STEP 1
# Make pseudobulk count matrix of cell types. You can skip this step if you have
# created such matrix by other means.
#   Input to this step is UMI count matrix per donor. Rows are genes, columns
# are cells, and elements are UMI counts of the gene of the cell. Examples can
# be found in "example/mat-per-donor" folder. Use readRDS() function of R to
# read it.
#   Output is (1) pseudobulk count matrix of cell types, and (2) a table of the
# number of cells per donor. In (1), rows are genes, columns are donors, and 
# elements are total UMI counts per gene per cell type per donor. An example for
# endothelial cells can be found at "example/End.rds". The table (2) shows the
# number of cells of the cell type for each donor. An example can be found at
# "example/End.tsv".

# Output directory
step1_outdir="1_pseudobulk_counts"

mkdir -p ${step1_outdir}
tail -n+2 ${celltype_file} | while read line; do
    a=(${line})
    prefix=${a[0]}

    ${run} scripts/compute_pseudobulk_counts.R ${annotation_file} \
      ${celltype_file} ${donor_matrix_dir} ${step1_outdir} ${prefix}
done

################################################################################
# STEP 2
# Exclude donors who have less than 10 cells. Also exclude cell types from
# analysis if they were found in less than 10 donors. Input to this step is
# pseudobulk count matrix of cell types. Output is pseudobulk count matrix where
# donors with few cells are excluded.

# Output directory
step2_outdir="2_filter_donors"

mkdir -p ${step2_outdir}
for input_rds in ${step1_outdir}/*.rds; do
    input_tsv=${input_rds/.rds/.tsv}
    ${run} scripts/filter-donors.R ${input_rds} ${input_tsv} ${step2_outdir}
done

################################################################################
# STEP 3
# Make pseudobulk gene expression matrix by computing log2 of counts per million 
# reads mapped (CPM). Also exclude genes whose mean log2CPM is less than 2.0.
# Input to this step is filtered pseudobulk count matrix of cell types. Output
# is pseudobulk gene expression matrix in log2CPM.

# Output directory
step3_outdir="3_log2cpm"

mkdir -p ${step3_outdir}
for input_rds in ${step2_outdir}/*.rds; do
    ${run} scripts/cpm.R ${input_rds} ${step3_outdir}
done

################################################################################
# STEP 4
# Adjust batch effects in the log2CPM matrix using ComBat. Input to this step is
# pseudobulk gene expression matrix in log2CPM. Output is adjusted log2CPM 
# matrix.

# Output directory
step4_outdir="4_adjust_batch"

mkdir -p ${step4_outdir}
for input_tsv in ${step3_outdir}/*.log2cpm.tsv; do
    ${run} scripts/adjust-batch.R ${batch_file} ${input_tsv} ${step4_outdir}
done

################################################################################
# STEP 5
# Quantile normalize the log2CPM expression. Input to this step is batch-
# adjusted pseudobulk gene expression matrix in log2CPM. Output is quantile-
# normalized matrix.

# Output directory
step5_outdir="5_quant_norm"

mkdir -p ${step5_outdir}
for input_tsv in ${step4_outdir}/*.log2cpm.tsv; do
    ${run} scripts/quant-norm.R ${input_tsv} ${step5_outdir}
done
