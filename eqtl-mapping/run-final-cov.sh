#!/bin/bash

set -eu

base_dir="/home/mf3362/snuc-eqtl/v20211109.celltypes"

script="/mnt/mfs/ctcn/team/masashi/snuc-eqtl/script/run-matrix-eql-with-genotype-covariates-subset.2021-11-17.R"

mkdir -p log
for celltype_dir in ${base_dir}/*; do
    celltype=$(basename ${celltype_dir})
    echo "${celltype}"
    output="${celltype_dir}/matrix-eqtl/covariates-20211118"
    covfile="${celltype_dir}/covariates-20211118.tsv"
    qsub -o log/ ${script} -i ${base_dir} -o ${output} -t ${celltype} -c ${covfile}
done
