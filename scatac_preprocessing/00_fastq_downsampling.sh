#!/bin/bash

n_lines=$1

well_ids="$(cat well_ids.txt)"
fastq_samples="$(echo ${well_ids} | sed -E 's/[BX][0-9]{3}-AP0C1W//g')"
batch_ids="$(echo ${well_ids} | sed -E 's/-AP0C1W[0-9]//g')"
in_dir=/mnt/atac-methods-in/raw
out_dir=/mnt/atac-methods-in/downsampled_125M

for r in I1 R1 R2 R3
do
parallel -j 8 \
  --bar \
  --link \
  "mkdir -p ${out_dir}/{1}/{2}/; \
  zcat ${in_dir}/{1}/{2}/{2}_S{3}_L001_${r}_001.fastq.gz \
  | head -${n_lines} \
  > ${out_dir}/{1}/{2}/{2}_S{3}_L001_${r}_001.fastq; \
  gzip ${out_dir}/{1}/{2}/{2}_S{3}_L001_${r}_001.fastq" \
  ::: ${batch_ids} ::: ${well_ids} ::: ${fastq_samples}
done
