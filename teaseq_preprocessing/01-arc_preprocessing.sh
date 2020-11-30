#!/bin/bash

# Trap exit 1 to allow termination within the check_param() function.
trap "exit 1" TERM
export TOP_PID=$$

# Time statement function

stm() {
	local ts=$(date +"%Y-%m-%d %H:%M:%S")
	echo "["$ts"] "$1
}

# Elapsed time function

format_time() {
  local h=$(($1/60/60%24))
  local m=$(($1/60%60))
  local s=$(($1%60))
  
  printf "%02d:%02d:%02d" $h $m $s
}

elt() {
  local end_time="$(date -u +%s)"
  local split_diff="$(($end_time-$1))"
  local total_diff="$(($end_time-$2))"
  
  echo "Total time: " $(format_time $total_diff) "| Split time: " $(format_time $split_diff)
}

check_param() {
  local pflag=$1
  local pname=$2
  local pvar=$3
  
  if [ -z ${pvar} ]; then
    echo $(stm "ERROR ${pflag} ${pname}: parameter not set. Exiting.")
    kill -s TERM $TOP_PID 
  else
    echo  $(stm "PARAM ${pflag} ${pname}: ${pvar}")
  fi
}

pipeline_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Parse command-line arguments

while getopts "b:i:g:j:o:p:t:w:" opt; do
  case $opt in
    b) bedtools_loc="$OPTARG"
    ;;
    i) in_fragments="$OPTARG"
    ;;
    g) genome="$OPTARG"
    ;;
    j) num_jobs="$OPTARG"
    ;;
    o) output_dir="$OPTARG"
    ;;
    p) in_per_barcode_metrics="$OPTARG"
    ;;
    t) temp_dir="$OPTARG"
    ;;
    w) well_id="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

echo $(stm "START scATAC bedtools analysis")
echo $(check_param "-b" "BedTools Location" ${bedtools_loc})
echo $(check_param "-i" "Input fragments.tsv.gz" ${in_fragments})
echo $(check_param "-g" "Genome (hg19 or hg38)" ${genome})
echo $(check_param "-j" "N Parallel Jobs" ${num_jobs})
echo $(check_param "-o" "Output Directory" ${output_dir})
echo $(check_param "-p" "Input per_barcode_metrics.csv" ${in_per_barcode_metrics})
echo $(check_param "-t" "Temporary Directory" ${temp_dir})
echo $(check_param "-w" "Well ID" ${well_id})
total_start_time="$(date -u +%s)"

# Generate temp files

echo $(stm "Making Temporary Directory ${temp_dir}")
mkdir ${temp_dir}

echo $(stm "Indexing and sorting references")

split_start_time="$(date -u +%s)"

# Index and sort references
mkdir ${temp_dir}/ref

ref_sources="$(ls -1 ${pipeline_dir}/reference/${genome}*.bed.gz | tr '\n' ' ')"
ref_names="$(echo -n $ref_sources | sed 's!'"${pipeline_dir}/reference/${genome}_"'!!g' | sed 's/.bed.gz//g' )"
ref_pattern="$(echo -n $ref_names | sed 's/ /|/g' )"

parallel -j ${num_jobs} \
  --link \
  --bar \
  "zcat {1} \
  | dos2unix \
  | awk -F'\t' -v OFS='\t' \
    '{ print \$1,\$2,\$3,NR }' \
  | sort -k1,1 -k2,2n -k3,3n \
  > ${temp_dir}/ref/{2}.bed" \
  ::: $ref_sources ::: $ref_names

ref_files="$(ls -1 ${temp_dir}/ref/*.bed | tr '\n' ' ')"

# Generate whitelists for each prefix

mkdir ${temp_dir}/whitelist

if [ -f ${in_per_barcode_metrics} ]; then
  echo $(stm "Generating whitelist from in_per_barcode_metrics.csv")

  cat ${in_per_barcode_metrics} \
  | awk -F',' -v OFS='\t' \
    -v a=${temp_dir}/whitelist/ \
    -v c=.tsv \
    '{ if($28 > 1000 && NR > 2) { b=substr($1,0,4); print $1 > a b c } }'
else
  in_per_barcode_metrics=$( echo ${in_per_barcode_metrics} | sed 's/per_barcode_metrics.csv/barcode_counts.tsv/g' )
  
  echo $(stm "No singlecell.csv; Switching to barcode_counts.tsv")
  
  if [ ! -f ${in_per_barcode_metrics} ]; then
    echo $(stm "Generating barcode_counts.tsv from $in_fragments")

    zcat ${in_fragments} \
    | awk -F'\t' '{print $4}' \
    | sort -k1,1 \
    | uniq -c \
    | awk -v OFS='\t' '{print $2,$1}' \
    > ${in_per_barcode_metrics}
  fi
  
  echo $(stm "Generating whitelist from $in_per_barcode_metrics")

  cat ${in_per_barcode_metrics} \
  | awk -F'\t' -v OFS='\t' \
    -v a=${temp_dir}/whitelist/ \
    -v c=.tsv \
    '{ if($2 > 1000 && NR > 2) { b=substr($1,0,4); print $1 > a b c } }'
fi

echo $(stm "Splitting input file $in_fragments by barcode prefix")

# Split fragments by barcode prefix
mkdir ${temp_dir}/frag

zcat ${in_fragments} \
| parallel -j ${num_jobs} \
  --pipe --block 100M \
  "awk -F'\t' -v OFS='\t' \
    -v a=${temp_dir}/frag/ \
    -v b={#}_ \
    -v d=.bed \
    '{ c=substr(\$4,0,4); print >> a b c d}' "

echo $(stm "Filtering and Resorting fragments")

prefixes="$(ls ${temp_dir}/whitelist | sed 's/[0-9]\+_//g' | sed 's/.tsv//g' | sort | uniq)"

mkdir ${temp_dir}/frag_sort/
parallel -j ${num_jobs} \
  --bar \
  "cat ${temp_dir}/whitelist/{}.tsv \
    ${temp_dir}/frag/*_{}.bed \
  | awk -F'\t' -v OFS='\t' \
    '{ if(NF==1){w[\$1]=1;next} else if(\$4 in w) {print} }' \
  | sort -k1,1 -k2,2n -k3,3n \
  > ${temp_dir}/frag_sort/{/.}.bed" \
  ::: $prefixes

echo $(stm "$(elt $split_start_time $total_start_time)" )

split_start_time="$(date -u +%s)"

echo $(stm "Starting fragment analysis")

echo $(stm "Performing overlaps with references")

# Perform feature overlaps
mkdir ${temp_dir}/ref_sums
mkdir ${temp_dir}/ref_ol
parallel -j ${num_jobs} \
  --bar \
  "$bedtools_loc \
    intersect \
    -sorted \
    -a ${temp_dir}/frag_sort/{/.}.bed \
    -b ${ref_files} \
    -names ${ref_names} \
    -loj \
  | sort -k4,4 -k6,6 -k10,10n \
  | $bedtools_loc \
    groupby \
    -i stdin \
    -g 4,6,10 \
    -c 10 \
    -o count \
  > ${temp_dir}/ref_ol/{/.}.tsv" \
  ::: ${temp_dir}/frag_sort/*.bed

echo $(stm "$(elt $split_start_time $total_start_time)" )
split_start_time="$(date -u +%s)"

echo $(stm "Counting reads overlapping references")

parallel -j ${num_jobs} \
  --bar \
  "$bedtools_loc \
    intersect \
    -sorted \
    -a ${temp_dir}/frag_sort/{/.}.bed \
    -b ${ref_files} \
    -names ${ref_names} \
    -C \
  | awk -F'\t' -v OFS='\t' \
    '{ if(\$7 > 0) {print \$4,\$6} }' \
  | sort -k1,1 -k2,2 \
  | $bedtools_loc \
    groupby \
    -i stdin \
    -g 1,2 \
    -c 2 \
    -o count \
  > ${temp_dir}/ref_sums/{/.}.tsv" \
  ::: ${temp_dir}/frag_sort/*.bed

echo $(stm "$(elt $split_start_time $total_start_time)" )
split_start_time="$(date -u +%s)"

echo $(stm "Counting saturation")

mkdir ${temp_dir}/saturation
parallel -j ${num_jobs} \
  --bar \
  "cat ${temp_dir}/frag/{/.}.bed \
  | awk -F'\t' \
    '{print \$5}' \
  | sort -k1,1n \
  | uniq -c \
  | awk -v OFS='\t' \
    '{print \$2,\$1}' \
  > ${temp_dir}/saturation/{/.}.tsv" \
  ::: ${temp_dir}/frag/*.bed

echo $(stm "$(elt $split_start_time $total_start_time)" )
split_start_time="$(date -u +%s)"

echo $(stm "Counting fragment widths")

mkdir ${temp_dir}/widths
parallel -j ${num_jobs} \
  --bar \
  "cat ${temp_dir}/frag/*_{}.bed \
  | awk -F'\t' -v OFS='\t' \
    '{w=\$3-\$2; print \$4,w}' \
  | sort -k1,1 -k2,2n \
  | uniq -c \
  | awk -v OFS='\t' \
    '{print \$2,\$3,\$1}' \
  > ${temp_dir}/widths/{/.}.tsv" \
  ::: $prefixes

echo $(stm "$(elt $split_start_time $total_start_time)" )
split_start_time="$(date -u +%s)"

echo $(stm "Starting parallel window count with ${num_jobs} threads")

# Count chromosome windows
mkdir ${temp_dir}/window_5k
mkdir ${temp_dir}/window_20k
#mkdir ${temp_dir}/window_100k
parallel -j ${num_jobs} \
  --bar \
  "cat ${temp_dir}/frag_sort/{/.}.bed \
  | awk -F'\t' -v OFS='\t' \
    '{ c = (\$3 + \$2) / 2; x = int(c / 100000) + 1; y = int(c / 20000) + 1; z = int(c / 5000) + 1; print \$4,\$1,x,y,z }' \
  | sort -k1,1 -k2,2 -k5,5n \
  | $bedtools_loc \
    groupby \
    -i stdin \
    -g 1,2,3,4,5 \
    -c 5 \
    -o count \
  | tee ${temp_dir}/window_5k/{/.}.tsv \
  | $bedtools_loc \
    groupby \
    -i stdin \
    -g 1,2,3,4 \
    -c 6 \
    -o sum \
  > ${temp_dir}/window_20k/{/.}.tsv" \
  ::: ${temp_dir}/frag_sort/*.bed
  # | tee ${temp_dir}/window_20k/{/.}.tsv \
  # | $bedtools_loc \
  #   groupby \
  #   -i stdin \
  #   -g 1,2,3 \
  #   -c 5 \
  #   -o sum \
  # > ${temp_dir}/window_100k/{/.}.tsv" \

echo $(stm "$(elt $split_start_time $total_start_time)" )
split_start_time="$(date -u +%s)"

# Assemble results across all barcodes
mkdir ${output_dir}

echo $(stm "Assembling sparse matrices")

echo $(stm "Splitting overlap results by reference")
mkdir ${temp_dir}/ref_ol_split
parallel -j ${num_jobs} \
  --bar \
  "cat {} \
  | grep -E '${ref_pattern}' \
  | awk -F'\t' -v OFS='\t' \
    -v a=${temp_dir}/ref_ol_split/{/.}_ \
    -v c=.tsv \
    '{print \$1,\$3,\$4 > a \$2 c}'" \
  ::: ${temp_dir}/ref_ol/*.tsv

echo $(stm "Generating final overlap matrices")
parallel -j ${num_jobs} \
  "cat ${temp_dir}/ref_ol_split/*{}.tsv \
  > ${output_dir}/${well_id}_{}_sparse_matrix.tsv" \
  ::: ${ref_names}
  
echo $(stm "$(elt $split_start_time $total_start_time)" )
split_start_time="$(date -u +%s)"

echo $(stm "Assembling reference counts")

echo $(stm "Splitting results by reference")
mkdir ${temp_dir}/ref_sums_split
parallel -j ${num_jobs} \
  --bar \
  "cat {} \
  | grep -E '${ref_pattern}' \
  | awk -F'\t' -v OFS='\t' \
    -v a=${temp_dir}/ref_sums_split/{/.}_ \
    -v c=.tsv \
    '{print \$1,\$3 > a \$2 c}'" \
  ::: ${temp_dir}/ref_sums/*.tsv

echo $(stm "Generating final counts")
parallel -j ${num_jobs} \
  "cat ${temp_dir}/ref_sums_split/*{}.tsv \
  > ${output_dir}/${well_id}_{}_total_counts.tsv" \
  ::: $ref_names

echo $(stm "$(elt $split_start_time $total_start_time)" )
split_start_time="$(date -u +%s)"

echo $(stm "Assembling window counts")
cat ${temp_dir}/window_5k/*.tsv \
| awk -F'\t' -v OFS='\t' \
  '{print $1,$2,$5,$6}' \
> ${output_dir}/${well_id}_window_5k_counts.tsv

cat ${temp_dir}/window_20k/*.tsv \
  | awk -F'\t' -v OFS='\t' \
    '{print $1,$2,$4,$5}' \
  > ${output_dir}/${well_id}_window_20k_counts.tsv
#
# cat ${temp_dir}/window_100k/*.tsv \
# > ${output_dir}/${well_id}_window_100k_counts.tsv

echo $(stm "$(elt $split_start_time $total_start_time)" )
split_start_time="$(date -u +%s)"

echo $(stm "Assembling saturation")

cat ${temp_dir}/saturation/*.tsv \
| sort -k1,1n \
| $bedtools_loc \
  groupby \
  -i stdin \
  -g 1 \
  -c 2 \
  -o sum \
> ${output_dir}/${well_id}_saturation_curve.tsv

echo $(stm "$(elt $split_start_time $total_start_time)" )
split_start_time="$(date -u +%s)"

echo $(stm "Assembling widths")

cat ${temp_dir}/widths/*.tsv \
  > ${output_dir}/${well_id}_fragment_widths.tsv

echo $(stm "$(elt $split_start_time $total_start_time)" )
split_start_time="$(date -u +%s)"

echo $(stm "Compressing outputs")

parallel -j ${num_jobs} \
  --bar \
  "gzip {}" \
  ::: ${output_dir}/*

echo $(stm "$(elt $split_start_time $total_start_time)" )
split_start_time="$(date -u +%s)"

echo $(stm "Cleaning up temporary files")

rm -r ${temp_dir}

echo $(stm "$(elt $split_start_time $total_start_time)" )
echo $(stm "END scATAC bedtools analysis")
