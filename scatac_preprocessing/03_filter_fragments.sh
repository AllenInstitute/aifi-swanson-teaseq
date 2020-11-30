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

while getopts "i:j:m:o:t:s:" opt; do
  case $opt in
    i) in_fragments="$OPTARG"
    ;;
    j) num_jobs="$OPTARG"
    ;;
    m) in_meta="$OPTARG"
    ;;
    o) output_dir="$OPTARG"
    ;;
    t) temp_dir="$OPTARG"
    ;;
    s) sample_name="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

echo $(stm "START scATAC fragment filtering")
echo $(check_param "-i" "Input fragments.tsv.gz" ${in_fragments})
echo $(check_param "-j" "N Parallel Jobs" ${num_jobs})
echo $(check_param "-m" "Filtered metadata.csv.gz" ${in_meta})
echo $(check_param "-o" "Output Directory" ${output_dir})
echo $(check_param "-t" "Temporary Directory" ${temp_dir})
echo $(check_param "-s" "Sample Name" ${sample_name})
total_start_time="$(date -u +%s)"

# Generate temp files

echo $(stm "Making Temporary Directory ${temp_dir}")
mkdir ${temp_dir}

# Generate whitelists for each prefix

mkdir ${temp_dir}/whitelist

echo $(stm "Generating whitelist from $in_meta")

zcat ${in_meta} \
  | awk -F',' -v OFS='\t' \
  -v a=${temp_dir}/whitelist/ \
  -v c=.tsv \
  '{ if(NR > 1) { b=substr($1,0,4); print $1,$18 > a b c } }'

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

echo $(stm "Filtering fragments")

split_start_time="$(date -u +%s)"

prefixes="$(ls ${temp_dir}/whitelist | sed 's/[0-9]\+_//g' | sed 's/.tsv//g' | sort | uniq)"

mkdir ${temp_dir}/frag_filt/
mkdir ${temp_dir}/frag_fail/ 
parallel -j ${num_jobs} \
  --bar \
  "cat ${temp_dir}/whitelist/{}.tsv \
    ${temp_dir}/frag/*_{}.bed \
  | awk -F'\t' -v OFS='\t' \
    -v p=${temp_dir}/frag_filt/{/.}.bed \
    -v f=${temp_dir}/frag_fail/{/.}.bed \
    '{ if(NF==2){w[\$1]=\$2;next} else if(\$4 in w) {print \$1,\$2,\$3,w[\$4],\$5 > p} else {print \$0 > f} }' " \
  ::: $prefixes

echo $(stm "$(elt $split_start_time $total_start_time)" )

split_start_time="$(date -u +%s)"

echo $(stm "Reassembling passing fragments")

cat ${temp_dir}/frag_filt/* \
  | sort -k1,1 -k2,2n -k3,3n \
  > ${output_dir}/${sample_name}_filtered_fragments.tsv

echo $(stm "Reassembling failing fragments")

cat ${temp_dir}/frag_fail/* \
  | sort -k1,1 -k2,2n -k3,3n \
  > ${output_dir}/${sample_name}_failed_fragments.tsv

echo $(stm "$(elt $split_start_time $total_start_time)" )

split_start_time="$(date -u +%s)"

echo $(stm "Compressing outputs")

bgzip ${output_dir}/${sample_name}_filtered_fragments.tsv
bgzip ${output_dir}/${sample_name}_failed_fragments.tsv

echo $(stm "$(elt $split_start_time $total_start_time)" )
split_start_time="$(date -u +%s)"

echo $(stm "Cleaning up temporary files")

rm -r ${temp_dir}

echo $(stm "$(elt $split_start_time $total_start_time)" )
echo $(stm "END scATAC fragment filtering analysis")
