#!/bin/bash

USAGE="${0} -r=<run dir> -o=<output dir> -g=<gene expression samplesheet> -a=<atac samplesheet>

Command line options:
    -r [ --run-dir ]            Path of Illumina BCL run folder.
    -o [ --output-dir ]         Output path for fastq generation.
    -g [ --gex-samplesheet ]    Path to Illumina Experiment Manager compatible samplesheet for Gene Expression libraries.
    -a [ --atac-samplesheet ]   Path to Illumina Experiment Manager compatible samplesheet for Gene Expression libraries.
    -c [ --cite-samplesheet ]   Path to Illumina Experiment Manager compatible samplesheet for ADT libraries.
    "

# check number of arguments
if [ $# -ne 5 ]
then
    echo "$USAGE"
    exit 1
fi

# initialize values to NULL
FLOWCELL_DIR=''
OUTPUT_DIR=''
GEX_SAMPLE_SHEET_PATH=''
ATAC_SAMPLE_SHEET_PATH=''
ADT_SAMPLE_SHEET_PATH=''

# parse arguments
for arg in "$@"
do
    case $arg in
        -r=*|--run-dir=*)
        FLOWCELL_DIR="${arg#*=}"
        shift # remove --bam from processing
        ;;
        -o=*|--output-dir=*)
        OUTPUT_DIR="${arg#*=}"
        shift # remove --barcodes from processing
        ;;
        -g=*|--gex-samplesheet=*)
        GEX_SAMPLE_SHEET_PATH="${arg#*=}"
        shift # remove --threads from processing
        ;;
        -a=*|--atac-samplesheet=*)
        ATAC_SAMPLE_SHEET_PATH="${arg#*=}"
        shift # remove --sample-name from processing
        ;;
        -c=*|--cite-samplesheet=*)
        ADT_SAMPLE_SHEET_PATH="${arg#*=}"
        shift # remove --sample-name from processing
        ;;
        *)
        echo "Unknown argument ${arg} provided. Exiting..."
        echo ""
        echo "$USAGE"
        exit 1
        ;;
    esac
done

# code from LG for checking user provided script parameters
# Trap exit 1 to allow termination within the check_param() function.
trap "exit 1" TERM
export TOP_PID=$$

# Time statement function
stm() {
    local ts=$(date +"%Y-%m-%d %H:%M:%S")
    echo "["$ts"] "$1
}

# check parameters
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

echo
echo $(stm "START Multiome bcl2fastq demultiplexing")
echo $(check_param "-r | --run-dir" "Run Directory" ${FLOWCELL_DIR})
echo $(check_param "-o | --output-dir" "Output Directory" ${OUTPUT_DIR})
echo $(check_param "-g | --gex-samplesheet" "Gene Expression Illumina Samplesheet" ${GEX_SAMPLE_SHEET_PATH})
echo $(check_param "-a | --atac-samplesheet" "ATAC Illumina Samplesheet" ${ATAC_SAMPLE_SHEET_PATH})
echo $(check_param "-c | --cite-samplesheet" "ADT Illumina Samplesheet" ${ADT_SAMPLE_SHEET_PATH})
echo

INTEROP_DIR="${FLOWCELL_DIR}/InterOp"

# demultiplex gene expression data
bcl2fastq --use-bases-mask=Y28n*,I10,I10n*,Y90n* \
        --create-fastq-for-index-reads \
        --minimum-trimmed-read-length=8 \
        --mask-short-adapter-reads=8 \
        --ignore-missing-positions \
        --ignore-missing-filter \
        --ignore-missing-bcls \
        -r 24 -w 24 -p 80 \
        -R ${FLOWCELL_DIR} \
        --output-dir="${OUTPUT_DIR}/GEX" \
        --interop-dir=${INTEROP_DIR} \
        --sample-sheet=${GEX_SAMPLE_SHEET_PATH}

# demultiplex ATAC data
bcl2fastq --use-bases-mask=Y50n*,I8n*,Y16,Y50n* \
        --create-fastq-for-index-reads \
        --minimum-trimmed-read-length=8 \
        --mask-short-adapter-reads=8 \
        --ignore-missing-positions \
        --ignore-missing-filter \
        --ignore-missing-bcls \
        -r 24 -w 24 -p 80 \
        -R ${FLOWCELL_DIR} \
        --output-dir="${OUTPUT_DIR}/ATAC" \
        --interop-dir=${INTEROP_DIR} \
        --sample-sheet=${ATAC_SAMPLE_SHEET_PATH}

# demultiplex ADT data
bcl2fastq --use-bases-mask=Y28n*,I8n*,n*,Y90n* \
        --create-fastq-for-index-reads \
        --minimum-trimmed-read-length=8 \
        --mask-short-adapter-reads=8 \
        --ignore-missing-positions \
        --ignore-missing-controls \
        --ignore-missing-filter \
        --ignore-missing-bcls \
        -r 24 -w 24 -p 80 \
        -R ${FLOWCELL_DIR} \
        --output-dir="${OUTPUT_DIR}/ADT" \
        --interop-dir=${INTEROP_DIR} \
        --sample-sheet=${ADT_SAMPLE_SHEET_PATH}