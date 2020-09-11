"""
Script written by Elliott Swanson of the Allen Institute for Immunology (elliott.swanson@alleninstitute.org). Free for academic use only.
See LICENSE for code reuse permissions.
"""

import os 
import sys
from glob import glob
import subprocess
from datetime import datetime
import pytz
import logging
import pandas as pd
import numpy as np
import csv
import concurrent.futures
import gzip
import argparse


def main():
    """ Take path to a data directory as input and create a list of all read 2 fastq files in that directory. FastQ files must follow the standard
        Illumina naming comvention and be .gz zipped.

        For each read pair:
        1. Trim 3' end of read 2 using fastp
        2. Align each read 2 fastq file to hg38 using bowtie2 and wite stdout to a SAM file.
        3. Add cell barcode and UMI fields to SAM alignments and filter out short or low quality alignments. Perform barcode correction.
        4. Create sorted, corresponding BAM file from SAM file using GATK.
        5. Index BAM file using GATK.
        6. Create a sorted BED file from BAM file using Bedtools.
        7. Create fragments.tsv file to be used for downstream 10X processing

        Stdout and Stderr will be tracked in a log file.

        Output subdirectories will be created for each file type and the respective files will be relocated.


        Dependencies:
        All programs executed within this script must be added to the users path before execution.

        fastp version 0.21.0
        python3.7.3
        bowtie2 version 2.3.0
        bedtools v2.26.0
        The Genome Analysis Toolkit (GATK) v4.1.4.0:
            HTSJDK Version: 2.20.3
            Picard Version: 2.21.1
        openjdk version "1.8.0_222" (java)
            OpenJDK Runtime Environment (build 1.8.0_222-8u222-b10-1~deb9u1-b10)
            OpenJDK 64-Bit Server VM (build 25.222-b10, mixed mode)

    """

    # parse command line arguments
    parser = argparse.ArgumentParser(description = "ICICLE-seq preprocessing pipeline for fastq sequencing reads to fragments.tsv. Performs cell barcode corection of at most one low quality (Q score < 20) basecall based on the user provided barcode whitelist.",
        epilog = "Script written by Elliott Swanson of the Allen Institute for Immunology (elliott.swanson@alleninstitute.org). Free for academic use only. See LICENSE for code reuse permissions.")
    parser.add_argument("-1", "--read1", required = True, metavar = '', nargs = "+", help = "list of read 1 files in fastq.gz format")
    parser.add_argument("-2", "--read2", required = True, metavar = '', nargs = "+", help = "list of read 2 files in fastq.gz format")
    parser.add_argument("-o", "--output", required = True, metavar = '', help = "output directory")
    parser.add_argument("-w", "--whitelist", required = True, metavar = '', help = "barcode whitelist in .txt or .gz format with one barcode sequence per line")
    parser.add_argument("-r", "--reference", required = True, metavar = '', help = "bowtie2 reference genome index")
    args = parser.parse_args()

    # sort fastq filenames
    r1s = sorted(args.read1)
    r2s = sorted(args.read2)

    if len(r1s) != len(r2s):
        sys.exit(f"An equal number of fastq files must be supplied for read1 and read2. {len(r1s)} read 1 files and {len(r2s)} read 2 files were supplied. Exiting...")

    sample_name = os.path.basename(r1s[0]).split('_')[0]
    output_dir = os.path.abspath(args.output)
    whitelist = os.path.abspath(args.whitelist)
    ref_index = os.path.abspath(args.reference)

    # ensure path to user provided bowtie2 index is valid
    if not os.path.isdir(os.path.dirname(ref_index)):
        sys.exit(f'Cannot access bowtie2 index {ref_index}')

    # Assumes all fastq files exist and are in standard Illumina .fastq.gz format with standard naming convention: SampleName_S1_L001_R2_001.fastq.gz
    for i in range(len(r1s)):
        # ensure full path is included for each fastq file
        r1s[i] = os.path.abspath(r1s[i])
        r2s[i] = os.path.abspath(r2s[i])
       
        # Ensure all read 1 and read 2 fastq file pairs have the same sample name
        if os.path.basename(r1s[i]).split('_')[0] != sample_name or os.path.basename(r2s[i]).split('_')[0] != sample_name:
            base_names = ''.join([os.path.basename(file) + '\n' for file in r1s] + [os.path.basename(file) + '\n' for file in r2s])
            sys.exit(f'All fastq files must have the same sample name. FastQs provided: \n{base_names}\n')

        # ensure that each read is formatted with correct Illumina read designation
        if '_R1_' not in os.path.basename(r1s[i]) or '_R2_' not in os.path.basename(r2s[i]):
            sys.exit(f'FastQ files are not in standard Illumina naming format. Read number must be specified (R1 or R2)\n{r1s[i]}\n{r2s[i]}\n')

        # ensure that all user provided FastQ files exit
        if not os.path.isfile(r1s[i]):
            sys.exit(f'Cannot access read1 path {r1s[i]}')
        if not os.path.isfile(r2s[i]):
            sys.exit(f'Cannot access read2 path {r2s[i]}')

    # ensure user provided barcode whitelist exists and has a valid extension
    if whitelist.endswith('.gz'):
        gzipped = True
    elif whitelist.endswith('.txt'):
        gzipped = False
    else:
        sys.exit(f"Barcode whitelist must be in either .txt or .gz format.\nUnrecognized file extension {os.path.splitext(wl)[1]} from user provided whitelist {whitelist}. Exiting...")

    # if output directory directory doesn't exist, create it
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # ensure directory was successfully created (failure would likely be due to lack of permissions), exit if not
    if not os.path.isdir(output_dir):
        sys.exit(f'Failed to create output directory {output_dir}. Exiting...')

    # format log file
    logging.basicConfig(filename='{}/{}_ICICLE_analysis_{}.log'.format(output_dir, sample_name, datetime.now(pytz.timezone('US/Pacific')).strftime('%Y-%m-%d')), level=logging.INFO, format='%(levelname)s:%(asctime)s:%(message)s')

    # Create log header
    CPUs = os.cpu_count()
    USER = os.environ['USER']
    logging.info(f'Analysis script "ICICLE-seq_preprocessing.py" is being run by {USER} using {CPUs} CPUs.\n{output_dir} will be the output directory.\n')
    logging.info('The following reads will be aligned:')
    for i in range(len(r1s)):
        logging.info(f'{r1s[i]}')
        logging.info(f'{r2s[i]}\n')
    logging.info(f'User provided barcode whitelist: {whitelist}\n')

    # ensure current directory is output_dir
    os.chdir(output_dir)

    min_qual = 30
    sam_file = sample_name + '.sam'
    tagged_sam = sam_file.replace('.sam','_TAGGED.sam')
    bam_file = tagged_sam.replace('.sam','.bam')
    bed_name = bam_file.replace('.bam','.bed')
    trimmed_read2_files = []

    # trim 3' adapter from each read2 fastq utilizing multiple threads
    with concurrent.futures.ThreadPoolExecutor() as executor:
        for read2 in sorted(r2s):
            trimmed_fastq = f"{output_dir}/{os.path.basename(read2).replace('.fastq.gz','_TRIMMED.fastq.gz')}"
            trimmed_read2_files.append(trimmed_fastq)
            executor.submit(trim_read2,read2,trimmed_fastq,output_dir)

    # make SAM file
    exit_code = fastq_to_sam(ref_index, trimmed_read2_files, output_dir)
    if exit_code != 0:
        logging.error(f'fastq_to_sam for sample {sample_name} returned error. Skipping further processing')
        return(1)

    # create tagged, filtered SAM file
    clusters = tag_filter_sam(sam_file, r1s, min_qual, whitelist, gzipped)

    # make BAM file
    exit_code = sam_to_bam(tagged_sam)
    if exit_code != 0:
        logging.error('function sam_to_bam with arg {} returned error. Skipping further processing'.format(sam_file))
        return(1)

    # index BAM file
    exit_code =index_bam(bam_file)
    if exit_code != 0:
        logging.error('function index_bam with arg {} returned error. Skipping further processing'.format(bam_file))
        return(1)

    # make BED file
    exit_code = bam_to_bed(bam_file, bed_name)
    if exit_code != 0:
        logging.error('function bam_to_bed with arg {} returned error. Skipping further processing'.format(bam_file))
        return(1)

    # make fragments.tsv file
    exit_code = make_fragments_file(clusters, bed_name)
    if exit_code != 0:
        logging.error(f'function make_fragments_file with args {clusters} and {bed_name} returned error. Skipping further processing')
        return(1)

    # make fragments.tsv file containing UMI counts per position and cell barcodes combination
    exit_code = make_fragments_UMIcount_file(clusters, bed_name)
    if exit_code != 0:
        logging.error(f'function make_fragments_UMIcount_file with args {clusters} and {bed_name} returned error. Skipping further processing')
        return(1)

    return(0)

def trim_read2(read2,trimmed_fastq,output_dir):
    # trim ME/Poly-A adapter off 3' end of read 2 fastq file

    # ensure trimmed fastq file doesn't already exist. If it does, return an error
    if os.path.isfile(trimmed_fastq):
        logging.error(f'Trimmed read2 fastq file {trimmed_fastq} already exists. {read2} will not be processed further.')
        return(2)

    fastp_html = os.path.basename(read2).replace('.fastq.gz','.html')
    shell_command = f'fastp --in1={read2} --adapter_sequence=CTGTCTCTTATACACATCT --cut_tail --trim_poly_x --html={output_dir}/{fastp_html} --out1={trimmed_fastq}'
    logging.info(f"\nTimming 3' end of {read2} fastq file. Resulting output file will be {trimmed_fastq}\nshell command: \n{shell_command}\n")
    exit_code = call_command(shell_command)
    if exit_code != 0:
        logging.error(f'function trim_read2 using {read2} returned error an error. Skipping further processing')
        sys.exit(f'function trim_read2 using {read2} returned error an error. Skipping further processing')
    
    return(True)


def fastq_to_sam(ref_index, read2, output_dir):
    # Align list of read2 files against specified bowtie2 reference genomce index "ref_index"
    sam_file = '{}/{}.sam'.format(output_dir, os.path.basename(read2[0]).split('_')[0])

    # ensure SAM file doesn't already exist. If it does, return an error
    if os.path.isfile(sam_file):
        logging.error(f'SAM file {sam_file} already exists. {os.path.basename(read2[0]).split("_")[0]} will not be processed further.')
        return(2)

    shell_command = f'bowtie2 --local --sensitive --no-unal --phred33 --threads 32 -x {ref_index} -U {",".join(read2)} > {sam_file}'
    logging.info('\nAligning {} fastq file(s) with index {}\nshell command: \n{}\n'.format(' '.join(read2),ref_index,shell_command))

    # call shell command and log output
    try:
        output = subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)

        # format log message
        if len(output.stdout) == 0:
            log_command = 'Shell command {} returned code {}\n'.format(output.args, output.returncode)
        elif len(output.stdout) > 0:
            log_command = 'Shell command {} returned code {}. \nOutput: {}\n'.format(output.args, output.returncode, output.stdout)
        else:
            log_command = ''
        if len(log_command) > 0:
            logging.info(log_command)

        exit_code = 0
    except Exception as error:
        logging.exception(error)
        exit_code = 1
    
    return(exit_code)


def tag_filter_sam(sam, read1, min_qual, whitelist, gzipped):
    # tag input sam file with cell BC and UMI fields using read1 fastq
    # filter out low quality alignments with alignment score below "min_qual"
    # returns 'clusters' dictionary of all cluster IDs with cell barcode and UMI associations

    new_sam = sam.replace('.sam','_TAGGED.sam')

    # ensure tagged SAM file doesn't already exist. If it does, return an error
    if os.path.isfile(new_sam):
        logging.error(f'Tagged SAM file {new_sam} already exists. {os.path.basename(read1[0]).split("_")[0]} will not be processed further.')
        return(2)

    # read all high quality alignments into dict
    clusters = dict()
    with open(sam,'r') as sam_file:
        sam_reader = csv.reader(sam_file, delimiter = '\t')
        for line in sam_reader:
            low_qual = False
            if line[0].startswith('@') == False:
                # filter out low quality alignments
                if int(line[4]) < min_qual:
                    low_qual = True
                if low_qual == False:
                    clusters[line[0]] = dict()

    # read barcode whitelist into dictionary for fast lookups
    barcodes = dict()
    if gzipped == True:
        with gzip.open(whitelist, 'rt') as fr:
            for line in fr:
                barcodes[line.rstrip()] = ''
    elif gzipped == False:
        with open(whitelist, 'r') as fr:
            for line in fr:
                barcodes[line.rstrip()] = ''
    # dictionary for lookups of bases for barcode correction
    bases = {'A':('C','G','T'),'C':('G','T','A'),'G':('T','A','C'),'T':('A','C','G'),'N':['A','C','G','T']}

    # assign BC and UMI seq and quals for each aligned cluster
    for file in read1:
        with gzip.open(file,'rt') as read1_file:
            while True:
                # parse each read. If reading a line returns an empty string the end of file has been reached, break.
                header = read1_file.readline().lstrip('@')
                # if readline() == EOF, break
                if header == '':
                    break

                cluster = header.split()[0]
                seq = read1_file.readline().rstrip('\n')
                read1_file.readline()
                qual = read1_file.readline().rstrip('\n')
                barcode = seq[0:16]

                # barcode correction: if barcode seq is NOT in barcode whitelist, 
                # change base at each basecall and choose the match resulting from change at lowest quality basecall
                if barcode not in barcodes.keys():
                    # create list of barcode match / basecall quality scores
                    matches = []
                    s = list(barcode)
                    q = list(qual[0:16])
                    corrected = None
                    for i in range(16):
                        temp = s[i]
                        for base in bases[temp]:
                            s[i] = base
                            if ''.join(s) in barcodes.keys():
                                matches.append((''.join(s),ord(q[i])))
                            # reset s[i] to original base
                            s[i] = temp
                    # if at least one match has been found, assign barcode by the lowest quality score
                    if len(matches) > 0:
                        corrected = sorted(matches, key = lambda x: x[1])[0][0]
                            
                else:
                    corrected = barcode

                # if cluster ID is in dictionary, split sequence into cell BC and UMI and add to dictionary
                if cluster in clusters.keys():
                    clusters[cluster]['cbc'] = corrected
                    clusters[cluster]['bcs'] = seq[0:16]
                    clusters[cluster]['umis'] = seq[16:28]
                    clusters[cluster]['bcq'] = qual[0:16]
                    clusters[cluster]['umiq'] = qual[16:28]

    # CB Z Chromium cellular barcode sequence that is error-corrected and confirmed against a barcode whitelist..
    # CR Z Chromium cellular barcode sequence as reported by the sequencer.
    # CY Z Chromium cellular barcode read quality. Phred scores as reported by sequencer.
    # UR Z Chromium molecular barcode sequence as reported by the sequencer.
    # UY Z Chromium molecular barcode read quality. Phred scores as reported by sequencer.

    # read SAM into memory and add BC and UMI tags. Write output to "new_sam"
    with open(sam,'r') as sam_file:
        with open(new_sam,'w') as out_sam:
            sam_reader = csv.reader(sam_file, delimiter = '\t')
            sam_writer = csv.writer(out_sam, delimiter = '\t')
            for line in sam_reader:
                if line[0].startswith('@') == False:
                    cluster = line[0]
                    if cluster in clusters.keys() and clusters[cluster]['cbc'] != None:
                        new = [f"CB:Z:{clusters[cluster]['cbc']}",f"CR:Z:{clusters[cluster]['bcs']}",f"CY:Z:{clusters[cluster]['bcq']}",f"UR:Z:{clusters[cluster]['umis']}",f"UY:Z:{clusters[cluster]['umiq']}"]
                        line = line + new
                        sam_writer.writerow(line)
                else:
                    sam_writer.writerow(line)

    return(clusters)


def sam_to_bam(sam_file):
    # Convert input SAM file to a BAM file sorted by coordinate using picard via gatk

    # Ensure input SAM file exists
    if not os.path.isfile(sam_file):
        logging.error('SAM file {} not found'.format(sam_file))
        return(1)

    bam_file = sam_file.replace('.sam', '.bam')
    # ensure BAM file doesn't already exist. If it does, return an error
    if os.path.isfile(bam_file):
        logging.error('BAM file {} already exists. {} will not be processed further.'.format(bam_file, sam_file.rstrip('.sam')))
        return(2)

    # use gatk to sort SAM file by coordinate into BAM file
    shell_command = f'gatk SortSam --INPUT={sam_file} --OUTPUT={bam_file} --SORT_ORDER=coordinate'
    logging.info('Creating sorted BAM file {} from SAM file {}. Shell command:\n{}\n'.format(bam_file, sam_file, shell_command))
    exit_code = call_command(shell_command)

    return(exit_code)

def index_bam(bam_file):
    # use gatk Picard to build ".bai" BAM index

    # Ensure input BAM file exists
    if not os.path.isfile(bam_file):
        logging.error('BAM file {} not found'.format(bam_file))
        return(1)

    index_name = bam_file.replace('.bam','.bam.bai')
    # ensure BAM index file doesn't already exist. If it does, return an error
    if os.path.isfile(index_name):
        logging.error('BAM index file {} already exists. {} will not be processed further.'.format(index_name, bam_file.rstrip('.bam')))
        return(2)

    # use BAM file as input to create .bai BAM index
    logging.info('BAM .bai index file {} from BAM file {}'.format(index_name, bam_file))
    shell_command = f'gatk BuildBamIndex --INPUT={bam_file} --OUTPUT={index_name}'
    exit_code = call_command(shell_command)

    return(exit_code)

def bam_to_bed(bam_file, bed_name):
    # take sorted BAM file as input and output a BED file in single end format.

    # Ensure input BAM file exists
    if not os.path.isfile(bam_file):
        logging.error('BAM file {} not found'.format(bam_file))
        return(1)

    # ensure BED file doesn't already exist. If it does, return an error
    if os.path.isfile(bed_name):
        logging.error('BED file {} already exists. {} will not be processed further.'.format(bed_name, bam_file.rstrip('.bam')))
        return(2)

    # use BAM file as input to create sorted BED file
    logging.info('Creating BED file {} from BAM file {}'.format(bed_name, bam_file))
    # shell_command = 'bedtools bamtobed -i {} -tag MAPQ | bedtools sort -chrThenSizeA > {}'.format(bam_file, bed_name)
    shell_command = 'bedtools bamtobed -i {} > {}'.format(bam_file, bed_name)
    exit_code = call_command(shell_command)

    return(exit_code)


def make_fragments_file(cluster_dict, bed_file):
    """ uses dictionary of cluster ID, cell barcode, umi associations (generated by tag_filter_sam) and input BED file to format data into 
    a new 'fragments.tsv' file for downstream analysis.
    
    Format: chr start end cell_barcode umi strand """

    # Ensure input BED file exists
    if not os.path.isfile(bed_file):
        logging.error('BED file {} not found'.format(bed_file))
        return(1)

    fragments = bed_file.replace('.bed','_fragments.tsv')
    # ensure fragments.tsv file doesn't already exist. If it does, return an error
    if os.path.isfile(fragments):
        logging.error('Fragments file {} already exists. {} will not be processed further.'.format(fragments, bed_file.rstrip('.bed')))
        return(2)

    # read BED file one line at a time
    with open(bed_file,'r') as bedr:
        with open(fragments,'w') as frag:
            frag_writer = csv.writer(frag, delimiter = '\t')
            bed_reader = csv.reader(bedr, delimiter = '\t')
            for line in bed_reader:
                chrom,start,end,cluster,score,strand = line[0:6]
                frag_writer.writerow([chrom,start,end,cluster_dict[cluster]['cbc'],cluster_dict[cluster]['umis'],strand])
            
    return(0)

def make_fragments_UMIcount_file(cluster_dict, bed_file):
    """ uses dictionary of cluster ID, cell barcode, umi associations (generated by tag_filter_sam) and input BED file to format data into 
    a new 'fragments.tsv' file for downstream analysis.

    # make dict of pos/cellBC from BED and read1 fastq --> set of UMIs
    # write fragments file using pos/cellBC key and len(UMI set) as UMI count
    
    Format: chr start end cell_barcode umi_count strand """

    # Ensure input BED file exists
    if not os.path.isfile(bed_file):
        logging.error('BED file {} not found'.format(bed_file))
        return(1)

    fragments = bed_file.replace('.bed','_UMIcount_fragments.tsv')
    # ensure fragments.tsv file doesn't already exist. If it does, return an error
    if os.path.isfile(fragments):
        logging.error('Fragments file {} already exists. {} will not be processed further.'.format(fragments, bed_file.rstrip('.bed')))
        return(2)

    # make dict of pos/cellBC from BED and read1 fastq --> set of UMIs
    # read BED file one line at a time
    pos_bc_counts = dict()
    with open(bed_file,'r') as bedr:
        bed_reader = csv.reader(bedr, delimiter = '\t')
        for line in bed_reader:
            chrom,start,end,cluster,score,strand = line[0:6]
            pos_bc = chrom + start + end + cluster_dict[cluster]['cbc']
            if pos_bc in pos_bc_counts.keys():
                pos_bc_counts[pos_bc]['umi_set'].add(cluster_dict[cluster]['umis'])
            else:
                pos_bc_counts[pos_bc] = {'umi_set':{cluster_dict[cluster]['umis']}, 'flag': False}

        # read BED file one line at a time
    with open(bed_file,'r') as bedr:
        with open(fragments,'w') as frag:
            frag_writer = csv.writer(frag, delimiter = '\t')
            bed_reader = csv.reader(bedr, delimiter = '\t')
            for line in bed_reader:
                chrom,start,end,cluster,score,strand = line[0:6]
                pos_bc = chrom + start + end + cluster_dict[cluster]['cbc']

                # ensure each position / cell barcode region is writen only once
                if pos_bc_counts[pos_bc]['flag'] == False:
                    pos_bc_counts[pos_bc]['flag'] = True
                    frag_writer.writerow([chrom,start,end,cluster_dict[cluster]['cbc'],len(pos_bc_counts[pos_bc]['umi_set']),strand])
            
    return(0)

def call_command(shell_command):
    # call shell command and log output
    try:
        output = subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)

        # format log message
        if len(output.stdout) == 0 and len(output.stderr) == 0:
            log_command = 'Shell command {} returned code {}\n'.format(output.args, output.returncode)
        elif len(output.stdout) > 0 and len(output.stderr) == 0:
            log_command = 'Shell command {} returned code {}. \nOutput: {}\n'.format(output.args, output.returncode, output.stdout)
        elif len(output.stderr) > 0 and len(output.stdout) == 0:
            log_command = 'Shell command {} returned code {}. \nStderr: {}\n'.format(output.args, output.returncode, output.stderr)
        else:
            log_command = 'Shell command {} returned code {}. \nOutput: {}\nStderr: {}\n'.format(output.args, output.returncode, output.stdout, output.stderr)
        
        logging.info(log_command)
        return(0)
    except Exception as error:
        logging.exception(error)
        return(1)


if __name__ == '__main__':
    main()

    