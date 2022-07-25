#!/bin/bash

########### File is designed to detect pathogens and extra non-human organisms in NGS Sammir's RNA-Seq DLBCL samples
########### File uses Bracken/Kraken to identify number of aligned reads to non-human species and puts them in separate fastq file. 
########### Make sure to run RNA-Seq or WGS pipeline before non-human species detection.
########## Set required variables 

export NTHREADS=16  ### number of threads
export path_to_fastq=/home/cluster/abalay/scratch/Sammirs_RNA_Seq/Cutadapt ### path to gzipped FASTQ files preferably trimmed, run rna_seq_sample_preprocessing.sh

########### Do not change these variables
export script_dir=/home/cluster/abalay/scratch/Pipelines_and_Files  ### directory to packages and extra scripts
export kraken2=$script_dir/Packages/kraken2   ### path to kraken2
export bracken=$script_dir/Packages/Bracken/bracken   ### path to bracken
export bracken_build=$script_dir/Packages/Bracken/bracken-build   ### path to bracken-build
export kreport2mpa=$script_dir/Packages/KrakenTools/kreport2mpa.py ### path to MPA style converter from Bracken report  
export combine_mpa=$script_dir/Packages/KrakenTools/combine_mpa.py ### path to combine MPA TEXT output
export extract_kraken_reads=$script_dir/Packages/KrakenTools/extract_kraken_reads.py  ### path to read  extractor from Kraken2 output

########### Part 1: Bracken/Kraken
########### Input: processed FASTQ files in Cutadapt
########### Output: Bracken-generated reports in Kraken2 manner with read percentages classified to each species

## Step 0: Build Kraken2.0 database: includes archaea, bacteria, fungi, plant, human, protozoa, virus genome collection + plasmids + vectors, adapters, linkers, primer sequences from NCBI 
[ ! -d $script_dir/kraken2_db ] && $script_dir/Files_for_scripts/db_kraken2_build.sh

### Set extra variables 
export KRAKEN2_DB=$script_dir/kraken2_db ### path to kraken2 database
export KMER_LEN=35 ### length of kmer used to build the Kraken2 database
export READ_LEN=35 ### minimum ideal read length in samples from NGS sequencing

## Step 1: Generate Bracken database file (database${READ_LEN}mers.kraken)
if [ ! -e $script_dir/kraken2_db/database35mers.kmer_distrib ] && [ ! -e $script_dir/kraken2_db/database35mers.kraken ]; then
    $bracken_build -d $KRAKEN2_DB  -t $NTHREADS -k $KMER_LEN -l $READ_LEN -x $kraken2
fi


## Step 2: Generate Kraken2.0 report file

### Create directory for Kraken2 reports and move all of them to single subfolder Kraken2_output
mkdir -p $path_to_fastq/Kraken2_output


for file in $(find $path_to_fastq -name "*_cutadapt.fastq.gz")
do
    curr=$(basename $file _cutadapt.fastq.gz)
        if [ "${curr: -2}" == _1 ]; then
            echo "Detected paired-end read file"
            currp=$(basename $file _1_cutadapt.fastq.gz)
            if [ -e $path_to_fastq/"$currp"_2_cutadapt.fastq.gz ] && [ ! -e $path_to_fastq/Kraken2_output/"$currp".kraken2 ]; then     
                        echo "Generating kraken2 report for the paired FASTQ files"
                        $kraken2/kraken2 --db $KRAKEN2_DB --threads $NTHREADS --report $path_to_fastq/Kraken2_output/"$currp".k2report --use-names --paired --gzip-compressed $path_to_fastq/"$currp"_1_cutadapt.fastq.gz
                        $path_to_fastq/"$currp"_2_cutadapt.fastq.gz > $path_to_fastq/Kraken2_output/"$currp".kraken2 
            fi
        fi
        if [ ! "${curr: -2}" = _1 ] && [ ! "${curr: -2}" = _2 ] && [ ! -e $path_to_fastq/Kraken2_output/"$curr".kraken2 ]; then
            echo "Detected single-end read file"
            echo "Generating kraken2 report for the single FASTQ file"
            $kraken2/kraken2 --db $KRAKEN2_DB --threads $NTHREADS --report $path_to_fastq/Kraken2_output/"$curr".k2report --use-names --gzip-compressed $path_to_fastq/"$curr"_cutadapt.fastq.gz > $path_to_fastq/Kraken2_output/"$curr".kraken2
        fi
done


## Step 3: Run Bracken for abundance estimation
for file in $(find $path_to_fastq/Kraken2_output -name "*.k2report")
do
    curr=$(basename $file .k2report)
    if [ ! -e $path_to_fastq/Kraken2_output/"$curr".breport ];then
        $bracken -d $KRAKEN2_DB -i $path_to_fastq/Kraken2_output/"$curr".k2report -o $path_to_fastq/Kraken2_output/"$curr".bracken -r $READ_LEN -w $path_to_fastq/Kraken2_output/"$curr".breport -l 'S' -t 10
    fi
done

## Step 4: Generate MPA style files from Bracken reports and at the end combine them
### Create a string to use it for combining MPA style TEXT format documents
string=""
for file in $(find $path_to_fastq/Kraken2_output -name "*.breport")
do
    curr=$(basename $file .breport)
    if [ ! -e $path_to_fastq/Kraken2_output/"$curr".mpa.txt ];then  
        $kreport2mpa -r $path_to_fastq/Kraken2_output/"$curr".breport -o $path_to_fastq/Kraken2_output/"$curr".mpa.txt --no-intermediate-ranks --read_count  --display-header
    fi
    string="${string} "$path_to_fastq"/Kraken2_output/"$curr".mpa.txt "
done


arr=( $string )

### Combine MPA TEXT files
cd $path_to_fastq/Kraken2_output
if [ ! -e $path_to_fastq/Kraken2_output/Sammirs_RNA_Seq.patho.combined.mpa.txt ];then
    $combine_mpa -i "${arr[@]}" -o $path_to_fastq/Kraken2_output/Sammirs_RNA_Seq.patho.combined.mpa.txt
fi

## Step 5: Extract reads unaligned to human genome

### Create directory tos tore unmapped fastq files
mkdir -p $path_to_fastq/Unmapped


for file in $(find $path_to_fastq -name "*_cutadapt.fastq.gz")
do
    curr=$(basename $file _cutadapt.fastq.gz)
        if [ "${curr: -2}" == _1 ]; then
            echo "Detected paired-end read file"
            currp=$(basename $file _1_cutadapt.fastq.gz)
            if [ -e $path_to_fastq/"$currp"_2_cutadapt.fastq.gz ] && [ ! -e $path_to_fastq/Unmapped/"$curr"_r1_unmapped.fastq ] || [ ! -e $path_to_fastq/Unmapped/"$curr"_r2_unmapped.fastq ]; then
                python3 $extract_kraken_reads -k $path_to_fastq/Kraken2_output/"$currp".kraken2 -1 $path_to_fastq/"$currp"_1_cutadapt.fastq.gz -2 $path_to_fastq/"$currp"_2_cutadapt.fastq.gz \
                -o $path_to_fastq/Unmapped/"$curr"_r1_unmapped.fastq -o2 $path_to_fastq/Unmapped/"$curr"_r2_unmapped.fastq --taxid 9606 --exclude
            fi
        fi
        if [ ! "${curr: -2}" = _1 ] && [ ! "${curr: -2}" = _2 ] && [ ! -e $path_to_fastq/Unmapped/"$curr"_unmapped.fastq ]; then
            echo "Detected single-end read file"
            python3 $extract_kraken_reads -k $path_to_fastq/Kraken2_output/"$curr".kraken2 -1 $path_to_fastq/"$curr"_cutadapt.fastq.gz -o $path_to_fastq/Unmapped/"$curr"_unmapped.fastq --taxid 9606 --exclude
        fi
done





