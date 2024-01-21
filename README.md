# PathogenTracker is designed to detect known human pathogens
# Pipeline takes raw fastq files from whole-genome sequencing (WGS) or bulk RNA sequencing data types and outputs list of inferred pathogens per sample. It consists of 2 main steps:
# First, run
./Geninput4pathodetect.sh

# Then, if data is WGS:
./Pathodetect4WGS.sh

# If RNA-seq:
./Pathodetect4bRNAseq.sh
