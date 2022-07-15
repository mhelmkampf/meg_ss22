### ================================= ###
### 06 Genome sequencing and assembly ###
### ================================= ###


## Preparations

# Change into project directory
cd ~/meg_ss22   # change to path (directory location) on your system

# Update Git repository
git pull

# On error, remove old project directory and download again
# git clone https://github.com/mhelmkampf/meg_ss22

# Software installation (macOS only)
# python3 -m pip install --user --upgrade cutadapt

# Additional software is provided in meg_ss22/apps



## Navigating the file system
ls -l
cd data/genome
pwd
cd ..   # .. stands for parent directory


## Raw data processing

# Check data integrity using md5 hash
md5sum Hpue_raw300_F.fastq.gz
cat md5/Hpue_raw300_F.fastq.gz.md5

# Decompress Fastq file and print header to screen
gzip -d Hpue_raw300_F.fastq.gz
head -n 4 Hpue_raw300_F.fastq

# Count number of reads in file
grep '@HW' Hpue_raw300_F.fastq | wc -l   # print all lines containing '@HW', then count

# Assess read quality with Fastqc
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc
cd ../../apps
unzip fastqc_v0.11.9.zip

# Run Fastqc
cd FastQC
./fastqc -h   # . stands for current directory
./fastqc ../../data/genome/Hpue_raw300_F.fastq

# Remove adapters with cutadapt
cd ..
./cutadapt-3.4.exe -help
./cutadapt-3.4.exe -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -q 10 -o ../data/genome/Hpue_raw300_F_trim.fastq ../data/genome/Hpue_raw300_F.fastq

# Rerun Fastqc
cd FastQC
./fastqc ../../data/genome/Hpue_raw300_F_trim.fastq
cd ../../data/genome && ls -l   # change directory and list contents in one command
mkdir qc
mv *fastqc* qc                  # move all files containing 'fastqc' in the file name to qc directory


## Genome assembly
zcat < Hpue_genome_LG12Mc.fas.gz | head -n 100   # read in compressed fasta file, then print first 100 lines
zcat < Hpue_genome_LG12Mc.fas.gz | tail          # read in compressed fasta file, then print last 10 (default) lines
zcat < Hpue_genome_LG12Mc.fas.gz | grep '>'      # print all fasta headers (sequence ids)

# Quality assessment (QUAST)
# Example data: http://cab.cc.spbu.ru/quast/reports/02_Jun_2022_10:18:17_569501/report.html
# H. puella reference genome: http://cab.cc.spbu.ru/quast/reports/02_Jun_2022_10:19:29_143981/report.html
