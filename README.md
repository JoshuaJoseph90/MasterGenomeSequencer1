# MasterGenomeSequencer1

## Overview
This repository provides a comprehensive guide to installing and setting up a host gene removal pipeline using the following tools:
- Bowtie2
- BWA
- Hocort
- Kraken2
- Scrubby 
<br /> 
The pipeline is specifically tailored for removing human (host) sequences from FASTQ files, optimizing for systems with limited storage (32 GB).
## Prerequisites
Before proceeding, ensure you have the following: 
<br /> 
- Operating System: A Linux distribution (e.g., Ubuntu 20.04 or later).
- User Permissions: Ability to install software and write to necessary directories (use sudo where required).
- Internet Connection: Required for downloading software and reference genomes.
- Terminal Access: Familiarity with basic Linux commands.

## Installing required programs
### Installing Conda
#### Download Miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/Downloads/miniconda.sh

#### Run the installer
bash ~/Downloads/miniconda.sh -b -p $HOME/miniconda

#### Initialize Conda
eval "$($HOME/miniconda/bin/conda shell.bash hook)"
conda init

#### Reload shell configuration
source ~/.bashrc

#### Verify Conda installation
conda --version

### Install bwa 
wget https://sourceforge.net/projects/bwa/files/bwa-0.7.17.tar.bz2 <br /> 
tar -xjf bwa-0.7.17.tar.bz2 <br /> 
cd bwa-0.7.17 <br /> 
make <br /> 
sudo mv bwa /usr/local/bin/ <br /> 
cd .. <br /> 
rm -rf bwa-0.7.17 bwa-0.7.17.tar.bz2 <br /> 

### Install Bowtie2
cd ~/Downloads <br /> 
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.5/bowtie2-2.4.5-linux-x86_64.zip <br /> 
unzip bowtie2-2.4.5-linux-x86_64.zip <br /> 
sudo mv bowtie2-2.4.5-linux-x86_64/* /usr/local/bin/ <br /> 
rm -rf bowtie2-2.4.5-linux-x86_64 bowtie2-2.4.5-linux-x86_64.zip <br /> 

### Install Hocort 
conda create -n hocort -c conda-forge -c bioconda hocort -y

### Install Scrubby 
conda create -n scrubby_env -c conda-forge -c bioconda scrubby -y

## Setting up Human Genome Reference and Indexes 

### Create directory for genomes
mkdir -p ~/genomes/human
cd ~/genomes/human

#### Download the human genome FASTA (GRCh38.p13)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz

#### Decompress the FASTA file
gunzip GCF_000001405.39_GRCh38.p13_genomic.fna.gz


### Bowtie2 Indexing
mkdir -p ~/genome_indexes/bowtie2<br /> 
bowtie2-build ~/genomes/human/GCF_000001405.39_GRCh38.p13_genomic.fna ~/genome_indexes/bowtie2/bowtie2_index

### BWA Indexing
mkdir -p ~/genome_indexes/bwa<br /> 
bwa index ~/genomes/human/GCF_000001405.39_GRCh38.p13_genomic.fna<br /> 
mv ~/genomes/human/GCF_000001405.39_GRCh38.p13_genomic.fna.* ~/genome_indexes/bwa/

## HMP sample in order to test benchmarks
#### Download the latest SRA Toolkit (Ubuntu 64-bit version)
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz <br /> 
tar -xvzf sratoolkit.current-ubuntu64.tar.gz
#### Change directory to the extracted folder (replace 'version' with actual folder name)
cd sratoolkit.*-ubuntu64
#### Add the binaries to your PATH (this is temporary, valid only for the current session)
export PATH=$(pwd)/bin:$PATH
#### Download paired-end FASTQ files
fastq-dump --split-files SRR12148457



## Executing Host removal script 
chmod +x ./host_removal.sh
./host_removal.sh ~/path/to/sample_reads.fastq




