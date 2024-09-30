#!/bin/bash

# =============================================================================
# Script: host_removal.sh
# Description: Removes human (host) sequences from a FASTQ file using Bowtie2,
#              BWA, Hocort, and Scrubby. Outputs cleaned FASTQ files and captures
#              execution metrics for comparison.
# Usage: ./host_removal.sh input.fastq
# =============================================================================

# Exit immediately if a command exits with a non-zero status,
# if an undefined variable is used, and if any command in a pipeline fails
set -euo pipefail

# -----------------------------------------------------------------------------
# Function: print_usage
# Description: Prints usage instructions.
# -----------------------------------------------------------------------------
print_usage() {
    echo "Usage: $0 input.fastq"
    echo "Example: $0 sample_reads.fastq"
}

# -----------------------------------------------------------------------------
# Function: check_dependencies
# Description: Checks if required tools are installed.
# -----------------------------------------------------------------------------
check_dependencies() {
    local dependencies=("bowtie2" "bwa" "samtools" "conda")
    for cmd in "${dependencies[@]}"; do
        if ! command -v "$cmd" &> /dev/null; then
            echo "Error: $cmd is not installed or not in PATH."
            exit 1
        fi
    done
}

# -----------------------------------------------------------------------------
# Function: run_command
# Description: Runs a command, captures execution metrics, and logs them.
# Arguments:
#   $1 - Tool name
#   $2 - Command to execute
#   $3 - Output file path
# -----------------------------------------------------------------------------
run_command() {
    local TOOL=$1
    local CMD=$2
    local OUTFILE=$3
    local METRIC_LOG="$OUTFILE.metrics.log"
    
    echo "[$(date)] Starting $TOOL..."
    
    # Run the command with /usr/bin/time to capture metrics
    /usr/bin/time -f "Time:%e,User:%U,Sys:%S,MaxMem:%M" bash -c "$CMD" 2> "$METRIC_LOG"
    
    # Extract metrics
    local TIME_REAL=$(grep "Time:" "$METRIC_LOG" | cut -d',' -f1 | cut -d':' -f2)
    local TIME_USER=$(grep "User:" "$METRIC_LOG" | cut -d',' -f2 | cut -d':' -f2)
    local TIME_SYS=$(grep "Sys:" "$METRIC_LOG" | cut -d',' -f3 | cut -d':' -f2)
    local MAX_MEM=$(grep "MaxMem:" "$METRIC_LOG" | cut -d',' -f4 | cut -d':' -f2)
    
    # Append to metrics summary
    echo "$TOOL,$TIME_REAL,$TIME_USER,$TIME_SYS,$MAX_MEM" >> "$METRICS_FILE"
    
    echo "[$(date)] $TOOL completed. Metrics: Time=${TIME_REAL}s, User=${TIME_USER}s, Sys=${TIME_SYS}s, MaxMem=${MAX_MEM}KB"
    
    # Remove metric log
    rm "$METRIC_LOG"
}

# -----------------------------------------------------------------------------
# Main Script Execution Starts Here
# -----------------------------------------------------------------------------

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Error: Incorrect number of arguments."
    print_usage
    exit 1
fi

INPUT_FASTQ_1=$1
INPUT_FASTQ_2=$2

# Check if the input FASTQ file exists
if [ ! -f "$INPUT_FASTQ_1" ] || [ ! -f "$INPUT_FASTQ_2" ]; then
    echo "Error: Input FASTQ files '$INPUT_FASTQ_1' or '$INPUT_FASTQ_2' not found."
    exit 1
fi

# Check for required dependencies
check_dependencies

# Source conda to allow environment activation
# Adjust the path to conda.sh based on your Conda installation
if [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
    source "/opt/conda/etc/profile.d/conda.sh"
else
    echo "Error: conda.sh not found. Please adjust the path in the script."
    exit 1
fi


# Define paths to indexes and reference
# ============================== IMPORTANT ==============================
# Replace the following placeholder paths with the actual paths on your system.
# Ensure that Bowtie2 and BWA indexes are built for the human genome.
# =============================================================================
BOWTIE2_INDEX="/home/codespace/genome_indexes/bowtie2/bowtie2_index.fa"  # e.g., ~/kraken2_custom_db/bowtie2_index
BWA_INDEX="/home/codespace/genomes/GCF_000001405.39_GRCh38.p13_genomic.fna"          # e.g., ~/kraken2_custom_db/bwa_index
HOST_REFERENCE="/home/codespace/genomes/GCF_000001405.39_GRCh38.p13_genomic.fna"  # e.g., ~/kraken2_custom_db/human_reference/GCF_000001405.39_GRCh38.p13_genomic.fna
SCRUBBY_DB_DIR="/home/codespace/genomes/GCF_000001405.39_GRCh38.p13_genomic.fna"    # e.g., ~/scrubby_db
# =============================================================================

# Validate that the defined paths exist
#if [ ! -f "$BOWTIE2_INDEX.1.bt2" ]; then
#    echo "Error: Bowtie2 index files not found at '$BOWTIE2_INDEX'."
#    exit 1
#fi

if [ ! -f "$BWA_INDEX.bwt" ]; then
    echo "Error: BWA index files not found at '$BWA_INDEX'."
    exit 1
fi

if [ ! -f "$HOST_REFERENCE" ]; then
    echo "Error: Human reference FASTA file not found at '$HOST_REFERENCE'."
    exit 1
fi

#if [ ! -d "$SCRUBBY_DB_DIR" ]; then
#    echo "Error: Scrubby database directory not found at '$SCRUBBY_DB_DIR'."
#    exit 1
#fi

# Create output directory
OUTPUT_DIR="host_removal_outputs"
mkdir -p "$OUTPUT_DIR"

# Initialize metrics summary CSV
METRICS_FILE="$OUTPUT_DIR/metrics_summary.csv"
echo "Tool,Time(s),User(s),Sys(s),MaxMem(KB)" > "$METRICS_FILE"

# -----------------------------------------------------------------------------
# Host Removal Using Bowtie2
# -----------------------------------------------------------------------------
#BOWTIE2_OUTPUT="$OUTPUT_DIR/bowtie2_removed.fastq"
#BOWTIE2_CMD="bowtie2 -x $BOWTIE2_INDEX -U $INPUT_FASTQ --un $BOWTIE2_OUTPUT -S /dev/null"

#run_command "Bowtie2" "$BOWTIE2_CMD" "$BOWTIE2_OUTPUT"

# -----------------------------------------------------------------------------
# Host Removal Using BWA
# -----------------------------------------------------------------------------
BWA_OUTPUT="$OUTPUT_DIR/bwa_removed.fastq"
BWA_CMD="bwa mem $BWA_INDEX $INPUT_FASTQ_1 $INPUT_FASTQ_2 | samtools view -f 4 -b | samtools fastq - > $BWA_OUTPUT"

run_command "BWA" "$BWA_CMD" "$BWA_OUTPUT"

# -----------------------------------------------------------------------------
# Host Removal Using Hocort
# -----------------------------------------------------------------------------
# Define output file names
HOCORT_OUTPUT_1="$OUTPUT_DIR/hocort_rm_1.fastq"
HOCORT_OUTPUT_2="$OUTPUT_DIR/hocort_rm_2.fastq"

# Activate the 'hocort' Conda environment
conda activate hocort

# Adjust the HOCORT_CMD to use the correct 'map' subcommand and bowtie2 index
HOCORT_CMD="hocort map bowtie2 -x $HOST_REFERENCE \
-i $INPUT_FASTQ_1 $INPUT_FASTQ_2 \
-o $HOCORT_OUTPUT_1 $HOCORT_OUTPUT_2"

# Run the hocort command
run_command "Hocort" "$HOCORT_CMD" "$HOCORT_OUTPUT_1 $HOCORT_OUTPUT_2"


# -----------------------------------------------------------------------------
# Host Removal Using Scrubby
# -----------------------------------------------------------------------------
#SCRUBBY_OUTPUT="$OUTPUT_DIR/scrubby_removed.fastq"

#SCRUBBY_CMD="scrubby filter --database $SCRUBBY_DB_DIR --input $INPUT_FASTQ_1 --input2 $INPUT_FASTQ_2 --output $SCRUBBY_OUTPUT --mode fastp"

#run_command "Scrubby" "$SCRUBBY_CMD" "$SCRUBBY_OUTPUT"

# Deactivate the 'hocort' Conda environment
#conda deactivate

# -----------------------------------------------------------------------------
# Final Metrics Summary
# -----------------------------------------------------------------------------
echo ""
echo "======================== Metrics Summary ========================"
cat "$METRICS_FILE"
echo "================================================================="
echo "Host removal using all tools completed successfully."
echo "Outputs are saved in the '$OUTPUT_DIR' directory."
