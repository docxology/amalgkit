#!/bin/bash

# Set up logging
LOG_DIR="logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/analysis_$(date +%Y%m%d_%H%M%S).log"
exec 1> >(tee -a "$LOG_FILE")
exec 2>&1

# Create timestamp for this run
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Default settings
FORCE_RERUN=false
SKIP_CONFIRMATION=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --force)
      FORCE_RERUN=true
      shift
      ;;
    --yes)
      SKIP_CONFIRMATION=true
      shift
      ;;
    *)
      echo "Unknown option: $1"
      echo "Usage: $0 [--force] [--yes]"
      echo "  --force: Force rerun of all steps even if output files exist"
      echo "  --yes: Skip confirmation prompts"
      exit 1
      ;;
  esac
done

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Function to check if a command succeeded
check_status() {
    if [ $? -eq 0 ]; then
        log_message "SUCCESS: $1"
    else
        log_message "ERROR: $1"
        exit 1
    fi
}

# Function to check and install necessary tools
check_and_install_dependencies() {
    # Check if seqkit is installed
    if ! command -v seqkit &> /dev/null; then
        log_message "SeqKit not found. Attempting to install via Bioconda..."
        
        # Check if conda is installed
        if command -v conda &> /dev/null; then
            # Install seqkit using conda
            log_message "Installing SeqKit via conda..."
            conda install -y -c bioconda seqkit
            check_status "Installing SeqKit via conda"
        elif command -v mamba &> /dev/null; then
            # Install using mamba if available (faster than conda)
            log_message "Installing SeqKit via mamba..."
            mamba install -y -c bioconda seqkit
            check_status "Installing SeqKit via mamba"
        else
            # If conda is not available, try to download the binary directly
            log_message "Conda not found. Attempting to download SeqKit binary directly..."
            mkdir -p bin
            
            # Determine system architecture
            ARCH=$(uname -m)
            OS=$(uname -s | tr '[:upper:]' '[:lower:]')
            
            if [[ "$ARCH" == "x86_64" ]]; then
                SEQKIT_URL="https://github.com/shenwei356/seqkit/releases/download/v2.5.1/seqkit_linux_amd64.tar.gz"
                if [[ "$OS" == "darwin" ]]; then
                    SEQKIT_URL="https://github.com/shenwei356/seqkit/releases/download/v2.5.1/seqkit_darwin_amd64.tar.gz"
                fi
            elif [[ "$ARCH" == "arm64" ]] || [[ "$ARCH" == "aarch64" ]]; then
                SEQKIT_URL="https://github.com/shenwei356/seqkit/releases/download/v2.5.1/seqkit_linux_arm64.tar.gz"
                if [[ "$OS" == "darwin" ]]; then
                    SEQKIT_URL="https://github.com/shenwei356/seqkit/releases/download/v2.5.1/seqkit_darwin_arm64.tar.gz"
                fi
            else
                log_message "ERROR: Unsupported architecture: $ARCH"
                exit 1
            fi
            
            # Download and extract SeqKit
            log_message "Downloading SeqKit from $SEQKIT_URL"
            wget -O seqkit.tar.gz "$SEQKIT_URL"
            tar -xzf seqkit.tar.gz -C bin
            rm seqkit.tar.gz
            chmod +x bin/seqkit
            check_status "Downloading SeqKit binary"
            
            # Add to PATH
            export PATH="$PATH:$(readlink -f bin)"
            log_message "SeqKit installed at $(which seqkit)"
        fi
        
        # Verify installation
        if command -v seqkit &> /dev/null; then
            log_message "SeqKit installed successfully: $(seqkit version)"
        else
            log_message "ERROR: Failed to install SeqKit"
            exit 1
        fi
    else
        log_message "SeqKit already installed: $(which seqkit)"
    fi
}

# Function to display disk usage
show_disk_usage() {
    log_message "Current disk usage:"
    du -sh data/* 2>/dev/null || true
    du -sh results/* 2>/dev/null || true
}

# Function to create directory and confirm it exists
create_dir() {
    local dir_path="$1"
    local dir_name="$2"
    
    mkdir -p "$dir_path"
    
    if [ -d "$dir_path" ]; then
        log_message "Directory created/exists: $dir_name ($dir_path)"
    else
        log_message "ERROR: Failed to create directory: $dir_name ($dir_path)"
        exit 1
    fi
}

# Function to check if output files already exist
should_skip_step() {
    local step_name="$1"
    local output_file="$2"
    
    if [ "$FORCE_RERUN" = true ]; then
        log_message "Force rerun enabled, processing $step_name regardless of existing files"
        return 1
    fi
    
    if [ -f "$output_file" ]; then
        log_message "Output file for $step_name already exists: $output_file"
        log_message "Skipping $step_name step"
        return 0
    else
        return 1
    fi
}

# Create all required directories
log_message "Creating required directories..."
create_dir "data/metadata" "metadata"
create_dir "data/selected" "selected data"
create_dir "data/selected/metadata" "selected metadata"
create_dir "data/fastq" "fastq files"
create_dir "data/reference" "reference genomes and indices"
create_dir "bin" "binary tools"
create_dir "lib" "library dependencies"
create_dir "results" "results"
create_dir "results/quant" "quantification results"
create_dir "results/correlation" "correlation results"
create_dir "results/sanity" "sanity check results"

# Add bin directory to PATH
export PATH="$PATH:$(readlink -f bin)"
log_message "Updated PATH to include local bin directory: $(readlink -f bin)"

echo "==========================================="
log_message "Starting Amalgkit analysis pipeline for Pogonomyrmex transcriptomic data"
log_message "Timestamp: $TIMESTAMP"
if [ "$FORCE_RERUN" = true ]; then
    log_message "FORCE RERUN: All steps will be executed regardless of existing files"
fi
echo "==========================================="

# Step 0: Prepare reference genomes and Kallisto indices
echo "-------------------------------------------"
log_message "Step 0: Preparing transcriptomes and Kallisto indices..."
if [ -f "data/reference/references_prepared.txt" ] && [ "$FORCE_RERUN" = false ]; then
    log_message "Reference preparation marker file exists, skipping transcriptome preparation"
else
    if [ -f "scripts/fetch_transcriptomes.sh" ]; then
        cd scripts
        ./fetch_transcriptomes.sh
        cd ..
        check_status "Transcriptome preparation"
        
        # Check if reference preparation was successful
        if [ ! -f "data/reference/references_prepared.txt" ]; then
            log_message "ERROR: Transcriptome preparation did not complete successfully"
            exit 1
        fi
    else
        log_message "ERROR: Transcriptome preparation script not found at scripts/fetch_transcriptomes.sh"
        exit 1
    fi
fi
show_disk_usage

# Step 1: Retrieve metadata from NCBI SRA
echo "-------------------------------------------"
log_message "Step 1: Retrieving metadata from NCBI SRA..."
if should_skip_step "metadata retrieval" "data/metadata/metadata.tsv"; then
    log_message "Using existing metadata"
else
    log_message "Running amalgkit metadata command..."
    amalgkit metadata \
        --search_string "(Pogonomyrmex[Organism]) AND (RNA-seq[Strategy])" \
        --out_dir data \
        --entrez_email "Daniel@ActiveInference.Institute" \
        --redo yes | grep -v "Converting 0th sample from XML to DataFrame" | grep -v "Processing 0th sample"
    check_status "Metadata retrieval"
fi
show_disk_usage

# Step 2: Select SRA entries
echo "-------------------------------------------"
log_message "Step 2: Selecting SRA entries based on criteria..."
# Ensure metadata file exists before proceeding
if [ ! -f "data/metadata/metadata.tsv" ]; then
    log_message "ERROR: Metadata file not found at data/metadata/metadata.tsv"
    exit 1
fi

if should_skip_step "SRA selection" "data/selected/selected_metadata.tsv"; then
    log_message "Using existing SRA selection"
else
    # Copy metadata file to selected directory
    cp data/metadata/metadata.tsv data/selected/metadata/metadata.tsv
    check_status "Copying metadata to selected directory"

    log_message "Running amalgkit select command..."
    amalgkit select \
        --metadata data/selected/metadata/metadata.tsv \
        --config_dir config \
        --out_dir data/selected \
        --min_nspots 1000000 \
        --max_sample 10 \
        --mark_redundant_biosamples yes | grep -v "Converting 0th sample from XML to DataFrame" | grep -v "Processing 0th sample"
    check_status "SRA selection"

    # Create a selected_metadata.tsv file by filtering for is_sampled=yes
    log_message "Creating selected metadata file..."
    awk -F'\t' 'NR==1 || $9=="yes"' data/selected/metadata/metadata.tsv > data/selected/selected_metadata.tsv
    check_status "Creating selected metadata file"

    # Create a file with just the SRR IDs for getfastq
    log_message "Creating SRR ID list for download..."
    awk -F'\t' 'NR>1 && $9=="yes" {print $16}' data/selected/metadata/metadata.tsv > data/selected/srr_list.txt
    check_status "Creating SRR ID list"
fi

# Display selection summary
echo "-------------------------------------------"
log_message "Selection Summary:"
if [ -f "data/selected/selected_metadata.tsv" ]; then
    selected_count=$(grep -v "^scientific_name" data/selected/selected_metadata.tsv | wc -l)
    log_message "Selected entries: $selected_count"
    log_message "Species distribution:"
    cut -f1 data/selected/selected_metadata.tsv | grep -v "^scientific_name" | sort | uniq -c | sort -nr
    log_message "Tissue distribution:"
    cut -f3 data/selected/selected_metadata.tsv | grep -v "^tissue" | sort | uniq -c | sort -nr
else
    log_message "ERROR: Selected metadata file not created"
    exit 1
fi
show_disk_usage

# Ask for confirmation before proceeding
if [ "$SKIP_CONFIRMATION" = false ]; then
    echo "-------------------------------------------"
    read -p "Do you want to proceed with downloading and analyzing these selected entries? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        log_message "Analysis stopped by user"
        exit 1
    fi
else
    log_message "Skipping confirmation as --yes flag was provided"
fi

# Step 3: Download fastq files
echo "-------------------------------------------"
log_message "Step 3: Downloading fastq files..."
# Check for and install SeqKit if needed
check_and_install_dependencies

# Count existing fastq files
FASTQ_COUNT=$(find data/fastq -name "*.fastq.gz" 2>/dev/null | wc -l)
SRR_COUNT=$(wc -l < data/selected/srr_list.txt 2>/dev/null || echo 0)

if [ $FASTQ_COUNT -ge $SRR_COUNT ] && [ "$FORCE_RERUN" = false ]; then
    log_message "All FASTQ files appear to be already downloaded ($FASTQ_COUNT files found for $SRR_COUNT SRR entries)"
    log_message "Skipping fastq download step"
else
    # Make the smart downloader script executable
    chmod +x scripts/smart_downloader.py
    
    # Use the smart downloader
    log_message "Downloading and processing $SRR_COUNT SRA entries with smart downloader..."
    
    # Force redownload if --force flag was used
    FORCE_FLAG=""
    if [ "$FORCE_RERUN" = true ]; then
        FORCE_FLAG="--force"
    fi
    
    scripts/smart_downloader.py \
        --id_list data/selected/srr_list.txt \
        --out_dir data/fastq \
        --threads 8 \
        --metadata data/selected/selected_metadata.tsv \
        --amalgkit_path $(which amalgkit) \
        $FORCE_FLAG
        
    check_status "Fastq download"
fi
show_disk_usage

# Step 4: Quantify transcript abundance
echo "-------------------------------------------"
log_message "Step 4: Quantifying transcript abundance..."
# Check if any fastq files were downloaded
if [ -z "$(ls -A data/fastq 2>/dev/null)" ]; then
    log_message "ERROR: No fastq files found in data/fastq directory. Cannot proceed with quantification."
    exit 1
fi

# Check if transcriptomes were prepared
if [ ! -f "data/reference/references_prepared.txt" ]; then
    log_message "ERROR: Transcriptomes were not properly prepared. Cannot proceed with quantification."
    exit 1
fi

# Check if Kallisto indices were created
if [ -z "$(ls -A data/reference/kallisto_indices/*.idx 2>/dev/null)" ]; then
    log_message "ERROR: No Kallisto indices found. Cannot proceed with quantification."
    exit 1
fi

# Check if quantification has already been done
QUANT_FILES=$(find results/quant -name "abundance.tsv" 2>/dev/null | wc -l)
if [ $QUANT_FILES -gt 0 ] && [ "$FORCE_RERUN" = false ]; then
    log_message "Quantification results already exist ($QUANT_FILES abundance files found)"
    log_message "Skipping quantification step"
else
    log_message "Running amalgkit quant command..."
    amalgkit quant \
        --metadata data/selected/selected_metadata.tsv \
        --out_dir results/quant \
        --threads 8 \
        --config config/species_config.yaml \
        --use_kallisto_index yes \
        --kallisto_binary bin/kallisto | grep -v "Processing 0th sample" | grep -v "Converting 0th"
    check_status "Transcript quantification"
fi
show_disk_usage

# Step 5: Merge abundance tables
echo "-------------------------------------------"
log_message "Step 5: Merging abundance tables..."
# Check if quant results were generated
if [ -z "$(ls -A results/quant 2>/dev/null)" ]; then
    log_message "ERROR: No quantification results found in results/quant directory. Cannot proceed with merging."
    exit 1
fi

MERGED_FILE="results/merged_counts_${TIMESTAMP}.csv"
PREVIOUS_MERGED=$(find results -name "merged_counts_*.csv" 2>/dev/null | sort -r | head -n 1)

if [ -n "$PREVIOUS_MERGED" ] && [ "$FORCE_RERUN" = false ]; then
    log_message "Previous merged file exists: $PREVIOUS_MERGED"
    log_message "Using existing merged file instead of creating a new one"
    MERGED_FILE=$PREVIOUS_MERGED
else
    log_message "Running amalgkit merge command..."
    amalgkit merge \
        --input results/quant \
        --out_dir results \
        --prefix merged_counts_${TIMESTAMP} | grep -v "Processing 0th sample" | grep -v "Converting 0th"
    check_status "Table merging"
fi
show_disk_usage

# Step 6: Perform cross-species TMM normalization
echo "-------------------------------------------"
log_message "Step 6: Performing cross-species TMM normalization..."
# Check if merged results file exists
if [ ! -f "$MERGED_FILE" ]; then
    log_message "ERROR: Merged counts file not found. Cannot proceed with normalization."
    exit 1
fi

NORMALIZED_FILE="results/normalized_${TIMESTAMP}.csv"
PREVIOUS_NORMALIZED=$(find results -name "normalized_*.csv" 2>/dev/null | sort -r | head -n 1)

if [ -n "$PREVIOUS_NORMALIZED" ] && [ "$FORCE_RERUN" = false ]; then
    log_message "Previous normalized file exists: $PREVIOUS_NORMALIZED"
    log_message "Using existing normalized file instead of creating a new one"
    NORMALIZED_FILE=$PREVIOUS_NORMALIZED
else
    log_message "Running amalgkit cstmm command..."
    amalgkit cstmm \
        --input "$MERGED_FILE" \
        --out_dir results \
        --config config/species_config.yaml \
        --prefix normalized_${TIMESTAMP} | grep -v "Processing 0th sample" | grep -v "Converting 0th"
    check_status "TMM normalization"
fi
show_disk_usage

# Step 7: Curate data
echo "-------------------------------------------"
log_message "Step 7: Curating data..."
# Check if normalized results file exists
if [ ! -f "$NORMALIZED_FILE" ]; then
    log_message "ERROR: Normalized file not found. Cannot proceed with curation."
    exit 1
fi

CURATED_FILE="results/curated_${TIMESTAMP}.csv"
PREVIOUS_CURATED=$(find results -name "curated_*.csv" 2>/dev/null | sort -r | head -n 1)

if [ -n "$PREVIOUS_CURATED" ] && [ "$FORCE_RERUN" = false ]; then
    log_message "Previous curated file exists: $PREVIOUS_CURATED"
    log_message "Using existing curated file instead of creating a new one"
    CURATED_FILE=$PREVIOUS_CURATED
else
    log_message "Running amalgkit curate command..."
    amalgkit curate \
        --input "$NORMALIZED_FILE" \
        --out_dir results \
        --config config/curate_config.yaml \
        --prefix curated_${TIMESTAMP} | grep -v "Processing 0th sample" | grep -v "Converting 0th"
    check_status "Data curation"
fi
show_disk_usage

# Step 8: Perform cross-species correlation analysis
echo "-------------------------------------------"
log_message "Step 8: Performing cross-species correlation analysis..."
# Check if curated results file exists
if [ ! -f "$CURATED_FILE" ]; then
    log_message "ERROR: Curated file not found. Cannot proceed with correlation analysis."
    exit 1
fi

CORRELATION_FILES=$(find results/correlation -type f 2>/dev/null | wc -l)
if [ $CORRELATION_FILES -gt 0 ] && [ "$FORCE_RERUN" = false ]; then
    log_message "Correlation results already exist ($CORRELATION_FILES files found)"
    log_message "Skipping correlation analysis step"
else
    log_message "Running amalgkit csca command..."
    amalgkit csca \
        --input "$CURATED_FILE" \
        --out_dir results/correlation \
        --config config/species_config.yaml | grep -v "Processing 0th sample" | grep -v "Converting 0th"
    check_status "Correlation analysis"
fi
show_disk_usage

# Step 9: Run sanity checks
echo "-------------------------------------------"
log_message "Step 9: Running sanity checks..."
# Check if any results were generated
if [ -z "$(ls -A results 2>/dev/null)" ]; then
    log_message "ERROR: No results found in results directory. Cannot proceed with sanity checks."
    exit 1
fi

SANITY_FILES=$(find results/sanity -type f 2>/dev/null | wc -l)
if [ $SANITY_FILES -gt 0 ] && [ "$FORCE_RERUN" = false ]; then
    log_message "Sanity check results already exist ($SANITY_FILES files found)"
    log_message "Skipping sanity checks step"
else
    log_message "Running amalgkit sanity command..."
    amalgkit sanity \
        --input results \
        --out_dir results/sanity \
        --prefix check_${TIMESTAMP} | grep -v "Processing 0th sample" | grep -v "Converting 0th"
    check_status "Sanity checks"
fi
show_disk_usage

echo "==========================================="
log_message "Analysis pipeline completed successfully!"
log_message "Results are available in the results directory"
log_message "Log file: $LOG_FILE"
echo "===========================================" 