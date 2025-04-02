#!/bin/bash

# Script to download transcriptomes from NCBI and build Kallisto indices

# Set up logging
LOG_DIR="../logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/transcriptome_fetch_$(date +%Y%m%d_%H%M%S).log"
exec 1> >(tee -a "$LOG_FILE")
exec 2>&1

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

# Create required directories
mkdir -p ../data/reference/transcriptomes
mkdir -p ../data/reference/kallisto_indices

# Species to process - add all target species here with their RefSeq accessions
declare -A SPECIES
SPECIES["Pogonomyrmex barbatus"]="GCF_000187915.1"
SPECIES["Pogonomyrmex californicus"]="GCF_023869295.1"
SPECIES["Camponotus floridanus"]="GCF_003227725.1"
SPECIES["Solenopsis invicta"]="GCF_000188075.2"
SPECIES["Harpegnathos saltator"]="GCF_003227715.1"

# Make sure we have required tools
log_message "Checking for required tools..."

# Check for kallisto
if ! command -v kallisto &> /dev/null; then
    log_message "Installing kallisto..."
    # Download and install latest kallisto release to local bin
    mkdir -p ../bin
    wget -O kallisto.tar.gz 'https://github.com/pachterlab/kallisto/releases/download/v0.48.0/kallisto_linux-v0.48.0.tar.gz'
    tar -xzf kallisto.tar.gz
    cp kallisto*/kallisto ../bin/
    chmod +x ../bin/kallisto
    rm -rf kallisto.tar.gz kallisto*
    export PATH="$PATH:$(readlink -f ../bin)"
    check_status "Install kallisto"
fi

# Function to process a species using direct RefSeq download
fetch_transcriptome() {
    local species_name="$1"
    local accession="$2"
    local safe_name=$(echo "$species_name" | tr ' ' '_')
    
    # Check if files already exist
    if [ -f "../data/reference/transcriptomes/${safe_name}_transcripts.fa" ] && 
       [ -f "../data/reference/kallisto_indices/${safe_name}.idx" ]; then
        log_message "Transcriptome and Kallisto index already exist for $species_name. Skipping."
        
        # Check if entry exists in species_config.yaml
        local entry_exists=$(python3 -c "
import yaml
try:
    with open('../config/species_config.yaml', 'r') as f:
        config = yaml.safe_load(f)
    for species in config.get('species', []):
        if species.get('name') == '$species_name':
            print('true')
            break
except:
    pass
")
        if [ "$entry_exists" = "true" ]; then
            log_message "Configuration entry already exists for $species_name."
            return 0
        else
            log_message "Updating species configuration for $species_name..."
            # Add entry to species_config.yaml
            python3 -c "
import yaml
import sys
try:
    with open('../config/species_config.yaml', 'r') as f:
        config = yaml.safe_load(f)
    
    if 'species' not in config:
        config['species'] = []
    
    new_species = {
        'name': '$species_name',
        'transcriptome': 'data/reference/transcriptomes/${safe_name}_transcripts.fa',
        'kallisto_index': 'data/reference/kallisto_indices/${safe_name}.idx',
        'accession': '$accession'
    }
    
    # Check if species already exists
    exists = False
    for species in config['species']:
        if species.get('name') == '$species_name':
            exists = True
            break
    
    if not exists:
        config['species'].append(new_species)
    
    with open('../config/species_config.yaml', 'w') as f:
        yaml.dump(config, f, default_flow_style=False)
        
except Exception as e:
    print(f'Error updating YAML: {e}', file=sys.stderr)
    sys.exit(1)
"
            check_status "Update species configuration for $species_name"
            return 0
        fi
    fi
    
    log_message "Processing transcriptome data for $species_name (accession: $accession)..."
    
    # Build the FTP URL for the transcriptome
    local refseq_dir=$(echo "$accession" | cut -d '.' -f 1)
    local gc_prefix=$(echo "$refseq_dir" | cut -c 1-3)
    local ftp_path="https://ftp.ncbi.nlm.nih.gov/genomes/all/${gc_prefix}/${refseq_dir}/${accession}_*"
    
    # Use wget to list the directory
    log_message "Looking up assembly directory..."
    local assembly_dir=$(wget -q -O - "$ftp_path" | grep -o "href=\"[^\"]*\"" | cut -d'"' -f2 | grep "$accession")
    
    if [ -z "$assembly_dir" ]; then
        log_message "Trying alternate directory structure..."
        ftp_path="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/${safe_name}/all_assembly_versions/${accession}"
        assembly_dir=$(wget -q -O - "$ftp_path" | grep -o "href=\"[^\"]*\"" | cut -d'"' -f2 | grep "$accession")
    fi
    
    if [ -z "$assembly_dir" ]; then
        log_message "Could not find assembly directory. Trying direct path..."
        # Try a more direct approach with a predictable path pattern
        assembly_dir="${accession}_latest_rna.fna.gz"
        ftp_path="https://ftp.ncbi.nlm.nih.gov/genomes/all/${gc_prefix}/${refseq_dir}/${accession}"
    fi
    
    # Download the RNA FASTA file
    local rna_url="${ftp_path}/${assembly_dir}/${accession}_rna.fna.gz"
    log_message "Attempting to download from: $rna_url"
    
    if ! wget -O "../data/reference/transcriptomes/${safe_name}_transcripts.fa.gz" "$rna_url"; then
        log_message "Failed with direct URL. Trying alternative URL patterns..."
        
        # Try alternative URL patterns
        rna_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/${gc_prefix}/${refseq_dir}/${accession}/${accession}_rna.fna.gz"
        if ! wget -O "../data/reference/transcriptomes/${safe_name}_transcripts.fa.gz" "$rna_url"; then
            
            # Try another pattern
            rna_url="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/$(echo "$species_name" | tr ' ' '_' | tr '[:upper:]' '[:lower:]')/latest_assembly_versions/${accession}/${accession}_rna.fna.gz"
            if ! wget -O "../data/reference/transcriptomes/${safe_name}_transcripts.fa.gz" "$rna_url"; then
                log_message "Failed to download transcriptome. Using hardcoded URLs as fallback."
                
                # Hardcoded URLs for specific species as fallback
                case "$species_name" in
                    "Pogonomyrmex barbatus")
                        wget -O "../data/reference/transcriptomes/${safe_name}_transcripts.fa.gz" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/187/915/GCF_000187915.1_Pbar_UMD_V03/GCF_000187915.1_Pbar_UMD_V03_rna.fna.gz"
                        ;;
                    "Pogonomyrmex californicus")
                        wget -O "../data/reference/transcriptomes/${safe_name}_transcripts.fa.gz" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/869/295/GCF_023869295.1_ASM2386929v1/GCF_023869295.1_ASM2386929v1_rna.fna.gz"
                        ;;
                    "Camponotus floridanus")
                        wget -O "../data/reference/transcriptomes/${safe_name}_transcripts.fa.gz" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/227/725/GCF_003227725.1_Cflo_v7.5/GCF_003227725.1_Cflo_v7.5_rna.fna.gz"
                        ;;
                    "Solenopsis invicta")
                        wget -O "../data/reference/transcriptomes/${safe_name}_transcripts.fa.gz" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/075/GCF_000188075.2_Si_gnH/GCF_000188075.2_Si_gnH_rna.fna.gz"
                        ;;
                    "Harpegnathos saltator")
                        wget -O "../data/reference/transcriptomes/${safe_name}_transcripts.fa.gz" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/227/715/GCF_003227715.1_Hsal_v8.5/GCF_003227715.1_Hsal_v8.5_rna.fna.gz"
                        ;;
                    *)
                        log_message "No fallback URL for $species_name"
                        return 1
                        ;;
                esac
            fi
        fi
    fi
    
    # Decompress the downloaded file
    gunzip -f "../data/reference/transcriptomes/${safe_name}_transcripts.fa.gz"
    
    # Check if download was successful
    if [ ! -f "../data/reference/transcriptomes/${safe_name}_transcripts.fa" ]; then
        log_message "Failed to download transcriptome for $species_name"
        return 1
    fi
    
    # Build Kallisto index
    log_message "Building Kallisto index for $species_name..."
    ../bin/kallisto index -i "../data/reference/kallisto_indices/${safe_name}.idx" "../data/reference/transcriptomes/${safe_name}_transcripts.fa"
    check_status "Build Kallisto index for $species_name"
    
    # Update the species_config.yaml with correct paths
    log_message "Updating species configuration for $species_name..."
    python3 -c "
import yaml
import sys

try:
    with open('../config/species_config.yaml', 'r') as f:
        config = yaml.safe_load(f)
    
    # Find the species entry
    found = False
    for species in config.get('species', []):
        if species.get('name') == '$species_name':
            species['transcriptome'] = 'data/reference/transcriptomes/${safe_name}_transcripts.fa'
            species['kallisto_index'] = 'data/reference/kallisto_indices/${safe_name}.idx'
            species['accession'] = '$accession'
            found = True
            break
    
    # If species not found, add it
    if not found:
        if 'species' not in config:
            config['species'] = []
        
        new_species = {
            'name': '$species_name',
            'transcriptome': 'data/reference/transcriptomes/${safe_name}_transcripts.fa',
            'kallisto_index': 'data/reference/kallisto_indices/${safe_name}.idx',
            'accession': '$accession'
        }
        config['species'].append(new_species)
    
    # Write updated config
    with open('../config/species_config.yaml', 'w') as f:
        yaml.dump(config, f, default_flow_style=False)
        
except Exception as e:
    print(f'Error updating YAML: {e}', file=sys.stderr)
    sys.exit(1)
"
    check_status "Update species configuration for $species_name"
    
    log_message "Completed processing for $species_name"
    return 0
}

# Process each species
log_message "Starting transcriptome fetch and Kallisto index building..."
processed_count=0
skipped_count=0

for species in "${!SPECIES[@]}"; do
    # Check if files already exist before processing
    safe_name=$(echo "$species" | tr ' ' '_')
    if [ -f "../data/reference/transcriptomes/${safe_name}_transcripts.fa" ] && 
       [ -f "../data/reference/kallisto_indices/${safe_name}.idx" ]; then
        log_message "Skipping $species - files already exist"
        skipped_count=$((skipped_count + 1))
        # Still run fetch_transcriptome to ensure config is updated
        fetch_transcriptome "$species" "${SPECIES[$species]}"
    else
        log_message "Processing $species"
        fetch_transcriptome "$species" "${SPECIES[$species]}"
        processed_count=$((processed_count + 1))
    fi
done

# Create/update file indicating successful reference preparation
echo "$(date '+%Y-%m-%d %H:%M:%S')" > ../data/reference/references_prepared.txt

log_message "Transcriptome processing complete: $processed_count species processed, $skipped_count species skipped"
log_message "All transcriptomes downloaded and Kallisto indices built successfully!"
exit 0 