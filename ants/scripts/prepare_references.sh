#!/bin/bash

# Script to download reference genomes from NCBI and build Kallisto indices

# Set up logging
LOG_DIR="../logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/reference_prep_$(date +%Y%m%d_%H%M%S).log"
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
mkdir -p ../data/reference/genomes
mkdir -p ../data/reference/transcriptomes
mkdir -p ../data/reference/kallisto_indices

# Parse species_config.yaml to get species and accession numbers
log_message "Parsing species configuration..."
SPECIES_CONFIG="../config/species_config.yaml"

if [ ! -f "$SPECIES_CONFIG" ]; then
    log_message "ERROR: Species configuration file not found at $SPECIES_CONFIG"
    exit 1
fi

# Extract species names and NCBI accession numbers
python3 -c "
import yaml
import sys

try:
    with open('$SPECIES_CONFIG', 'r') as f:
        config = yaml.safe_load(f)
    
    for species in config.get('species', []):
        name = species.get('name', '')
        genome_acc = species.get('genome', '')
        if name and genome_acc:
            print(f\"{name}|{genome_acc}\")
except Exception as e:
    print(f'Error parsing YAML: {e}', file=sys.stderr)
    sys.exit(1)
" > ../data/reference/species_list.txt

check_status "Parsing species configuration"

# Add Pogonomyrmex species if not already in config
echo "Pogonomyrmex barbatus|GCF_000187915.1" >> ../data/reference/species_list.txt
echo "Pogonomyrmex californicus|GCF_023869295.1" >> ../data/reference/species_list.txt

# Function to download and process genome for a species
process_species() {
    local species_info="$1"
    local species_name=$(echo "$species_info" | cut -d'|' -f1)
    local accession=$(echo "$species_info" | cut -d'|' -f2)
    
    # Make species name filesystem-friendly
    local safe_name=$(echo "$species_name" | tr ' ' '_')
    
    log_message "Processing $species_name (Accession: $accession)"
    
    # Get genome and annotation files using datasets
    log_message "Downloading genome data for $species_name..."
    datasets download genome accession $accession --filename ../data/reference/genomes/${safe_name}_${accession}.zip
    check_status "Download genome for $species_name"
    
    # Extract the downloaded ZIP file
    unzip -o ../data/reference/genomes/${safe_name}_${accession}.zip -d ../data/reference/genomes/${safe_name}_${accession}
    check_status "Extract genome for $species_name"
    
    # Find the transcriptome (RNA) FASTA file
    RNA_FASTA=$(find ../data/reference/genomes/${safe_name}_${accession} -name "*_rna.fna" | head -1)
    
    if [ -z "$RNA_FASTA" ]; then
        log_message "No RNA FASTA found for $species_name. Trying to create transcriptome from genomic and GFF..."
        
        # Find genomic FASTA and GFF files
        GENOMIC_FASTA=$(find ../data/reference/genomes/${safe_name}_${accession} -name "*_genomic.fna" | head -1)
        GFF_FILE=$(find ../data/reference/genomes/${safe_name}_${accession} -name "*.gff" | head -1)
        
        if [ -n "$GENOMIC_FASTA" ] && [ -n "$GFF_FILE" ]; then
            log_message "Found genomic FASTA and GFF. Using gffread to extract transcripts..."
            
            # Use gffread to extract transcripts
            gffread -w ../data/reference/transcriptomes/${safe_name}_${accession}_transcripts.fa -g $GENOMIC_FASTA $GFF_FILE
            check_status "Extract transcripts using gffread for $species_name"
            
            RNA_FASTA="../data/reference/transcriptomes/${safe_name}_${accession}_transcripts.fa"
        else
            log_message "ERROR: Could not find necessary files to create transcriptome for $species_name"
            return 1
        fi
    else
        # Copy the RNA FASTA to transcriptomes directory
        cp "$RNA_FASTA" "../data/reference/transcriptomes/${safe_name}_${accession}_transcripts.fa"
        check_status "Copy RNA FASTA for $species_name"
        
        RNA_FASTA="../data/reference/transcriptomes/${safe_name}_${accession}_transcripts.fa"
    fi
    
    # Build Kallisto index
    log_message "Building Kallisto index for $species_name..."
    kallisto index -i ../data/reference/kallisto_indices/${safe_name}_${accession}.idx "$RNA_FASTA"
    check_status "Build Kallisto index for $species_name"
    
    # Update the species_config.yaml with correct paths
    log_message "Updating species configuration for $species_name..."
    python3 -c "
import yaml
import sys

try:
    with open('$SPECIES_CONFIG', 'r') as f:
        config = yaml.safe_load(f)
    
    # Find the species entry
    found = False
    for species in config.get('species', []):
        if species.get('name') == '$species_name':
            species['genome'] = '$accession'
            species['annotation'] = '$accession'
            species['transcriptome'] = 'data/reference/transcriptomes/${safe_name}_${accession}_transcripts.fa'
            species['kallisto_index'] = 'data/reference/kallisto_indices/${safe_name}_${accession}.idx'
            found = True
            break
    
    # If species not found, add it
    if not found:
        if 'species' not in config:
            config['species'] = []
        
        config['species'].append({
            'name': '$species_name',
            'genome': '$accession',
            'annotation': '$accession',
            'transcriptome': 'data/reference/transcriptomes/${safe_name}_${accession}_transcripts.fa',
            'kallisto_index': 'data/reference/kallisto_indices/${safe_name}_${accession}.idx'
        })
    
    # Write updated config
    with open('$SPECIES_CONFIG', 'w') as f:
        yaml.dump(config, f, default_flow_style=False)
        
except Exception as e:
    print(f'Error updating YAML: {e}', file=sys.stderr)
    sys.exit(1)
"
    check_status "Update species configuration for $species_name"
    
    log_message "Completed processing for $species_name"
}

# Make sure we have required tools
log_message "Checking for required tools..."

# Check for datasets command (NCBI datasets)
if ! command -v datasets &> /dev/null; then
    log_message "Installing NCBI datasets command-line tools..."
    # Install datasets command to local bin directory
    mkdir -p ../bin
    wget -O ../bin/datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
    chmod +x ../bin/datasets
    export PATH="$PATH:$(readlink -f ../bin)"
    check_status "Install NCBI datasets tool"
fi

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

# Check for gffread
if ! command -v gffread &> /dev/null; then
    log_message "Installing gffread..."
    mkdir -p ../bin
    # Try using pip to install gffread
    if command -v pip &> /dev/null; then
        pip install --target=../lib gffread
        # Create wrapper script for gffread
        cat > ../bin/gffread << 'EOF'
#!/bin/bash
PYTHONPATH=$(readlink -f ../lib):$PYTHONPATH python -m gffread "$@"
EOF
        chmod +x ../bin/gffread
    else
        log_message "ERROR: pip not found, cannot install gffread"
        exit 1
    fi
    check_status "Install gffread"
fi

# Process each species
while IFS= read -r species_info; do
    process_species "$species_info"
done < ../data/reference/species_list.txt

log_message "All reference genomes downloaded and Kallisto indices built successfully!"

# Create/update file indicating successful reference preparation
echo "$(date '+%Y-%m-%d %H:%M:%S')" > ../data/reference/references_prepared.txt

exit 0 