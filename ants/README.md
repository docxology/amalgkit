# Ant Transcriptomics Analysis Pipeline

This repository contains an analysis pipeline for transcriptomic data from Pogonomyrmex ant species using the [amalgkit](https://github.com/kfuku52/amalgkit) toolkit.

## Directory Structure

```
ants/
├── config/                # Configuration files for the analysis
│   ├── control_term.config       # Control terms for sample classification
│   ├── curate_config.yaml        # Configuration for data curation
│   ├── exclude_keyword.config    # Keywords to exclude samples
│   ├── group_attribute.config    # Sample grouping attributes
│   ├── select_config.yaml        # Configuration for sample selection
│   └── species_config.yaml       # Species configuration
├── data/                  # Data directory
│   ├── fastq/                    # Downloaded FASTQ files 
│   ├── metadata/                 # Original metadata from NCBI SRA
│   └── selected/                 # Selected metadata and analysis files
│       ├── metadata/             # Processed metadata
│       │   ├── metadata.tsv              # Complete metadata
│       │   ├── metadata_original.tsv     # Original copy of metadata
│       │   ├── pivot_qualified.tsv       # Qualified samples pivot table
│       │   └── pivot_selected.tsv        # Selected samples pivot table
│       ├── selected_metadata.tsv         # Filtered metadata with only selected samples
│       └── srr_list.txt                  # List of SRR IDs for download
├── logs/                  # Analysis logs
├── results/               # Analysis results
│   ├── correlation/               # Cross-species correlation results
│   ├── quant/                     # Transcript quantification results
│   └── sanity/                    # Sanity check results
└── run_analysis.sh        # Main analysis script
```

## Key Features

1. **Robust Directory Management**:
   - The script automatically creates all required directories
   - Each directory creation is verified and logged
   - Errors during directory creation will halt the pipeline

2. **Enhanced Error Handling**:
   - Validation checks before each major processing step
   - File existence verification between pipeline stages
   - Detailed error messages with troubleshooting hints

3. **Step-by-Step Data Processing**:
   - The `data/selected` directory contains metadata for selected samples
   - The `data/selected/selected_metadata.tsv` file contains only samples marked for analysis (`is_sampled=yes`)
   - The `data/selected/srr_list.txt` file contains just the SRR IDs needed for downloading
   - The `data/fastq` directory stores downloaded fastq files
   - The `results` directory stores all analysis outputs

4. **Pipeline Workflow**:
   - Retrieving metadata from NCBI SRA
   - Selecting SRA entries based on criteria
   - Generating a list of SRR IDs for download
   - Downloading fastq files (using standard SRA tools only, no extra dependencies)
   - Quantifying transcript abundance
   - Merging abundance tables
   - Performing cross-species TMM normalization
   - Curating data
   - Performing cross-species correlation analysis
   - Running sanity checks

## Usage

To run the full analysis pipeline:

```bash
./run_analysis.sh
```

The script will:
1. Create all required directories and verify they exist
2. Retrieve metadata from NCBI SRA for Pogonomyrmex RNA-seq data
3. Select entries based on criteria in the configuration files
4. Create a filtered list of entries to analyze (`selected_metadata.tsv`)
5. Generate a list of SRR IDs for download (`srr_list.txt`)
6. Prompt for confirmation before proceeding with downloads
7. Download and process the fastq files (using minimal dependencies)
8. Perform transcript quantification and analysis

## Requirements

### Core Dependencies
- [amalgkit](https://github.com/kfuku52/amalgkit) toolkit
- SRA Toolkit for downloading fastq files
- [SeqKit](https://bioinf.shenwei.me/seqkit/) for sequence processing (required for getfastq)
- [Kallisto](https://pachterlab.github.io/kallisto/) for transcript quantification
- AWK for text processing
- Bash shell

### Installation Instructions

1. **amalgkit**:
   ```bash
   pip install amalgkit
   ```

2. **SRA Toolkit**:
   - Follow installation instructions at: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

3. **SeqKit** (required):
   - Using conda: `conda install -c bioconda seqkit`
   - Using homebrew: `brew install seqkit`
   - Download binaries from: https://github.com/shenwei356/seqkit/releases
   - See full installation instructions: https://bioinf.shenwei.me/seqkit/download/

4. **Kallisto** (required):
   - Follow installation instructions at: https://pachterlab.github.io/kallisto/download

### Optional/Bypassed Dependencies
The script has been configured to work with minimal dependencies:

- **parallel-fastq-dump**: For parallel downloading of fastq files
  - **Status**: DISABLED in current configuration (using `--pfd no`)
  - To enable: Install `parallel-fastq-dump` and remove the `--pfd no` option
  
- **fastp**: For quality filtering of FASTQ files
  - **Status**: DISABLED in current configuration (using `--fastp no`)
  - To enable: Install `fastp` and remove the `--fastp no` option

## Configuration

You can modify the analysis parameters by editing the configuration files in the `config/` directory:

- **control_term.config**: Define control terms for sample classification
- **curate_config.yaml**: Configure data curation parameters
- **exclude_keyword.config**: Define keywords to exclude samples
- **group_attribute.config**: Configure sample grouping attributes
- **select_config.yaml**: Set sample selection criteria
- **species_config.yaml**: Configure species information for analysis

## Log Files

Each run of the analysis creates a timestamped log file in the `logs/` directory, which includes:
- Detailed information about each step of the process
- Success/failure statuses for each operation
- Directory creation confirmation
- File existence verification
- Disk usage statistics at each major step 