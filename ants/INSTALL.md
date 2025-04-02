# Installation Guide

This document provides detailed instructions for installing all dependencies required to run the ant transcriptomics analysis pipeline.

## Required Dependencies

### 1. amalgkit

The main toolkit for transcriptome analysis:

```bash
# Using pip
pip install amalgkit

# For development version
git clone https://github.com/kfuku52/amalgkit.git
cd amalgkit
pip install -e .
```

### 2. SRA Toolkit

Required for downloading SRA data:

#### Linux
```bash
# Using conda
conda install -c bioconda sra-tools

# Using apt (Ubuntu/Debian)
sudo apt-get install sra-tools
```

#### macOS
```bash
# Using homebrew
brew install sra-tools

# Using conda
conda install -c bioconda sra-tools
```

For detailed installation instructions, see: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

### 3. SeqKit

Required for sequence processing during getfastq step:

```bash
# Using conda (recommended)
conda install -c bioconda seqkit

# Using homebrew
brew install seqkit

# Using Go
go get -u github.com/shenwei356/seqkit/seqkit
```

Alternatively, download pre-compiled binaries from: https://github.com/shenwei356/seqkit/releases

For detailed installation instructions, see: https://bioinf.shenwei.me/seqkit/download/

### 4. Kallisto

Required for transcript quantification:

#### Linux
```bash
# Using conda
conda install -c bioconda kallisto

# Using apt (Ubuntu/Debian)
sudo apt-get install kallisto
```

#### macOS
```bash
# Using homebrew
brew install kallisto

# Using conda
conda install -c bioconda kallisto
```

For detailed installation instructions, see: https://pachterlab.github.io/kallisto/download

## Optional Dependencies

These dependencies are optional and disabled by default in the pipeline:

### 1. parallel-fastq-dump

For faster downloading of FASTQ files:

```bash
# Using pip
pip install parallel-fastq-dump

# Using conda
conda install -c bioconda parallel-fastq-dump
```

### 2. fastp

For FASTQ quality filtering:

```bash
# Using conda
conda install -c bioconda fastp

# Using homebrew
brew install fastp
```

For detailed installation instructions, see: https://github.com/OpenGene/fastp

## Quick Setup with Conda

The easiest way to install all dependencies is using Conda:

```bash
# Create a new environment for ant transcriptomics
conda create -n ant_transcriptomics python=3.9

# Activate the environment
conda activate ant_transcriptomics

# Install all required dependencies
conda install -c bioconda -c conda-forge amalgkit sra-tools seqkit kallisto

# Optional: Install additional tools
conda install -c bioconda -c conda-forge parallel-fastq-dump fastp
```

## Verifying Installation

To verify that all dependencies are correctly installed:

```bash
# Check amalgkit
amalgkit --version

# Check SRA Toolkit
prefetch --version

# Check SeqKit
seqkit version

# Check Kallisto
kallisto version
```

If any of these commands fail, revisit the installation instructions for the specific tool. 