# Troubleshooting Guide

This document provides solutions for common issues you might encounter when running the ant transcriptomics analysis pipeline.

## Common Issues and Solutions

### Missing Dependencies

#### Problem: Error about `parallel-fastq-dump` not found
```
FileNotFoundError: [Errno 2] No such file or directory: 'parallel-fastq-dump'
```

**Solution:** This error is expected and handled. The script automatically disables the use of `parallel-fastq-dump` with the `--pfd no` option. No action is required.

#### Problem: Error about `fastp` not found
```
FileNotFoundError: [Errno 2] No such file or directory: 'fastp'
```

**Solution:** This error is expected and handled. The script automatically disables the use of `fastp` with the `--fastp no` option. No action is required.

#### Problem: Error about `seqkit` not found
```
Exception: SeqKit not found. Please make sure SeqKit is installed properly.
```

**Solution:** SeqKit is required for the getfastq step. Install it with one of these methods:

- Using conda: `conda install -c bioconda seqkit`
- Using homebrew: `brew install seqkit`
- Download binaries from: https://github.com/shenwei356/seqkit/releases
- See full installation instructions: https://bioinf.shenwei.me/seqkit/download/

#### Problem: Error about `kallisto` not found
```
FileNotFoundError: [Errno 2] No such file or directory: 'kallisto'
```

**Solution:** Kallisto is required for the transcript quantification step. Install it following the instructions at [https://pachterlab.github.io/kallisto/download](https://pachterlab.github.io/kallisto/download).

### File Not Found Errors

#### Problem: Metadata file not found
```
ERROR: Metadata file not found at data/metadata/metadata.tsv
```

**Solution:** Ensure the metadata retrieval step completed successfully. Check the log file for any errors during the metadata retrieval process.

#### Problem: No fastq files found
```
ERROR: No fastq files found in data/fastq directory. Cannot proceed with quantification.
```

**Solution:** 
1. Verify that the SRA toolkit is installed and functioning
2. Check that your internet connection is stable
3. Verify that the selected SRA entries are accessible
4. Try manually downloading one of the SRR entries using the SRA toolkit

#### Problem: Merged counts file not found
```
ERROR: Merged counts file not found. Cannot proceed with normalization.
```

**Solution:** This indicates that the merging step failed or was not completed. Check the log file for errors during the quantification and merging steps.

### Data Analysis Issues

#### Problem: Error during normalization or correlation analysis
```
ERROR: TMM normalization
```

**Solution:**
1. Verify that the input files exist and have data
2. Check that the species config file contains the correct orthology information
3. Ensure there are enough common genes across samples to perform normalization

#### Problem: Empty results or unexpected output
```
ERROR: No results found in results directory
```

**Solution:**
1. Check each step in the log file to identify where the process failed
2. Verify that the selection criteria didn't filter out all samples
3. Ensure the species configuration is correctly set up

## Running Individual Steps

If the full pipeline fails, you can run individual steps manually:

### Metadata Retrieval Only
```bash
amalgkit metadata --search_string "(Pogonomyrmex[Organism]) AND (RNA-seq[Strategy])" --out_dir data --entrez_email "your.email@example.com"
```

### Selection Only
```bash
amalgkit select --metadata data/metadata/metadata.tsv --config_dir config --out_dir data/selected --min_nspots 1000000 --max_sample 10
```

### Download Only
```bash
amalgkit getfastq --id_list data/selected/srr_list.txt --out_dir data/fastq --pfd no --fastp no
```

## Checking Log Files

The log files in the `logs` directory contain detailed information about each step of the process. To check for specific errors:

```bash
grep "ERROR" logs/analysis_TIMESTAMP.log
```

Replace `TIMESTAMP` with the actual timestamp of your analysis run.

## Getting Help

If you continue to experience issues:

1. Check the [amalgkit GitHub repository](https://github.com/kfuku52/amalgkit) for known issues
2. Review the error messages in the log file
3. Verify all dependencies are installed correctly
4. Consider running steps manually to isolate the problem 