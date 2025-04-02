#!/usr/bin/env python3
"""
Smart SRA Downloader for Amalgkit
--------------------------------
This script improves upon the default amalgkit getfastq process by:
1. Checking for existing files and truly skipping them
2. Providing clear, concise progress reporting
3. Handling errors gracefully
4. Avoiding redundant processing
5. Supporting parallel downloads
"""

import os
import sys
import argparse
import subprocess
import time
import concurrent.futures
from queue import Queue
from threading import Lock
try:
    import pyarrow
except ImportError:
    # Install pyarrow for future pandas compatibility
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pyarrow"], 
                             stdout=subprocess.DEVNULL, 
                             stderr=subprocess.DEVNULL)
        import pyarrow
    except:
        pass  # Silently continue if pyarrow installation fails
import pandas as pd
import logging
from datetime import datetime
import shutil
import gzip
import re

# ===== DEFAULT CONFIGURATION (MODIFY THESE VALUES) =====
# Path to ID list file (one SRR ID per line)
DEFAULT_ID_LIST = "data/selected/srr_list.txt"
# Output directory for FASTQ files
DEFAULT_OUT_DIR = "data/fastq"
# Number of threads per download
DEFAULT_THREADS = 2
# Maximum concurrent downloads
DEFAULT_MAX_CONCURRENT = 6
# Metadata TSV file with SRR information
DEFAULT_METADATA = "data/selected/selected_metadata.tsv"
# Path to amalgkit executable
DEFAULT_AMALGKIT_PATH = "amalgkit"
# Whether to force redownload of existing files
DEFAULT_FORCE = False
# Terminal width for progress bars (0 for auto-detect)
DEFAULT_TERM_WIDTH = 0
# ===== END OF DEFAULT CONFIGURATION =====

# Suppress warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=DeprecationWarning)
pd.set_option('future.no_silent_downcasting', True)  # Opt-in to future behavior

# Set up logging with thread-safety
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("smart_downloader")

# Create a lock for thread-safe logging
log_lock = Lock()
# Progress display lock
progress_lock = Lock()

# Terminal width for progress bars
if DEFAULT_TERM_WIDTH <= 0:
    try:
        terminal_width = shutil.get_terminal_size().columns
    except:
        terminal_width = 80
else:
    terminal_width = DEFAULT_TERM_WIDTH

def safe_log(level, message):
    """Thread-safe logging function."""
    with log_lock:
        if level == "info":
            logger.info(message)
        elif level == "error":
            logger.error(message)
        elif level == "warning":
            logger.warning(message)
        elif level == "debug":
            logger.debug(message)

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Smart SRA Downloader for Amalgkit")
    parser.add_argument("--id_list", default=DEFAULT_ID_LIST, help=f"File containing SRR IDs, one per line (default: {DEFAULT_ID_LIST})")
    parser.add_argument("--out_dir", default=DEFAULT_OUT_DIR, help=f"Output directory for FASTQ files (default: {DEFAULT_OUT_DIR})")
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS, help=f"Number of threads per download (default: {DEFAULT_THREADS})")
    parser.add_argument("--max_concurrent", type=int, default=DEFAULT_MAX_CONCURRENT, help=f"Maximum concurrent downloads (default: {DEFAULT_MAX_CONCURRENT})")
    parser.add_argument("--metadata", default=DEFAULT_METADATA, help=f"Metadata TSV file with SRR information (default: {DEFAULT_METADATA})")
    parser.add_argument("--force", action="store_true", default=DEFAULT_FORCE, help="Force redownload of existing files")
    parser.add_argument("--amalgkit_path", default=DEFAULT_AMALGKIT_PATH, help=f"Path to amalgkit executable (default: {DEFAULT_AMALGKIT_PATH})")
    parser.add_argument("--test", action="store_true", help="Only check connectivity and SRA toolkit installation without downloading")
    return parser.parse_args()

def is_valid_fastq_file(file_path):
    """Check if a file is a valid FASTQ file and not an HTML error page."""
    try:
        # Check file size first (less than 1KB is suspicious)
        file_size = os.path.getsize(file_path)
        if file_size < 1024:
            # Files less than 1KB are probably not FASTQ files
            # Check if it's an HTML error page
            try:
                with open(file_path, 'rb') as f:
                    header = f.read(200).decode('utf-8', errors='ignore').lower()
                    if any(error_marker in header for error_marker in ['<!doctype html', '<html', 'error', 'not found', '404']):
                        safe_log("error", f"  File {file_path} is an HTML error page, not a FASTQ file")
                        return False
            except Exception as e:
                safe_log("error", f"  Error reading file header {file_path}: {e}")
                return False
                
        # If gzipped, try to check the first few bytes
        if file_path.endswith('.gz'):
            try:
                with gzip.open(file_path, 'rt', errors='ignore') as f:
                    # Try to read the first 4 lines of the FASTQ file
                    lines = []
                    for _ in range(4):
                        try:
                            line = f.readline()
                            if not line:  # End of file
                                break
                            lines.append(line)
                        except Exception as e:
                            safe_log("error", f"  Error reading line from gzipped file {file_path}: {e}")
                            return False
                    
                    # Not enough lines
                    if len(lines) < 4:
                        safe_log("error", f"  File {file_path} has fewer than 4 lines, not a valid FASTQ")
                        return False
                    
                    # Valid FASTQ starts with @, then sequence, then +, then quality
                    if not lines[0].startswith('@'):
                        safe_log("error", f"  File {file_path} does not start with @ character expected in FASTQ format")
                        return False
                    
                    if not lines[2].startswith('+'):
                        safe_log("error", f"  File {file_path} does not have the + character separating sequence and quality")
                        return False
                    
                    return True
            except (gzip.BadGzipFile, IOError, UnicodeDecodeError) as e:
                safe_log("error", f"  File {file_path} is not a valid gzip file or has encoding issues: {e}")
                return False
        # If not gzipped, check directly
        else:
            try:
                with open(file_path, 'rt', errors='ignore') as f:
                    # Try to read the first 4 lines of the FASTQ file
                    lines = []
                    for _ in range(4):
                        try:
                            line = f.readline()
                            if not line:  # End of file
                                break
                            lines.append(line)
                        except Exception as e:
                            safe_log("error", f"  Error reading line from file {file_path}: {e}")
                            return False
                    
                    # Not enough lines
                    if len(lines) < 4:
                        safe_log("error", f"  File {file_path} has fewer than 4 lines, not a valid FASTQ")
                        return False
                    
                    # Valid FASTQ starts with @, then sequence, then +, then quality
                    if not lines[0].startswith('@'):
                        safe_log("error", f"  File {file_path} does not start with @ character expected in FASTQ format")
                        return False
                    
                    if not lines[2].startswith('+'):
                        safe_log("error", f"  File {file_path} does not have the + character separating sequence and quality")
                        return False
                    
                    return True
            except (IOError, UnicodeDecodeError) as e:
                safe_log("error", f"  Could not read file {file_path} or has encoding issues: {e}")
                return False
    except Exception as e:
        safe_log("error", f"  Error validating FASTQ file {file_path}: {e}")
        return False
    
    # If we get here, the file passed basic checks
    return True

def verify_and_clean_downloads(fastq_dir, srr_id, is_paired=False):
    """Verify downloaded files are valid FASTQ and remove invalid ones."""
    # Check for all possible file patterns
    valid_files_found = False
    
    # Patterns for FASTQ files
    if is_paired:
        file_pairs = [
            (f"{srr_id}_1.fastq.gz", f"{srr_id}_2.fastq.gz"),
            (f"{srr_id}_1.fastq", f"{srr_id}_2.fastq"),
            (f"{srr_id}.1.fastq.gz", f"{srr_id}.2.fastq.gz"),
            (f"{srr_id}.1.fastq", f"{srr_id}.2.fastq"),
            (f"{srr_id}_1.fq.gz", f"{srr_id}_2.fq.gz"),
            (f"{srr_id}_1.fq", f"{srr_id}_2.fq")
        ]
        
        for file_pair in file_pairs:
            file1_path = os.path.join(fastq_dir, file_pair[0])
            file2_path = os.path.join(fastq_dir, file_pair[1])
            
            if os.path.exists(file1_path) and os.path.exists(file2_path):
                valid1 = is_valid_fastq_file(file1_path)
                valid2 = is_valid_fastq_file(file2_path)
                
                if valid1 and valid2:
                    safe_log("info", f"  Verified valid FASTQ pair: {file_pair[0]} and {file_pair[1]}")
                    valid_files_found = True
                    break
                elif not valid1:
                    safe_log("warning", f"  Invalid FASTQ file: {file_pair[0]} - Removing")
                    try:
                        os.remove(file1_path)
                    except Exception as e:
                        safe_log("error", f"  Failed to remove invalid file {file1_path}: {e}")
                elif not valid2:
                    safe_log("warning", f"  Invalid FASTQ file: {file_pair[1]} - Removing")
                    try:
                        os.remove(file2_path)
                    except Exception as e:
                        safe_log("error", f"  Failed to remove invalid file {file2_path}: {e}")
    else:
        # Single-end patterns
        single_files = [
            f"{srr_id}.fastq.gz",
            f"{srr_id}.fastq",
            f"{srr_id}.fq.gz",
            f"{srr_id}.fq"
        ]
        
        for file_name in single_files:
            file_path = os.path.join(fastq_dir, file_name)
            
            if os.path.exists(file_path):
                valid = is_valid_fastq_file(file_path)
                
                if valid:
                    safe_log("info", f"  Verified valid FASTQ file: {file_name}")
                    valid_files_found = True
                    break
                else:
                    safe_log("warning", f"  Invalid FASTQ file: {file_name} - Removing")
                    try:
                        os.remove(file_path)
                    except Exception as e:
                        safe_log("error", f"  Failed to remove invalid file {file_path}: {e}")
    
    # Remove the completion marker if no valid files were found
    if not valid_files_found:
        marker_file = os.path.join(fastq_dir, f"{srr_id}.completed")
        if os.path.exists(marker_file):
            try:
                os.remove(marker_file)
                safe_log("warning", f"  Removed invalid completion marker for {srr_id}")
            except Exception as e:
                safe_log("error", f"  Failed to remove completion marker {marker_file}: {e}")
    
    return valid_files_found

def check_existing_files(srr_id, out_dir, paired=False):
    """Check if files for an SRR ID already exist."""
    fastq_dir = os.path.join(out_dir, srr_id)
    
    # Check for completed marker file
    marker_file = os.path.join(fastq_dir, f"{srr_id}.completed")
    if os.path.exists(marker_file):
        # Verify the downloads are valid even if the marker exists
        if verify_and_clean_downloads(fastq_dir, srr_id, paired):
            return True
        else:
            # If verification failed, remove the marker
            try:
                os.remove(marker_file)
            except:
                pass
            return False
    
    # Check all possible filename patterns
    patterns = []
    
    # Common patterns produced by different tools
    if paired:
        patterns.extend([
            # SRA toolkit patterns
            (f"{srr_id}_1.fastq.gz", f"{srr_id}_2.fastq.gz"),
            (f"{srr_id}_1.fastq", f"{srr_id}_2.fastq"),
            (f"{srr_id}.1.fastq.gz", f"{srr_id}.2.fastq.gz"),
            (f"{srr_id}.1.fastq", f"{srr_id}.2.fastq"),
            # Amalgkit patterns
            (f"{srr_id}_1.fq.gz", f"{srr_id}_2.fq.gz"),
            (f"{srr_id}_1.fq", f"{srr_id}_2.fq"),
        ])
    else:
        patterns.extend([
            # Single-end patterns
            f"{srr_id}.fastq.gz",
            f"{srr_id}.fastq",
            f"{srr_id}.fq.gz",
            f"{srr_id}.fq"
        ])
    
    # Debug existing files
    if os.path.exists(fastq_dir):
        existing_files = os.listdir(fastq_dir)
        safe_log("debug", f"Existing files in {fastq_dir}: {existing_files}")
        
        # Check paired patterns
        if paired:
            for pattern_pair in patterns:
                file1_path = os.path.join(fastq_dir, pattern_pair[0])
                file2_path = os.path.join(fastq_dir, pattern_pair[1])
                if (os.path.exists(file1_path) and os.path.getsize(file1_path) > 0 and
                    os.path.exists(file2_path) and os.path.getsize(file2_path) > 0):
                    
                    # Verify these are valid FASTQ files
                    if is_valid_fastq_file(file1_path) and is_valid_fastq_file(file2_path):
                        # Files exist and are valid, create marker file
                        with open(marker_file, 'w') as f:
                            f.write(f"Downloaded and verified on {datetime.now().isoformat()}")
                        return True
                    else:
                        # Invalid files, remove them
                        if not is_valid_fastq_file(file1_path):
                            try:
                                os.remove(file1_path)
                            except:
                                pass
                        if not is_valid_fastq_file(file2_path):
                            try:
                                os.remove(file2_path)
                            except:
                                pass
        else:
            # Check single-end patterns
            for pattern in patterns:
                file_path = os.path.join(fastq_dir, pattern)
                if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
                    # Verify it's a valid FASTQ file
                    if is_valid_fastq_file(file_path):
                        # File exists and is valid, create marker file
                        with open(marker_file, 'w') as f:
                            f.write(f"Downloaded and verified on {datetime.now().isoformat()}")
                        return True
                    else:
                        # Invalid file, remove it
                        try:
                            os.remove(file_path)
                        except:
                            pass
    
    return False

def get_layout_from_metadata(metadata_file, srr_id):
    """Determine if an SRR entry is single or paired-end from metadata."""
    try:
        metadata = pd.read_csv(metadata_file, sep='\t')
        
        # Debug available columns
        if 'run_accession' not in metadata.columns:
            # Use 'run' column instead of 'run_accession'
            if 'run' in metadata.columns:
                entry = metadata[metadata['run'] == srr_id]
                if not entry.empty and 'lib_layout' in metadata.columns:
                    layout = entry['lib_layout'].iloc[0]
                    return layout.lower() == 'paired'
            return None  # Unknown
        
        entry = metadata[metadata['run_accession'] == srr_id]
        if not entry.empty and 'library_layout' in entry.columns:
            layout = entry['library_layout'].iloc[0]
            return layout.lower() == 'paired'
        return None  # Unknown
    except Exception as e:
        safe_log("warning", f"Error reading metadata for {srr_id}: {e}")
        return None

def fix_directory_permissions(directory):
    """Ensure directory has proper permissions."""
    try:
        # Make sure directory exists
        os.makedirs(directory, exist_ok=True)
        
        # Set permissions to writable
        os.chmod(directory, 0o755)  # rwxr-xr-x
        
        return True
    except Exception as e:
        safe_log("error", f"Error fixing permissions for {directory}: {e}")
        return False

def download_sra_directly(srr_id, out_dir, threads):
    """Use SRA toolkit directly to download faster."""
    try:
        # Create output directory
        fastq_dir = os.path.join(out_dir, srr_id)
        os.makedirs(fastq_dir, exist_ok=True)
        
        # Ensure directory is writable
        fix_directory_permissions(fastq_dir)
        
        # First use prefetch to get the SRA file
        safe_log("info", f"  Running prefetch for {srr_id}...")
        prefetch_cmd = ["prefetch", srr_id, "--progress", "--max-size", "50G"]
        prefetch_proc = subprocess.run(prefetch_cmd, 
                                    stdout=subprocess.PIPE, 
                                    stderr=subprocess.STDOUT,
                                    text=True, check=False)
        
        if prefetch_proc.returncode != 0:
            safe_log("error", f"  Prefetch failed: {prefetch_proc.stdout}")
            return False
            
        # Then convert to fastq using fasterq-dump
        safe_log("info", f"  Running fasterq-dump for {srr_id}...")
        fasterq_cmd = ["fasterq-dump", srr_id, 
                      "--outdir", fastq_dir,
                      "--threads", str(threads),
                      "--progress",
                      "--split-files"]
        
        fasterq_proc = subprocess.run(fasterq_cmd, 
                                    stdout=subprocess.PIPE, 
                                    stderr=subprocess.STDOUT,
                                    text=True, check=False)
        
        if fasterq_proc.returncode != 0:
            safe_log("error", f"  Fasterq-dump failed: {fasterq_proc.stdout}")
            return False
            
        # Compress the files
        safe_log("info", f"  Compressing fastq files for {srr_id}...")
        for file in os.listdir(fastq_dir):
            if file.endswith('.fastq'):
                fastq_path = os.path.join(fastq_dir, file)
                gzip_cmd = ["gzip", fastq_path]
                subprocess.run(gzip_cmd, check=False)
        
        # Verify the output files are valid FASTQ files
        is_paired = len([f for f in os.listdir(fastq_dir) if f.endswith('_1.fastq.gz') or f.endswith('_2.fastq.gz')]) > 1
        if verify_and_clean_downloads(fastq_dir, srr_id, is_paired):
            return True
        else:
            safe_log("error", f"  SRA download failed validation for {srr_id}")
            return False
        
    except Exception as e:
        safe_log("error", f"Error in direct SRA download: {e}")
        return False

def configure_proxy():
    """
    Configure proxy settings if needed.
    Some institutions require proxy settings for external connections.
    """
    print("[INFO] Checking for proxy configuration...")
    
    # Check if proxy variables are already set in environment
    http_proxy = os.environ.get('http_proxy') or os.environ.get('HTTP_PROXY')
    https_proxy = os.environ.get('https_proxy') or os.environ.get('HTTPS_PROXY')
    
    if http_proxy or https_proxy:
        print(f"[INFO] Proxy already configured: HTTP={http_proxy}, HTTPS={https_proxy}")
        return
    
    # Ask if user wants to configure a proxy
    print("[INFO] No proxy detected. If you're behind a institutional firewall, you may need to configure a proxy.")
    print("[INFO] To set proxy later, run: export http_proxy=http://proxy.example.com:8080 export https_proxy=http://proxy.example.com:8080")
    
    # Note: We don't actually prompt here since it would block the script
    # In a real interactive script, you might want to add a prompt

def try_download_with_curl(srr_id, out_dir, threads, is_paired):
    """Try to download FASTQ files directly from EBI with curl."""
    try:
        fastq_dir = os.path.join(out_dir, srr_id)
        os.makedirs(fastq_dir, exist_ok=True)
        
        # Ensure directory is writable
        fix_directory_permissions(fastq_dir)
        
        # Try different URL patterns for EBI
        success = False
        
        # Pattern 1: Common EBI SRA pattern with fastq subdirectory hierarchy
        try:
            # European Nucleotide Archive URL pattern 1
            base_url = "https://ftp.sra.ebi.ac.uk/vol1/fastq"
            
            # Get the first 6 characters of the SRR ID
            srr_prefix = srr_id[:6]
            
            # Calculate additional path components
            remaining_digits = len(srr_id) - 6
            if remaining_digits > 0:
                additional_path = "00" + srr_id[-remaining_digits:] if remaining_digits <= 2 else "0" + srr_id[-3:]
            else:
                additional_path = "00"
            
            # Construct the URL
            url_path = f"{base_url}/{srr_prefix}/{additional_path}/{srr_id}"
            
            # Download files based on layout
            if is_paired:
                file_names = [f"{srr_id}_1.fastq.gz", f"{srr_id}_2.fastq.gz"]
            else:
                file_names = [f"{srr_id}.fastq.gz"]
            
            # Try to download each file
            curl_success = True
            downloaded_files = []
            
            for file_name in file_names:
                url = f"{url_path}/{file_name}"
                output_path = os.path.join(fastq_dir, file_name)
                downloaded_files.append(output_path)
                
                safe_log("info", f"  Trying to download {file_name} from EBI (pattern 1)...")
                
                # Use curl to download the file with timeout and user agent
                curl_cmd = [
                    "curl", "-L", "-o", output_path, url, "-f", 
                    "--connect-timeout", "30", 
                    "--max-time", "300",
                    "-A", "Mozilla/5.0 Amalgkit/1.0 (SRA Download Tool; https://github.com/amalgkit; Please contact your@email.com if download issues)"
                ]
                curl_process = subprocess.run(
                    curl_cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=False
                )
                
                # Check if download was successful
                if curl_process.returncode != 0 or not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
                    safe_log("warning", f"  Curl download failed for {file_name} (pattern 1): {curl_process.returncode}")
                    curl_success = False
                    break
                
                # Verify it's a valid FASTQ file, not an HTML error page
                if not is_valid_fastq_file(output_path):
                    safe_log("warning", f"  Downloaded file is not a valid FASTQ file: {file_name}")
                    curl_success = False
                    break
            
            if curl_success:
                # Create marker file only if all files are valid
                marker_file = os.path.join(fastq_dir, f"{srr_id}.completed")
                with open(marker_file, 'w') as f:
                    f.write(f"Downloaded and verified on {datetime.now().isoformat()}")
                safe_log("info", f"  Successfully downloaded from EBI (pattern 1)")
                return True
            else:
                # Clean up any partial downloads
                for file_path in downloaded_files:
                    if os.path.exists(file_path):
                        try:
                            os.remove(file_path)
                        except:
                            pass
        except Exception as e:
            safe_log("warning", f"  Error in EBI pattern 1 download: {e}")
        
        # Pattern 2: Alternative EBI pattern (directly in fastq directory)
        try:
            # European Nucleotide Archive URL pattern 2
            base_url = "https://ftp.sra.ebi.ac.uk/vol1/fastq"
            url_path = f"{base_url}/{srr_id}"
            
            # Download files based on layout
            if is_paired:
                file_names = [f"{srr_id}_1.fastq.gz", f"{srr_id}_2.fastq.gz"]
            else:
                file_names = [f"{srr_id}.fastq.gz"]
            
            # Try to download each file
            curl_success = True
            downloaded_files = []
            
            for file_name in file_names:
                url = f"{url_path}/{file_name}"
                output_path = os.path.join(fastq_dir, file_name)
                downloaded_files.append(output_path)
                
                safe_log("info", f"  Trying to download {file_name} from EBI (pattern 2)...")
                
                # Use curl to download the file with timeout and user agent
                curl_cmd = [
                    "curl", "-L", "-o", output_path, url, "-f", 
                    "--connect-timeout", "30", 
                    "--max-time", "300",
                    "-A", "Mozilla/5.0 Amalgkit/1.0 (SRA Download Tool; https://github.com/amalgkit; Please contact your@email.com if download issues)"
                ]
                curl_process = subprocess.run(
                    curl_cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=False
                )
                
                # Check if download was successful
                if curl_process.returncode != 0 or not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
                    safe_log("warning", f"  Curl download failed for {file_name} (pattern 2): {curl_process.returncode}")
                    curl_success = False
                    break
                
                # Verify it's a valid FASTQ file, not an HTML error page
                if not is_valid_fastq_file(output_path):
                    safe_log("warning", f"  Downloaded file is not a valid FASTQ file: {file_name}")
                    curl_success = False
                    break
            
            if curl_success:
                # Create marker file only if all files are valid
                marker_file = os.path.join(fastq_dir, f"{srr_id}.completed")
                with open(marker_file, 'w') as f:
                    f.write(f"Downloaded and verified on {datetime.now().isoformat()}")
                safe_log("info", f"  Successfully downloaded from EBI (pattern 2)")
                return True
            else:
                # Clean up any partial downloads
                for file_path in downloaded_files:
                    if os.path.exists(file_path):
                        try:
                            os.remove(file_path)
                        except:
                            pass
        except Exception as e:
            safe_log("warning", f"  Error in EBI pattern 2 download: {e}")
        
        # Pattern 3: ERA pattern for European submissions
        try:
            # European Nucleotide Archive URL pattern 3 (ERA)
            base_url = "https://ftp.sra.ebi.ac.uk/vol1/ERA"
            # For ERA, we need different path construction
            if srr_id.startswith("ERR"):
                era_prefix = srr_id[:6]
                era_path = f"{base_url}/{era_prefix}/{srr_id}"
                
                # Download files based on layout
                if is_paired:
                    file_names = [f"{srr_id}_1.fastq.gz", f"{srr_id}_2.fastq.gz"]
                else:
                    file_names = [f"{srr_id}.fastq.gz"]
                
                # Try to download each file
                curl_success = True
                downloaded_files = []
                
                for file_name in file_names:
                    url = f"{era_path}/{file_name}"
                    output_path = os.path.join(fastq_dir, file_name)
                    downloaded_files.append(output_path)
                    
                    safe_log("info", f"  Trying to download {file_name} from EBI (ERA pattern)...")
                    
                    # Use curl to download the file with timeout and user agent
                    curl_cmd = [
                        "curl", "-L", "-o", output_path, url, "-f", 
                        "--connect-timeout", "30", 
                        "--max-time", "300",
                        "-A", "Mozilla/5.0 Amalgkit/1.0 (SRA Download Tool; https://github.com/amalgkit; Please contact your@email.com if download issues)"
                    ]
                    curl_process = subprocess.run(
                        curl_cmd,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                        check=False
                    )
                    
                    # Check if download was successful
                    if curl_process.returncode != 0 or not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
                        safe_log("warning", f"  Curl download failed for {file_name} (ERA pattern): {curl_process.returncode}")
                        curl_success = False
                        break
                    
                    # Verify it's a valid FASTQ file, not an HTML error page
                    if not is_valid_fastq_file(output_path):
                        safe_log("warning", f"  Downloaded file is not a valid FASTQ file: {file_name}")
                        curl_success = False
                        break
                
                if curl_success:
                    # Create marker file only if all files are valid
                    marker_file = os.path.join(fastq_dir, f"{srr_id}.completed")
                    with open(marker_file, 'w') as f:
                        f.write(f"Downloaded and verified on {datetime.now().isoformat()}")
                    safe_log("info", f"  Successfully downloaded from EBI (ERA pattern)")
                    return True
                else:
                    # Clean up any partial downloads
                    for file_path in downloaded_files:
                        if os.path.exists(file_path):
                            try:
                                os.remove(file_path)
                            except:
                                pass
        except Exception as e:
            safe_log("warning", f"  Error in EBI ERA pattern download: {e}")
        
        # If we got here, all patterns failed
        safe_log("warning", f"  All EBI download patterns failed for {srr_id}")
        return False
        
    except Exception as e:
        safe_log("warning", f"Error in curl download: {e}")
        return False

def download_with_fastq_dump(srr_id, out_dir, threads):
    """Use fastq-dump as a last resort."""
    try:
        fastq_dir = os.path.join(out_dir, srr_id)
        os.makedirs(fastq_dir, exist_ok=True)
        
        # Ensure directory is writable
        fix_directory_permissions(fastq_dir)
        
        # Download with fastq-dump
        safe_log("info", f"  Trying fastq-dump for {srr_id}...")
        
        # First, try with --split-files for paired-end data
        fastq_dump_cmd = [
            "fastq-dump", srr_id,
            "--outdir", fastq_dir,
            "--gzip",
            "--split-files"
        ]
        
        process = subprocess.run(
            fastq_dump_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False
        )
        
        # If --split-files failed, try without it for single-end data
        if process.returncode != 0:
            safe_log("info", f"  Split-files fastq-dump failed, trying single-end...")
            fastq_dump_cmd = [
                "fastq-dump", srr_id,
                "--outdir", fastq_dir,
                "--gzip"
            ]
            
            process = subprocess.run(
                fastq_dump_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=False
            )
        
        # Check if download was successful by looking for files
        if process.returncode == 0:
            # Check if any files were created
            files = [f for f in os.listdir(fastq_dir) if f.endswith('.fastq.gz') or f.endswith('.fastq')]
            if files:
                safe_log("info", f"  fastq-dump successful, created files: {files}")
                
                # Create marker file
                marker_file = os.path.join(fastq_dir, f"{srr_id}.completed")
                with open(marker_file, 'w') as f:
                    f.write(f"Downloaded and verified on {datetime.now().isoformat()}")
                
                return True
        
        safe_log("warning", f"  fastq-dump failed: {process.stdout}")
        return False
        
    except Exception as e:
        safe_log("warning", f"Error in fastq-dump: {e}")
        return False

def try_ncbi_download(srr_id, out_dir, threads, is_paired):
    """Try to download from NCBI directly."""
    try:
        fastq_dir = os.path.join(out_dir, srr_id)
        os.makedirs(fastq_dir, exist_ok=True)
        
        # Ensure directory is writable
        fix_directory_permissions(fastq_dir)
        
        # Try several NCBI URL patterns
        
        # Pattern 1: SRA direct fastq API
        try:
            base_url = "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq"
            url = f"{base_url}?acc={srr_id}"
            
            safe_log("info", f"  Trying to download {srr_id} from NCBI (pattern 1)...")
            
            # Download to a temporary file
            temp_file = os.path.join(fastq_dir, f"{srr_id}.temp")
            curl_cmd = [
                "curl", "-L", "-o", temp_file, url, "-f", 
                "--connect-timeout", "30", 
                "--max-time", "600",
                "-A", "Mozilla/5.0 Amalgkit/1.0 (SRA Download Tool; https://github.com/amalgkit; Please contact your@email.com if download issues)"
            ]
            curl_process = subprocess.run(
                curl_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=False
            )
            
            # Check if download was successful
            if curl_process.returncode == 0 and os.path.exists(temp_file) and os.path.getsize(temp_file) > 1024:
                # Process the downloaded file (might be a zip or tar archive)
                # Try to extract if it's an archive
                extract_success = False
                
                # Check if it's a zip file
                try:
                    import zipfile
                    if zipfile.is_zipfile(temp_file):
                        with zipfile.ZipFile(temp_file, 'r') as zip_ref:
                            zip_ref.extractall(fastq_dir)
                        safe_log("info", f"  Successfully extracted ZIP archive")
                        extract_success = True
                except Exception as e:
                    safe_log("warning", f"  Error extracting zip: {e}")
                
                # If not a zip file, try as tar
                if not extract_success:
                    try:
                        import tarfile
                        if tarfile.is_tarfile(temp_file):
                            with tarfile.open(temp_file) as tar:
                                tar.extractall(path=fastq_dir)
                            safe_log("info", f"  Successfully extracted TAR archive")
                            extract_success = True
                    except Exception as e:
                        safe_log("warning", f"  Error extracting tar: {e}")
                
                # If extraction failed or the file is not an archive, check if it's a direct fastq file
                if not extract_success:
                    if is_valid_fastq_file(temp_file):
                        if is_paired:
                            # For paired data, not sure what to do with a single file
                            safe_log("warning", f"  NCBI returned a single file for paired data")
                            os.remove(temp_file)
                            return False
                        else:
                            # For single-end, rename to standard format
                            fastq_file = os.path.join(fastq_dir, f"{srr_id}.fastq")
                            os.rename(temp_file, fastq_file)
                            
                            # Compress it
                            gzip_cmd = ["gzip", fastq_file]
                            subprocess.run(gzip_cmd, check=False)
                            safe_log("info", f"  Successfully processed direct FASTQ file")
                            extract_success = True
                    else:
                        # Not valid - clean up
                        os.remove(temp_file)
                
                # Clean up temporary file if it still exists
                if os.path.exists(temp_file):
                    os.remove(temp_file)
                
                # Verify files are valid FASTQ
                if extract_success and verify_and_clean_downloads(fastq_dir, srr_id, is_paired):
                    marker_file = os.path.join(fastq_dir, f"{srr_id}.completed")
                    with open(marker_file, 'w') as f:
                        f.write(f"Downloaded and verified on {datetime.now().isoformat()}")
                    safe_log("info", f"  Successfully downloaded from NCBI (pattern 1)")
                    return True
            else:
                # Clean up failed download
                if os.path.exists(temp_file):
                    os.remove(temp_file)
                safe_log("warning", f"  NCBI pattern 1 download failed: {curl_process.returncode}")
        except Exception as e:
            safe_log("warning", f"  Error in NCBI pattern 1 download: {e}")
        
        # Pattern 2: NCBI sra-pub-run direct download
        try:
            base_url = "https://sra-pub-run-odp.s3.amazonaws.com/sra"
            url = f"{base_url}/{srr_id}/{srr_id}"
            
            safe_log("info", f"  Trying to download {srr_id} from NCBI S3 (pattern 2)...")
            
            # Download to a temporary file
            temp_file = os.path.join(fastq_dir, f"{srr_id}.sra")
            curl_cmd = [
                "curl", "-L", "-o", temp_file, url, "-f", 
                "--connect-timeout", "30", 
                "--max-time", "600",
                "-A", "Mozilla/5.0 Amalgkit/1.0 (SRA Download Tool; https://github.com/amalgkit; Please contact your@email.com if download issues)"
            ]
            curl_process = subprocess.run(
                curl_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=False
            )
            
            # Check if download was successful
            if curl_process.returncode == 0 and os.path.exists(temp_file) and os.path.getsize(temp_file) > 1024:
                # Now use fastq-dump to extract the SRA file if available
                try:
                    # Check if fastq-dump is available
                    fastq_dump_cmd = ["which", "fastq-dump"]
                    subprocess.run(fastq_dump_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    
                    # Extract with fastq-dump
                    safe_log("info", f"  Converting SRA to FASTQ with fastq-dump...")
                    extract_cmd = [
                        "fastq-dump",
                        "--outdir", fastq_dir,
                        "--gzip",
                        "--split-files" if is_paired else "",
                        temp_file
                    ]
                    extract_cmd = [cmd for cmd in extract_cmd if cmd]  # Remove empty elements
                    
                    extract_process = subprocess.run(
                        extract_cmd,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                        check=False
                    )
                    
                    if extract_process.returncode == 0:
                        # Verify the extracted files
                        if verify_and_clean_downloads(fastq_dir, srr_id, is_paired):
                            marker_file = os.path.join(fastq_dir, f"{srr_id}.completed")
                            with open(marker_file, 'w') as f:
                                f.write(f"Downloaded and verified on {datetime.now().isoformat()}")
                            safe_log("info", f"  Successfully downloaded from NCBI S3 (pattern 2)")
                            # Clean up SRA file
                            if os.path.exists(temp_file):
                                os.remove(temp_file)
                            return True
                    else:
                        safe_log("warning", f"  fastq-dump extraction failed: {extract_process.returncode}")
                except Exception as e:
                    safe_log("warning", f"  fastq-dump not available or failed: {e}")
                
                # If extraction failed or fastq-dump not available, keep the SRA file
                # but report as failure since we couldn't convert to FASTQ
                if os.path.exists(temp_file):
                    safe_log("warning", f"  Downloaded SRA file but couldn't convert to FASTQ")
                    os.remove(temp_file)
            else:
                # Clean up failed download
                if os.path.exists(temp_file):
                    os.remove(temp_file)
                safe_log("warning", f"  NCBI S3 pattern 2 download failed: {curl_process.returncode}")
        except Exception as e:
            safe_log("warning", f"  Error in NCBI S3 pattern 2 download: {e}")
        
        # If we got here, all patterns failed
        safe_log("warning", f"  All NCBI download patterns failed for {srr_id}")
        return False
    
    except Exception as e:
        safe_log("warning", f"Error in NCBI download: {e}")
        temp_file = os.path.join(fastq_dir, f"{srr_id}.temp")
        if os.path.exists(temp_file):
            os.remove(temp_file)
        return False

def download_single_sra(srr_id, out_dir, threads, amalgkit_path, metadata_file):
    """Download a single SRR entry using amalgkit, with improved logging."""
    fastq_dir = os.path.join(out_dir, srr_id)
    os.makedirs(fastq_dir, exist_ok=True)
    
    # Ensure directory is writable
    fix_directory_permissions(fastq_dir)
    
    # Get workspace root directory
    workspace_root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
    bin_dir = os.path.join(workspace_root, 'bin')
    
    # Check layout from metadata
    is_paired = get_layout_from_metadata(metadata_file, srr_id)
    layout_str = 'paired' if is_paired else 'single' if is_paired is not None else 'unknown'
    safe_log("info", f"Processing {srr_id} (Layout: {layout_str})")
    
    # Check if files already exist
    if check_existing_files(srr_id, out_dir, paired=is_paired):
        safe_log("info", f"✓ {srr_id}: FASTQ files already exist, skipping download")
        return True
    
    # Try direct SRA toolkit download if prefetch and fasterq-dump are available
    try:
        # Check if prefetch is available
        subprocess.run(["which", "prefetch"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        safe_log("info", f"⬇️  Downloading {srr_id} with direct SRA toolkit...")
        start_time = time.time()
        success = download_sra_directly(srr_id, out_dir, threads)
        
        if success:
            elapsed_time = time.time() - start_time
            safe_log("info", f"✅ {srr_id}: Direct download successful in {elapsed_time:.1f} seconds")
            return True
    except Exception as e:
        safe_log("info", f"Direct SRA toolkit not available: {e}")
    
    # Try downloading with curl from EBI
    try:
        safe_log("info", f"⬇️  Trying to download {srr_id} from EBI...")
        start_time = time.time()
        success = try_download_with_curl(srr_id, out_dir, threads, is_paired)
        
        if success:
            elapsed_time = time.time() - start_time
            safe_log("info", f"✅ {srr_id}: EBI download successful in {elapsed_time:.1f} seconds")
            return True
    except Exception as e:
        safe_log("info", f"EBI download failed: {e}")
    
    # Try downloading from NCBI directly
    try:
        safe_log("info", f"⬇️  Trying to download {srr_id} from NCBI...")
        start_time = time.time()
        success = try_ncbi_download(srr_id, out_dir, threads, is_paired)
        
        if success:
            elapsed_time = time.time() - start_time
            safe_log("info", f"✅ {srr_id}: NCBI download successful in {elapsed_time:.1f} seconds")
            return True
    except Exception as e:
        safe_log("info", f"NCBI download failed: {e}")
    
    # If all else fails, try amalgkit
    safe_log("info", f"⬇️  Downloading {srr_id} with amalgkit...")
    
    # Prepare command
    temp_id_file = os.path.join(fastq_dir, f"{srr_id}.id")
    with open(temp_id_file, 'w') as f:
        f.write(f"{srr_id}\n")
    
    # Add bin directory to PATH for seqkit
    env = os.environ.copy()
    if os.path.exists(bin_dir):
        env["PATH"] = bin_dir + os.pathsep + env.get("PATH", "")
        safe_log("info", f"  Added {bin_dir} to PATH for seqkit")
    
    # Use amalgkit getfastq with different parameters
    cmd = [
        amalgkit_path, "getfastq",
        "--id_list", temp_id_file,
        "--out_dir", out_dir,
        "--threads", str(threads),
        "--metadata", metadata_file,
        "--pfd", "no",  # Try without prefetch/fasterq-dump
        "--fastp", "no"
    ]
    
    try:
        # Run command with filtered output
        start_time = time.time()
        process = subprocess.Popen(
            cmd, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            env=env  # Use the modified environment with bin in PATH
        )
        
        # Filter and display important messages only
        for line in process.stdout:
            line = line.strip()
            if any(pattern in line for pattern in [
                "Total bases:", 
                "Downloading SRA", 
                "Time elapsed", 
                "Library layout",
                "ERROR", 
                "Warning",
                "Exception"
            ]):
                safe_log("info", f"  {line}")
        
        process.wait()
        elapsed_time = time.time() - start_time
        
        # Check if download was successful
        success = process.returncode == 0
        
        if success:
            # Check both standard amalgkit output location and the specific fastq directory
            potential_directories = [
                fastq_dir,
                os.path.join(out_dir, "getfastq", srr_id),  # Some amalgkit versions put files here
                os.path.join(workspace_root, "data", "fastq", srr_id),
                os.path.join(workspace_root, "data", "getfastq", srr_id),
            ]
            
            for check_dir in potential_directories:
                if os.path.exists(check_dir):
                    safe_log("info", f"  Checking directory: {check_dir}")
                    
                    # List files in the output directory
                    files = os.listdir(check_dir)
                    safe_log("info", f"  Files in {check_dir}: {files}")
                    
                    # Look for fastq files with any extension/pattern
                    fastq_files = [f for f in files if f.endswith(('.fastq.gz', '.fastq', '.fq.gz', '.fq'))]
                    if fastq_files:
                        safe_log("info", f"  Found FASTQ files: {fastq_files}")
                        
                        # If files are in a different directory, move them to the expected location
                        if check_dir != fastq_dir:
                            os.makedirs(fastq_dir, exist_ok=True)
                            for file in fastq_files:
                                src = os.path.join(check_dir, file)
                                dst = os.path.join(fastq_dir, file)
                                safe_log("info", f"  Moving {src} to {dst}")
                                shutil.move(src, dst)
                        
                        # Create marker file
                        marker_file = os.path.join(fastq_dir, f"{srr_id}.completed")
                        with open(marker_file, 'w') as f:
                            f.write(f"Downloaded and verified on {datetime.now().isoformat()}")
                        
                        # Verify files are valid FASTQ
                        if verify_and_clean_downloads(fastq_dir, srr_id, is_paired):
                            safe_log("info", f"✅ {srr_id}: Download successful in {elapsed_time:.1f} seconds")
                            return True
            
            # If we get here, no valid files were found in any directory
            safe_log("warning", f"❌ {srr_id}: Download completed but no valid FASTQ files found")
            
            # Try fastq-dump as a last resort
            safe_log("info", f"  Amalgkit failed, trying fastq-dump as last resort...")
            if download_with_fastq_dump(srr_id, out_dir, threads):
                return True
            
            # List directory contents for debugging
            safe_log("error", f"  Directory contents of {fastq_dir}: {os.listdir(fastq_dir)}")
            return False
        else:
            safe_log("error", f"❌ {srr_id}: Download failed with return code {process.returncode}")
            
            # If amalgkit failed, try fastq-dump as a last resort
            safe_log("info", f"  Amalgkit failed, trying fastq-dump as last resort...")
            if download_with_fastq_dump(srr_id, out_dir, threads):
                return True
            
            return False
            
    except Exception as e:
        safe_log("error", f"❌ {srr_id}: Error during download: {e}")
        return False
    finally:
        # Clean up temp file
        if os.path.exists(temp_id_file):
            try:
                os.remove(temp_id_file)
            except:
                pass

def download_worker(srr_id, out_dir, threads, amalgkit_path, metadata_file, force=False):
    """Worker function for threaded downloads."""
    # Check if files already exist (without force flag)
    if not force and check_existing_files(srr_id, out_dir):
        safe_log("info", f"✓ {srr_id}: FASTQ files already exist, skipping")
        return (srr_id, "skipped")
    
    # Download the SRR entry
    success = download_single_sra(
        srr_id, 
        out_dir, 
        threads, 
        amalgkit_path,
        metadata_file
    )
    
    if success:
        return (srr_id, "success")
    else:
        return (srr_id, "failed")

def display_emoji_progress_bar(completed, total, successful, skipped, failed, elapsed_time):
    """Display a beautiful emoji progress bar."""
    with progress_lock:
        # Calculate percentage and remaining time
        progress_pct = (completed / total) * 100 if total > 0 else 0
        estimated_total = (elapsed_time / completed) * total if completed > 0 else 0
        remaining_time = max(0, estimated_total - elapsed_time) if estimated_total > 0 else 0
        
        # Determine bar width based on terminal size (leave room for text)
        text_space = 40  # Space for text indicators
        bar_width = max(10, terminal_width - text_space)
        
        # Create the progress bar
        filled_length = int(bar_width * completed // total)
        empty_length = bar_width - filled_length
        
        # Choose progress bar characters based on success/fail ratio
        if completed == 0:
            bar_char = "🔷"
        elif failed > (completed * 0.5):
            bar_char = "🔴"  # Mostly failures
        elif failed > 0:
            bar_char = "🟠"  # Some failures
        else:
            bar_char = "🟢"  # All successful/skipped
            
        # Build the bar
        bar = bar_char * filled_length + "⬜" * empty_length
        
        # Status icons
        success_icon = "✅"
        skipped_icon = "⏭️"
        failed_icon = "❌"
        time_icon = "⏱️"
        
        # Print the complete progress bar
        print(f"\r{bar} {progress_pct:5.1f}% | " +
              f"{success_icon}{successful} {skipped_icon}{skipped} {failed_icon}{failed} | " +
              f"{time_icon} {remaining_time/60:.1f}m left", end="")
        
        # Add a newline if complete
        if completed >= total:
            print()

def ensure_sra_toolkit():
    """Checks if SRA toolkit is installed, and tries to install it if not."""
    # Define local logging function since safe_log might not be available yet
    def local_log(level, message):
        """Local logging function when safe_log is not yet defined."""
        if level == "info":
            print(f"[INFO] {message}")
        elif level == "warning":
            print(f"[WARNING] {message}")
        elif level == "error":
            print(f"[ERROR] {message}")
    
    try:
        # Check if prefetch is available
        subprocess.run(["which", "prefetch"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        local_log("info", "SRA toolkit is already installed")
        return True
    except subprocess.CalledProcessError:
        local_log("warning", "SRA toolkit not found, attempting to install...")
        
        # Try using conda if available
        try:
            subprocess.run(["which", "conda"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            local_log("info", "Conda found, installing SRA toolkit with conda...")
            
            install_cmd = ["conda", "install", "-c", "bioconda", "-y", "sra-tools"]
            process = subprocess.run(
                install_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=False
            )
            
            if process.returncode == 0:
                local_log("info", "Successfully installed SRA toolkit with conda")
                return True
            else:
                local_log("warning", f"Failed to install SRA toolkit with conda: {process.stdout}")
        except subprocess.CalledProcessError:
            local_log("warning", "Conda not available")
        
        # Try apt-get if on a Debian/Ubuntu system
        try:
            subprocess.run(["which", "apt-get"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            local_log("info", "apt-get found, trying to install SRA toolkit...")
            
            local_log("info", "This may require sudo privileges. Please enter your password if prompted.")
            install_cmd = ["sudo", "apt-get", "update"]
            subprocess.run(install_cmd, check=False)
            
            install_cmd = ["sudo", "apt-get", "install", "-y", "sra-toolkit"]
            process = subprocess.run(
                install_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=False
            )
            
            if process.returncode == 0:
                local_log("info", "Successfully installed SRA toolkit with apt-get")
                return True
            else:
                local_log("warning", f"Failed to install SRA toolkit with apt-get: {process.stdout}")
        except subprocess.CalledProcessError:
            local_log("warning", "apt-get not available")
        
        local_log("error", "Could not install SRA toolkit automatically. Please install it manually.")
        local_log("info", "Instructions: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit")
        return False

def check_network_connectivity():
    """
    Check network connectivity to EBI and NCBI servers.
    Returns a dictionary with connectivity status for each server.
    """
    print("[INFO] Checking network connectivity to sequence repositories...")
    
    # List of URLs to check for connectivity
    urls_to_check = {
        "EBI_1": "https://ftp.sra.ebi.ac.uk/vol1/fastq",
        "EBI_2": "https://era-fasp.sra.ebi.ac.uk",
        "NCBI_1": "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq",
        "NCBI_2": "https://sra-pub-run-odp.s3.amazonaws.com/sra"
    }
    
    results = {}
    
    for name, url in urls_to_check.items():
        try:
            # Use curl with a 10-second timeout, only getting headers
            cmd = ["curl", "-I", "-s", "--connect-timeout", "10", url]
            process = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False
            )
            
            # Check for successful response (including redirects and HTTP/2)
            success = False
            if process.returncode == 0:
                status_line = process.stdout.splitlines()[0] if process.stdout and process.stdout.splitlines() else ""
                
                # Check for HTTP 1.1 statuses
                if any(status in process.stdout for status in [
                    "200 OK", "301 Moved", "302 Found", "303 See Other", 
                    "307 Temporary Redirect", "308 Permanent Redirect"
                ]):
                    success = True
                
                # Check for HTTP/2 success
                elif "HTTP/2" in status_line and "200" in status_line:
                    success = True
            
            if success:
                print(f"[INFO] ✅ Connection to {name} ({url}) successful")
                results[name] = True
            else:
                print(f"[WARNING] ⚠️ Connection to {name} ({url}) failed or returned unexpected status")
                print(f"[DEBUG] Response status: {process.stdout.splitlines()[0] if process.stdout and process.stdout.splitlines() else 'No response'}")
                results[name] = False
                
        except Exception as e:
            print(f"[ERROR] ❌ Error checking connection to {name} ({url}): {e}")
            results[name] = False
    
    # Check if any server is reachable
    if not any(results.values()):
        print("[ERROR] ❌ No sequence repositories are reachable. Please check your internet connection.")
        print("[INFO] Possible solutions:")
        print("  1. Check if you need to use a proxy server")
        print("  2. Verify firewall settings aren't blocking bioinformatics domains")
        print("  3. Try connecting to a different network")
        print("  4. Try using a VPN to bypass potential geographical restrictions")
    else:
        print("[INFO] At least one repository is reachable. Download attempts will be made.")
    
    return results

def main():
    """Main entry point for the script."""
    args = parse_args()
    
    # Try to install SRA toolkit if not found
    sra_installed = ensure_sra_toolkit()
    
    # Configure proxy if needed
    configure_proxy()
    
    # Check network connectivity to download servers
    connectivity = check_network_connectivity()
    
    # If just testing, exit here
    if args.test:
        print("\n[INFO] Test mode - skipping downloads")
        print(f"[INFO] SRA toolkit installed: {sra_installed}")
        print(f"[INFO] Network connectivity test results:")
        for server, status in connectivity.items():
            print(f"  - {server}: {'✅ Connected' if status else '❌ Failed'}")
        print("\n[INFO] Test completed. Fix any issues before running the actual download.")
        return
    
    # Create missing paths if appropriate
    workspace_root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
    
    # Convert relative paths to absolute if they're relative to workspace
    if not os.path.isabs(args.id_list):
        args.id_list = os.path.join(workspace_root, args.id_list)
    
    if not os.path.isabs(args.out_dir):
        args.out_dir = os.path.join(workspace_root, args.out_dir)
    
    if not os.path.isabs(args.metadata):
        args.metadata = os.path.join(workspace_root, args.metadata)
    
    # Validate inputs
    if not os.path.exists(args.id_list):
        safe_log("error", f"ID list file does not exist: {args.id_list}")
        sys.exit(1)
    
    if not os.path.exists(args.metadata):
        safe_log("error", f"Metadata file does not exist: {args.metadata}")
        sys.exit(1)
    
    # Ensure output directory exists and is writable
    os.makedirs(args.out_dir, exist_ok=True)
    fix_directory_permissions(args.out_dir)
    
    # Read SRR IDs
    with open(args.id_list, 'r') as f:
        srr_ids = [line.strip() for line in f if line.strip()]
    
    total_ids = len(srr_ids)
    
    # Display startup message
    print("\n" + "=" * terminal_width)
    print(f"🧬 Smart SRA Downloader for Amalgkit")
    print(f"📊 Starting download of {total_ids} SRA entries with up to {args.max_concurrent} concurrent downloads")
    print("=" * terminal_width)
    
    safe_log("info", f"Starting smart download of {total_ids} SRA entries with up to {args.max_concurrent} concurrent downloads")
    
    # Debug metadata column mapping
    try:
        metadata = pd.read_csv(args.metadata, sep='\t')
        metadata_columns = list(metadata.columns)
        # Only print the first 5 columns and "..." if there are more
        columns_to_show = metadata_columns[:5]
        columns_str = ", ".join(columns_to_show)
        if len(metadata_columns) > 5:
            columns_str += ", ..."
        safe_log("info", f"Metadata columns available: {columns_str}")
        if 'run' in metadata.columns:
            safe_log("info", f"Using 'run' column instead of 'run_accession'")
    except Exception as e:
        safe_log("warning", f"Could not read metadata for debugging: {e}")
    
    # Start timing
    start_time_total = time.time()
    
    # Results counters
    successful = 0
    skipped = 0
    failed = 0
    completed = 0
    
    # Use ThreadPoolExecutor for parallel downloads
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_concurrent) as executor:
        # Submit all tasks
        future_to_srr = {
            executor.submit(
                download_worker, 
                srr_id, 
                args.out_dir, 
                args.threads, 
                args.amalgkit_path, 
                args.metadata,
                args.force
            ): srr_id for srr_id in srr_ids
        }
        
        # Initial empty progress bar
        display_emoji_progress_bar(0, total_ids, 0, 0, 0, 0)
        
        # Process results as they complete
        for future in concurrent.futures.as_completed(future_to_srr):
            srr_id = future_to_srr[future]
            try:
                srr_id, result = future.result()
                completed += 1
                
                if result == "success":
                    successful += 1
                elif result == "skipped":
                    skipped += 1
                else:
                    failed += 1
                
                # Show progress
                elapsed_time = time.time() - start_time_total
                display_emoji_progress_bar(completed, total_ids, successful, skipped, failed, elapsed_time)
                
                # Also log this in the log file
                progress_pct = (completed / total_ids) * 100
                estimated_total = (elapsed_time / completed) * total_ids if completed > 0 else 0
                remaining_time = estimated_total - elapsed_time if estimated_total > 0 else 0
                
                safe_log("info", f"Progress: {progress_pct:.1f}% ({completed}/{total_ids}) - "
                        f"Est. remaining: {remaining_time/60:.1f} min - "
                        f"Success: {successful}, Skipped: {skipped}, Failed: {failed}")
                
            except Exception as e:
                safe_log("error", f"Error processing {srr_id}: {e}")
                failed += 1
                completed += 1
                
                # Update progress bar after error
                elapsed_time = time.time() - start_time_total
                display_emoji_progress_bar(completed, total_ids, successful, skipped, failed, elapsed_time)
    
    # Final summary
    total_time = time.time() - start_time_total
    
    # Print beautiful summary
    print("\n" + "=" * terminal_width)
    print(f"📊 Download Summary:")
    print(f"🔢 Total SRA entries: {total_ids}")
    print(f"✅ Successfully downloaded: {successful}")
    print(f"⏭️  Skipped (already exists): {skipped}")
    print(f"❌ Failed: {failed}")
    print(f"⏱️  Total time: {total_time/60:.1f} minutes")
    print("=" * terminal_width + "\n")
    
    safe_log("info", "=" * 60)
    safe_log("info", f"Download Summary:")
    safe_log("info", f"  Total SRA entries: {total_ids}")
    safe_log("info", f"  Successfully downloaded: {successful}")
    safe_log("info", f"  Skipped (already exists): {skipped}")
    safe_log("info", f"  Failed: {failed}")
    safe_log("info", f"  Total time: {total_time/60:.1f} minutes")
    safe_log("info", "=" * 60)
    
    # Create a marker file to indicate successful completion
    summary_file = os.path.join(args.out_dir, "download_summary.txt")
    with open(summary_file, 'w') as f:
        f.write(f"Download completed on {datetime.now().isoformat()}\n")
        f.write(f"Total SRA entries: {total_ids}\n")
        f.write(f"Successfully downloaded: {successful}\n")
        f.write(f"Skipped (already exists): {skipped}\n")
        f.write(f"Failed: {failed}\n")
        f.write(f"Total time: {total_time/60:.1f} minutes\n")
    
    # Return appropriate exit code
    if failed > 0:
        sys.exit(1)
    sys.exit(0)

if __name__ == "__main__":
    main() 