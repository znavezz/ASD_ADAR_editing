import pandas as pd  
import numpy as np   
import os
import math
from pyfaidx import Fasta, FetchError
from pandarallel import pandarallel
import json
import requests
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional    

try:
    from tqdm.notebook import tqdm
    _HAS_TQDM = True
    print("tqdm imported successfully")
except Exception:
    _HAS_TQDM = False
    raise RuntimeError("tqdm could not be imported")

import threading
_thread_local = threading.local()

# Initialize pandarallel with a progress bar
pandarallel.initialize(progress_bar=True, nb_workers=20)

def _print(msg: str) -> None:
    print(msg, flush=True)

# Display settings
pd.set_option('display.max_columns', None)  # Show all columns in outputs




def load_vcf(vcf_path: str, first_col_name: str = "CHROM") -> pd.DataFrame:
    """
    Optimized VCF loader - extracts header then reads data efficiently.
    
    Parameters:
    ----------
    vcf_path : str
        Path to the VCF file.
    first_col_name : str
        Name of the first column in the VCF (default: 'CHROM').
    
    Returns:
    -------
    pd.DataFrame
        DataFrame containing the VCF data.
    """
    # Extract header line efficiently
    header = None
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith(f'#{first_col_name}'):
                header = line.strip().split('\t')
                break
    
    if header is None:
        raise ValueError(f"Header line starting with '#{first_col_name}' not found")
    
    # Read data efficiently, skipping all comment lines
    return pd.read_csv(
        vcf_path,
        sep='\t',
        comment='#',  # Skip all lines starting with #
        names=header,  # Use extracted header
        engine='c',
        low_memory=False
    )

def parse_info_column(df: pd.DataFrame, col_to_expand: str="INFO") -> pd.DataFrame:
    """Parse INFO column and expand into separate columns WITHOUT duplicating rows."""
    
    if col_to_expand not in df.columns:
        print(f"⚠️  Column '{col_to_expand}' not found")
        return df
    
    # CRITICAL FIX: Reset index before concatenation
    df_reset = df.reset_index(drop=True)
    
    def create_info_dict(info_str: str) -> dict:
        if pd.isna(info_str) or info_str == '.':
            return {}
        return {k: v for k, v in 
                [item.split('=', 1) for item in str(info_str).split(';') if '=' in item]}
    
    # Parse INFO column
    info_dicts = df_reset[col_to_expand].apply(create_info_dict)
    info_df = pd.DataFrame(info_dicts.tolist())
    
    # Concatenate with matching indices
    df_expanded = pd.concat([df_reset, info_df], axis=1)
    
    print(f"Added {len(info_df.columns)} columns from {col_to_expand} field")
    print(f"New columns: {list(info_df.columns)}")
    
    return df_expanded

def check_columns(df: pd.DataFrame, columns: list[str]) -> None:
    """Raise ValueError if any required column is missing."""
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {', '.join(missing)}")
    
def _get_fasta_for_thread(fasta_path: str) -> Fasta:
    fa: Optional[Fasta] = getattr(_thread_local, "fasta", None)
    if fa is None or getattr(_thread_local, "fasta_path", None) != fasta_path:
        # Open read-only; pyfaidx handles its own indexing if .fai exists
        _thread_local.fasta = Fasta(fasta_path)
        _thread_local.fasta_path = fasta_path
    return _thread_local.fasta

def _lookup_one(fasta_path: str, genome_version: str, chrom: str, pos: int, ref: str) -> str:
    fa = _get_fasta_for_thread(fasta_path)
    seq = fa[chrom][pos - 1 : pos - 1 + len(ref)].seq.upper()
    return seq


def add_genome_ref_column(
    df: pd.DataFrame,
    fasta_path: str,
    genome_version: str = "hg19",
    workers: Optional[int] = None,
    chunk_size: int = 10_000
) -> None:
    """
    Create column hg38 or hg19 in *df* by looking up reference sequence in FASTA.
    Parallelized with a thread pool and per-thread Fasta handles.
    """
    if genome_version not in ["hg19", "hg38"]:
        raise ValueError(f"Invalid genome_version: {genome_version}.\n Expected 'hg19' or 'hg38'.")
    check_columns(df, ["#CHROM", "POS", "REF"])

    # Validate FASTA first (open once in main thread)
    try:
        fa_test = Fasta(fasta_path)
    except Exception as e:
        raise RuntimeError(f"Could not open FASTA: {e}") from None

    # Auto-detect chr prefix: if FASTA keys start with 'chr', prepend it to VCF chroms
    fasta_keys = set(fa_test.keys())
    needs_chr_prefix = any(k.startswith("chr") for k in fasta_keys)

    n = len(df)
    _print(f"→ FASTA lookups: {n:,} rows | genome={genome_version} | fasta={os.path.basename(fasta_path)}")
    workers = workers or os.cpu_count() or 4
    _print(f"→ Using {workers} worker threads")

    # Prepare inputs
    chroms = df["#CHROM"].astype(str).to_numpy()
    if needs_chr_prefix:
        chroms = np.array(["chr" + c if not c.startswith("chr") else c for c in chroms])
    poss   = df["POS"].astype(int).to_numpy()
    refs   = df["REF"].astype(str).to_numpy()

    # Work in chunks to keep memory bounded and keep progress meaningful
    results = [None] * n
    num_chunks = math.ceil(n / chunk_size)
    iterator = range(num_chunks)
    if _HAS_TQDM:
        iterator = tqdm(iterator, desc="FASTA lookups (chunks)", unit="chunk")

    for ci in iterator:
        start = ci * chunk_size
        end = min(start + chunk_size, n)
        batch_len = end - start

        # inner progress for each chunk (optional)
        if _HAS_TQDM:
            pbar = tqdm(total=batch_len, leave=False, unit="row", desc=f"Chunk {ci+1}/{num_chunks}")
        else:
            pbar = None

        with ThreadPoolExecutor(max_workers=workers) as ex:
            futs = []
            for i in range(start, end):
                fut = ex.submit(_lookup_one, fasta_path, genome_version, chroms[i], poss[i], refs[i])
                futs.append((i, fut))

            for i, fut in futs:
                try:
                    seq = fut.result()
                except (KeyError, FetchError) as e:
                    raise RuntimeError(f"FASTA lookup failed for {chroms[i]}:{poss[i]}: {e}") from None
                results[i] = seq
                if pbar: pbar.update(1)

        if pbar: pbar.close()

    df[genome_version] = results  # type: ignore[list-item]




# def add_off_targets(df: pd.DataFrame, blat_path: str, genome_reference: str) -> pd.DataFrame:
#     """
#     Runs BLAT using parameters that mimic the UCSC Genome Browser web tool settings, 
#     processes the results in PSL format, and aggregates hit counts.

#     Note: The 'cds_reference' argument was removed as it wasn't used in the original code snippet.
#     """
#     df_copy = df.copy()

#     def create_fasta_file_for_all_seq(kind_of_seq, path_to_write):
#         # create list of fasta format to all sequence
#         fasta_list = [">" + str(x) + "\n" + str(y) + "\n" for x, y in df_copy[kind_of_seq].items()]
#         # open file and save all fasta seq 
#         with open(path_to_write, "w") as ofile:
#             for fasta_entry in fasta_list:
#                 ofile.write(fasta_entry)

#     def sum_amount_of_hits_for_each_query_and_assign_results(blat_results_df, original_df, col_name_min_identity):
#         # Create dataframe that sums hits per query
#         amount_of_hits_df = blat_results_df.groupby(['Query_id']).size().reset_index(name='Amount')
#         amount_of_hits_df = amount_of_hits_df.set_index('Query_id')
        
#         # Initialize column with NaN
#         original_df[col_name_min_identity] = np.nan
        
#         # Map only to existing indices
#         original_df[col_name_min_identity] = original_df.index.map(amount_of_hits_df['Amount'])

#     def read_blat_results_psl(blat_file, min_alignment_length, min_identity_precent):
#         # The PSL format has 21 columns. The first 5 lines are headers to be skipped.
#         col_names = [
#             'matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 
#             'strand', 'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 
#             'blockSizes', 'qStarts', 'tStarts'
#         ]
        
#         # Read the PSL file, skipping the standard 5-line header
#         blat_results = pd.read_table(blat_file, names=col_names, skiprows=5)
        
#         # Rename qName to Query_id for consistency with your existing functions
#         blat_results = blat_results.rename(columns={'qName': 'Query_id'})

#         # Calculate identity percentage manually for PSL format: 
#         # Identity = (matches + repMatches) / (matches + misMatches + repMatches + qNumInsert) * 100
#         # A simpler approximation often used in UCSC tools is just (matches / (matches + mismatches))
#         # We will use the common implementation for identity percentage calculation:
#         blat_results['identity_precent'] = (
#             (blat_results['matches'] / (blat_results['matches'] + blat_results['misMatches'] + blat_results['repMatches'])) * 100
#         )
        
#         # Calculate alignment length as sum of blockSizes (approximate length check)
#         # We will use qEnd - qStart as a simpler alignment length approximation for filtering
#         blat_results['alignment_length'] = blat_results['qEnd'] - blat_results['qStart']

#         # filter results that the alignment length is shorter than 21
#         blat_results = blat_results[blat_results['alignment_length'] > min_alignment_length]
#         # filter results that the minimal identity is lower than
#         blat_results = blat_results[blat_results['identity_precent'] >= min_identity_precent]
#         return blat_results[['Query_id', 'identity_precent', 'alignment_length']] # Return only necessary columns

#     # Create a temporary directory
#     tmp_dir = tempfile.mkdtemp()
#     try:
#         queries_path = os.path.join(tmp_dir, "variants_dna_context.fa")
#         create_fasta_file_for_all_seq("DNA_context", queries_path)

#         results_path = os.path.join(tmp_dir, "blat_result_variants_dna_context.psl")

#         # --- MODIFIED BLAT COMMAND WITH WEB PARAMETERS AND PSL OUTPUT ---
#         # References: genome.ucsc.edu
#         # Reference for web parameters: https://genome.ucsc.edu/FAQ/FAQblat.html#:~:text=The%20default%20setting%20for%20gfServer%20dna%20matches,repMatch%20must%20be%20specified%20when%20using%20BLAT.
#         subprocess.run([
#             blat_path, 
#             genome_reference, 
#             queries_path,
#             '-out=psl',           # Use standard PSL output format
#             '-stepSize=5',        # Web default setting
#             '-repMatch=2253',     # Web default setting calculation (1024 * 11/5)
#             '-minScore=20',       # Web default setting
#             '-minIdentity=0',     # Web default setting
#             results_path
#         ], check=True) # Use check=True to raise an error if blat fails

#         # A list of the percentage of identity we want to get in the results of the ballot
#         list_of_percent = [85, 90, 95, 100]
#         for percent in list_of_percent:
#             # We now use the psl reader function
#             blat_results_df = read_blat_results_psl(results_path, 20, percent)
#             col_results_name = 'DNA_off_target_hits_minIdentity=' + str(percent)
#             sum_amount_of_hits_for_each_query_and_assign_results(blat_results_df, df_copy, col_results_name)

#         for percent in list_of_percent:
#             # create column name
#             col_results_name = 'DNA_off_target_hits_minIdentity=' + str(percent)
#             df_copy[col_results_name] = df_copy[col_results_name].fillna(0)
#             #convert the column DNA off-target hits from float to intger
#             df_copy[col_results_name] = df_copy[col_results_name].astype('int64')
#             # Subtract 1 from the amount of hits to lower the hits that are actually against the given sequence
#             # NOTE: This assumes your query sequence exactly matches a location in your reference file.
#             df_copy[col_results_name] = df_copy[col_results_name] - 1
            
#     finally:
#         # Clean up temporary directory at the end
#         shutil.rmtree(tmp_dir, ignore_errors=True)

#     return df_copy


# def add_off_targets(df: pd.DataFrame, blat_path: str, genome_reference: str) -> pd.DataFrame:
#     """
#     Runs BLAT with different identity thresholds and counts off-target hits.
#     Uses BLAT's -minIdentity parameter for accurate filtering.
#     """
#     df_copy = df.copy()

#     def create_fasta_file(sequences_dict, output_path):
#         """Create FASTA file from sequence dictionary."""
#         with open(output_path, "w") as f:
#             for idx, seq in sequences_dict.items():
#                 f.write(f">{idx}\n{seq}\n")

#     def count_hits_from_psl(psl_file):
#         """Read PSL file and count hits per query (skipping header)."""
#         try:
#             # PSL has 5-line header, then tab-delimited data
#             df_psl = pd.read_csv(psl_file, sep='\t', skiprows=5, usecols=[9], names=['qName'])
#             return df_psl['qName'].value_counts().to_dict()
#         except pd.errors.EmptyDataError:
#             return {}

#     # Create temporary directory
#     tmp_dir = tempfile.mkdtemp()
    
#     try:
#         # Write sequences to FASTA
#         queries_path = os.path.join(tmp_dir, "queries.fa")
#         create_fasta_file(df_copy["DNA_context"], queries_path)

#         # Run BLAT for each identity threshold
#         identity_thresholds = [85, 90, 95, 100]
        
#         for min_identity in identity_thresholds:
#             print(f"→ Running BLAT with minIdentity={min_identity}%...")
            
#             results_path = os.path.join(tmp_dir, f"blat_{min_identity}.psl")
#             # --- RUN BLAT COMMAND WITH WEB PARAMETERS AND PSL OUTPUT ---
#             # Reference for web parameters: https://genome.ucsc.edu/FAQ/FAQblat.html#:~:text=The%20default%20setting%20for%20gfServer%20dna%20matches,repMatch%20must%20be%20specified%20when%20using%20BLAT.
#             # Run BLAT with identity threshold
#             subprocess.run([
#                 blat_path,
#                 genome_reference,
#                 queries_path,
#                 results_path,
#                 '-out=psl', # Use standard PSL output format
#                 '-stepSize=5', # Web default setting
#                 '-repMatch=2253', # Web default setting calculation (1024 * 11/5)
#                 '-minScore=20', # Web default setting
#                 f'-minIdentity={min_identity}'
#             ], check=True, capture_output=True)
            
#             # Count hits per sequence
#             hit_counts = count_hits_from_psl(results_path)
            
#             # Create column and assign counts
#             col_name = f'DNA_off_target_hits_minIdentity={min_identity}'
#             df_copy[col_name] = df_copy.index.map(hit_counts).fillna(0).astype('int64')
            
#             # Subtract 1 to exclude self-hit
#             df_copy[col_name] = (df_copy[col_name] - 1).clip(lower=0)
            
#             print(f"✓ Found {df_copy[col_name].sum():,} total off-targets at {min_identity}%")

#     finally:
#         shutil.rmtree(tmp_dir, ignore_errors=True)

#     print("✅ BLAT off-target analysis complete!")
#     return df_copy


# def add_off_targets(df: pd.DataFrame, blat_path: str, genome_reference: str, workers: int = 4) -> pd.DataFrame:
#     """
#     Runs BLAT with different identity thresholds in parallel and counts off-target hits.
    
#     Parameters:
#     -----------
#     workers : int or None
#         Number of parallel BLAT processes (default: None = use all available cores)
#     """
#     df_copy = df.copy()

#     def create_fasta_file(sequences_dict, output_path):
#         """Create FASTA file from sequence dictionary."""
#         with open(output_path, "w") as f:
#             for idx, seq in sequences_dict.items():
#                 f.write(f">{idx}\n{seq}\n")

#     def count_hits_from_psl(psl_file):
#         """Read PSL file and count hits per query (skipping header)."""
#         try:
#             df_psl = pd.read_csv(psl_file, sep='\t', skiprows=5, usecols=[9], names=['qName'])
#             return df_psl['qName'].value_counts().to_dict()
#         except pd.errors.EmptyDataError:
#             return {}

#     def run_blat_for_identity(args):
#         """Worker function to run BLAT for a single identity threshold."""
#         min_identity, queries_path, genome_reference, blat_path, tmp_dir = args
        
#         results_path = os.path.join(tmp_dir, f"blat_{min_identity}.psl")
        
#         print(f"→ Running BLAT with minIdentity={min_identity}%...")
        
#         subprocess.run([
#             blat_path,
#             genome_reference,
#             queries_path,
#             results_path,
#             '-out=psl',
#             '-stepSize=5',
#             '-repMatch=2253',
#             '-minScore=20',
#             f'-minIdentity={min_identity}'
#         ], check=True, capture_output=True, stderr=subprocess.DEVNULL)
        
#         # Count hits
#         hit_counts = count_hits_from_psl(results_path)
        
#         print(f"✓ Completed BLAT for minIdentity={min_identity}%")
        
#         return (min_identity, hit_counts)

#     # Create temporary directory
#     tmp_dir = tempfile.mkdtemp()
    
#     try:
#         # Write sequences to FASTA once
#         queries_path = os.path.join(tmp_dir, "queries.fa")
#         create_fasta_file(df_copy["DNA_context"], queries_path)

#         # Prepare arguments for parallel execution
#         identity_thresholds = [85, 90, 95, 100]
#         args_list = [
#             (min_id, queries_path, genome_reference, blat_path, tmp_dir)
#             for min_id in identity_thresholds
#         ]
        
#         # Run BLAT in parallel (Pool auto-limits to number of tasks if workers > tasks)
#         num_processes = workers if workers is not None else len(identity_thresholds)
#         print(f"→ Running {len(identity_thresholds)} BLAT processes using {num_processes} workers...")
#         with Pool(processes=num_processes) as pool:
#             results = pool.map(run_blat_for_identity, args_list)
        
#         # Process results and create columns
#         for min_identity, hit_counts in results:
#             col_name = f'DNA_off_target_hits_minIdentity={min_identity}'
#             df_copy[col_name] = df_copy.index.map(hit_counts).fillna(0).astype('int64')
            
#             # Subtract 1 to exclude self-hit
#             df_copy[col_name] = (df_copy[col_name] - 1).clip(lower=0)
            
#             print(f"✓ Processed {df_copy[col_name].sum():,} total off-targets at {min_identity}%")

#     finally:
#         shutil.rmtree(tmp_dir, ignore_errors=True)

#     print("✅ BLAT off-target analysis complete!")
#     return df_copy




# def _run_blat_for_identity(args):
#     """Worker function to run BLAT for a single identity threshold."""
#     min_identity, queries_path, genome_reference, blat_path, tmp_dir = args
    
#     results_path = os.path.join(tmp_dir, f"blat_{min_identity}.psl")
    
#     print(f"→ Running BLAT with minIdentity={min_identity}%...")
#     # --- RUN BLAT COMMAND WITH WEB PARAMETERS AND PSL OUTPUT ---
#     # Reference for web parameters: https://genome.ucsc.edu/FAQ/FAQblat.html#:~:text=The%20default%20setting%20for%20gfServer%20dna%20matches,repMatch%20must%20be%20specified%20when%20using%20BLAT.

#     subprocess.run([
#         blat_path,
#         genome_reference,
#         queries_path,
#         results_path,
#         '-out=psl', # Use standard PSL output format
#         '-stepSize=5', # Web default setting
#         '-repMatch=2253', # Web default setting calculation (1024 * 11/5)
#         '-minScore=20', # Web default setting
#         f'-minIdentity={min_identity}'
#     ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)  
    
#     # Count hits (reuse count_hits_from_psl logic inline)
#     try:
#         df_psl = pd.read_csv(results_path, sep='\t', skiprows=5, usecols=[9], names=['qName'])
#         hit_counts = df_psl['qName'].value_counts().to_dict()
#     except pd.errors.EmptyDataError:
#         hit_counts = {}
    
#     print(f"✓ Completed BLAT for minIdentity={min_identity}%")
    
#     return (min_identity, hit_counts)


# def add_off_targets(df: pd.DataFrame, blat_path: str, genome_reference: str, workers: int = 4) -> pd.DataFrame:
#     """Runs BLAT with different identity thresholds in parallel."""
#     df_copy = df.copy()

#     def create_fasta_file(sequences_dict, output_path):
#         with open(output_path, "w") as f:
#             for idx, seq in sequences_dict.items():
#                 f.write(f">{idx}\n{seq}\n")

#     tmp_dir = tempfile.mkdtemp()
    
#     try:
#         queries_path = os.path.join(tmp_dir, "queries.fa")
#         create_fasta_file(df_copy["DNA_context"], queries_path)

#         identity_thresholds = [85, 90, 95, 100]
#         args_list = [
#             (min_id, queries_path, genome_reference, blat_path, tmp_dir)
#             for min_id in identity_thresholds
#         ]
        
#         # Run BLAT in parallel (Pool auto-limits to number of tasks if workers > tasks)
#         num_processes = workers if workers is not None else len(identity_thresholds)
#         print(f"→ Running {len(identity_thresholds)} BLAT processes using {num_processes} workers...")
#         with Pool(processes=num_processes) as pool:
#             results = pool.map(_run_blat_for_identity, args_list)

#         for min_identity, hit_counts in results:
#             col_name = f'DNA_off_target_hits_minIdentity={min_identity}'
#             df_copy[col_name] = df_copy.index.map(hit_counts).fillna(0).astype('int64')
#             df_copy[col_name] = (df_copy[col_name] - 1).clip(lower=0)
#             print(f"✓ Processed {df_copy[col_name].sum():,} total off-targets at {min_identity}%")

#     finally:
#         shutil.rmtree(tmp_dir, ignore_errors=True)

#     print("✅ BLAT off-target analysis complete!")
#     return df_copy



def calculate_psl_percent_identity(matches, misMatches, repMatches, qNumInsert, tNumInsert, 
                                   qStart, qEnd, tStart, tEnd, sizeMul=1):
    """
    Calculate percent identity exactly as UCSC BLAT web interface does.
    Based on pslCalcMilliBad from UCSC's `https://raw.githubusercontent.com/ucscGenomeBrowser/kent/master/src/utils/pslScore/pslScore.pl` file.
    """
    qAliSize = sizeMul * (qEnd - qStart)
    tAliSize = tEnd - tStart
    aliSize = min(qAliSize, tAliSize)
    
    if aliSize <= 0:
        return 0.0
    
    sizeDif = abs(qAliSize - tAliSize)
    insertFactor = qNumInsert + tNumInsert
    
    total = sizeMul * (matches + repMatches + misMatches)
    if total == 0:
        return 0.0
    
    # UCSC's log-scaled penalty for size differences
    roundAwayFromZero = 3 * math.log(1 + sizeDif)
    roundAwayFromZero = int(roundAwayFromZero + 0.5 if roundAwayFromZero >= 0 else roundAwayFromZero - 0.5)
    
    milliBad = (1000 * (misMatches * sizeMul + insertFactor + roundAwayFromZero)) / total
    
    return 100.0 - (milliBad * 0.1)





# ═══════════════════════════════════════════════════════════════════
# GENETIC CODE CONSTANTS
# ═══════════════════════════════════════════════════════════════════
GENETIC_CODE = {
    'UUU':'F','UUC':'F','UUA':'L','UUG':'L',
    'UCU':'S','UCC':'S','UCA':'S','UCG':'S',
    'UAU':'Y','UAC':'Y','UAA':'*','UAG':'*',
    'UGU':'C','UGC':'C','UGA':'*','UGG':'W',
    'CUU':'L','CUC':'L','CUA':'L','CUG':'L',
    'CCU':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAU':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'CGU':'R','CGC':'R','CGA':'R','CGG':'R',
    'AUU':'I','AUC':'I','AUA':'I','AUG':'M',
    'ACU':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAU':'N','AAC':'N','AAA':'K','AAG':'K',
    'AGU':'S','AGC':'S','AGA':'R','AGG':'R',
    'GUU':'V','GUC':'V','GUA':'V','GUG':'V',
    'GCU':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAU':'D','GAC':'D','GAA':'E','GAG':'E',
    'GGU':'G','GGC':'G','GGA':'G','GGG':'G'
}

GENETIC_CODE_3LETTER = {
    'UUU':'Phe','UUC':'Phe','UUA':'Leu','UUG':'Leu',
    'UCU':'Ser','UCC':'Ser','UCA':'Ser','UCG':'Ser',
    'UAU':'Tyr','UAC':'Tyr','UAA':'Ter','UAG':'Ter',
    'UGU':'Cys','UGC':'Cys','UGA':'Ter','UGG':'Trp',
    'CUU':'Leu','CUC':'Leu','CUA':'Leu','CUG':'Leu',
    'CCU':'Pro','CCC':'Pro','CCA':'Pro','CCG':'Pro',
    'CAU':'His','CAC':'His','CAA':'Gln','CAG':'Gln',
    'CGU':'Arg','CGC':'Arg','CGA':'Arg','CGG':'Arg',
    'AUU':'Ile','AUC':'Ile','AUA':'Ile','AUG':'Met',
    'ACU':'Thr','ACC':'Thr','ACA':'Thr','ACG':'Thr',
    'AAU':'Asn','AAC':'Asn','AAA':'Lys','AAG':'Lys',
    'AGU':'Ser','AGC':'Ser','AGA':'Arg','AGG':'Arg',
    'GUU':'Val','GUC':'Val','GUA':'Val','GUG':'Val',
    'GCU':'Ala','GCC':'Ala','GCA':'Ala','GCG':'Ala',
    'GAU':'Asp','GAC':'Asp','GAA':'Glu','GAG':'Glu',
    'GGU':'Gly','GGC':'Gly','GGA':'Gly','GGG':'Gly'
}

STOP_CODONS = ['UAA', 'UAG', 'UGA']

# Ensembl REST API configuration
ENSEMBL_SERVER_38 = "https://rest.ensembl.org"
ENSEMBL_SERVER_37 = "https://grch37.rest.ensembl.org"

# VEP API configuration
VEP_SERVER = ENSEMBL_SERVER_37
VEP_HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}
VEP_EXT = "/vep/human/hgvs"

# ═══════════════════════════════════════════════════════════════════
# PROTEIN SEQUENCE FETCHING FUNCTIONS
# ═══════════════════════════════════════════════════════════════════



# ── Consequence term → region mapping ──
# Single source of truth shared by phase2 (variants_consequences) and
# phase3 (bystanders_consequences) so the two can't drift. The `consequences`
# table stores the region *name* (e.g. "Missense"), not the raw VEP term
# (e.g. "missense_variant").
GENOMIC_CATEGORIES = {
    "Intron": ["intron_variant"],
    "StopGained": ["stop_gained"],
    "StopLost": ["stop_lost"],
    "Missense": ["missense_variant"],
    "Frameshift": ["frameshift_variant"],
    "Synonymous": ["synonymous_variant"],
    "3UTR": ["3_prime_UTR_variant"],
    "5UTR": ["5_prime_UTR_variant"],
    "Splice": [
        "splice_donor_variant",
        "splice_acceptor_variant",
        "splice_donor_region_variant",
        "splice_donor_5th_base_variant",
        "splice_region_variant",
        "splice_polypyrimidine_tract_variant",
    ],
    "Upstream": ["upstream_gene_variant"],
    "Downstream": ["downstream_gene_variant"],
    "StartLost": ["start_lost"],
    "StopRetained": ["stop_retained_variant"],
    "Regulatory": ["regulatory_region_variant"],
}

CONSEQUENCE_TO_REGION = {
    term: region
    for region, terms in GENOMIC_CATEGORIES.items()
    for term in terms
}

ALL_CONSEQUENCES = [term for terms in GENOMIC_CATEGORIES.values() for term in terms]


def consequence_terms_to_regions(terms) -> list:
    """Map a list of raw VEP consequence terms to deduped region names.

    Unknown terms map to "Other" (matching phase2 behaviour). Order is
    preserved and duplicates removed so the result can feed a nested
    join-table insert without (bystander_id, consequence_id) conflicts.
    """
    regions = []
    for term in terms or []:
        if term is None:
            continue
        region = CONSEQUENCE_TO_REGION.get(term, "Other") if term in CONSEQUENCE_TO_REGION else "Other"
        if region not in regions:
            regions.append(region)
    return regions




def normalize_missing(
    df: pd.DataFrame,
    null_tokens: set[str],
    strip_whitespace: bool = True
) -> pd.DataFrame:
    """
    Normalize missing values so Postgres COPY loads them as NULL.

    - whitespace-only strings -> NaN
    - empty strings -> NaN
    - tokens in null_tokens -> NaN
    """
    df = df.copy()

    # Operate only on object/string columns
    obj_cols = df.select_dtypes(include="object").columns

    for col in obj_cols:
        s = df[col]

        if strip_whitespace:
            s = s.str.strip()

        # Replace empty strings with NaN
        s = s.replace("", np.nan)

        # Replace explicit null tokens (e.g. '.', '-')
        if null_tokens:
            s = s.replace(list(null_tokens), np.nan)

        df[col] = s

    return df


def gql(query: str, variables: dict) -> dict:
    hasura_url = str(os.environ.get("HASURA_URL"))
    hasura_admin_secret = str(os.environ.get("HASURA_ADMIN_SECRET"))
    headers = {"x-hasura-admin-secret": hasura_admin_secret}
    for attempt in range(5):
        r = requests.post(
            hasura_url,
            json={"query": query, "variables": variables},
            headers=headers,
            timeout=180,
        )
        if r.status_code == 200:
            payload = r.json()
            if "errors" not in payload:
                # print(f"Hasura URL: {HASURA_URL}, \n Query: {query}, \nVariables: {variables}, \n Headers: {HEADERS}, \nResponse: {payload}")
                return payload["data"]
            raise RuntimeError(payload["errors"])
        if r.status_code in (502, 503, 504):
            time.sleep(0.5 * (2 ** attempt))
            continue
        raise RuntimeError(f"HTTP {r.status_code}: {r.text[:800]}")
    raise RuntimeError("Hasura request failed after retries")

