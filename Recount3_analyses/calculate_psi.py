import pandas as pd
import numpy as np
from collections import defaultdict
import time
import os

def parse_samples(sample_str):
    if pd.isna(sample_str):
        return {}
    try:
        return {int(float(s.split(':')[0])): int(s.split(':')[1]) for s in str(sample_str).strip(',').split(',') if s}
    except (ValueError, IndexError):
        print(f"Warning: Unable to parse sample string: {sample_str}")
        return {}

def calculate_psi(group, min_junction_reads):
    samples_dict = defaultdict(lambda: {'inclusion_left': 0, 'inclusion_right': 0, 'exclusion': 0})
    
    for _, row in group.iterrows():
        parsed_samples = parse_samples(row['samples'])
        for sample_id, count in parsed_samples.items():
            samples_dict[sample_id][row['junction_type']] += count
    
    psi_dict = {}
    for sample_id, counts in samples_dict.items():
        inclusion = (counts['inclusion_left'] + counts['inclusion_right']) / 2
        total = inclusion + counts['exclusion']
        psi = inclusion / total if total >= min_junction_reads else np.nan # change total number of reads considered for the analysis
        psi_dict[sample_id] = psi
    
    return pd.Series(psi_dict)

def process_chunk(chunk, min_junction_reads):
    return chunk.groupby('cryptic_ID').apply(calculate_psi, min_junction_reads).reset_index()

def main(input_file, min_junction_reads, chunksize):

    output_file = 'PSI_results_'+input_file
    print(f"Processing file: {input_file}")
    print(f"File size: {os.path.getsize(input_file) / (1024 * 1024):.2f} MB")

    # Read the CSV file with 'samples' column as string
    reader = pd.read_csv(input_file, chunksize=chunksize, dtype={'samples': str}, low_memory = False)

    results = []
    for i, chunk in enumerate(reader):
        chunk_start_time = time.time()
        print(f"Processing chunk {i+1}")
        result = process_chunk(chunk, min_junction_reads)
        results.append(result)
        chunk_end_time = time.time()
        print(f"Chunk {i+1} took {chunk_end_time - chunk_start_time:.2f} seconds")
    
    print("Concatenating results...")
    final_result = pd.concat(results, ignore_index=True)

    print(f"Saving results to {output_file}")
    final_result.to_csv(output_file, index=False)
    print("PSI calculation complete.")
    print(f"Output file size: {os.path.getsize(output_file) / (1024 * 1024):.2f} MB")

if __name__ == "__main__":
    start_time = time.time()

    #input file, minimum junction reads and chunksize # Adjust this based on your available memory
    main('CEs_total_cry_minus_C_40S_FP_wSJs_wide_WithSamples.csv', 10, 10000)
    main('CEs_total_cry_minus_C_CLIP_ind_FP_wSJs_wide_WithSamples.csv', 10, 10000)

    end_time = time.time()
    print(f"Total runtime: {end_time - start_time:.2f} seconds")