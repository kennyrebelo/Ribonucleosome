import os
import pandas as pd
import numpy as np

# Directory containing your '_tpm_results.csv' files
folder_path = './'

# Output file
output_file = 'tpm_results_log2.csv'

# Collect all files with '_tpm_results.csv' suffix
files = [f for f in os.listdir(folder_path) if f.endswith('_tpm_results.csv')]

# Placeholder for dataframes
dataframes = []

# Process files
for i, file in enumerate(files):
    file_path = os.path.join(folder_path, file)
    # Read the file
    df = pd.read_csv(file_path)
    if i > 0:
        # Skip the header for all files except the first
        df = df[1:]
    dataframes.append(df)

# Concatenate all dataframes
final_df = pd.concat(dataframes, ignore_index=True)

# Convert the TPM column to numeric
final_df['TPM'] = pd.to_numeric(final_df['TPM'], errors='coerce')
    
# Apply log2 transformation
final_df['TPM'] = np.log2(final_df['TPM'] + 1)

# Write to the output file with a single header
final_df.to_csv(output_file, index=False)
