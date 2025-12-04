import pandas as pd
import sqlite3
import time
import re

def time_function(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"{func.__name__} took {end_time - start_time:.2f} seconds to run.")
        return result
    return wrapper

@time_function
def load_data(file_path):
    return pd.read_csv(file_path, delimiter='\t') # KR: delimiter='\t'

@time_function
def process_data(new_exons, db_path):
    conn = sqlite3.connect(db_path)
    
    new_exons.to_sql('temp_junctions', conn, if_exists='replace', index=False)
    
    # KR: changed column names to be matched below
    query = """
    SELECT t.*, i.samples
    FROM temp_junctions t
    LEFT JOIN intron i ON 
        i.chrom = t.chrom_jx AND
        i.start = t.start_jx AND 
        i.end = t.end_jx AND 
        i.strand = t.strand
    """
    
    new_exons_with_samples = pd.read_sql_query(query, conn)
    conn.close()
    
    return new_exons_with_samples

@time_function
def main(input_sample, DB):
    # Load the data
    new_exons = load_data(input_sample)
    
    # Process the data
    new_exons_with_samples = process_data(new_exons, DB) #KR: filepath
    
    # Check the results
    print(new_exons_with_samples.shape)
    print(new_exons_with_samples.head())
    
    # Check how many junctions were found
    found_junctions = new_exons_with_samples['samples'].notna().sum()
    print(f"Junctions found in the database: {found_junctions}")

    # Replace anything from the '.' to the end of the string with '_WithSamples.csv'
    new_file_name = re.sub(r'\..*$', '_WithSamples.csv', input_sample)
    
    # Save the results
    new_exons_with_samples.to_csv(new_file_name, index=False)

if __name__ == "__main__":
    start_time = time.time()


    main('CEs_total_cry_minus_C_40S_FP_wSJs_wide.bed','junctions.sqlite')
    main('CEs_total_cry_minus_C_CLIP_ind_FP_wSJs_wide.bed','junctions.sqlite')


    end_time = time.time()
    print(f"Total runtime: {end_time - start_time:.2f} seconds")