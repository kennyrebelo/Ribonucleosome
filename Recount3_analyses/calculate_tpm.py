import sqlite3
import pandas as pd
import time
import argparse

def read_sample_ids(file_path):
    with open(file_path, 'r') as file:
        return [int(line.strip()) for line in file if line.strip()]

def create_temp_tables(conn, cursor):
    print("Creating temporary tables...")
    cursor.execute("""
    CREATE TEMPORARY TABLE gene_list (
        gene_name TEXT,
        gene_id TEXT,
        total_exonic_length FLOAT
    )
    """)
    
    cursor.execute("""
    CREATE TEMPORARY TABLE samples (
        sample_id INTEGER PRIMARY KEY,
        gene_fc_count_unique_total FLOAT
    )
    """)
    print("Temporary tables created.")

def load_data(conn, cursor, gene_name, sample_ids):
    print("Loading gene list data...")
    gene_df = pd.read_csv('table_transcript_exonic_lengths.tab', delimiter='\t')
    gene_df = gene_df[gene_df['gene_name'] == gene_name]
    gene_df.to_sql('gene_list', conn, if_exists='append', index=False)
    print(f"Loaded gene {gene_name} into temporary table.")

    print("Loading samples data...")
    samples_df = pd.read_csv('samples_higher_10000000_column_star.uniquely_mapped_reads_number.tsv', sep='\t', low_memory=False)
    samples_df = samples_df.rename(columns={'rail_id': 'sample_id', 'gene_fc_count_unique.total': 'gene_fc_count_unique_total'})
    samples_df = samples_df[samples_df['sample_id'].isin(sample_ids)]
    samples_df[['sample_id', 'gene_fc_count_unique_total']].to_sql('samples', conn, if_exists='append', index=False)
    print(f"Loaded {len(samples_df)} samples into temporary table.")

def calculate_tpm(conn, cursor):
    print("Executing main query to calculate TPM...")
    query = """
    WITH parsed_samples AS (
        SELECT 
            i.right_annotated,
            CAST(SUBSTR(s.value, 1, INSTR(s.value, ':') - 1) AS INTEGER) AS sample_id,
            CAST(SUBSTR(s.value, INSTR(s.value, ':') + 1) AS INTEGER) AS count
        FROM intron i
        CROSS JOIN json_each('["' || REPLACE(i.samples, ',', '","') || '"]') AS s
        WHERE EXISTS (
            SELECT 1 FROM gene_list gl
            WHERE i.right_annotated LIKE '%' || gl.gene_name || '%'
        )
    ),
    rpk_calc AS (
        SELECT 
            gl.gene_name,
            ps.sample_id,
            SUM(ps.count / (gl.total_exonic_length / 1000.0)) AS RPK
        FROM parsed_samples ps
        JOIN gene_list gl ON ps.right_annotated LIKE '%' || gl.gene_name || '%'
        GROUP BY gl.gene_name, ps.sample_id
    )
    SELECT 
        rc.gene_name,
        rc.sample_id,
        rc.RPK / (s.gene_fc_count_unique_total / 1e6) AS TPM
    FROM rpk_calc rc
    JOIN samples s ON rc.sample_id = s.sample_id
    """
    
    results = pd.read_sql_query(query, conn)
    print(f"Query executed. Retrieved {len(results)} rows.")
    return results

def main():
    parser = argparse.ArgumentParser(description='Calculate TPM for a specific gene and list of sample IDs.')
    parser.add_argument('gene_name', type=str, help='Name of the gene to calculate TPM for')
    parser.add_argument('sample_file', type=str, help='Path to the file containing sample IDs (one per line)')
    args = parser.parse_args()

    sample_ids = read_sample_ids(args.sample_file)

    start_time = time.time()

    print("Connecting to SQLite database...")
    conn = sqlite3.connect('genes.sqlite')
    cursor = conn.cursor()
    print("Connected to database.")

    create_temp_tables(conn, cursor)
    load_data(conn, cursor, args.gene_name, sample_ids)
    
    results = calculate_tpm(conn, cursor)
    
    print("Saving results to CSV...")
    results.to_csv(f'{args.gene_name}_tpm_results.csv', index=False)
    print(f"Saved TPM data to '{args.gene_name}_tpm_results.csv'")

    conn.close()
    print("Database connection closed.")

    end_time = time.time()
    execution_time = end_time - start_time
    print(f"\nScript execution time: {execution_time:.2f} seconds")

if __name__ == "__main__":
    main()
