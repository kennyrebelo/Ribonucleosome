
The scripts and workflows outlined below for calculating splicing quantifications across Recount3 were adapted from work originally done by [Charlotte Capitanchik](https://github.com/CharlotteAnne)

## Workflow

1. Junction query, collect spliced reads that overlap junctions of interest

1. 1. define BED file containing the two inclusion and one exclusion junction for each exon of interest.

Example tables:
	- CEs_total_cry_minus_C_40S_FP_wSJs_wide.bed
	- CEs_total_cry_minus_C_CLIP_ind_FP_wSJs_wide.bed

1. 2. Spliced read junction database from Recount3

Download from:
	- https://snaptron.cs.jhu.edu/data/tcgav2/junctions.sqlite

Documentation:
https://snaptron.cs.jhu.edu/data.html

1. 3. Run junction query script: junction_query.py

inputs:
	- defined in 1.1. and 1.2.
	- hard-coded in the main function at the end of the script

output:
	- Each will be termed {input_bed_file}_WithSamples.csv
	- Examples: CEs_total_cry_minus_C_40S_FP_wSJs_wide_WithSamples.csv, CEs_total_cry_minus_C_CLIP_ind_FP_wSJs_wide_WithSamples.csv


2. Calculate PSI for each exon of interest: calculate_psi.py

input:
	- output file(s) from 1.3.
	- minimum threshold of total reads to quantify PSI for specific junction. Set at 10.
	- number of rows processed at a time. Useful to read the data into smaller more manageable segments when dealing with large inputs. Set at 10000.
	- hard-coded in the main function at the end of the script

output:
	- PSI_results_{input_bed_file}
	- Examples: PSI_results_CEs_total_cry_minus_C_40S_FP_wSJs_wide_WithSamples.csv, PSI_results_CEs_total_cry_minus_C_CLIP_ind_FP_wSJs_wide_WithSamples.csv


3. calculate TPM expression for samples and genes of interest

3. 1. Obtain list of rail_IDs for each sample of the current Recount3 database

Download sample table from:
	- https://snaptron.cs.jhu.edu/data/tcgav2/samples.tsv

Extract rail_ids for all samples:
	- example bash command: cut -f1 samples.tsv | sort | uniq > id_all_samples.txt
	- remove row with 'rail_id' from the id_all_samples.txt file


3. 2. Calculate TPMs: calculate_tpm.py

Download Recount3 gene database:
	- https://snaptron.cs.jhu.edu/data/tcgav2/genes.sqlite

inputs:
	- tab delimited file that contains gene names and transcript exonic lengths. Example: table_transcript_exonic_lengths.tab. (Definition of this input is currently hard-coded inside the script)
	- samples.tsv table. Example: samples_higher_10000000_column_star.uniquely_mapped_reads_number.tsv (Definition of this input is currently hard-coded inside the script)
	- https://snaptron.cs.jhu.edu/data/tcgav2/genes.sqlite (Definition of this input is currently hard-coded inside the script)
	- Name of the gene of interest. Example: HNRNPC
	- file that contains rail_ids for all samples. Example: id_all_samples.txt

Example run to calculate TPM for HNRNPC and HNRNPA1:
~~~~
> python calculate_tpm.py HNRNPC id_all_samples.txt
> python calculate_tpm.py HNRNPA1 id_all_samples.txt
~~~~

outputs:
	- HNRNPC_tpm_results.csv
	- HNRNPA1_tpm_results.csv


3. 3. Merge TPMs into a single file and Log2 transform the expression data: concat_tpm_results.py

inputs:
	- directory that contains the results from 3.2.

output:
	- tpm_results.csv



4. Plot inclusion levels of cryptic splicing events in function of hnRNP expression: cryptic_inclusion_levels_by_hnRNP_exp.ipynb

inputs:
	- output from 2. (PSI_results_CEs_total_cry_minus_C_40S_FP_wSJs_wide_WithSamples.csv, PSI_results_CEs_total_cry_minus_C_CLIP_ind_FP_wSJs_wide_WithSamples.csv)
	- https://snaptron.cs.jhu.edu/data/tcgav2/samples.tsv (downloaded in 3.1.)
	- tpm_results.csv

outputs:
	- Corr_C-40S.pdf
	- Corr_C-iCLIP.pdf
