# read AFDB foldseek output table and convert AFDB ids to uniprot identifiers in the 
# proper format for the pipeline

import csv
import sys

# Input and output file paths
input_csv_file = sys.argv[1]  # Input CSV file provided as a command-line argument
output_file = sys.argv[2] # Output file to store extracted UniProt IDs

# Function to extract UniProt ID from AlphaFold DB ID
def extract_uniprot_id(alphafold_id):
    # Split the ID by '-' and take the middle part
    parts = alphafold_id.split('-')
    if len(parts) == 3:  # Ensure the format is correct
        return parts[1]
    return None

# Read the CSV file and extract UniProt IDs
uniprot_ids = []
with open(input_csv_file, mode='r') as csv_file:
    csv_reader = csv.reader(csv_file)  # Use csv.reader to access columns by index
    header = next(csv_reader)  # Skip the header row

    for row in csv_reader:
        # Skip rows that are empty or do not have at least two columns
        if len(row) < 2:
            continue
        
        # Access the second column  where the AF ID is usually located(index 1)
        alphafold_id = row[1]
        if alphafold_id:  # Ensure it's not empty
            uniprot_id = extract_uniprot_id(alphafold_id)
            if uniprot_id:
                uniprot_ids.append(uniprot_id)

# Write the UniProt IDs to the output file
with open(output_file, mode='w') as file:
    for uniprot_id in uniprot_ids:
        file.write(uniprot_id + '\n')

print(f"Extracted {len(uniprot_ids)} UniProt IDs and saved them to {output_file}.")

