import csv

# Input and output file paths
input_csv_file = "JASSY_Q9CAL6_AFDB_FS.csv"  # Replace with your input CSV file name
output_file = "AFDB_ids.txt"      # Replace with your desired output file name

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
    csv_reader = csv.DictReader(csv_file)  # Use DictReader to access columns by name
    for row in csv_reader:
        # Skip empty rows
        if not row:  # Check if the row is empty
            continue
        # Access the 'AFDB accession' column
        alphafold_id = row.get('AFDB accession')
        if alphafold_id:  # Ensure the column is not empty
            uniprot_id = extract_uniprot_id(alphafold_id)
            if uniprot_id:
                uniprot_ids.append(uniprot_id)

# Write the UniProt IDs to the output file
with open(output_file, mode='w') as file:
    for uniprot_id in uniprot_ids:
        file.write(alphafold_id + '\n')

print(f"Extracted {len(uniprot_ids)} AlphaFold IDs and saved them to {output_file}.")
