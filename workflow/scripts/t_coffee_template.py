from Bio import SeqIO
import os
import sys

# Snakemake inputs
fasta_file = snakemake.input.fasta  # FASTA file path
pdb_folder = snakemake.input.pdb_dir  # Folder containing PDB files
output_template = snakemake.output.templates  # Output file path

# Read sequence names from FASTA file
sequence_names = [record.id for record in SeqIO.parse(fasta_file, "fasta")]

# Write all template entries to a single file
with open(output_template, "w") as out_file:
    for seq_name in sequence_names:
        pdb_file = os.path.join(pdb_folder, f"{seq_name}.pdb")  # Expected PDB file path
        if os.path.exists(pdb_file):
            out_file.write(f">{seq_name} _P_ {pdb_file}\n")
        else:
            print(f"Warning: No PDB file found for {seq_name}, skipping.")

print(f"Template file '{output_template}' created successfully!")




