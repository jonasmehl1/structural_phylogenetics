# IMPORTED FROM FOLDTREE

import glob
import os
from Bio import PDB
from Bio.Data.PDBData import protein_letters_3to1

# Snakemake input/output
pdb_dir = snakemake.input[0]
output_fasta = snakemake.output[0]


# Initialize PDB parser
parser = PDB.PDBParser(QUIET=True)

with open(output_fasta, 'w') as f:
    pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))

    for pdb in pdb_files:
        try:
            # Extract filename without extension
            pdb_basename = os.path.basename(pdb).split('.')[0]
            structure = parser.get_structure(pdb_basename, pdb)

            # Default: Take the first model and chain A (if it exists)
            model = structure[0]  # Model 0
            if 'A' in model:
                chain = model['A']
            else:
                print(f"Skipping {pdb}: No chain A found.", flush=True)
                continue

            # Extract sequence
            sequence = ""
            for residue in chain.get_residues():
                resname = residue.get_resname().upper()
                try:
                    sequence += protein_letters_3to1[resname]
                except KeyError:
                    print(f"Skipping unknown residue: {residue.get_resname()} in {pdb_basename}", flush=True)
                    continue  # Skip unknown residues

            # Write FASTA entry
            f.write(f">{pdb_basename}\n{sequence}\n")
        
        except Exception as e:
            print(f"Error processing {pdb}: {e}", flush=True)

