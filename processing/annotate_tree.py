# Match the final pdb structures to their uniprot species annotation and format for upload in itol

import pandas as pd

# Read your mapping table
df = pd.read_csv('finalset_JASSY.csv', sep=',')

# Create iTOL annotation file
with open('itol_species_labels.txt', 'w') as f:  
    f.write('LABELS\n')
    f.write('SEPARATOR COMMA\n')
    f.write('DATA\n')
    f.write('NODE_ID,LABEL,CLASS\n')
    for idx, row in df.iterrows():
        f.write(f"{row['Entry']},{row['Organism']}\n")

