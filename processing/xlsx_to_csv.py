# Convert .xlsx file to .csv

#!/usr/bin/env python3

import os
import glob
import sys
import pandas as pd

def convert_excel_to_csv(folder):
    # Get list of Excel files (.xls and .xlsx) in the folder
    excel_files = glob.glob(os.path.join(folder, '*.xls')) + glob.glob(os.path.join(folder, '*.xlsx'))
    
    for excel_file in excel_files:
        print(f"Processing {excel_file}...")
        try:
            # Load the Excel file and get its first sheet name
            xls = pd.ExcelFile(excel_file)
            first_sheet = xls.sheet_names[0]
            
            # Read only the first sheet
            df = pd.read_excel(excel_file, sheet_name=first_sheet)
            base_name = os.path.splitext(os.path.basename(excel_file))[0]
            csv_filename = f"{base_name}.csv"
            csv_path = os.path.join(folder, csv_filename)
            
            # Save the data from the first sheet as CSV
            df.to_csv(csv_path, index=False)
            print(f"Saved CSV: {csv_path}")
        except Exception as e:
            print(f"Error processing {excel_file}: {e}")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python convert_excel_to_csv.py /path/to/folder")
        sys.exit(1)
    
    folder = sys.argv[1]
    if not os.path.isdir(folder):
        print(f"Error: {folder} is not a directory.")
        sys.exit(1)
    
    convert_excel_to_csv(folder)
