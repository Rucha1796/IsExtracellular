import time
from google.colab import files
import io
import pandas as pd
import requests
from functools import lru_cache
from concurrent.futures import ThreadPoolExecutor, as_completed
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

# Load the peptide data
file_path = '/content/combined_modified_peptide.tsv'
df_peptide = pd.read_csv(file_path, sep='\t')

# Load the Wollscheid list
with open('/content/surfaceome_accession_wollschied (2).txt', 'r') as file:
    wollscheid_set = set(file.read().splitlines())

# Filter out rows without necessary data in the peptide dataframe
df_peptide = df_peptide[df_peptide['Peptide Sequence'].notna() & df_peptide['Assigned Modifications'].notna() & df_peptide['Assigned Modifications'].str.strip() != '']

# Check the number of rows after filtering
#num_rows = len(df_peptide)
#print(f"The number of rows in the dataset after filtering: {num_rows}")

# Retry configuration
retry_strategy = Retry(
    total=5,
    backoff_factor=1,
    status_forcelist=[429, 500, 502, 503, 504],
    allowed_methods=["HEAD", "GET", "OPTIONS"]
)
adapter = HTTPAdapter(max_retries=retry_strategy)
http = requests.Session()
http.mount("https://", adapter)
http.mount("http://", adapter)

@lru_cache(maxsize=5000)
def get_protein_sequence(uniprot_id):
    """ Cached version to avoid repeated network requests for the same UniProt ID. """
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = http.get(url, timeout=10)
    if response.status_code == 200:
        return ''.join(response.text.split('\n')[1:])
    return None

# function to retrieve and check extracellular domains
def get_extracellular_domains(uniprot_id):
    """ Retrieves and checks for extracellular domains. """
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.txt"
    response = http.get(url, timeout=10)
    domains = []
    found_topo_dom = False
    if response.status_code == 200:
      lines = response.text.split('\n')
      for i in range(len(lines)):
        line = lines[i]
        if 'TOPO_DOM' in line:
            found_topo_dom = True
            if 'Extracellular' in next_line(lines, i+1):
                parts = line.split()
                if len(parts) > 2:
                    range_str = parts[2]
                    if '..' in range_str:
                      start, end = map(int, range_str.split('..'))
                      domains.append((start, end))
                    else:
                        print(f"Unexpected range format in line: {line}")
                else:
                        print(f"Unexpected format, not enough parts in line: {line}")
      if not found_topo_dom:
            return None
    return domains

def next_line(lines, index):
    """ Helper function to safely return the next line if it exists. """
    if index < len(lines):
        return lines[index]
    return ""

# Functions to find positions and check for modifications
def find_position_in_protein(full_sequence, peptide_sequence, amino_acid_position):
  index = full_sequence.find(peptide_sequence)
  return index + amino_acid_position if index != -1 else -1

def is_position_extracellular(position, domains):
  return any(start <= position <= end for start, end in domains)

def parse_modifications(mod_str):
  if isinstance(mod_str, float):
    mod_str = str(mod_str)
  return [int(m.split('(')[0][:-1]) for m in mod_str.split(', ') if '(' in m and m.split('(')[0][:-1].isdigit()]

# Function to process each row with additional Wollscheid check
def process_row(index, row, results):
  start_time = time.time()
  protein_sequence = get_protein_sequence(row['Protein ID'])
  is_wollscheid = 'Yes' if row['Protein ID'] in wollscheid_set else 'No'
  if protein_sequence:
    mod_positions = parse_modifications(row['Assigned Modifications'])
    domains = get_extracellular_domains(row['Protein ID'])
    row_results = []
    for mod_pos in mod_positions:
      pos_in_protein = find_position_in_protein(protein_sequence, row['Peptide Sequence'], mod_pos)
      if pos_in_protein != -1:
        if domains is None:
           extracellular_status = 'Not Available'
        else:
           extracellular_status = 'Yes' if is_position_extracellular(pos_in_protein, domains) else 'No'
        row_results.append((pos_in_protein, extracellular_status, is_wollscheid))
    results[index] = row_results
  else:
    results[index] = [('Protein sequence not found', 'Protein sequence not found', is_wollscheid)]
  print(f"Processed row {index} in {time.time() - start_time:.2f} seconds")

# Add new columns for results
df_peptide['Position in Protein'] = None
df_peptide['Is Extracellular'] = None
df_peptide['IsWollschied'] = None

# Use ThreadPoolExecutor to process rows concurrently
start_time = time.time()
results = {}
with ThreadPoolExecutor(max_workers=20) as executor:
  futures = [executor.submit(process_row, index, row, results) for index, row in df_peptide.iterrows()]
  for future in as_completed(futures):
    future.result()  # Ensure any exceptions are raised

# Update DataFrame with results
for index, row_results in results.items():
  df_peptide.at[index, 'Position in Protein'] = ', '.join(str(res[0]) for res in row_results if res[0] is not None)
  df_peptide.at[index, 'Is Extracellular'] = ', '.join(res[1] for res in row_results)
  df_peptide.at[index, 'IsWollschied'] = ', '.join(res[2] for res in row_results)

print(f"Processed all rows in {time.time() - start_time:.2f} seconds")

# Save the updated DataFrame to a new Excel file
updated_file_path = 'updated_peptide.xlsx'
df_peptide.to_excel(updated_file_path, index=False)
files.download(updated_file_path)
