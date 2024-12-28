# IsExtracellular
The IsExtracellular Python script reads a TSV file of peptide data, maps post-translational modifications to their positions in protein sequences, and determines their subcellular location (extracellular or not) using UniProt topology data. It also flags proteins in the Wollscheid surfaceome list and saves annotated results as an Excel file.

# Post Translational Modification Data Analysis Script version 1

## Overview
This script processes peptide data from an Excel file, checks extracellular domain data via
UniProt, and compares protein IDs against a provided Wollscheid list.

## Requirements
- Python 3
- Pandas, Requests libraries
- Access to UniProt API
  
## Files
- ‘peptide.xlsx’: Source Excel file containing peptide data.
- ‘surfaceome_accession_wollschied.txt’: Text file containing the list of Wollscheid protein IDs
(make sure these are Uniprot Accession IDs)

## Usage
1. Open the peptide.tsv file with Excel and save it as peptide.xlsx.
2. Place ‘peptide.xlsx’ and ‘surfaceome_accession_wollschied.txt’ in your working
directory.
3. Update the file path to your working directory. The first file path is the path of the
source Excel file containing the peptide and the modification data. The second file path
is the path where the woolshed list text file is placed.
4. Run the script in a Python environment that supports concurrent processing. If using
Google Colab, use the .ipynb extension python file).
5. The script outputs ‘updated_peptide.xlsx’ with new columns: ‘Position in Protein,’ ‘Is
Extracellular,’ and ‘IsWollschied.’

## Output
- ‘updated_peptide.xlsx’: Excel file with additional data on peptide positions, extracellular
status, and Wollscheid list comparison.
