###
# NOTE TO SELF: This script takes the .csv from the paper and parses it out into three datasets,
# each with their own wildtypes and scores. The wildtype sequences are obtained by going to accession number
# at UniProt and extracting the FASTA from there. The PDBs come from the UniProt's link.
# In the for loop, each gene is iterated through to generate the data in the right format needed for the
# function "get_design_sequence".


import os
import pandas as pd
from data_loader import Dataset
from biotite.sequence import ProteinSequence, NucleotideSequence
from biotite.sequence.io import fasta
from biotite.application.blast import BlastWebApp
from biotite.database import entrez
from tempfile import gettempdir

geneNames = ['fabZ', 'lpxC', 'murA']
pdbIDs = ['6N3P', '4MQY', '1UAE']
authorName = 'dewachter2023'
datasetNames = []

# Create dataset names
for gene in geneNames:
    datasetNames.append(f'{authorName}_{gene}')

# Create Dataset objects and store them in a new list
datasets = [Dataset(name) for name in datasetNames]

urls = [
    'https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-35940-3/MediaObjects/41467_2023_35940_MOESM7_ESM.csv',
    'https://rest.uniprot.org/uniprotkb/P0A6Q6.fasta',  # fabZ
    'https://rest.uniprot.org/uniprotkb/P0A725.fasta',  # lpxC
    'https://rest.uniprot.org/uniprotkb/P0A749.fasta'  # murA
]

url_files = [
    'score.txt',
    f'{geneNames[0]}_fasta.fa',
    f'{geneNames[1]}_fasta.fa',
    f'{geneNames[2]}_fasta.fa'
]

# Download files
for i, dataset in enumerate(datasets):
    # Download the score file
    dataset.download_from_url(urls[0], url_files[0])

    # Download the corresponding FASTA file based on the gene name
    fasta_index = i + 1  # Index 1 corresponds to the first FASTA file (fabZ)
    dataset.download_from_url(urls[fasta_index], url_files[fasta_index])

# Read wild-type sequences
wt_fabZ = fasta.get_sequence(fasta.FastaFile.read(os.path.join(datasets[0].dir, url_files[1])))
wt_lpxC = fasta.get_sequence(fasta.FastaFile.read(os.path.join(datasets[1].dir, url_files[2])))
wt_murA = fasta.get_sequence(fasta.FastaFile.read(os.path.join(datasets[2].dir, url_files[3])))

for i, gene in enumerate(geneNames):
    # Load the score data
    data = pd.read_csv(os.path.join(datasets[i].dir, url_files[0]), sep=';')

    # Replace comma with dot and convert to float
    data['CompetitionCoefficient'] = data['CompetitionCoefficient'].str.replace(',', '.').astype(float)

    # Filter the data to keep only rows with the relevant GeneName
    data = data[data['GeneName'] == gene]

    # Save the filtered data back to a new score file (or overwrite)
    data.to_csv(os.path.join(datasets[i].dir, url_files[0]), sep=';', index=False)

    # Modify the mutate_sequence function to dynamically pick the correct WT sequence based on gene name
    def mutate_sequence(wt_sequence, position, mutated_aa, target_aa):
        # Convert the biotite sequence object to a string
        wt_sequence_str = str(wt_sequence)

        # Convert the 1-based position to 0-based index for Python string manipulation
        pos_index = position - 1

        if pos_index < 0 or pos_index >= len(wt_sequence_str):
            return None  # Index out of bounds

        if wt_sequence_str[pos_index] == target_aa:
            # Perform the mutation
            return wt_sequence_str[:pos_index] + mutated_aa + wt_sequence_str[pos_index + 1:]
        else:
            # If the expected amino acid does not match the wild-type sequence at the position
            return None

    # Apply the mutation function to each row for the current dataset
    data['mutated_seq'] = data.apply(lambda row: mutate_sequence(
        wt_fabZ if gene == 'fabZ' else (wt_lpxC if gene == 'lpxC' else wt_murA),
        row['TargetAaPosition'],
        row['MutatedAa'],
        row['TargetAa']
    ), axis=1)

    # Create the ID column that combines GeneName and cumcount with an underscore
    data['ID'] = data['GeneName'] + '_' + data.groupby('GeneName').cumcount().astype(str)

    # Save the modified DataFrame to a new CSV file
    data.to_csv(os.path.join(datasets[i].dir, f'filtered_{gene}_output_data.csv'), index=False)
    datasets[i].get_design_sequence(list(data['mutated_seq']), list(data['ID']), list(data['CompetitionCoefficient']))
    datasets[i].get_pdb_sequence(pdbIDs[i])