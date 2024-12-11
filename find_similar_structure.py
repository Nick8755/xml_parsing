import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

# Load the CSV file
df = pd.read_csv('structure links.csv')
#print(df.info())

# Check if the SMILES record is present for each molecule
#print(df['SMILES'].isnull().sum())
if 'SMILES' in df.columns:
    df_with_smiles = df.dropna(subset=['SMILES'])
    df_without_smiles = df[df['SMILES'].isna()]

    # Save the two dataframes to separate CSV files
    df_with_smiles.to_csv('structures_with_smiles.csv', index=False)
    df_without_smiles.to_csv('structures_without_smiles.csv', index=False)

# Use the data with SMILES records for similarity calculations
df_new = df_with_smiles[['DrugBank ID','Name', 'SMILES']]
print(df_new.head())

# Ask the user for an external SMILES input
external_smiles = input("Enter a SMILES structure for comparison: ")

# Generate RDKit molecule object for the input SMILES
external_mol = Chem.MolFromSmiles(external_smiles)
if external_mol is None:
    print("Invalid SMILES format entered. Exiting.")
    exit()

# Use MorganGenerator for fingerprint generation
morgan_gen = GetMorganGenerator(radius=2, fpSize=2048)
external_fp = morgan_gen.GetFingerprint(external_mol)

# Calculate the similarity of each molecule in the CSV file to the external SMILES
similarity_scores = []
for _, row in df_new.iterrows():
    smiles = row['SMILES']
    name = row['Name']
    drugbank_id = row['DrugBank ID']
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = morgan_gen.GetFingerprint(mol)
            similarity = DataStructs.TanimotoSimilarity(external_fp, fp)  # Using Tanimoto similarity
            similarity_scores.append({'DrugBank ID': drugbank_id, 'Name': name, 'SMILES': smiles, 'Similarity': similarity})
        else:
            print(f"Invalid SMILES in dataset: {smiles}")
    except Exception as e:
        print(f"Error processing molecule {name} with SMILES {smiles}: {e}")

# Create a DataFrame for results and sort by similarity
similarity_df = pd.DataFrame(similarity_scores)
similarity_df = similarity_df.sort_values(by='Similarity', ascending=False)

# Save the results to a CSV file
#similarity_df.to_csv('external_similarity_results.csv', index=False)

# Display the top results
print("\nTop similar structures:")
print(similarity_df.head(40))