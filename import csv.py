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

# Find the row for the template, e.g., Bazedoxifene
template_row = df_new[df_new['Name'] == 'Bazedoxifene']
if not template_row.empty:
    print(template_row)
    template_smiles = template_row['SMILES'].values[0]
    print(f"Template SMILES: {template_smiles}")
else:
    print('No rows found for Bazedoxifene')
    exit()

# Generate RDKit molecule objects for the template
template_mol = Chem.MolFromSmiles(template_smiles)
if template_mol is None:
    print('Could not generate template molecule')
    exit()

# Use MorganGenerator for fingerprint generation
morgan_gen = GetMorganGenerator(radius=2, fpSize=2048)
template_fp = morgan_gen.GetFingerprint(template_mol)

# Choose the similarity metric
similarity_metric = 'Tanimoto' # Options: Tanimoto, Dice, Cosine, etc.

# Calculate the similarity of each molecule to the template with 'SMILES' record
similarity_scores = []
for _, row in df_new.iterrows():
    smiles = row['SMILES']
    name = row['Name']
    drugbank_id = row['DrugBank ID']
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = morgan_gen.GetFingerprint(mol)
            if similarity_metric == "Tanimoto":
                similarity = DataStructs.TanimotoSimilarity(template_fp, fp)
            elif similarity_metric == "Dice":
                similarity = DataStructs.DiceSimilarity(template_fp, fp)
            elif similarity_metric == "Cosine":
                similarity = DataStructs.BulkCosineSimilarity(template_fp, [fp])[0]
            else:
                raise ValueError(f"Unknown similarity metric: {similarity_metric}")

            if similarity >=0.3: # Set the similarity threshold
                similarity_scores.append({'DrugBank ID': drugbank_id,'Name': name, 'Similarity': similarity}) # add 'SMILES': smiles, optionaly
        else:
            print(f"Invalid SMILES: {smiles}")
    except Exception as e:
        print(f"Error processing molecule {name} with SMILES {smiles}: {e}")

# Create a DataFrame for results and sort by similarity
similarity_df = pd.DataFrame(similarity_scores)
similarity_df = similarity_df.sort_values(by='Similarity', ascending=False)

# Display the top results
print(similarity_df.head(20))
similarity_df.to_csv('similarity_results.csv', index=False)


'''
Heat Map
# Create a pairwise similarity matrix
unique_names = similarity_df['DrugBank ID'].unique()
similarity_matrix = pd.DataFrame(0, index=unique_names, columns=unique_names)

# Populate the matrix
for _, row in similarity_df.iterrows():
    name = row['DrugBank ID']
    similarity_matrix.loc[name, name] = row['Similarity']

# Fill diagonal with similarity 1.0 (self-similarity)
np.fill_diagonal(similarity_matrix.values, 1.0)

# Plot the heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(similarity_matrix, annot=True, fmt=".2f", cmap="coolwarm", cbar=True)
plt.title("Heatmap of Similarity Scores")
plt.xlabel("DrugBank ID")
plt.ylabel("DrugBank ID")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()
'''