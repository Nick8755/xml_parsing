# create calculation of similarity between two SMILES
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

# input SMILES string of the 1. molecule
first_SMILES = input("Enter a SMILES structure for 1. molecule: ")
# input SMILES string of the 2. molecule
second_SMILES = input("Enter a SMILES structure for 2. molecule: ")

# Generate RDKit molecule object for the input SMILES
first_mol = Chem.MolFromSmiles(first_SMILES)
second_mol = Chem.MolFromSmiles(second_SMILES)

# Use MorganGenerator for fingerprint generation
morgan_gen = GetMorganGenerator(radius=2, fpSize=2048)
first_fp = morgan_gen.GetFingerprint(first_mol)
second_fp = morgan_gen.GetFingerprint(second_mol)

# Calculate the similarity of the two molecules
similarity = DataStructs.TanimotoSimilarity(first_fp, second_fp)  # Using Tanimoto similarity
print(f"Similarity between the two molecules: {similarity}")

