import pandas as pd

# Load the data
path = r"C:\Users\Nikita Abramenko\Documents\SynologyDrive\abramenko\uniprot_data.csv"
df = pd.read_csv(path, low_memory=False)
dhodh_rec = "dhodh.csv"
""" 
obtain the column names from the file and save them to a text file 
df = pd.read_csv(path, low_memory=False, nrows=1, index_col=0)
columns = df.columns.tolist()
columns_file = 'columns_list.txt'
with open(columns_file, 'w') as file:
    for column in columns:
        file.write(f"{column}\n")
print(f"Column names have been saved to {columns_file}")
"""

# Select the records for the DHODH protein
column_name = 'Protein names'
filter_value = 'Dihydroorotate dehydrogenase'
filtered_df = df.loc[df[column_name].str.contains(filter_value, na=False, case=False)]
filtered_df.to_csv(dhodh_rec, index=False)