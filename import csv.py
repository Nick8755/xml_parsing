import pandas as pd
df = pd.read_csv('structure links.csv')
df_new = df[['DrugBank ID', 'Name', 'SMILES']]
print(df_new['Name'])

estr_rows = df_new[df_new['Name'] == 'Estradiol']

if not estr_rows.empty:
    print(estr_rows)