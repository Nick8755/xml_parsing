import biopython as Bio
import pandas as pd
from Bio import Entrez
import csv
import time

def fetch_pubmed_data(query, email, batch_size=200, retmax=100000, out_csv='dhodh_inh_pubmed.csv'):
    # Fetch the data
    Entrez.email = "Nikita.Abramenko@lf1.cuni.cz"

# 1. Search for the records
    handle = Entrez.esearch(
        db="pubmed",
        term=query,
        retmax=100000)
    results = Entrez.read(handle)
    handle.close()

    id_list = results['IdList']
    count = len(id_list)
    print(f"Found {count} records")

# 2. CSV output
    with open(out_csv, 'w', encoding='utf-8', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(['PMID', 'Title', 'DOI'])

        # 3. Fetch the records in batches
        for start in range(0, count, batch_size):
            end = min(count, start + batch_size)
            id_batch = id_list[start:end]

        # 4. esummary
            summary_handle = Entrez.esummary(
                db="pubmed",
                id=id_batch)
            summary_results = Entrez.read(summary_handle)
            summary_handle.close()

            for record in summary_results:
                pmid = record.get('Id', '')
                title = record.get('Title', '')
                doi = record.get('DOI', '')
                csv_writer.writerow([pmid, title, doi])
            time.sleep(1)

    print(f"Data saved to {out_csv}")

if __name__ == '__main__':
    fetch_pubmed_data("dihydroorotate dehydrogenase inhibitors",
                      email = "Nikita.Abramenko@lf1.cuni.cz",
                    batch_size=200, retmax=100000, out_csv='dhodh_inh_pubmed.csv')


