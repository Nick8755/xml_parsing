import pandas as pd
import json
import biopython as Bio
from Bio import Entrez
import openpyxl as op

"""
Email: 'nikita.abramenko@lf1.cuni.cz'
API key: '10bff70a7e2be3a20316562bd3bfe7d34009'
"""



# Set the email address to avoid any potential issues with Entrez
Entrez.email = 'nikita.abramenko@lf1.cuni.cz'

# Define lists of authors and topics
#authors = ['Bryan Holland', 'Mehmet Oz', 'Anthony Fauci']  # Example authors, adjust as needed
topics = ['photodynamic therapy AND immune checkpoint inhibitor']  # Example topics, adjust as needed

# Define date range
date_range = '("2000/03/01"[Date - Create] : "2024/12/31"[Date - Create])'

# Build the query dynamically based on the available authors and topics
queries = []

#if authors:
#    author_queries = ['{}[Author]'.format(author) for author in authors]
#    queries.append('(' + ' OR '.join(author_queries) + ')')

if topics:
    topic_queries = ['{}[Title/Abstract]'.format(topic) for topic in topics]
    queries.append('(' + ' OR '.join(topic_queries) + ')')

full_query = ' AND '.join(queries) + ' AND ' + date_range

# Search PubMed for relevant records
handle = Entrez.esearch(db='pubmed', retmax=4000, term=full_query)
record = Entrez.read(handle)
id_list = record['IdList']

# DataFrame to store the extracted data
df = pd.DataFrame(columns=['PMID', 'Title', 'Abstract', 'Authors', 'Journal', 'Keywords', 'URL', 'Affiliations'])

# Fetch information for each record in the id_list
for pmid in id_list:
    handle = Entrez.efetch(db='pubmed', id=pmid, retmode='xml')
    records = Entrez.read(handle)

    # Process each PubMed article in the response
    for record in records['PubmedArticle']:
        # Print the record in a formatted JSON style
        print(json.dumps(record, indent=4, default=str))  # default=str handles types JSON can't serialize like datetime
        title = record['MedlineCitation']['Article']['ArticleTitle']
        abstract = ' '.join(record['MedlineCitation']['Article']['Abstract']['AbstractText']) if 'Abstract' in record['MedlineCitation']['Article'] and 'AbstractText' in record['MedlineCitation']['Article']['Abstract'] else ''
        keywords = ', '.join(keyword['DescriptorName'] for keyword in record['MedlineCitation']['MeshHeadingList']) if 'MeshHeadingList' in record['MedlineCitation'] else ''
        #authors = ', '.join(author.get('LastName', '') + ' ' + author.get('ForeName', '') for author in record['MedlineCitation']['Article']['AuthorList'])
        
        '''
        affiliations = []
        for author in record['MedlineCitation']['Article']['AuthorList']:
            if 'AffiliationInfo' in author and author['AffiliationInfo']:
                affiliations.append(author['AffiliationInfo'][0]['Affiliation'])
        affiliations = '; '.join(set(affiliations))

        journal = record['MedlineCitation']['Article']['Journal']['Title']
        keywords = ', '.join(keyword['DescriptorName'] for keyword in record['MedlineCitation']['MeshHeadingList']) if 'MeshHeadingList' in record['MedlineCitation'] else ''
        url = f"https://www.ncbi.nlm.nih.gov/pubmed/{pmid}"
    

        new_row = pd.DataFrame({
            'PMID': [pmid],
            'Title': [title],
            'Abstract': [abstract],
            'Authors': [authors],
            'Journal': [journal],
            'Keywords': [keywords],
            'URL': [url],
            'Affiliations': [affiliations]
        })
        '''

        new_row_short = pd.DataFrame({
            'PMID': [pmid],
            'Title': [title],
            'Abstract': [abstract],
            'Keywords': [keywords]
        })

        df = pd.concat([df, new_row_short], ignore_index=True)

# Save DataFrame to an Excel file
df.to_csv('PDT and ICI.csv', index=False)
