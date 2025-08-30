#if the script doesn't work, make sure you have installed all 
packages: mysql-connector-python requests pandas pybiomart 
biopython 
 
#import all packages we need 
import mysql.connector 
from pybiomart import Server 
import pandas as pd 
from Bio import Entrez 
import requests 
 
#connect to my database and create a cursor 
db = mysql.connector.connect ( 
 host="localhost", 
 user="B271542", 
 password="c3xMLhYK", 
 database="B271542") 
 
cursor = db.cursor() 
 
#put all the ensembl id in gene_ids 
gene_ids = [ 
 "ENSMUSG00000036061", "ENSMUSG00000000555", 
"ENSMUSG00000023055",  "ENSMUSG00000075394", "ENSMUSG00000001655", 
"ENSMUSG00000022485", 
 "ENSMUSG00000001657", "ENSMUSG00000001661", 
"ENSMUSG00000076010", 
 "ENSMUSG00000023048"] 
 
##############################################################
##############################################################
########1 
 
#use pybiomart to get data from Ensembl(the first database) 
 
# Connect to Ensembl BioMart server 
server = Server(host='http://www.ensembl.org') 
 
# Access the Ensembl Genes database and the mouse gene dataset 
(Mus musculus) 
database = server['ENSEMBL_MART_ENSEMBL'] 
dataset = database['mmusculus_gene_ensembl'] 
 
# Query the dataset for specific gene IDs 
result=dataset.query( 
 attributes=['ensembl_gene_id','entrezgene_id','mgi_id','ext
ernal_gene_name','chromosome_name','description'], 
filters={'link_ensembl_gene_id': gene_ids} 
) 
 
# Convert result to DataFrame 
df_ensembl = pd.DataFrame(result)  
##############################################################
##############################################################
######2 
 
#using biopython to get data from NCBI(the second database) 
 
#Take NCBI ids from df_ensembl and format them 
NCBI_IDs = df_ensembl['NCBI gene (formerly Entrezgene) 
ID'].tolist() 
 
# Provide my email (required by NCBI) 
Entrez.email = "sxxxxxxxx@ed.ac.uk" 
 
# Function to fetch xml data and transfer it to a readable format 
 
fetch_handle = Entrez.efetch(db="gene", id=NCBI_IDs, 
retmode="xml") 
gene_record = Entrez.read(fetch_handle) 
fetch_handle.close() 
 
# Prepare list to hold information I need 
gene_data = [] 
 
# Extract information I need from gene_record 
for record in gene_record: 
 gene_info = { 
 "GeneID": record["Entrezgene_track-info"]["Gene-track"]["Gene-track_geneid"], 
 "Name": record["Entrezgene_gene"]["Gene-ref"]["Generef_locus"],

 "Description": record["Entrezgene_gene"]["Generef"].get("Gene-ref_desc",
 ""), 
 "Organism": 
record["Entrezgene_source"]["BioSource"]["BioSource_org"]["Org
-ref"]["Org-ref_taxname"], 
 "Synonyms": ", ".join(record["Entrezgene_gene"]["Generef"].get("Gene-ref_syn",
 [])), 
 "Chromosome": 
record["Entrezgene_source"]["BioSource"]["BioSource_subtype"][
0]["SubSource_name"], 
 "Genomic_Start": record["Entrezgene_locus"][0]["Genecommentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]["Seqinterval_from"],

 "Genomic_End": record["Entrezgene_locus"][0]["Genecommentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]["Seqinterval_to"],

 "Gene_Type": record["Entrezgene_type"] 
 } 
 gene_data.append(gene_info) 
 
# Convert gene_data to DataFrame 
df_NCBI = pd.DataFrame(gene_data) 
 
##############################################################
############################################################3 
 
#use request to get data from Uniprot(the third database)  
# Use pybioMart to get Uniprot/Swiss ID, don't combine this with 
the first pybiomart query!!! 
# Because these is a gene having no Uniprot/Swiss ID, it will 
disappear in result if you combined the queries 
 
result=dataset.query( 
attributes=['ensembl_gene_id','uniprotswissprot'], 
filters={'link_ensembl_gene_id': gene_ids} 
) 
 
# Convert result to DataFrame, and remove the meaningless lines 
having NaN(if you print result you will see them) 
df = pd.DataFrame(result) 
filtered_df = df.dropna(subset=['UniProtKB/Swiss-Prot ID']) 
 
# Get Uniprot ids from dataframe 
uniprot_ids = filtered_df['UniProtKB/Swiss-Prot ID'].tolist() 
 
# Set url for UniProt API 
url = "https://rest.uniprot.org/uniprotkb/search" 
 
# Prepare a list to store results 
protein_data = [] 
 
# Fetch data for each UniProt ID 
for uniprot_id in uniprot_ids: 
 #set parameters and get result in json format  params = { 
 "query": f"accession:{uniprot_id}", 
 "fields": 
"accession,protein_name,length,cc_function,cc_subcellular_loca
tion", 
 "format": "json" 
 } 
 r = requests.get(url, params=params) 
 data = r.json() 
 
 #extract the data I want and store them into protein_data 
 entry = data["results"][0] 
 protein_data.append({ 
 "UniProtKB/Swiss-Prot ID": uniprot_id, 
 "Protein_Name": entry.get("proteinDescription", 
{}).get("recommendedName", {}).get("fullName", {}).get("value", 
"N/A"), 
 "Protein_Length": entry.get("sequence", 
{}).get("length", "N/A"), 
 "Subcellular Location": 
next((item["subcellularLocations"][0]["location"]["value"] 
 for item in entry.get("comments", []) 
 if item.get("commentType") == "SUBCELLULAR 
LOCATION"), "N/A") 
 }) 
 
# Convert the list to a DataFrame 
df = pd.DataFrame(protein_data) 
 # Merge it with filtered_df on 'UniProtKB/Swiss-Prot ID' column 
result = pd.merge(filtered_df, df, on='UniProtKB/Swiss-Prot ID', 
how='inner') 
 
# Add a new row for that missing gene (it will not produce 
protein) 
new_row = {'Gene stable ID': 'ENSMUSG00000076010', 
 'UniProtKB/Swiss-Prot ID': 'N/A', 'Protein_Name': 
'N/A', 'Protein_Length': 'N/A', 'Subcellular Location': 'N/A'} 
df_UniProt = result._append(new_row, ignore_index=True) 
 
##############################################################
##############################################################
#### 
 
#create three table for in database 
 
#first, Delete tables that have the same name as the table we 
are going to create 
cursor.execute("DROP TABLE IF EXISTS ensembl") 
cursor.execute("DROP TABLE IF EXISTS NCBI") 
cursor.execute("DROP TABLE IF EXISTS UniProt") 
 
#create three tables 
create_table_ensembl = """ 
 CREATE TABLE ensembl ( 
 ensembl_gene_id VARCHAR(50) NOT NULL, 
 entrezgene_id INT, 
 mgi_id VARCHAR(50),  external_gene_name VARCHAR(50), 
 chromosome_name VARCHAR(50), 
 description VARCHAR(100), 
 PRIMARY KEY (ensembl_gene_id) 
 ) 
 """ 
 
create_table_NCBI = """ 
 CREATE TABLE NCBI ( 
 GeneID VARCHAR(50) NOT NULL, 
 Name VARCHAR(50), 
 Description VARCHAR(100), 
 Organism VARCHAR(50), 
 Synonyms VARCHAR(50), 
 Chromosome VARCHAR(50), 
 Genomic_Start VARCHAR(50), 
 Genomic_End VARCHAR(50), 
 Gene_Type VARCHAR(50) 
 ) 
 """ 
 
create_table_UniProt = """ 
 CREATE TABLE UniProt ( 
 Gene_stable_ID VARCHAR(50) NOT NULL, 
 UniProt_ID VARCHAR(50), 
 Protein_Name VARCHAR(100), 
 Protein_Length VARCHAR(50),  Subcellular_Location VARCHAR(50) 
 ) 
 """ 
 
cursor.execute(create_table_ensembl) 
cursor.execute(create_table_NCBI) 
cursor.execute(create_table_UniProt) 
 
##############################################################
##############################################################
### 
 
# Insert data into MySQL tables 
for _, row in df_ensembl.iterrows(): 
 cursor.execute(""" 
 INSERT INTO ensembl (ensembl_gene_id, 
entrezgene_id, mgi_id, external_gene_name, chromosome_name, 
description) 
 VALUES (%s, %s, %s, %s, %s, %s) 
 """, list(row)) 
 
for _, row in df_NCBI.iterrows(): 
 cursor.execute(""" 
 INSERT INTO NCBI (GeneID, Name, Description, 
Organism, Synonyms, Chromosome, Genomic_Start, Genomic_End, 
Gene_Type) 
 VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s) 
 """, list(row)) 
 for _, row in df_UniProt.iterrows(): 
 cursor.execute(""" 
 INSERT INTO UniProt (Gene_stable_ID, UniProt_ID, 
Protein_Name, Protein_Length, Subcellular_Location) 
 VALUES (%s, %s, %s, %s, %s) 
 """, list(row)) 
 
# Commit changes 
db.commit() 
 
##############################################################
##############################################################
### 
 
# Define the SQL query to join the tables and fetch the data 
query = """ 
SELECT 
 ensembl.ensembl_gene_id, 
 ensembl.entrezgene_id, 
 ensembl.mgi_id, 
 ensembl.external_gene_name, 
 ensembl.chromosome_name, 
 NCBI.Organism, 
 NCBI.Synonyms, 
 NCBI.Genomic_Start, 
 NCBI.Genomic_End, 
 NCBI.Gene_Type, 
 UniProt.UniProt_ID,  UniProt.Protein_Length, 
 UniProt.Subcellular_Location, 
 ensembl.description, 
 UniProt.Protein_Name 
FROM ensembl 
LEFT JOIN NCBI ON ensembl.entrezgene_id = NCBI.GeneID 
LEFT JOIN UniProt ON ensembl.ensembl_gene_id = 
UniProt.Gene_stable_ID 
""" 
 
# Execute the query and fetch the result 
cursor.execute(query) 
results = cursor.fetchall() 
 
#Make result a table(dataframe) 
columns = ['ensembl_gene_id', 'entrezgene_id', 'mgi_id', 
'external_gene_name', 'chromosome_name', 
 'Organism', 'Synonyms', 'Genomic_Start', 
'Genomic_End', 'Gene_Type', 'UniProt_ID', 
 'Protein_Length', 'Subcellular_Location', 
'description', 'Protein_Name'] 
summary = pd.DataFrame(results, columns=columns) 
 
#Change the option to make all data visable, output the table 
pd.set_option('display.max_columns', None) 
pd.set_option('display.max_colwidth', 100) 
print(summary) 
