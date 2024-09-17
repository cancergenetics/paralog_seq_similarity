import os
import sys
import warnings
from itertools import groupby
from multiprocessing import Pool
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.spatial.distance import cdist, squareform
from sklearn.manifold import trustworthiness
import torch
import argparse

sys.dont_write_bytecode = True
np.set_printoptions(precision=6, suppress=True)

from my_library import Database, ESM_Model, Metrics, ESM_Representations, read_fasta, gzip_tensor, validate_fasta, validate_database, neighbor_joining, cophenetic_distmat, silhouette

#==============================================================================
# Parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("threads", help = "Number of threads to allocate")
    parser.add_argument("fasta_file", help = "Fasta file of sequences used to generate embeddings")
    parser.add_argument("hgnc", help = "Gene pairs informations (eg gene_with_protein_product.txt)")
    parser.add_argument("paralog_pairs", help = "Paralog pairs of interest (need a sorted_gene_pair column)")
    parser.add_argument("device", help = "Device to allocate (eg cuda or cpu)")
    return parser.parse_args()

args = parser()

# Allocate ressources
THREADS = f'{args.threads}'
DEVICE  = f'{args.device}'
THREADS = 30
DEVICE  = 'cpu'

# Input files
FASTA_FILE = args.fasta_file #'paralog_sequences.fasta'
DB_FILE = 'esm2.db'
HGNC_FILE = args.hgnc #'/Users/olivierdennler/Documents/data/SLI_2023/gene_with_protein_product.txt'
LABEL_CSV  = ''
OUTPUT_DIR = 'paralog_models_esm2'
OUTPUT_FILE = 'esm2_embeddings_pair_distances.csv'
PAIRS_FILE = args.paralog_pairs

# ====================================================================
# Load Protein Language Model
# ====================================================================

try:
    encoder # if the model has already been loaded, do not reload
except:
    encoder = ESM_Model('esm2_t48_15B_UR50D') 

# ====================================================================
# Generate sequence embeddings
# ====================================================================

# ensure that an existing database does not exist
if os.path.exists(DB_FILE):
    print('Use existing DB ...')

else:
    # create a SQLite database to store the embeddings
    db = Database(DB_FILE)
    db.create_table(columns=[
        ('header',    'TEXT'),
        ('sequence',  'TEXT'),
        ('embedding', 'BLOB'),
    ])

    # define function for generating, then gzipping the embeddings for storage
    func_encode = lambda s: gzip_tensor(encoder.encode(s, device=DEVICE, threads=THREADS).type(torch.float16))

    # Safe function, in case of crash for one embeddigs, the whole execution does not crash (added)
    def safe_func_encode(h, s):
        try:
            return func_encode(s)
        except Exception as e:
            print(f"Error encoding sequence: {h}\n{s}. Error: {e}")
            return None

    # iteratively run the function on each sequence
    queue = ((h, s, safe_func_encode(h, s)) for h, s in read_fasta(FASTA_FILE))
    db.add_rows(('header', 'sequence', 'embedding'), queue)

# ====================================================================
# Distances
# ====================================================================

# load the database
db = Database(DB_FILE)

# Modified lambda function to handle NoneType
_format = lambda x: (x['header'], x['sequence'], gzip_tensor(x['embedding']).numpy()) if x['embedding'] is not None else None
# Use list comprehension with a condition to filter out None values
formatted_data = [_format(i) for i in db.retrieve() if i['embedding'] is not None]
# Unpack the headers, sequences, and embeddings, skipping entries that are None
headers, sequences, embeddings = zip(*[item for item in formatted_data if item is not None])

headers    = np.array(headers   , dtype=object)
accessions = np.array([i.split()[0] for i in headers], dtype=object)
sequences  = np.array(sequences , dtype=object)
embeddings = np.array(embeddings, dtype=object)

# if labels are available load them from csv
if os.path.exists(LABEL_CSV):
    labels = dict(pd.read_csv(LABEL_CSV)[['accession','label']].fillna('').values)
    labels = np.array([labels[i] for i in accessions], dtype=object)
else:
    labels = []

# representations and metrics can be added and removed from this list
representations = {
    'beginning_of_sequence'  : ESM_Representations.beginning_of_sequence,
    'end_of_sequence'        : ESM_Representations.end_of_sequence,
    'mean_of_special_tokens' : ESM_Representations.mean_special_tokens,
    'mean_of_residue_tokens' : ESM_Representations.mean_residue_tokens,
    }

metrics = {
    'cosine'        : Metrics.cosine,
    'euclidean'     : Metrics.euclidean,
    'manhattan'     : Metrics.manhattan,
    'ts_ss'         : Metrics.ts_ss,
    # 'jensenshannon' : Metrics.jensenshannon, #not me
    }

# Init distance metrics
# create base directory
if not os.path.exists(OUTPUT_DIR):
    os.mkdir(OUTPUT_DIR)
if not os.path.exists(f'{OUTPUT_DIR}/fixedsize'):
    os.mkdir(f'{OUTPUT_DIR}/fixedsize')

models = []
for _rep in representations:
    rep = np.array([representations[_rep](i) for i in embeddings])
    np.savez_compressed(f'{OUTPUT_DIR}/fixedsize/{_rep}.npz', 
                        **{'headers':headers, 'embedding':rep})
    for _met in metrics:
        distmat = metrics[_met](rep,rep)
        models += [{
            'representation'          : _rep,
            'metric'                  : _met,
            'distmat(repr)'           : distmat}]
        if len(labels) == distmat.shape[0]:
            models[-1]['silhouette(repr)'] = silhouette(distmat, labels)
        
        sys.stderr.write(f'Calculating distances using "{_rep}" with "{_met}"\n')

# ====================================================================
# Format distance as pair metrics
# ====================================================================

# Optimized by ChatGPT

# Use uniprot gene mapping
# Efficient data loading
cols_to_use_uniprot = ['gene', 'uniprot']  # Specify the columns you need
df_uniprot_infos = pd.read_table(HGNC_FILE, usecols=cols_to_use_uniprot, sep=',', low_memory=False)

# Preprocessing to speed up lookups
uniprot_to_symbol = df_uniprot_infos.set_index('uniprot')['gene'].to_dict()

# Load paralog gene pair to only write line for them
gene_pairs_full = pd.read_csv(PAIRS_FILE, usecols=['sorted_gene_pair'])
# Convert gene_pairs_full to a set for faster lookup
gene_pairs_set = set(gene_pairs_full['sorted_gene_pair'])

# Initialize DataFrame only once
df = pd.DataFrame()

# Pre-define a list to collect all rows across all types of distances
all_rows = []

# Placeholder for the new column names to ensure they match with each model's representation-metric
new_column_names = []

for type_dist in models:
    representation = type_dist['representation']
    metric = type_dist['metric']
    matrice = type_dist['distmat(repr)']  # Assuming this returns the distance matrix
    
    # Add the new column name for the current representation-metric to the list
    new_column_name = f'{representation}-{metric}'
    new_column_names.append(new_column_name)

    for i in range(len(headers)):
        uniprot_i = headers[i].split('|')[1]
        locus_i = uniprot_to_symbol.get(uniprot_i, '')

        for j in range(len(headers)):
            uniprot_j = headers[j].split('|')[1]
            locus_j = uniprot_to_symbol.get(uniprot_j, '')
            key = f"{locus_i}_{locus_j}"
            
            # Only look at known paralog pairs
            if key in gene_pairs_set:
                all_rows.append([locus_i, locus_j, matrice[i][j], new_column_name])

# Convert all collected data into a DataFrame
columns=['locus_i', 'locus_j', 'distance', 'distance_type']
temp_df = pd.DataFrame(all_rows, columns=columns)

# Pivot the DataFrame to have each distance type as a column
pivot_df = temp_df.pivot_table(index=['locus_i', 'locus_j'], columns='distance_type', values='distance', aggfunc='first').reset_index()

# Merge 'locus_i' and 'locus_j' into a single column for the sorted gene pair
pivot_df['sorted_gene_pair'] = pivot_df['locus_i'] + '_' + pivot_df['locus_j']
pivot_df.drop(['locus_i', 'locus_j'], axis=1, inplace=True)

# Reorder the DataFrame to have 'sorted_gene_pair' as the first column
col_order = ['sorted_gene_pair'] + [col for col in pivot_df.columns if col != 'sorted_gene_pair']
df = pivot_df[col_order]

# At this point, `df` contains your optimized DataFrame


# Write dataset in csv file
df.to_csv(OUTPUT_FILE, sep=',')
