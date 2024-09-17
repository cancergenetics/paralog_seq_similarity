#!/bin/python3

"""Compare and compute TM-score between structures of paralog pairs"""


import os
import argparse
import subprocess
import pandas as pd
import gzip
from pathlib import Path
from tqdm.auto import tqdm
from bs4 import BeautifulSoup

#==============================================================================
# Align all pairs
#==============================================================================

def align_pairs(paralog_pairs, hgnc_infos, structures_directory, aln_type) -> None:
    """
    For all paralog pairs, align their protein structures
    """
    # Building dataframes
    df_paralogs_pairs = pd.read_csv(paralog_pairs)
    df_hgnc_infos = pd.read_table(hgnc_infos, low_memory=False)

    # Iterate on all rows ie, paralog pairs, paralogs being genes (maybe change it later by an vectoriel operation...)
    for index, row in tqdm(df_paralogs_pairs.iterrows()):

        # Entrez_id of the 2 paralog genes
        entrez_id_A1 = row['A1_entrez']
        entrez_id_A2 = row['A2_entrez']
        # Get the corresponding uniprot id
        uniprot_id_A1 = df_hgnc_infos.loc[df_hgnc_infos['entrez_id'] == entrez_id_A1, 'uniprot_ids'].squeeze()
        uniprot_id_A2 = df_hgnc_infos.loc[df_hgnc_infos['entrez_id'] == entrez_id_A2, 'uniprot_ids'].squeeze()
        
        # Get the corresponding AlphaFold filenames
        # TODO : fix or clarify these exceptions
        try:
            fn_struct_A1 = uniprot_to_AFfn(uniprot_id_A1, structures_directory)
        except:
            print(f'Issue with entrez_id {entrez_id_A1}, or uniprot {uniprot_id_A1}. Check for HGNC or AlphaFold available informations !')
        try:
            fn_struct_A2 = uniprot_to_AFfn(uniprot_id_A2, structures_directory)
        except:
            print(f'Issue with entrez_id {entrez_id_A2}, or uniprot {uniprot_id_A2}. Check for HGNC or AlphaFold available informations !')
        # Align the structures
        try:
            if aln_type == 'foldseek':
                pairwise_foldseek(fn_struct_A1, fn_struct_A2, row['sorted_gene_pair']) #Foldseek, 3Di local
            elif aln_type == 'foldseek_tmalign':
                pairwise_foldseek_tmalign_tab(fn_struct_A1, fn_struct_A2, row['sorted_gene_pair']) #Foldseek, tmalign global, tab output
        except:
            print(f'Issue with AlphaFold files, check if there is a pdb file for : {uniprot_id_A1}({fn_struct_A1}) and {uniprot_id_A2}({fn_struct_A1})')

def align_vsDB(paralog_pairs, fn_db, hgnc_infos, gene_uniq_uniprot, structures_directory, aln_type) -> None:
    """
    For all paralog pairs, align their protein structures against a structural proteom database
    1 foldseek aln for each paralog of the pair
    """
    # Building dataframes
    df_paralogs_pairs = pd.read_csv(paralog_pairs)
    #df_hgnc_infos = pd.read_table(hgnc_infos, low_memory=False)
    df_gene_uniprot = pd.read_csv(gene_uniq_uniprot, usecols=['gene', 'uniprot'])

    # Output file building (will be completed during iteration/computation)
    r_line = '\n'
    with open('vsAll_scores.csv', 'w') as csv_file:
        csv_file.write(f'sorted_gene_pair,A1,A2,uniprot_id_A1,uniprot_id_A2,A1_A2_rank,A2_A1_rank,A1_nb_human,A2_nb_human,A1_nb_taxid,A2_nb_taxid,A1_bits,A2_bits,A1_fident,A2_fident,A1_alnlen,A2_alnlen,A1_evalue,A2_evalue,A1_alntmscore,A2_alntmscore,A1_qtmscore,A2_qtmscore,A1_ttmscore,A2_ttmscore,A1_prob,A2_prob,A1_lddt,A2_lddt,A1_lddtfull,A2_lddtfull{r_line}')

        # Iterate on all rows ie, paralog pairs, paralogs being genes (maybe change it later by an vectoriel operation...)
        for index, row in tqdm(df_paralogs_pairs.iterrows()):

            # New version, use the file with mapping - gene : uniprot
            # We can directly use gene pair
            pair = row['sorted_gene_pair']
            A1, A2 = pair.split('_')
            # Iterate on gene and get their uniprots
            for gene_id_i in [A1,A2]:
                uniprot_id_i = str(df_gene_uniprot.loc[df_gene_uniprot['gene'] == gene_id_i, 'uniprot'].squeeze())
                if gene_id_i == A1: 
                    paralog = 'A1'
                    uniprot_id_A1 = uniprot_id_i
                elif gene_id_i == A2: 
                    paralog = 'A2'
                    uniprot_id_A2 = uniprot_id_i

                # Get the corresponding AlphaFold filenames
                # TODO : fix or clarify these exceptions
                try:
                    fn_struct_i = uniprot_to_AFfn(uniprot_id_i, structures_directory)
                except:
                    print(f'Issue with entrez_id, or uniprot {uniprot_id_i}. Check for HGNC or AlphaFold available informations !')
                # Align against database
                try:
                    vsDB_foldseek(fn_struct_i, fn_db, row[paralog], aln_type)
                except:
                    print(f'Issue with AlphaFold files, check if there is a pdb file for : {uniprot_id_i}({fn_struct_i})')

            # Define a pair based on both paralogs vs all
            pair = row['sorted_gene_pair']
            A1, A2 = pair.split('_')
            """
            A distance btw both paralogs A1 and A2 of a pair is described with :
            For A1 vs all (db) :
                - A2 rank
                - Number of human proteins better than A2
                - Number of different taxids with proteins better than A2
            For A2 vs all (db) :
                - A1 rank
                - Number of human proteins better than A1
                - Number of different taxids with proteins better than A1
            """
            #try:
            # Check if there is a aln file, and a uniprot for the paralog to search for
            if Path(f"foldseek_{aln_type}_{A1}.tsv").is_file() and len(uniprot_id_A2) > 0:
                # A1 vs all information
                df_A1_all = pd.read_table(f"foldseek_{aln_type}_{A1}.tsv")
                # Start by searching for A2 hit
                A2_row = df_A1_all[df_A1_all['target'].str.contains(uniprot_id_A2)]
                # If A2 is in the hits of A1 vs all
                if len(A2_row.index.tolist()) > 0:
                    # The rank is the index of the line (note that the query is supposed to be the index 0)
                    A1_A2_rank = int(A2_row.index.tolist()[0])
                    # Select the row corresponding better hits than the paralog
                    df_better_than_A2 = df_A1_all.loc[range(1,A1_A2_rank,1),]
                    # Number of human protein that are better hits
                    A1_nb_human = len(df_better_than_A2.loc[df_better_than_A2['taxid'] == 9606])
                    # Number of different taxid that have proteins being better hits
                    A1_nb_taxid = len(set(df_better_than_A2.loc[df_better_than_A2['taxid'] != 9606]['taxid'].tolist()))
                    # Other aln features (pairwise)
                    A1_bits = A2_row['bits'].values[0]
                    A1_fident = A2_row['fident'].values[0]
                    A1_alnlen = A2_row['alnlen'].values[0]
                    A1_evalue = A2_row['evalue'].values[0]
                    A1_alntmscore = A2_row['alntmscore'].values[0]
                    A1_qtmscore = A2_row['qtmscore'].values[0]
                    A1_ttmscore = A2_row['ttmscore'].values[0]
                    A1_prob = A2_row['prob'].values[0]
                    A1_lddt = A2_row['lddt'].values[0]
                    A1_lddtfull = A2_row['lddtfull'].values[0]
                # If A2 is not in the hits of A1 vs all
                else:
                    print(f'{A2} is not a hit of {A1} vs all -> introduce NA_hit')
                    A1_A2_rank = 'NA_hit'
                    A1_nb_human = 'NA_hit'
                    A1_nb_taxid = 'NA_hit'
                    A1_bits = 'NA_hit'
                    A1_fident = 'NA_hit'
                    A1_alnlen = 'NA_hit'
                    A1_evalue = 'NA_hit'
                    A1_alntmscore = 'NA_hit'
                    A1_qtmscore = 'NA_hit'
                    A1_ttmscore = 'NA_hit'
                    A1_prob = 'NA_hit'
                    A1_lddt = 'NA_hit'
                    A1_lddtfull = 'NA_hit'
            # If there is no aln file
            else:
                print(f'No aln file for {A1} vs all -> introduce NA_aln')
                A1_A2_rank = 'NA_aln'
                A1_nb_human = 'NA_aln'
                A1_nb_taxid = 'NA_aln'
                A1_bits = 'NA_aln'
                A1_fident = 'NA_aln'
                A1_alnlen = 'NA_aln'
                A1_evalue = 'NA_aln'
                A1_alntmscore = 'NA_aln'
                A1_qtmscore = 'NA_aln'
                A1_ttmscore = 'NA_aln'
                A1_prob = 'NA_aln'
                A1_lddt = 'NA_aln'
                A1_lddtfull = 'NA_aln'

            # Check if there is a aln file, and a uniprot for the paralog to search for
            if Path(f"foldseek_{aln_type}_{A2}.tsv").is_file() and len(uniprot_id_A1) > 0:
                # A2 vs all information
                df_A2_all = pd.read_table(f"foldseek_{aln_type}_{A2}.tsv")
                # Start by searching for A1 hit
                A1_row = df_A2_all[df_A2_all['target'].str.contains(uniprot_id_A1)]
                # If A1 is in the hits of A2 vs all
                if len(A1_row.index.tolist()) > 0:
                    # The rank is the index of the line (note that the query is supposed to be the index 0)
                    A2_A1_rank = int(A1_row.index.tolist()[0])
                    # Select the row corresponding better hits than the paralog
                    df_better_than_A1 = df_A2_all.loc[range(1,A2_A1_rank,1),]
                    # Number of human protein that are better hits
                    A2_nb_human = len(df_better_than_A1.loc[df_better_than_A1['taxid'] == 9606])
                    # Number of different taxid that have proteins being better hits
                    A2_nb_taxid = len(set(df_better_than_A1.loc[df_better_than_A1['taxid'] != 9606]['taxid'].tolist()))
                    # Other aln features (pairwise)
                    A2_bits = A1_row['bits'].values[0]
                    A2_fident = A1_row['fident'].values[0]
                    A2_alnlen = A1_row['alnlen'].values[0]
                    A2_evalue = A1_row['evalue'].values[0]
                    A2_alntmscore = A1_row['alntmscore'].values[0]
                    A2_qtmscore = A1_row['qtmscore'].values[0]
                    A2_ttmscore = A1_row['ttmscore'].values[0]
                    A2_prob = A1_row['prob'].values[0]
                    A2_lddt = A1_row['lddt'].values[0]
                    A2_lddtfull = A1_row['lddtfull'].values[0]

                # If A1 is not in the hits of A2 vs all
                else:
                    print(f'{A1} is not a hit of {A2} vs all -> introduce NA_hit')
                    A2_A1_rank = 'NA_hit'
                    A2_nb_human = 'NA_hit'
                    A2_nb_taxid = 'NA_hit'
                    A2_bits = 'NA_hit'
                    A2_fident = 'NA_hit'
                    A2_alnlen = 'NA_hit'
                    A2_evalue = 'NA_hit'
                    A2_alntmscore = 'NA_hit'
                    A2_qtmscore = 'NA_hit'
                    A2_ttmscore = 'NA_hit'
                    A2_prob = 'NA_hit'
                    A2_lddt = 'NA_hit'
                    A2_lddtfull = 'NA_hit'
            # If there is no aln file
            else:
                print(f'No aln file for {A2} vs all -> introduce NA_aln')
                A2_A1_rank = 'NA_aln'
                A2_nb_human = 'NA_aln'
                A2_nb_taxid = 'NA_aln'
                A2_bits = 'NA_aln'
                A2_fident = 'NA_aln'
                A2_alnlen = 'NA_aln'
                A2_evalue = 'NA_aln'
                A2_alntmscore = 'NA_aln'
                A2_qtmscore = 'NA_aln'
                A2_ttmscore = 'NA_aln'
                A2_prob = 'NA_aln'
                A2_lddt = 'NA_aln'
                A2_lddtfull = 'NA_aln'

            # Write a line in output, 1 line : 1 pair
            csv_file.write(f'{pair},{A1},{A2},{uniprot_id_A1},{uniprot_id_A2},{A1_A2_rank},{A2_A1_rank},{A1_nb_human},{A2_nb_human},{A1_nb_taxid},{A2_nb_taxid},{A1_bits},{A2_bits},{A1_fident},{A2_fident},{A1_alnlen},{A2_alnlen},{A1_evalue},{A2_evalue},{A1_alntmscore},{A2_alntmscore},{A1_qtmscore},{A2_qtmscore},{A1_ttmscore},{A2_ttmscore},{A1_prob},{A2_prob},{A1_lddt},{A2_lddt},{A1_lddtfull},{A2_lddtfull}{r_line}')

#==============================================================================
# Foldseek call
#==============================================================================

def pairwise_foldseek(fn_struct_A1, fn_struct_A2, fn_out) -> None:
    """
    Run foldseek to pairwise align 2 structures
    """
    # Check if aln exists, if not, align
    if Path(f"foldseek_3DiAA_{fn_out}.tsv").is_file() != True:
        # Get the length of the structured proteins
        len_A1 = length_from_pdbgz(fn_struct_A1)
        len_A2 = length_from_pdbgz(fn_struct_A2)
        # Align the shortest (query) on the longest (target)
        if int(len_A1) > int(len_A2):
            query = fn_struct_A2
            target = fn_struct_A1
        else:
            query = fn_struct_A1
            target = fn_struct_A2
        # Align with foldseek
        subprocess.run([
            "foldseek", "easy-search",
            f"{query}",
            f"{target}",
            f"foldseek_3DiAA_{fn_out}.tsv",
            "tmp",
            "-v", f"3",
            "--format-mode", "4",
            '--format-output', "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore,qtmscore,ttmscore,prob",
            '--alignment-type', '2'
            ],
            stdout=open(os.devnull, 'wb')
                    )

def pairwise_foldseek_tmalign_tab(fn_struct_A1, fn_struct_A2, fn_out) -> None:
    """
    Run foldseek to pairwise align 2 structures
    """
    # Check if aln exists, if not, align
    if Path(f"foldseek_tmalign_{fn_out}.tsv").is_file() != True:
        # Get the length of the structured proteins
        len_A1 = length_from_pdbgz(fn_struct_A1)
        len_A2 = length_from_pdbgz(fn_struct_A2)
        # Align the shortest (query) on the longest (target)
        if int(len_A1) > int(len_A2):
            query = fn_struct_A2
            target = fn_struct_A1
        else:
            query = fn_struct_A1
            target = fn_struct_A2
        # Align with foldseek
        subprocess.run([
            "foldseek", "easy-search",
            f"{query}",
            f"{target}",
            f"foldseek_tmalign_{fn_out}.tsv",
            "tmp",
            "-v", f"3",
            "--format-mode", "4",
            '--format-output', "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore,qtmscore,ttmscore,prob",
            '--alignment-type', '1'
                ],
            stdout=open(os.devnull, 'wb')
                    )

def vsDB_foldseek(fn_struct_A, fn_db, fn_out, aln_type='foldseek') -> None:
    """
    Run foldseek to align a protein structure against a structural database
    => Mainly used, cause context searches also include the pairwise features
    """
    # Check if aln exists, if not, align
    if Path(f"foldseek_{aln_type}_{fn_out}.tsv").is_file() != True:
        query = fn_struct_A
        # Aln type in foldseek format
        if aln_type == 'foldseek_tmalign':
            fa_type = '1'
        elif aln_type == 'foldseek':
            fa_type = '2'
        # Align with foldseek
        subprocess.run([
            "foldseek", "easy-search",
            f"{query}",
            f"{fn_db}",
            f"foldseek_{aln_type}_{fn_out}.tsv",
            "tmp",
            "-v", f"3",
            "--format-mode", "4",
            '--format-output', "query,target,fident,alnlen,taxid,taxname,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore,qtmscore,ttmscore,prob,lddt,lddtfull",
            '--alignment-type', f'{fa_type}',
            '--remove-tmp-files',
	        '--exhaustive-search'
            # Non inf version, for speed up
            #'-e', 'inf'
                ],
            stdout=open(os.devnull, 'wb')
                    )

#==============================================================================
# Tools
#==============================================================================

def length_from_pdbgz(fn_struct) -> int:
    """
    Read the protein lengths in the pdb compressed files
    """
    with gzip.open(fn_struct,'rb') as file:
        for line in file:
            l = line.decode('utf-8')
            l = l.replace('\n', '')
            if 'DBREF' in l:
                length = int([i for i in l.split(' ') if len(i) > 0][-1])
                return length
                break

def plddt_from_pdbgz(fn_struct) -> int:
    """
    Read the protein plddt in the pdb compressed files
    """
    # TODO
    plddt_list = []
    with gzip.open(fn_struct,'rb') as file:
        for line in file:
            l = line.decode('utf-8')
            l = l.replace('\n', '')
            # All atoms of a residue have the same plddt (still need to finish)
            if 'ATOM' in l:
                res_plddt = float([i for i in l.split(' ') if len(i) > 0][-2])
                continue
                break

def uniprot_to_AFfn(uniprot_id, struct_directory) -> str:
    """
    Get the filename of the AlphaFold structure (present in the structures directory) corresponding to the uniprot id
    """
    files = struct_directory.glob('*.pdb.gz')
    for file in files:
        if uniprot_id in file.stem:
            return file
            break

#==============================================================================
# Output files
#==============================================================================

def parse_foldseekHTML(fn_foldseek_html) -> dict:
    """
    Parse an foldseek output file in html format
    For a pairwise aln
    And return the different aln informations
    """
    with open(fn_foldseek_html, 'r') as fn_html:
        soup = BeautifulSoup(fn_html, 'html.parser')
        pair = '_'.join(fn_foldseek_html.stem.split('_')[1:])
        # Not elegant scraping, but will do what I need
        try:
            seq_identity = float(soup.prettify().split('"seqId": ')[1].split(',')[0])
            score = float(soup.prettify().split('"score": ')[1].split(',')[0])
            e_value = float(soup.prettify().split('"eval": ')[1].split(',')[0])
            return {
                'pair' : pair,
                'seq_identity' : seq_identity,
                'score' : score,
                'e_value' : e_value,
                }
        except:
            print(f'Delete {fn_foldseek_html} file (no scores in it...)')
            fn_foldseek_html.unlink(missing_ok=False)
            return 'Failed'

def cat_html_in_csv() -> None:
    """
    Extract all html informations into one csv file
    Plus, add the pair information in it (from filename)
    """
    html_files = Path.cwd().glob('*.html')
    scores_list = []
    for file in html_files:
        scores_list.append(parse_foldseekHTML(file))
    scores_list = [score for score in scores_list if score != 'Failed']
    with open('html_scores.csv', 'w') as csv_file:
        csv_file.write('sorted_gene_pair,seq_identity,score,e_value\n')
        for p in scores_list:
            csv_file.write(f'{p["pair"]},{p["seq_identity"]},{p["score"]},{p["e_value"]}\n')

def cat_tsv_in_csv() -> None:
    """
    Regroup all tsv informations into one csv file
    Plus, add the pair information in it (from filename)
    """
    tsv_files = Path.cwd().glob('foldseek*.tsv')
    header = False
    tab = '\t'
    r_line = '\n'
    with open('tsv_scores.csv', 'w') as csv_file:
        for file in tsv_files:
            with open(file, 'r') as tsv_str:
                for line in tsv_str:
                    # Header case (extract if not already done, else juste skip it)
                    if 'query' in line:
                        if header == False:
                            csv_file.write(f'sorted_gene_pair,{",".join(line.split(tab))}')
                            header = True
                    # Content lines
                    else:
                        csv_file.write(f'{"_".join(file.stem.split("_")[2:])},{",".join(line.split(tab))}')

#==============================================================================
# Parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("paralog_pairs",
                        help = "Csv file that describe paralog pairs (eg SI 8)",
                        type=str)
    parser.add_argument("hgnc_infos",
                        help = "Tsv file that describe all genes",
                        type=str)
    parser.add_argument("structures_directory",
                        help = "Directory containing all pdb files (eg UP000005640_9606_HUMAN_V4)",
                        type=str)
    parser.add_argument("map_gene_uniprot",
                        help = "csv file containing gene : uniprot mapping to use (eg map_gene_uniprot_ens111.csv)",
                        type=str)
    parser.add_argument("output",
                        help = "Name to give to the output csv file (format : sorted_gene_pair,TM-score,????)",
                        type=str)
    parser.add_argument("-V",
                        help = "Verbose level of foldseek (NYI)",
                        nargs='?',
                        const="3",
                        type=str)
    parser.add_argument("--aln",
                        help = "Type of alignement [foldeseek,foldseek_tmalign], mandatory if aln are wanted",
                        type=str)
    parser.add_argument("--only_concat",
                        help = "Concat result file from already computed alignments, in the current directory (tsv or html from foldseek)",
                        action="store_true")
    parser.add_argument("--pairwise",
                        help = "Pairwise aln btw both paralogs of a pair",
                        action="store_true")
    parser.add_argument("--db",
                        help = "Align each paralog against this database",
                        type=str)
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

def main():

    args = parser()

    # If alignments already exists and we only want to regroup results
    if args.only_concat:
        cat_html_in_csv()
        cat_tsv_in_csv()

    # To compute alignments
    if args.aln:
        # Load all pair info
        paralog_pairs = Path(args.paralog_pairs).resolve()
        hgnc_infos = Path(args.hgnc_infos).resolve()
        structures_directory = Path(args.structures_directory).resolve()
        map_gene_uniprot = Path(args.map_gene_uniprot).resolve()

        # Align pairs (pairwise)
        if args.pairwise:
            align_pairs(paralog_pairs, hgnc_infos, structures_directory, args.aln)

        # Align each paralog aigainst a database (vs all, including other species)
        if args.db:
            align_vsDB(paralog_pairs, Path(args.db).resolve(), hgnc_infos, map_gene_uniprot, structures_directory, args.aln)

        # Concat results
        #cat_html_in_csv()
        cat_tsv_in_csv()




if __name__ == "__main__":
    main()
