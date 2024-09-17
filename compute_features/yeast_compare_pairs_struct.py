#!/bin/python3

"""Compare structures of gene pairs, YEAST adapted version based on paralog pairs analysis"""


import os
import argparse
import subprocess
import pandas as pd
import gzip
from pathlib import Path

#==============================================================================
# Align all pairs
#==============================================================================

def align_vsDB(paralog_pairs, fn_db, mapping_infos, structures_directory, aln_type) -> None:
    """
    For all pairs, align their protein structures against a structural proteom database
    1 foldseek aln for each gene of the pair
    """
    # Building dataframes
    df_paralogs_pairs = pd.read_csv(paralog_pairs, sep=',')
    df_mapping_infos = pd.read_csv(mapping_infos, sep=',')

    # Output file building (will be completed during iteration/computation)
    r_line = '\n'
    with open('vsAll_scores.csv', 'w') as csv_file:
        csv_file.write(f'sorted_gene_pair,A1,A2,uniprot_id_A1,uniprot_id_A2,A1_A2_rank,A2_A1_rank,A1_nb_human,A2_nb_human,A1_nb_taxid,A2_nb_taxid,A1_bits,A2_bits,A1_fident,A2_fident,A1_alnlen,A2_alnlen,A1_evalue,A2_evalue,A1_alntmscore,A2_alntmscore,A1_qtmscore,A2_qtmscore,A1_ttmscore,A2_ttmscore,A1_prob,A2_prob,A1_lddt,A2_lddt{r_line}')

        # Iterate on all rows ie, paralog pairs, paralogs being genes (maybe change it later by an vectoriel operation...)
        for index, row in df_paralogs_pairs.iterrows():

            # New version, use the file with mapping - gene : uniprot
            # We can directly use gene pair
            pair = row['sorted_gene_pair']
            A1, A2 = pair.split('_')
            # Iterate on gene and get their uniprots
            for gene_id_i in [A1,A2]:
                # Get the corresponding uniprot id
                uniprot_ids = df_mapping_infos.loc[df_mapping_infos['locus'] == gene_id_i, 'uniprot']
                # If there is multiple uniprot ids
                if len(uniprot_ids) >= 2:
                    uniprot_id_i = '|'.join(uniprot_ids)
                # If one, simply squeeze it
                else:
                    uniprot_id_i = uniprot_ids.squeeze()
                if gene_id_i == A1: 
                    paralog = 'A1'
                    if ',' in f'{uniprot_id_i}':
                        uniprot_id_A1 = 'None'
                    else:
                        uniprot_id_A1 = uniprot_id_i
                elif gene_id_i == A2: 
                    paralog = 'A2'
                    if ',' in f'{uniprot_id_i}':
                        uniprot_id_A2 = 'None'
                    else:
                        uniprot_id_A2 = uniprot_id_i
            
                # Get the corresponding AlphaFold filenames
                # TODO : fix or clarify these exceptions
                try:
                    fn_struct_i = uniprot_to_AFfn(uniprot_id_i, structures_directory)
                except:
                    print(f'Issue with entrez_id, or uniprot {uniprot_id_i}. Check for AlphaFold available informations !')
                # Align against database
                try:
                    vsDB_foldseek(fn_struct_i, fn_db, row[f'{paralog}_ensembl'], aln_type)
                except:
                    print(f'Issue with AlphaFold files, check if there is a pdb file for : {uniprot_id_i}({fn_struct_i})')

            # Define a pair based on both paralogs vs all
            pair = row['sorted_gene_pair']
            A1, A2 = pair.split('_')
            """
            A distance btw both paralogs A1 and A2 of a pair is described with :
            For A1 vs all (db) :
                - A2 rank
                - Number of same species proteins better than A2
                - Number of different taxids with proteins better than A2
            For A2 vs all (db) :
                - A1 rank
                - Number of same species proteins better than A1
                - Number of different taxids with proteins better than A1
            """
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
                    A1_nb_human = len(df_better_than_A2.loc[df_better_than_A2['taxid'] == 559292])
                    # Number of different taxid that have proteins being better hits
                    A1_nb_taxid = len(set(df_better_than_A2.loc[df_better_than_A2['taxid'] != 559292]['taxid'].tolist()))
                    # Other aln features
                    A1_bits = A2_row['bits'].values[0]
                    A1_fident = A2_row['fident'].values[0]
                    A1_alnlen = A2_row['alnlen'].values[0]
                    A1_evalue = A2_row['evalue'].values[0]
                    A1_alntmscore = A2_row['alntmscore'].values[0]
                    A1_qtmscore = A2_row['qtmscore'].values[0]
                    A1_ttmscore = A2_row['ttmscore'].values[0]
                    A1_prob = A2_row['prob'].values[0]
                    A1_lddt = A2_row['lddt'].values[0]
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
                    A2_nb_human = len(df_better_than_A1.loc[df_better_than_A1['taxid'] == 559292])
                    # Number of different taxid that have proteins being better hits
                    A2_nb_taxid = len(set(df_better_than_A1.loc[df_better_than_A1['taxid'] != 559292]['taxid'].tolist()))
                    # Other aln features
                    A2_bits = A1_row['bits'].values[0]
                    A2_fident = A1_row['fident'].values[0]
                    A2_alnlen = A1_row['alnlen'].values[0]
                    A2_evalue = A1_row['evalue'].values[0]
                    A2_alntmscore = A1_row['alntmscore'].values[0]
                    A2_qtmscore = A1_row['qtmscore'].values[0]
                    A2_ttmscore = A1_row['ttmscore'].values[0]
                    A2_prob = A1_row['prob'].values[0]
                    A2_lddt = A1_row['lddt'].values[0]
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

            # Write a line in output, 1 line : 1 pair
            csv_file.write(f'{pair},{A1},{A2},{uniprot_id_A1},{uniprot_id_A2},{A1_A2_rank},{A2_A1_rank},{A1_nb_human},{A2_nb_human},{A1_nb_taxid},{A2_nb_taxid},{A1_bits},{A2_bits},{A1_fident},{A2_fident},{A1_alnlen},{A2_alnlen},{A1_evalue},{A2_evalue},{A1_alntmscore},{A2_alntmscore},{A1_qtmscore},{A2_qtmscore},{A1_ttmscore},{A2_ttmscore},{A1_prob},{A2_prob},{A1_lddt},{A2_lddt}{r_line}')

#==============================================================================
# Foldseek call
#==============================================================================

def vsDB_foldseek(fn_struct_A, fn_db, fn_out, aln_type='foldseek') -> None:
    """
    Run foldseek to align a protein structure against a structural database
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
                ],
            stdout=open(os.devnull, 'wb')
                    )
        # For some reasons, there is always some tmp artefacts ... we want to clean them
        subprocess.run([
            'rm', '-rf', 'tmp'
            ])

#==============================================================================
# Tools
#==============================================================================

def uniprot_to_AFfn(uniprot_id, struct_directory) -> str:
    """
    Get the filename of the AlphaFold structure (present in the structures directory) corresponding to the uniprot id
    """
    # If there is multiple uniprot id, we take the longest available sequence
    if '|' in str(uniprot_id):
        ids = uniprot_id.split('|')
        # Init a trivial protein information in a list: uniprot_id, sequence lenght, fasta_str
        longest = ['NA', 0, 'NA']
        for i in ids:
            print(ids,i)
            # stop there, not always
            files = struct_directory.glob('*.pdb.gz')
            for file in files:
                if i in file.stem:
                    print(i)
                    len_i = length_from_pdbgz(file)
                    print(len_i)
                    if int(len_i) > int(longest[1]): 
                        longest = [i, len_i, file]
                        print(longest)
        # If still NA, report it
        if longest[0] == 'NA': print(f'Issue with {ids} (longest not found)')
        assert longest[0] != 'NA', f'Issue with {ids} (longest not found)'
        # After iteration and length comparaison, write the longest
        return longest[2]
    
    # If there is only one uniprot id
    else:
        files = struct_directory.glob('*.pdb.gz')
        for file in files:
            if uniprot_id in file.stem:
                return file
                break

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

#==============================================================================
# Output files
#==============================================================================

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
                        help = "Csv file that describe paralog pairs (eg SGA_NxN.txt)",
                        type=str)
    parser.add_argument("mapping_uniprot",
                        help = "Csv file with mapping locus/uniprot for each gene (build from build_locus_info.py)",
                        type=str)
    parser.add_argument("structures_directory",
                        help = "Directory containing all pdb files (eg UP000005640_9606_HUMAN_V4)",
                        type=str)
    parser.add_argument("-V",
                        help = "Verbose level of foldseek (NYI)",
                        nargs='?',
                        const="3",
                        type=str)
    parser.add_argument("--aln",
                        help = "Type of alignement [foldseek,foldseek_tmalign], mandatory if aln are wanted",
                        type=str)
    parser.add_argument("--only_concat",
                        help = "Concat result file from already computed alignments, in the current directory (tsv or html from foldseek)",
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
        cat_tsv_in_csv()

    # To compute alignments
    if args.aln:
        # Load all pair info
        paralog_pairs = Path(args.paralog_pairs).resolve()
        mapping_infos = Path(args.mapping_uniprot).resolve()
        structures_directory = Path(args.structures_directory).resolve()

        # Align each paralog aigainst a database (vs all, including other species)
        if args.db:
            align_vsDB(paralog_pairs, Path(args.db).resolve(), mapping_infos, structures_directory, args.aln)

        # Concat results
        cat_tsv_in_csv()


if __name__ == "__main__":
    main()
