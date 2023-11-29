import pandas as pd
import csv
import gzip

GENE_INFO_PATH = 'data/input/protein_coding_gene.csv'
EXPRESSION_INFO_PATH = 'data/input/Homo_sapiens_expr_simple.tsv.gz'
NODE_OUTPUT_PATH = 'data/output/csv/node_Anatomy.csv'
EDGE_OUTPUT_PATH = 'data/output/csv/edge_Gene_expressedIn_Anatomy.csv'
NODE_COLUMNS = {'Anatomical entity ID': 'identifier', 'Anatomical entity name': 'name'}
EDGE_COLUMNS = {'NCBI Gene ID': 'source_id', 'Anatomical entity ID': 'target_id', 'is_gold_quality': 'is_gold_quality'}

def load_gene_data_from_csv(file_path):
    ensembl_to_geneid = {}
    symbol_to_geneid = {}

    with open(file_path, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            gene_id = row['GeneID']
            gene_symbol = row['Symbol']
            ensembl = row['Ensembl']
            ensembl_to_geneid[ensembl] = gene_id
            symbol_to_geneid[gene_symbol] = gene_id

    return ensembl_to_geneid, symbol_to_geneid

def read_and_filter_expression_data(file_path):
    """Read expression data and add is_gold_quality column."""
    try:
        with gzip.open(file_path, 'rt') as file:
            df = pd.read_csv(file, delimiter='\t')
            df['is_gold_quality'] = df['Call quality'] == 'gold quality'
            return df[df['Expression'] == 'present']
    except IOError as e:
        print(f"Error reading file {file_path}: {e}")
        return pd.DataFrame()

def map_gene_ids(expression_df, ensembl_to_geneid, symbol_to_geneid):
    """Map Ensembl ID to NCBI Gene ID with fallback to Symbol."""
    expression_df['NCBI Gene ID'] = expression_df['Gene ID'].map(ensembl_to_geneid)
    expression_df['NCBI Gene ID'].fillna(expression_df['Gene name'].map(symbol_to_geneid), inplace=True)
    return expression_df.dropna(subset=['NCBI Gene ID']).astype({'NCBI Gene ID': int})

def save_to_csv(df, file_path, rename_columns):
    """Save specified columns of DataFrame to CSV with renamed columns."""
    df[list(rename_columns.keys())].drop_duplicates().rename(columns=rename_columns).to_csv(file_path, index=False)

def main():
    ensembl_to_geneid, symbol_to_geneid = load_gene_data_from_csv(GENE_INFO_PATH)
    expression_df = read_and_filter_expression_data(EXPRESSION_INFO_PATH)
    expression_df = map_gene_ids(expression_df, ensembl_to_geneid, symbol_to_geneid)

    save_to_csv(expression_df, NODE_OUTPUT_PATH, NODE_COLUMNS)
    save_to_csv(expression_df, EDGE_OUTPUT_PATH, EDGE_COLUMNS)

if __name__ == "__main__":
    main()