import pandas as pd
import gzip

GENE_INFO_PATH = 'data/input/Homo_sapiens.gene_info.gz'
EXPRESSION_INFO_PATH = 'data/input/Homo_sapiens_expr_simple.tsv.gz'
NODE_OUTPUT_PATH = 'data/output/node_Anatomy.csv'
EDGE_OUTPUT_PATH = 'data/output/edge_Gene_expressedIn_Anatomy.csv'
NODE_COLUMNS = {'Anatomical entity ID': 'identifier', 'Anatomical entity name': 'name'}
EDGE_COLUMNS = {'NCBI Gene ID': 'source_id', 'Anatomical entity ID': 'target_id', 'is_gold_quality': 'is_gold_quality'}

def read_gene_info(file_path):
    """Read gene information, filter for protein-coding genes, and extract Ensembl IDs."""
    try:
        with gzip.open(file_path, 'rt') as file:
            gene_df = pd.read_csv(file, delimiter='\t')
            protein_coding_genes = gene_df[gene_df['type_of_gene'] == 'protein-coding'].copy()
            protein_coding_genes.loc[:, 'Ensembl'] = protein_coding_genes['dbXrefs'].str.extract('Ensembl:(ENSG[0-9]+)')
            return protein_coding_genes
    except IOError as e:
        print(f"Error reading file {file_path}: {e}")
        return pd.DataFrame()

def create_mappings(gene_df):
    """Create mappings from Ensembl ID and Symbol to GeneID."""
    ensembl_to_geneid = pd.Series(gene_df['GeneID'].values, index=gene_df['Ensembl']).to_dict()
    symbol_to_geneid = pd.Series(gene_df['GeneID'].values, index=gene_df['Symbol']).to_dict()
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
    gene_df = read_gene_info(GENE_INFO_PATH)
    ensembl_to_geneid, symbol_to_geneid = create_mappings(gene_df)
    expression_df = read_and_filter_expression_data(EXPRESSION_INFO_PATH)
    expression_df = map_gene_ids(expression_df, ensembl_to_geneid, symbol_to_geneid)

    save_to_csv(expression_df, NODE_OUTPUT_PATH, NODE_COLUMNS)
    save_to_csv(expression_df, EDGE_OUTPUT_PATH, EDGE_COLUMNS)

if __name__ == "__main__":
    main()