import pandas as pd
import gzip

def read_gene_info(file_path):
    """Read gene information and extract Ensembl IDs."""
    with gzip.open(file_path, 'rt') as file:
        gene_df = pd.read_csv(file, delimiter='\t')
        gene_df['Ensembl'] = gene_df['dbXrefs'].str.extract('Ensembl:(ENSG[0-9]+)')
    return gene_df

def create_mappings(gene_df):
    """Create Ensembl ID to GeneID and Symbol to GeneID mappings."""
    ensembl_to_geneid = pd.Series(gene_df['GeneID'].values, index=gene_df['Ensembl']).to_dict()
    symbol_to_geneid = pd.Series(gene_df['GeneID'].values, index=gene_df['Symbol']).to_dict()
    return ensembl_to_geneid, symbol_to_geneid

def read_and_filter_expression_data(file_path):
    """Read expression data and apply filters."""
    with gzip.open(file_path, 'rt') as file:
        df = pd.read_csv(file, delimiter='\t')
        df = df[(df['Expression'] == 'present') & (df['Call quality'] == 'gold quality')]
    return df

def map_gene_ids(expression_df, ensembl_to_geneid, symbol_to_geneid):
    """Map Ensembl ID to NCBI Gene ID with fallback to Symbol."""
    expression_df['NCBI Gene ID'] = expression_df['Gene ID'].map(ensembl_to_geneid)
    missing_idx = expression_df['NCBI Gene ID'].isna()
    expression_df.loc[missing_idx, 'NCBI Gene ID'] = expression_df['Gene name'].map(symbol_to_geneid)
    return expression_df.dropna(subset=['NCBI Gene ID']).astype({'NCBI Gene ID': int})

def save_to_csv(df, file_path, columns):
    """Save specified columns of DataFrame to CSV."""
    df[columns].drop_duplicates().to_csv(file_path, index=False)

def main():
    gene_info_path = 'data/input/Homo_sapiens.gene_info.gz'
    expression_info_path = 'data/input/Homo_sapiens_expr_simple.tsv.gz'
    
    gene_df = read_gene_info(gene_info_path)
    ensembl_to_geneid, symbol_to_geneid = create_mappings(gene_df)
    expression_df = read_and_filter_expression_data(expression_info_path)
    expression_df = map_gene_ids(expression_df, ensembl_to_geneid, symbol_to_geneid)

    save_to_csv(expression_df, 'data/output/node_Anatomy.csv', ['Anatomical entity ID', 'Anatomical entity name'])
    save_to_csv(expression_df, 'data/output/edge_Gene_expressedIn_Anatomy.csv', ['Anatomical entity ID', 'NCBI Gene ID'])

if __name__ == "__main__":
    main()
