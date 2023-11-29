# Download the bgee gene expression calls
wget --timestamping --directory-prefix data/input/ https://www.bgee.org/ftp/current/download/calls/expr_calls/Homo_sapiens_expr_simple.tsv.gz

# Download NCBI Entrez protein-coding Genes
wget --timestamping --directory-prefix data/input/ https://raw.githubusercontent.com/nickzren/human-gene/main/data/output/protein_coding_gene.csv