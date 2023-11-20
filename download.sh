# Download the bgee gene expression calls
wget --timestamping --directory-prefix data/input/ https://www.bgee.org/ftp/current/download/calls/expr_calls/Homo_sapiens_expr_simple.tsv.gz

# Download the Human gene_info.gz file of Entrez Genes
wget --timestamping --directory-prefix data/input/ ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz