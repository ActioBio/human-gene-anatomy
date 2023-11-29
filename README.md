# human-gene-anatomy

The repository is dedicated to analyzing and mapping the relationships between genes and anatomical entities, highlighting where specific genes are expressed in the human body.

### Execution

```
conda env create -f environment.yml 

conda activate human-gene-anatomy

bash download.sh

python process.py
```

### Input

- protein_coding_gene.csv
  - The file is a CSV containing extracted data on protein-coding genes from the NCBI dataset.
- Homo_sapiens_expr_simple.tsv.gz
  - The file from [Bgee](https://www.bgee.org/) is a dataset providing summarized baseline presence or absence expression calls for human genes, aggregating information across various data types to indicate where each gene is actively expressed in the body.

### Output

- node_Anatomy.csv.gz
  - Contains unique anatomical entities, listing their IDs and names
- edge_Gene_expressedIn_Anatomy.csv.gz
  - Maps anatomical entity IDs to their corresponding expressed NCBI Gene IDs, illustrating the relationship between specific anatomical parts and genes.
