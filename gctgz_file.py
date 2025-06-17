import pandas as pd
import gzip

with gzip.open("CCLE_RNAseq_genes_rpkm_20180929.gct.gz", "rt") as f:
    df = pd.read_csv(f, sep="\t", skiprows=2)

print(df.head())
df.to_csv('CCLE_data.csv', index = False)