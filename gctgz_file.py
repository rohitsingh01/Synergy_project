import pandas as pd
import gzip

with gzip.open("Anderson_1433_kinome_4373aa-4373aa_1seqs_4373aa_15.tar.gz", "rt") as f:
    df = pd.read_csv(f, sep="\t", skiprows=2)

